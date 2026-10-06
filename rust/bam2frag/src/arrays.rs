use crate::{
    config::ResolvedConfig,
    model::{MaterializedFragment, ValidatedFragmentCandidate},
    reservoir::ReservoirStats,
};
use anyhow::{bail, Context, Result};
use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};
use ndarray::Array1;
use serde::Serialize;
use std::{
    collections::BTreeMap,
    fs::File,
    io::{BufReader, BufWriter, Cursor, Read, Seek, Write},
    path::Path,
};

pub const PACKED_BANK_FILENAME: &str = "fragment_bank.bin";
pub const PACKED_INDEX_FILENAME: &str = "fragment_bank_index.bin";

// These describe how a fragment is *represented*, never how it was chosen, so a
// reader can validate a bank without knowing the sampling strategy. They live
// here, beside the encoder that produces the planes they describe.
pub const ENDPOINT_DEFINITION: &str = "validated_outer_5p";
pub const CYCLE_PLANE_ENCODING: &str =
    "u8: 0=absent; 1..read_length=original sequencing cycle plus one";
pub const CYCLE_NORMALIZATION: &str = "stored_cycle_index/read_length; zero only where absent";
pub const STORED_STRUCTURAL_METADATA_FIELDS: [&str; 2] =
    ["fragment_span_bp", "fragment_orientation"];
// Every payload is stored as one zstd frame. The payload inside is byte-identical
// to an uncompressed fragment payload and the index is unchanged, with each row's
// `nbytes` the stored frame size: compression is how the bank is written, not a
// different representation. Most of a payload is zero padding outside the reads
// and monotone cycle ramps, which is why a whole-sample bank compresses about
// 12x (44 GB to about 3.6 GB on one colon sample). Level 3 already reaches that
// ratio; higher levels cost write time for little more.
pub const BANK_COMPRESSION_CODEC: &str = "zstd";
pub const BANK_COMPRESSION_LEVEL: i32 = 3;
const PAYLOAD_MAGIC: &[u8; 8] = b"FRAGBANK";
const INDEX_MAGIC: &[u8; 8] = b"FRAGIDX\0";
const FORMAT_VERSION: u16 = 1;
const PLANE_COUNT: usize = 9;
const STRUCTURAL_RECORD_BYTES: usize = 4;
const PAYLOAD_HEADER_BYTES: usize = 32 + PLANE_COUNT * 8;
pub const PAYLOAD_KINDS: usize = 1;

fn payload_magic() -> &'static [u8; 8] {
    PAYLOAD_MAGIC
}

fn index_magic() -> &'static [u8; 8] {
    INDEX_MAGIC
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
#[repr(u8)]
pub enum PayloadKind {
    Natural = 2,
}

#[derive(Clone, Copy, Debug, Default)]
pub struct PayloadIndex {
    pub count: u32,
    pub offset: u64,
    pub nbytes: u64,
}

#[derive(Clone, Debug)]
pub struct BankIndexRow {
    pub region_index: u32,
    pub region_id: i64,
    pub payloads: Vec<PayloadIndex>,
}

#[derive(Clone, Debug)]
pub struct PackedRegion {
    pub region_id: i64,
    pub payloads: Vec<CompressedPayload>,
    pub counts: Vec<u32>,
}

impl PackedRegion {
    pub fn fragment_count(&self) -> usize {
        self.counts.iter().map(|count| *count as usize).sum()
    }
}

/// One fragment payload as the bank stores it: a single zstd frame, and the
/// number of bytes it decompresses to.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CompressedPayload {
    pub frame: Vec<u8>,
    pub payload_nbytes: u64,
}

/// Compress one encoded fragment payload into one zstd frame.
///
/// Single-threaded and without a dictionary, so one payload always gives the
/// same bytes. The frame carries its content size and a checksum, so a reader
/// can verify both the length and the content of every stored payload.
pub fn compress_payload(payload: &[u8]) -> Result<CompressedPayload> {
    let mut compressor = zstd::bulk::Compressor::new(BANK_COMPRESSION_LEVEL)
        .with_context(|| format!("creating a zstd level-{BANK_COMPRESSION_LEVEL} compressor"))?;
    compressor
        .include_checksum(true)
        .context("enabling the zstd frame checksum")?;
    compressor
        .include_contentsize(true)
        .context("enabling the zstd frame content size")?;
    let frame = compressor
        .compress(payload)
        .context("compressing a fragment payload")?;
    Ok(CompressedPayload {
        frame,
        payload_nbytes: payload.len() as u64,
    })
}

/// Decompress one stored frame back to its fragment payload.
///
/// `expected_payload_nbytes` comes from the index row's fragment count, which
/// fixes a payload's exact length. The stored bytes must be exactly one zstd
/// frame, the frame must declare that length, and its checksum must match;
/// anything else is a corrupt or mis-indexed bank, not something to decode
/// around. Checking the declared length first also bounds the allocation.
pub fn decompress_payload(frame: &[u8], expected_payload_nbytes: usize) -> Result<Vec<u8>> {
    if frame.starts_with(payload_magic()) {
        bail!(
            "stored payload is an uncompressed fragment payload; banks are written as zstd frames, \
             so this sample predates compression and must be re-extracted"
        );
    }
    let frame_bytes = zstd::zstd_safe::find_frame_compressed_size(frame).map_err(|code| {
        anyhow::anyhow!(
            "stored payload is not a zstd frame: {}",
            zstd::zstd_safe::get_error_name(code)
        )
    })?;
    if frame_bytes != frame.len() {
        bail!(
            "stored payload holds a {frame_bytes}-byte zstd frame followed by {} further bytes",
            frame.len() - frame_bytes
        );
    }
    match zstd::zstd_safe::get_frame_content_size(frame) {
        Ok(Some(declared)) if declared == expected_payload_nbytes as u64 => {}
        Ok(Some(declared)) => bail!(
            "zstd frame declares {declared} payload bytes but its fragment count implies \
             {expected_payload_nbytes}"
        ),
        Ok(None) => bail!("zstd frame does not declare its payload size"),
        Err(_) => bail!("zstd frame content size cannot be read"),
    }
    let payload = zstd::bulk::decompress(frame, expected_payload_nbytes)
        .context("decompressing a fragment payload frame")?;
    if payload.len() != expected_payload_nbytes {
        bail!(
            "zstd frame decompressed to {} bytes but declares {expected_payload_nbytes}",
            payload.len()
        );
    }
    Ok(payload)
}

/// How a bank's payloads are stored, recorded identically by both commands.
#[derive(Clone, Debug, Serialize)]
pub struct BankCompressionMeta {
    pub codec: &'static str,
    pub level: i32,
    pub frame_content_size: bool,
    pub frame_checksum: bool,
}

pub fn bank_compression_meta() -> BankCompressionMeta {
    BankCompressionMeta {
        codec: BANK_COMPRESSION_CODEC,
        level: BANK_COMPRESSION_LEVEL,
        frame_content_size: true,
        frame_checksum: true,
    }
}

pub struct FragmentBankWriter {
    writer: BufWriter<File>,
    rows: Vec<BankIndexRow>,
    payload_bytes: u64,
    stored_bytes: u64,
}

impl FragmentBankWriter {
    pub fn create(output_dir: &Path, expected_regions: usize) -> Result<Self> {
        let bank_path = output_dir.join(PACKED_BANK_FILENAME);
        Ok(Self {
            writer: BufWriter::with_capacity(16 * 1024 * 1024, File::create(&bank_path)?),
            rows: Vec::with_capacity(expected_regions),
            payload_bytes: 0,
            stored_bytes: 0,
        })
    }

    /// Decompressed and stored bytes written so far.
    pub fn byte_totals(&self) -> (u64, u64) {
        (self.payload_bytes, self.stored_bytes)
    }

    pub fn write_region(&mut self, region: PackedRegion) -> Result<()> {
        let expected_payloads = 1;
        if region.payloads.len() != expected_payloads || region.counts.len() != expected_payloads {
            anyhow::bail!(
                "packed region payload count {} does not match the natural expected count {}",
                region.payloads.len(),
                expected_payloads
            );
        }
        let region_index = self.rows.len() as u32;
        let mut payloads = Vec::with_capacity(expected_payloads);
        for (kind_index, payload) in region.payloads.into_iter().enumerate() {
            if payload.frame.is_empty() {
                bail!(
                    "region {} has an empty compressed payload frame",
                    region.region_id
                );
            }
            let offset = self.writer.stream_position()?;
            self.writer.write_all(&payload.frame)?;
            self.payload_bytes += payload.payload_nbytes;
            self.stored_bytes += payload.frame.len() as u64;
            payloads.push(PayloadIndex {
                count: region.counts[kind_index],
                offset,
                nbytes: payload.frame.len() as u64,
            });
        }
        self.rows.push(BankIndexRow {
            region_index,
            region_id: region.region_id,
            payloads,
        });
        Ok(())
    }

    pub fn finish(mut self, output_dir: &Path) -> Result<Vec<BankIndexRow>> {
        self.writer.flush()?;
        write_binary_index(&output_dir.join(PACKED_INDEX_FILENAME), &self.rows)?;
        Ok(self.rows)
    }
}
pub struct ArrayMeta {
    pub n_regions: usize,
    pub kept: usize,
    pub dropped: usize,
    pub fragment_bank_file: String,
    pub fragment_bank_index_file: String,
}

pub fn write<T, D>(path: &Path, array: &ndarray::Array<T, D>) -> Result<()>
where
    T: ndarray_npy::WritableElement,
    D: ndarray::Dimension,
{
    ndarray_npy::write_npy(path, array).with_context(|| format!("writing {}", path.display()))?;
    Ok(())
}

pub fn pack_fixed_bits(values: &[u8], bits: u8) -> Result<Vec<u8>> {
    if bits == 0 || bits > 8 {
        anyhow::bail!("bits must be in 1..=8, got {bits}");
    }
    let max_value = if bits == 8 {
        u8::MAX
    } else {
        (1u8 << bits) - 1
    };
    if let Some(value) = values.iter().find(|value| **value > max_value) {
        anyhow::bail!("value {value} does not fit in {bits} bits");
    }
    let total_bits = values.len() * bits as usize;
    let mut output = vec![0u8; total_bits.div_ceil(8)];
    let mut bit_offset = 0usize;
    for value in values {
        let byte_offset = bit_offset / 8;
        let shift = bit_offset % 8;
        let expanded = (*value as u16) << shift;
        output[byte_offset] |= expanded as u8;
        if shift + bits as usize > 8 {
            output[byte_offset + 1] |= (expanded >> 8) as u8;
        }
        bit_offset += bits as usize;
    }
    Ok(output)
}

pub fn unpack_fixed_bits(bytes: &[u8], count: usize, bits: u8) -> Vec<u8> {
    let mask = if bits == 8 {
        0xffu16
    } else {
        (1u16 << bits) - 1
    };
    let mut values = Vec::with_capacity(count);
    let mut bit_offset = 0usize;
    for _ in 0..count {
        let byte_offset = bit_offset / 8;
        let shift = bit_offset % 8;
        let mut value = (bytes.get(byte_offset).copied().unwrap_or(0) as u16) >> shift;
        if shift + bits as usize > 8 {
            value |= (bytes.get(byte_offset + 1).copied().unwrap_or(0) as u16) << (8 - shift);
        }
        values.push((value & mask) as u8);
        bit_offset += bits as usize;
    }
    values
}

#[derive(Clone, Debug, PartialEq)]
pub struct DecodedPayload {
    pub region_id: i64,
    pub width: usize,
    pub ref_aln: Vec<u8>,
    pub r1_aln: Vec<u8>,
    pub r2_aln: Vec<u8>,
    pub r1_qual: Vec<u8>,
    pub r2_qual: Vec<u8>,
    pub r1_cycle: Vec<u8>,
    pub r2_cycle: Vec<u8>,
    pub atac_overlap_fractions: Vec<[f32; 3]>,
    pub fragment_span_bp: Vec<u16>,
    pub orientation_code: Vec<u8>,
}

/// The plane byte lengths of a natural payload holding `count` fragments.
fn expected_plane_lengths(count: usize, width: usize) -> Result<[usize; PLANE_COUNT]> {
    let cells = count.checked_mul(width).context("cell-count overflow")?;
    Ok([
        (cells * 4).div_ceil(8),
        (cells * 4).div_ceil(8),
        (cells * 4).div_ceil(8),
        (cells * 6).div_ceil(8),
        (cells * 6).div_ceil(8),
        cells,
        cells,
        count * 3 * std::mem::size_of::<f32>(),
        count * STRUCTURAL_RECORD_BYTES,
    ])
}

/// The exact byte length of a natural fragment payload holding `count`
/// fragments, which is what a stored frame must decompress to.
pub fn payload_nbytes(count: usize) -> Result<usize> {
    let planes = expected_plane_lengths(count, crate::config::MOLECULAR_WINDOW_BP)?;
    Ok(PAYLOAD_HEADER_BYTES + planes.iter().sum::<usize>())
}

pub fn decode_payload(bytes: &[u8]) -> Result<DecodedPayload> {
    let header_bytes = PAYLOAD_HEADER_BYTES;
    if bytes.len() < header_bytes {
        bail!("fragment payload is truncated before its header");
    }
    let mut cursor = Cursor::new(bytes);
    let mut magic = [0u8; 8];
    cursor.read_exact(&mut magic)?;
    if &magic != payload_magic() {
        bail!("invalid fragment payload magic");
    }
    let version = cursor.read_u16::<LittleEndian>()?;
    if version != FORMAT_VERSION {
        bail!("unsupported fragment payload version {version}");
    }
    if cursor.read_u8()? != PayloadKind::Natural as u8 {
        bail!("fragment payload kind is not natural");
    }
    if cursor.read_u8()? as usize != PLANE_COUNT {
        bail!("fragment payload must declare {PLANE_COUNT} planes");
    }
    let count = cursor.read_u32::<LittleEndian>()? as usize;
    let width = cursor.read_u16::<LittleEndian>()? as usize;
    let numeric_features = cursor.read_u16::<LittleEndian>()? as usize;
    let cycle_bits = cursor.read_u16::<LittleEndian>()?;
    let structural_bytes = cursor.read_u16::<LittleEndian>()? as usize;
    let region_id = cursor.read_i64::<LittleEndian>()?;
    if width != crate::config::MOLECULAR_WINDOW_BP
        || numeric_features != 3
        || cycle_bits != 8
        || structural_bytes != STRUCTURAL_RECORD_BYTES
    {
        bail!("invalid width/feature/record declaration");
    }
    let mut plane_lengths = [0usize; PLANE_COUNT];
    for length in &mut plane_lengths {
        *length = cursor.read_u64::<LittleEndian>()? as usize;
    }
    let cells = count.checked_mul(width).context("cell-count overflow")?;
    let expected = expected_plane_lengths(count, width)?;
    if plane_lengths != expected {
        bail!("plane sizes {plane_lengths:?} do not match {expected:?}");
    }
    if header_bytes + plane_lengths.iter().sum::<usize>() != bytes.len() {
        bail!("fragment payload byte length does not match its plane table");
    }
    let mut offset = header_bytes;
    let mut planes = Vec::with_capacity(PLANE_COUNT);
    for length in plane_lengths {
        planes.push(&bytes[offset..offset + length]);
        offset += length;
    }
    let ref_aln = unpack_fixed_bits(planes[0], cells, 4);
    let r1_aln = unpack_fixed_bits(planes[1], cells, 4);
    let r2_aln = unpack_fixed_bits(planes[2], cells, 4);
    let r1_qual = unpack_fixed_bits(planes[3], cells, 6);
    let r2_qual = unpack_fixed_bits(planes[4], cells, 6);
    let r1_cycle = planes[5].to_vec();
    let r2_cycle = planes[6].to_vec();
    let numeric = planes[7];
    let structural = planes[8];
    let mut atac_overlap_fractions = Vec::with_capacity(count);
    let mut numeric_cursor = Cursor::new(numeric);
    for _ in 0..count {
        let values = [
            numeric_cursor.read_f32::<LittleEndian>()?,
            numeric_cursor.read_f32::<LittleEndian>()?,
            numeric_cursor.read_f32::<LittleEndian>()?,
        ];
        if values
            .iter()
            .any(|value| !value.is_finite() || !(0.0..=1.0).contains(value))
        {
            bail!("ATAC overlap fraction is outside [0,1]");
        }
        atac_overlap_fractions.push(values);
    }
    let mut fragment_span_bp = Vec::with_capacity(count);
    let mut orientation_code = Vec::with_capacity(count);
    let mut structural_cursor = Cursor::new(structural);
    for _ in 0..count {
        let span = structural_cursor.read_u16::<LittleEndian>()?;
        let orientation = structural_cursor.read_u8()?;
        let reserved = structural_cursor.read_u8()?;
        if !(30..=299).contains(&span) || orientation > 1 || reserved != 0 {
            bail!("invalid structural record");
        }
        fragment_span_bp.push(span);
        orientation_code.push(orientation);
    }
    if ref_aln.iter().any(|base| *base > 4) {
        bail!("REF plane contains a non N/A/C/G/T code");
    }
    if r1_aln.iter().chain(&r2_aln).any(|base| *base > 4) {
        bail!("read plane contains a non N/A/C/G/T code");
    }
    for row in 0..count {
        let start = row * width;
        let end = start + width;
        let left_cut = (width - fragment_span_bp[row] as usize) / 2;
        let right_cut = left_cut + fragment_span_bp[row] as usize - 1;
        let (r1_endpoint, r2_endpoint, r1_step, r2_step) = if orientation_code[row] == 0 {
            (left_cut, right_cut, 1i16, -1i16)
        } else {
            (right_cut, left_cut, -1i16, 1i16)
        };
        if r1_cycle[start + r1_endpoint] != 1 || r2_cycle[start + r2_endpoint] != 1 {
            bail!("cycle 1 is not placed at both centred validated 5-prime endpoints");
        }
        for (name, bases, qualities, cycles, step) in [
            ("R1", &r1_aln, &r1_qual, &r1_cycle, r1_step),
            ("R2", &r2_aln, &r2_qual, &r2_cycle, r2_step),
        ] {
            let positions = cycles[start..end]
                .iter()
                .enumerate()
                .filter(|(_, cycle)| **cycle > 0)
                .collect::<Vec<_>>();
            for (position, cycle) in cycles[start..end].iter().enumerate() {
                if *cycle == 0 && (bases[start + position] != 0 || qualities[start + position] != 0)
                {
                    bail!("{name} has base/quality content where its cycle is absent");
                }
            }
            let has_cycle_one = positions.iter().any(|(_, cycle)| **cycle == 1);
            if positions.is_empty() || !has_cycle_one {
                bail!("{name} cycle plane has no cycle 1");
            }
            for pair in positions.windows(2) {
                let genomic_delta = pair[1].0 as i16 - pair[0].0 as i16;
                let cycle_delta = *pair[1].1 as i16 - *pair[0].1 as i16;
                if genomic_delta != 1 || cycle_delta != step {
                    bail!("{name} cycles are not contiguous in sequencing-cycle order");
                }
            }
        }
    }
    Ok(DecodedPayload {
        region_id,
        width,
        ref_aln,
        r1_aln,
        r2_aln,
        r1_qual,
        r2_qual,
        r1_cycle,
        r2_cycle,
        atac_overlap_fractions,
        fragment_span_bp,
        orientation_code,
    })
}

/// Validate a bank against its bank index, one stored frame at a time.
///
/// Frames are read sequentially rather than loading the bank, and each is
/// decompressed and decoded with the full fragment payload validation. Frames
/// must be contiguous from offset zero with nothing after the last one.
pub fn validate_bank_files(bank_path: &Path, index_path: &Path) -> Result<usize> {
    let bank_len = std::fs::metadata(bank_path)
        .with_context(|| format!("reading {}", bank_path.display()))?
        .len() as usize;
    let mut bank = BufReader::new(
        File::open(bank_path).with_context(|| format!("opening {}", bank_path.display()))?,
    );
    let index =
        std::fs::read(index_path).with_context(|| format!("reading {}", index_path.display()))?;
    if index.len() < 16 {
        bail!("bank index is truncated");
    }
    let mut cursor = Cursor::new(index.as_slice());
    let mut magic = [0u8; 8];
    cursor.read_exact(&mut magic)?;
    if &magic != index_magic() {
        bail!("invalid bank index magic");
    }
    let version = cursor.read_u16::<LittleEndian>()?;
    if version != FORMAT_VERSION {
        bail!("unsupported bank index version {version}");
    }
    if cursor.read_u16::<LittleEndian>()? != PAYLOAD_KINDS as u16 {
        bail!("bank index must declare one natural payload kind");
    }
    let rows = cursor.read_u32::<LittleEndian>()? as usize;
    let expected_index_bytes = 16 + rows * (4 + 8 + 4 + 8 + 8);
    if index.len() != expected_index_bytes {
        bail!("bank index length does not match its row count");
    }
    let mut previous_end = 0usize;
    let mut frame = Vec::new();
    for expected_row in 0..rows {
        if cursor.read_u32::<LittleEndian>()? as usize != expected_row {
            bail!("bank index region indices are not contiguous");
        }
        let region_id = cursor.read_i64::<LittleEndian>()?;
        let count = cursor.read_u32::<LittleEndian>()? as usize;
        let offset = cursor.read_u64::<LittleEndian>()? as usize;
        let nbytes = cursor.read_u64::<LittleEndian>()? as usize;
        if offset != previous_end
            || nbytes == 0
            || offset.checked_add(nbytes).is_none_or(|end| end > bank_len)
        {
            bail!("bank index payload offset/length is outside the packed bank");
        }
        frame.resize(nbytes, 0);
        bank.read_exact(&mut frame)
            .with_context(|| format!("reading bank index row {expected_row} payload frame"))?;
        let payload = decompress_payload(&frame, payload_nbytes(count)?)
            .with_context(|| format!("bank index row {expected_row}"))?;
        let decoded =
            decode_payload(&payload).with_context(|| format!("bank index row {expected_row}"))?;
        if decoded.region_id != region_id || decoded.fragment_span_bp.len() != count {
            bail!("bank index row does not match its fragment payload header");
        }
        previous_end = offset + nbytes;
    }
    if previous_end != bank_len {
        bail!("fragment bank has unindexed trailing bytes");
    }
    Ok(rows)
}

struct BitWriter {
    bytes: Vec<u8>,
    bit_offset: usize,
}

impl BitWriter {
    fn with_capacity(values: usize, bits: u8) -> Self {
        Self {
            bytes: Vec::with_capacity((values * bits as usize).div_ceil(8)),
            bit_offset: 0,
        }
    }

    fn push(&mut self, value: u8, bits: u8) -> Result<()> {
        if bits == 0 || bits > 8 {
            anyhow::bail!("bits must be in 1..=8, got {bits}");
        }
        let max_value = if bits == 8 {
            u8::MAX
        } else {
            (1u8 << bits) - 1
        };
        if value > max_value {
            anyhow::bail!("value {value} does not fit in {bits} bits");
        }
        let byte_offset = self.bit_offset / 8;
        let shift = self.bit_offset % 8;
        let extra_byte = if shift + bits as usize > 8 { 1 } else { 0 };
        while self.bytes.len() <= byte_offset + extra_byte {
            self.bytes.push(0);
        }
        let expanded = (value as u16) << shift;
        self.bytes[byte_offset] |= expanded as u8;
        if shift + bits as usize > 8 {
            self.bytes[byte_offset + 1] |= (expanded >> 8) as u8;
        }
        self.bit_offset += bits as usize;
        Ok(())
    }

    fn finish(self) -> Vec<u8> {
        self.bytes
    }
}

struct PayloadBuilder {
    kind: PayloadKind,
    count: u32,
    dna: [BitWriter; 3],
    quality: [BitWriter; 2],
    cycle: [Vec<u8>; 2],
    numeric: Vec<u8>,
    structural: Vec<u8>,
}

impl PayloadBuilder {
    fn new(kind: PayloadKind, expected_fragments: usize, width: usize) -> Self {
        Self {
            kind,
            count: 0,
            dna: [
                BitWriter::with_capacity(expected_fragments * width, 4),
                BitWriter::with_capacity(expected_fragments * width, 4),
                BitWriter::with_capacity(expected_fragments * width, 4),
            ],
            quality: [
                BitWriter::with_capacity(expected_fragments * width, 6),
                BitWriter::with_capacity(expected_fragments * width, 6),
            ],
            cycle: [
                Vec::with_capacity(expected_fragments * width),
                Vec::with_capacity(expected_fragments * width),
            ],
            numeric: Vec::with_capacity(expected_fragments * 12),
            structural: Vec::with_capacity(expected_fragments * STRUCTURAL_RECORD_BYTES),
        }
    }

    fn push_scratch(
        &mut self,
        scratch: &crate::features::MaterializeScratch,
        scalars: crate::features::MaterializedScalars,
        cfg: &ResolvedConfig,
    ) -> Result<()> {
        let width = cfg.typed.feature_sets.fragment_sequence.seq_window_w;
        for values in [&scratch.ref_aln, &scratch.r1_aln, &scratch.r2_aln] {
            if values.len() != width {
                anyhow::bail!("DNA width {} != {}", values.len(), width);
            }
        }
        if scratch.r1_qual.len() != width
            || scratch.r2_qual.len() != width
            || scratch.r1_cycle.len() != width
            || scratch.r2_cycle.len() != width
        {
            anyhow::bail!("quality/cycle width does not match configured width {width}");
        }
        for (writer, values) in
            self.dna
                .iter_mut()
                .zip([&scratch.ref_aln, &scratch.r1_aln, &scratch.r2_aln])
        {
            for &value in values {
                if value > 4 {
                    bail!("base code {value} is outside N/A/C/G/T");
                }
                writer.push(value, 4)?;
            }
        }
        for (writer, values) in self
            .quality
            .iter_mut()
            .zip([&scratch.r1_qual, &scratch.r2_qual])
        {
            for &value in values {
                if value > 63 {
                    bail!("quality value {value} exceeds the 6-bit range");
                }
                writer.push(value, 6)?;
            }
        }
        self.cycle[0].extend_from_slice(&scratch.r1_cycle);
        self.cycle[1].extend_from_slice(&scratch.r2_cycle);
        if !(30..=299).contains(&scalars.fragment_span_bp) || scalars.orientation_code > 1 {
            bail!("span/orientation is outside the validated contract");
        }
        for value in scalars.atac_overlap_fractions {
            if !value.is_finite() || !(0.0..=1.0).contains(&value) {
                bail!("ATAC overlap fraction {value} is outside [0,1]");
            }
            self.numeric.write_f32::<LittleEndian>(value)?;
        }

        self.structural
            .write_u16::<LittleEndian>(scalars.fragment_span_bp)?;
        self.structural.write_u8(scalars.orientation_code)?;
        self.structural.write_u8(0)?;
        self.count = self.count.saturating_add(1);
        Ok(())
    }

    fn finish(self, region_id: i64, cfg: &ResolvedConfig) -> Result<Vec<u8>> {
        let width = cfg.typed.feature_sets.fragment_sequence.seq_window_w;
        let mut planes: [Vec<u8>; PLANE_COUNT] = Default::default();
        let [dna0, dna1, dna2] = self.dna;
        let [qual0, qual1] = self.quality;
        planes[0] = dna0.finish();
        planes[1] = dna1.finish();
        planes[2] = dna2.finish();
        planes[3] = qual0.finish();
        planes[4] = qual1.finish();
        planes[5] = self.cycle[0].clone();
        planes[6] = self.cycle[1].clone();
        planes[7] = self.numeric;
        planes[8] = self.structural;

        let mut output = Vec::new();
        output.extend_from_slice(payload_magic());
        output.write_u16::<LittleEndian>(FORMAT_VERSION)?;
        output.write_u8(self.kind as u8)?;
        output.write_u8(PLANE_COUNT as u8)?;
        output.write_u32::<LittleEndian>(self.count)?;
        output.write_u16::<LittleEndian>(width as u16)?;
        output.write_u16::<LittleEndian>(3)?;
        output.write_u16::<LittleEndian>(8)?; // cycle-index bit width
        output.write_u16::<LittleEndian>(STRUCTURAL_RECORD_BYTES as u16)?;
        output.write_i64::<LittleEndian>(region_id)?;
        for plane in &planes {
            output.write_u64::<LittleEndian>(plane.len() as u64)?;
        }
        for plane in planes {
            output.extend_from_slice(&plane);
        }
        Ok(output)
    }
}

/// Direct region-level SoA packer used by the production hot path. It accepts
/// one reusable materialization scratch buffer at a time and immediately
/// appends the fragment to the final compact planes.
pub struct RegionPayloadBuilder {
    region_id: i64,
    payload: PayloadBuilder,
}

impl RegionPayloadBuilder {
    pub fn new(region_id: i64, expected_fragments: usize, cfg: &ResolvedConfig) -> Self {
        let width = cfg.typed.feature_sets.fragment_sequence.seq_window_w;
        Self {
            region_id,
            payload: PayloadBuilder::new(PayloadKind::Natural, expected_fragments, width),
        }
    }

    pub fn push(
        &mut self,
        scratch: &crate::features::MaterializeScratch,
        scalars: crate::features::MaterializedScalars,
        cfg: &ResolvedConfig,
    ) -> Result<()> {
        self.payload.push_scratch(scratch, scalars, cfg)
    }

    /// Encode and compress the region. This runs on the worker that built the
    /// region, so compression is spread over the extraction threads rather than
    /// serialized behind the single bank writer.
    pub fn finish(self, cfg: &ResolvedConfig) -> Result<PackedRegion> {
        let count = self.payload.count;
        let payload = self.payload.finish(self.region_id, cfg)?;
        Ok(PackedRegion {
            region_id: self.region_id,
            payloads: vec![compress_payload(&payload)?],
            counts: vec![count],
        })
    }
}

/// Encode one region/stream fragment payload directly from materialized
/// fragments. DNA and quality planes are bit-packed while iterating the
/// fragments, avoiding the previous count*width temporary byte arrays.
#[doc(hidden)]
pub fn encode_payload(
    region_id: i64,
    kind: PayloadKind,
    fragments: &[&MaterializedFragment],
    cfg: &ResolvedConfig,
) -> Result<Vec<u8>> {
    let count = fragments.len();
    let width = cfg.typed.feature_sets.fragment_sequence.seq_window_w;
    let mut planes: [Vec<u8>; PLANE_COUNT] = Default::default();

    for field in 0..3 {
        let mut writer = BitWriter::with_capacity(count * width, 4);
        for fragment in fragments {
            let values = match field {
                0 => &fragment.ref_aln,
                1 => &fragment.r1_aln,
                2 => &fragment.r2_aln,
                _ => unreachable!(),
            };
            if values.len() != width {
                anyhow::bail!("DNA width {} != {}", values.len(), width);
            }
            for &value in values {
                if value > 4 {
                    bail!("base code {value} is outside N/A/C/G/T");
                }
                writer.push(value, 4)?;
            }
        }
        planes[field] = writer.finish();
    }

    for field in 0..2 {
        let mut writer = BitWriter::with_capacity(count * width, 6);
        for fragment in fragments {
            let values = if field == 0 {
                &fragment.r1_qual
            } else {
                &fragment.r2_qual
            };
            if values.len() != width {
                anyhow::bail!("quality width does not match configured width {width}");
            }
            for &value in values {
                if value > 63 {
                    bail!("quality value {value} exceeds the 6-bit range");
                }
                writer.push(value, 6)?;
            }
        }
        planes[3 + field] = writer.finish();
    }

    for fragment in fragments {
        if fragment.r1_cycle.len() != width || fragment.r2_cycle.len() != width {
            bail!("cycle width does not match configured width {width}");
        }
        planes[5].extend_from_slice(&fragment.r1_cycle);
        planes[6].extend_from_slice(&fragment.r2_cycle);
    }

    for fragment in fragments {
        for value in fragment.atac_overlap_fractions {
            if !value.is_finite() || !(0.0..=1.0).contains(&value) {
                bail!("ATAC overlap fraction {value} is outside [0,1]");
            }
            planes[7].write_f32::<LittleEndian>(value)?;
        }
    }

    for fragment in fragments {
        if !(30..=299).contains(&fragment.fragment_span_bp) || fragment.orientation_code > 1 {
            bail!("span/orientation is outside the validated contract");
        }
        planes[8].write_u16::<LittleEndian>(fragment.fragment_span_bp)?;
        planes[8].write_u8(fragment.orientation_code)?;
        planes[8].write_u8(0)?;
    }

    let mut output = Vec::new();
    output.extend_from_slice(payload_magic());
    output.write_u16::<LittleEndian>(FORMAT_VERSION)?;
    output.write_u8(kind as u8)?;
    output.write_u8(PLANE_COUNT as u8)?;
    output.write_u32::<LittleEndian>(count as u32)?;
    output.write_u16::<LittleEndian>(width as u16)?;
    output.write_u16::<LittleEndian>(3)?;
    output.write_u16::<LittleEndian>(8)?;
    output.write_u16::<LittleEndian>(STRUCTURAL_RECORD_BYTES as u16)?;
    output.write_i64::<LittleEndian>(region_id)?;
    for plane in &planes {
        output.write_u64::<LittleEndian>(plane.len() as u64)?;
    }
    for plane in planes {
        output.extend_from_slice(&plane);
    }
    Ok(output)
}

/// Pack one complete region into its final compressed natural fragment payload.
/// The caller can drop all expanded `MaterializedFragment`s immediately after
/// this function returns.
pub fn pack_region(
    region_id: i64,
    fragments: &[MaterializedFragment],
    cfg: &ResolvedConfig,
) -> Result<PackedRegion> {
    for fragment in fragments {
        if fragment.key.region_id != region_id {
            anyhow::bail!(
                "fragment region_id {} does not match packed region {}",
                fragment.key.region_id,
                region_id
            );
        }
    }

    let mut natural: Vec<&MaterializedFragment> = fragments.iter().collect();
    natural.sort_by_key(|fragment| {
        (
            fragment.anchor_rpos0,
            fragment.fragment_span_bp,
            fragment.orientation_code,
        )
    });
    let payload = encode_payload(region_id, PayloadKind::Natural, &natural, cfg)?;
    Ok(PackedRegion {
        region_id,
        payloads: vec![compress_payload(&payload)?],
        counts: vec![natural.len() as u32],
    })
}

fn write_binary_index(path: &Path, rows: &[BankIndexRow]) -> Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    writer.write_all(index_magic())?;
    writer.write_u16::<LittleEndian>(FORMAT_VERSION)?;
    let payload_kinds = 1u16;
    writer.write_u16::<LittleEndian>(payload_kinds)?;
    writer.write_u32::<LittleEndian>(rows.len() as u32)?;
    for row in rows {
        if row.payloads.len() != payload_kinds as usize {
            anyhow::bail!(
                "index row has {} payloads but natural payloads require {}",
                row.payloads.len(),
                payload_kinds
            );
        }
        writer.write_u32::<LittleEndian>(row.region_index)?;
        writer.write_i64::<LittleEndian>(row.region_id)?;
        for payload in &row.payloads {
            writer.write_u32::<LittleEndian>(payload.count)?;
            writer.write_u64::<LittleEndian>(payload.offset)?;
            writer.write_u64::<LittleEndian>(payload.nbytes)?;
        }
    }
    writer.flush()?;
    Ok(())
}

/// Return active region keys in the exact deterministic order used by the
/// bank index and all region-level arrays.
pub fn region_keys_from_candidates(
    candidates: &[ValidatedFragmentCandidate],
    cfg: &ResolvedConfig,
) -> Vec<i64> {
    candidates
        .iter()
        .map(|candidate| candidate.key.region_id)
        .collect::<std::collections::BTreeSet<_>>()
        .into_iter()
        .collect()
}

/// Write region metadata and SSL targets directly from compact retained
/// candidates and reservoir summaries, before fragment sequence materialization.
/// This lets the subsequent materialize/pack stage consume candidates region by
/// region without retaining sample-wide expanded fragments.
/// Write the only region-level sidecar used downstream: the number of eligible
/// fragments seen before reservoir sampling. Rows exactly follow the bank index.
pub fn write_region_available_counts(
    output_dir: &Path,
    index_rows: &[BankIndexRow],
    stats: &[ReservoirStats],
) -> Result<ArrayMeta> {
    let lookup: BTreeMap<_, _> = stats.iter().map(|row| (row.region_id, row)).collect();
    let mut counts = Vec::with_capacity(index_rows.len());
    let mut kept = 0usize;
    for row in index_rows {
        let key = row.region_id;
        let stat = lookup.get(&key).ok_or_else(|| {
            anyhow::anyhow!("missing reservoir statistics for bank index key {key:?}")
        })?;
        if stat.available_count < stat.sampled_count {
            bail!("available count is smaller than sampled count for bank index key {key:?}");
        }
        counts.push(stat.available_count);
        kept += stat.sampled_count as usize;
    }
    write(
        &output_dir.join("region_available_count.npy"),
        &Array1::from_vec(counts),
    )?;

    Ok(ArrayMeta {
        n_regions: index_rows.len(),
        kept,
        dropped: 0,
        fragment_bank_file: PACKED_BANK_FILENAME.into(),
        fragment_bank_index_file: PACKED_INDEX_FILENAME.into(),
    })
}
