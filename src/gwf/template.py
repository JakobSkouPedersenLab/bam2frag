from pathlib import Path

from gwf import AnonymousTarget


def make_fragments(bam_file_path,
                   reference_path,
                   config_path,
                   chromatin_dir=None,
                   genomic_blacklist_file=None,
                   donor_blacklist_file=None,
                   snp_blacklist_path=None,
                   lowpower_blacklist_path=None,
                   problematic_blacklist_path=None,
                   sample_id=None,
                   patient_id=None,
                   account=None,
                   output_dir=None,
                   cores=12,
                   memory="50g",
                   walltime="24:00:00",
                   threads=None,
                   chromosome_batch_size=12,
                   seed=0):
    """
    Run the per-BAM cfDNA fragment extraction and write one sample directory.
    """
    if not sample_id:
        raise ValueError("sample_id is required")
    if not output_dir:
        raise ValueError("output_dir is required")

    bam = Path(bam_file_path)
    reference = Path(reference_path)
    sample_dir = Path(output_dir)

    inputs = [
        str(bam),
        f"{bam}.bai",
        str(reference),
        f"{reference}.fai",
        str(config_path),
    ]

    optional_paths = {
        "--chromatin_dir": chromatin_dir,
        "--genomic_blacklist_file": genomic_blacklist_file,
        "--donor_blacklist_file": donor_blacklist_file,
        "--snp_blacklist_path": snp_blacklist_path,
        "--lowpower_blacklist_path": lowpower_blacklist_path,
        "--problematic_blacklist_path": problematic_blacklist_path,
    }

    # Add each optional resource only if present
    for path in optional_paths.values():
        if path:
            inputs.append(str(path))

    # `sample_meta.json` is written only after every array and sidecar is in
    # place, so gwf treats a sample as done only when extraction fully finished.
    outputs = [
        str(sample_dir / "sample_meta.json"),
        str(sample_dir / "fragment_bank.bin"),
        str(sample_dir / "fragment_bank_index.bin"),
        str(sample_dir / "region_available_count.npy"),
        str(sample_dir / "extraction_stats.csv"),
    ]

    options = dict(cores=str(cores), memory=memory, walltime=walltime)

    # Cluster account, where the backend needs one. The caller supplies it; none is hardcoded.
    if account:
        options["account"] = account

    optional_args = dict(optional_paths)
    optional_args["--patient_id"] = patient_id

    # Only pass arguments that are set: an empty one would render as a bare flag and
    # argparse would take the following flag to be its value.
    argv = [str(bam), "$TEMP/ref.fa", str(config_path)]
    argv += [f"{flag} {value}" for flag, value in optional_args.items() if value]
    argv += [
        f"--sample_id {sample_id}",
        f"--output_dir {sample_dir}",
        f"--threads {threads if threads is not None else cores}",
        f"--chromosome_batch_size {chromosome_batch_size}",
        f"--seed {seed}",
        "--overwrite",
        "--verbose",
    ]
    args_block = " \\\n        ".join(argv)

    spec = f"""
    set -euo pipefail

    mkdir -p $TEMP/ {sample_dir.parent}

    # Every extraction thread reads the reference through faidx for the whole
    # run, so a cohort-wide submission would hammer one file on shared storage.
    # Copy it to node-local scratch once and read from there instead.
    cp {reference} $TEMP/ref.fa
    cp {reference}.fai $TEMP/ref.fa.fai

    # Clear any partial directory left behind by an interrupted run.
    rm -rf {sample_dir}.incomplete.*

    bam2frag \\
        {args_block}

    rm -f $TEMP/ref.fa $TEMP/ref.fa.fai
    """

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
