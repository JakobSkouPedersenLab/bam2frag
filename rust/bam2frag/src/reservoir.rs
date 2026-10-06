use crate::model::ValidatedFragmentCandidate;
use anyhow::Result;
use hashbrown::{HashMap, HashSet};
use rand::{Rng, SeedableRng};
use rand_pcg::Pcg64Mcg;
use serde::Serialize;

#[derive(Clone, Debug, Serialize)]
pub struct ReservoirStats {
    pub region_id: i64,
    pub local_pos: u16,
    pub available_count: u64,
    pub sampled_count: u64,
}

impl ReservoirStats {
    fn new(region_id: i64, local_pos: u16) -> Self {
        Self {
            region_id,
            local_pos,
            available_count: 0,
            sampled_count: 0,
        }
    }

    fn update(&mut self, _candidate: &ValidatedFragmentCandidate) {
        self.available_count += 1;
    }
}

/// One uniform natural reservoir per region; molecular differences never affect admission.
pub struct KeyedReservoir {
    capacity: usize,
    max_keys: Option<usize>,
    buckets: HashMap<(i64, u16), Vec<ValidatedFragmentCandidate>>,
    seen: HashMap<(i64, u16), u64>,
    region_keys: HashSet<(i64, u16)>,
    stats: HashMap<(i64, u16), ReservoirStats>,
    rng: Pcg64Mcg,
}

impl KeyedReservoir {
    pub fn new(stored_fragments_per_region: usize, seed: u64, max_keys: Option<usize>) -> Self {
        Self {
            capacity: stored_fragments_per_region,
            max_keys,
            buckets: HashMap::new(),
            seen: HashMap::new(),
            region_keys: HashSet::new(),
            stats: HashMap::new(),
            rng: Pcg64Mcg::seed_from_u64(seed),
        }
    }

    pub fn consider(&mut self, candidate: ValidatedFragmentCandidate) -> Result<()> {
        let region_key = (candidate.key.region_id, candidate.key.local_pos);
        if !self.stats.contains_key(&region_key) {
            if !self.region_keys.contains(&region_key)
                && self
                    .max_keys
                    .is_some_and(|max| self.region_keys.len() >= max)
            {
                return Ok(());
            }
            self.region_keys.insert(region_key);
            self.stats
                .insert(region_key, ReservoirStats::new(region_key.0, region_key.1));
        }
        self.stats.get_mut(&region_key).unwrap().update(&candidate);
        let seen = self.seen.entry(region_key).or_insert(0);
        *seen += 1;
        if self.capacity > 0 {
            let bucket = self.buckets.entry(region_key).or_default();
            if bucket.len() < self.capacity {
                bucket.push(candidate);
            } else {
                let replacement = self.rng.gen_range(0..*seen) as usize;
                if replacement < self.capacity {
                    bucket[replacement] = candidate;
                }
            }
            self.stats.get_mut(&region_key).unwrap().sampled_count = bucket.len() as u64;
        }
        Ok(())
    }

    pub fn into_parts(self) -> (Vec<ValidatedFragmentCandidate>, Vec<ReservoirStats>) {
        let mut buckets: Vec<_> = self.buckets.into_iter().collect();
        buckets.sort_by_key(|(key, _)| *key);
        let mut items = Vec::new();
        for (_, mut bucket) in buckets {
            bucket.sort_by_key(|candidate| {
                (
                    candidate.anchor_rpos0,
                    candidate.pair.left_5p_genomic,
                    candidate.pair.right_5p_genomic,
                    candidate.pair.fragment_hash,
                )
            });
            items.append(&mut bucket);
        }
        let mut stats: Vec<_> = self.stats.into_values().collect();
        stats.sort_by_key(|row| (row.region_id, row.local_pos));
        (items, stats)
    }
}
