mod common;

use bam2frag::scan_plan::scan_intervals;
use serde_json::json;

#[test]
fn a_chromosome_is_scanned_in_full() {
    // Regions tile the chromosome contiguously and each is scanned whole, so
    // there is one interval covering everything. The block/window subsampler
    // this replaced always merged its windows back into exactly this.
    let cfg = common::resolved_config();

    assert_eq!(
        scan_intervals("chr1", 100_000, &cfg).unwrap(),
        [(0, 100_000)]
    );
}

#[test]
fn a_chromosome_shorter_than_one_region_is_still_scanned() {
    let cfg = common::resolved_config();

    assert_eq!(scan_intervals("chr1", 10, &cfg).unwrap(), [(0, 10)]);
}

#[test]
fn region_ids_cannot_collide_across_chromosomes() {
    // `region_stride` partitions the region-id space per chromosome. A region
    // width small enough to overflow that partition must fail loudly: silently
    // merging two chromosomes' fragments into one reservoir is unrecoverable
    // and invisible in the output.
    let mut config = common::valid_config_json();
    config["coordinate_system"]["region_window_bp"] = json!(1_000);
    config["coordinate_system"]["region_stride"] = json!(10);
    let cfg = common::load_from_value(config).1.unwrap();

    // 100_000 bp at 1 kb regions is 100 bins, well past a stride of 10.
    let error = scan_intervals("chr1", 100_000, &cfg).unwrap_err();

    assert!(
        format!("{error:#}").contains("region ids would collide"),
        "unexpected error: {error:#}"
    );
}

#[test]
fn a_region_count_exactly_at_the_stride_is_allowed() {
    // The guard is an overflow check, not an off-by-one: ids 0..=stride-1 fit.
    let mut config = common::valid_config_json();
    config["coordinate_system"]["region_window_bp"] = json!(1_000);
    config["coordinate_system"]["region_stride"] = json!(100);
    let cfg = common::load_from_value(config).1.unwrap();

    assert!(scan_intervals("chr1", 100_000, &cfg).is_ok());
}
