//! The symmetric-difference distance between two alignments.
//!
//! Given two [`HomologyView`]s built against one shared [`Registry`], the distance is
//!
//! ```text
//!            sum over residues r of  |A(r) Δ B(r)|
//!   d(a,b) = ---------------------------------------
//!            sum over residues r of  |A(r)| + |B(r)|
//! ```
//!
//! where `A(r)` is the set of elements sharing `r`'s column in alignment `a`, and `B(r)`
//! the same in `b`. Both sums range over the residues that exist in *both* views — a
//! `None` slot means "this sequence has no residue at that position", i.e. past the end
//! of its residues, and contributes nothing to either sum.
//!
//! # What changed in Stage 4
//!
//! Two things, both from `CODE_REVIEW.md` §0:
//!
//! - **The iteration bound is `registry.len()`**, not the first alignment's sequence
//!   count. `main.rs` used to pass `msa_a.num_seqs`, so a comparison against an
//!   alignment with more sequences was silently computed over a truncated overlap and
//!   `d(a, b) != d(b, a)`. [`Registry::for_pair`] now rejects mismatched name sets
//!   outright, so the sequence dimension is shared by construction.
//! - **The denominator is `|A| + |B|`**, not `2·|A|`. Taking both counts from the A side
//!   made the ratio depend on which alignment was passed first whenever the two sides'
//!   set sizes could differ.
//!
//! Worth writing down, because it is not obvious and it is the kind of thing that looks
//! like a bug later: on any pair that gets this far, **the two denominators are equal**.
//! A residue's homology set is its whole column minus itself, and every element in a
//! column carries a distinct `seq`, so no two collapse in the hash set and
//! `|A(r)| = |B(r)| = num_seqs - 1` for *every* residue. The registry has already
//! forced both alignments to the same sequence count. So the fix is numerically inert
//! today. `|A| + |B|` is the definition of the metric; `2·|A|` merely happened to agree
//! with it, and would stop agreeing under any move to residues-only or gap-filtered
//! sets. `every_homology_set_has_size_num_seqs_minus_one` pins the invariant.
//!
//! # What changed in the `CODE_REVIEW.md` §2 fix
//!
//! The sets this reads are now **columns**, shared by every residue in them, rather
//! than one materialised homology set per residue. That took live memory from
//! O(num_seqs² × width) per alignment to O(num_seqs × width). See [`HomologyView`].
//!
//! # And what changed after it
//!
//! The columns are no longer `HashSet`s at all, but runs of `u32` codes, and
//! `|A Δ B|` is counted by scanning two of them in step — see
//! [`symmetric_difference_size`] for why that is the same number. The per-slot
//! contributions are also reduced in place instead of being collected into a `Vec` and
//! folded afterwards. Neither changes a result: on a 500 x 5000 measurement pair the
//! distance is identical to the last digit, while peak RSS goes from 330 MB to 49 MB and
//! the run from 3.5 s to 0.7 s.
//!
//! It is exact, not an approximation. A column includes the residue itself, and the
//! homology set is the column minus that one element `x = {sequence, position, gap:
//! false}`. Because a residue's identity does not depend on the gaps around it, `x` is
//! the same value in both alignments, so it is in both columns and therefore in neither
//! symmetric difference: `(Ca \ {x}) Δ (Cb \ {x}) == Ca Δ Cb`. The numerator needs no
//! adjustment at all; only the denominator subtracts the one element per side.
//! `sharing_the_column_sets_computes_exactly_the_materialised_distance` checks this
//! against a literal implementation over every 2x3 and 3x2 alignment.

use crate::homology::{HomologyView, Registry};
use anyhow::{Result, bail};
use colored::Colorize;
use itertools::Itertools;
use rayon::prelude::*;

/// The per-residue contribution to the distance: one symmetric-difference count and the
/// two set sizes that go under the line.
#[derive(Clone, Copy)]
struct ResidueContribution {
    /// `|A Δ B|` — the numerator's share.
    symmetric_difference: usize,
    /// `|A| + |B|` — the denominator's share.
    set_sizes: usize,
}

impl ResidueContribution {
    /// What a slot contributes when at least one side has no residue there: nothing.
    /// Also the identity for [`ResidueContribution::combine`].
    const NONE: ResidueContribution = ResidueContribution {
        symmetric_difference: 0,
        set_sizes: 0,
    };

    /// Adds two contributions. Both fields are integers, so this is associative and
    /// commutative *exactly* — which is why the distance does not depend on how rayon
    /// happens to split the work. A `f64` accumulator would have made the result
    /// thread-count dependent in the last bits.
    fn combine(self, other: ResidueContribution) -> ResidueContribution {
        ResidueContribution {
            symmetric_difference: self.symmetric_difference + other.symmetric_difference,
            set_sizes: self.set_sizes + other.set_sizes,
        }
    }
}

/// `|A Δ B|` for two columns of the same pair of views.
///
/// # Why this is a scan and not a set operation
///
/// Both slices hold one element per sequence, indexed by registry index, and an
/// element's identity is `(seq, position, gap)`. Two elements from *different*
/// sequences can therefore never be equal, so the intersection can only pair `a[s]`
/// with `b[s]`:
///
/// ```text
///   |A ∩ B| = #{s : a[s] = b[s]}
///   |A Δ B| = |A| + |B| - 2|A ∩ B| = 2n - 2·#{s : a[s] = b[s]}
///           = 2·#{s : a[s] ≠ b[s]}
/// ```
///
/// With `s` fixed on both sides the `seq` component is common, which is exactly why the
/// stored codes can omit it (see `homology::residue_code`). So the whole set machinery
/// collapses to counting disagreements down two slices — no hashing, no allocation, and
/// no `HashSet` to have built in the first place.
///
/// `sharing_the_column_sets_computes_exactly_the_materialised_distance` checks this
/// against literal `HashSet` symmetric differences over every 2x3 and 3x2 alignment.
fn symmetric_difference_size(column_a: &[u32], column_b: &[u32]) -> usize {
    debug_assert_eq!(
        column_a.len(),
        column_b.len(),
        "both views are sized by the shared registry"
    );
    2 * column_a
        .iter()
        .zip(column_b)
        .filter(|(a, b)| a != b)
        .count()
}

/// The symmetric-difference distance between two homology views built against
/// `registry`.
///
/// Returns `Err` if the two views have no residues in common to compare — a pair of
/// single-sequence alignments (where every homology set is empty), or alignments made
/// entirely of gaps. The ratio is `0/0` there; the old code returned `NaN`, which would
/// be written into the CSV as a confident-looking result.
pub fn compute_symmetric_difference(
    view_a: &HomologyView,
    view_b: &HomologyView,
    registry: &Registry,
) -> Result<f64> {
    log::info!("{}", "Computing Symmetric Difference".bold().purple());

    // Both views are indexed by registry index and sized by `registry.len()`, so the
    // sequence dimension needs no reconciling. The residue dimension is each
    // alignment's own width, which may differ between the two; take the larger and let
    // `column_of` return `None` past the end of the shorter.
    let width = view_a.width().max(view_b.width());

    // `ResidueContribution::NONE` means the slot pair contributed nothing: at least one
    // of the two sequences has no residue at this position. That is the ordinary case
    // for every position past a sequence's residue count, not an error.
    //
    // Reduced as it goes rather than collected first. The previous shape built a
    // `Vec<Option<ResidueContribution>>` of one 24-byte entry per (sequence, position)
    // slot — 60 MB on a 500 x 5000 pair — and then immediately folded it to two numbers.
    let totals = (0..registry.len())
        .cartesian_product(0..width)
        .par_bridge()
        .map(|(sequence, position)| {
            // Note there is no fallback arm reporting a missing sequence here. The old
            // code logged a warning in that position claiming the sequence was absent
            // from the hash set; it could not fire for the stated reason (the index was
            // bounded by A's own length), it fired once per column when it did fire, and
            // `Registry::for_pair` now makes a genuinely missing sequence impossible.
            //
            // These are the residue's *column*, which includes the residue itself; the
            // homology set proper is the column minus that one element. The two are
            // reconciled differently in the two halves of the fraction:
            //
            // - Numerator: the excluded element is `{sequence, position, gap: false}`
            //   in both alignments, so it lies in both columns and therefore in
            //   neither symmetric difference. `Ca Δ Cb` is already exactly `A Δ B`,
            //   with nothing to subtract. See `HomologyView::column_of`.
            // - Denominator: `|A| = |Ca| - 1`, so each side loses one. A column always
            //   holds one element per sequence, so this is `num_seqs - 1` on both
            //   sides, and a single-sequence alignment correctly contributes 0 —
            //   leaving the 0/0 check below to report it.
            let (Some(column_a), Some(column_b)) = (
                view_a.column_of(sequence, position),
                view_b.column_of(sequence, position),
            ) else {
                return ResidueContribution::NONE;
            };

            ResidueContribution {
                symmetric_difference: symmetric_difference_size(column_a, column_b),
                set_sizes: (column_a.len() - 1) + (column_b.len() - 1),
            }
        })
        .reduce(|| ResidueContribution::NONE, ResidueContribution::combine);

    let (numerator, denominator) = (totals.symmetric_difference, totals.set_sizes);

    if denominator == 0 {
        bail!(
            "The two alignments have no comparable residues: every homology set is empty, \
             so the distance would be 0/0. This happens when the alignments hold a single \
             sequence each (a residue's homology set is its column minus itself), or when \
             they contain nothing but gaps."
        );
    }

    let distance = (numerator as f64) / (denominator as f64);
    log::debug!("Symmetric difference: {numerator}/{denominator} = {distance}");

    Ok(distance)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::homology::homology_view;
    use crate::msa::{Msa, read_msa};

    fn msa(records: &[(&str, &str)]) -> Msa {
        Msa::new(
            records.iter().map(|(n, _)| n.to_string()).collect(),
            records.iter().map(|(_, r)| r.as_bytes().to_vec()).collect(),
        )
        .expect("test fixture should be a valid Msa")
    }

    /// The whole pipeline for one pair of in-memory alignments, as `main.rs` runs it.
    fn distance(a: &Msa, b: &Msa) -> Result<f64> {
        let registry = Registry::for_pair(a, b)?;
        let view_a = homology_view(a, &registry)?;
        let view_b = homology_view(b, &registry)?;
        compute_symmetric_difference(&view_a, &view_b, &registry)
    }

    /// The same, from two FASTA paths.
    fn distance_from_files(path_a: &str, path_b: &str) -> Result<f64> {
        let a = read_msa(path_a)?;
        let b = read_msa(path_b)?;
        distance(&a, &b)
    }

    // -----------------------------------------------------------------------------
    // The property the denominator fix exists for
    // -----------------------------------------------------------------------------

    #[test]
    fn the_distance_is_symmetric_on_the_fixture_pair() {
        // `d(a, b) == d(b, a)` — exactly, as `f64` bits, not to within a tolerance.
        // Before Stage 4 this could not be relied on: the iteration ran over A's
        // sequence count and the denominator counted A's sets twice, so swapping the
        // arguments swapped which alignment defined both bounds (`CODE_REVIEW.md` §0).
        let forwards = distance_from_files("test/test.fasta", "test/test2.fasta")
            .expect("fixtures should compare");
        let backwards = distance_from_files("test/test2.fasta", "test/test.fasta")
            .expect("fixtures should compare");

        assert_eq!(
            forwards, backwards,
            "the distance must not depend on argument order"
        );
    }

    #[test]
    fn the_distance_is_symmetric_for_differing_gap_patterns() {
        // A second pair, chosen so the two alignments disagree about gap placement in
        // every row rather than in one, and so the two rows have different residue
        // counts (2 and 3). Same requirement: exact equality both ways round.
        let a = msa(&[("s1", "A-C--"), ("s2", "AC-G-"), ("s3", "-A-CG")]);
        let b = msa(&[("s1", "--A-C"), ("s2", "-ACG-"), ("s3", "AC--G")]);

        let forwards = distance(&a, &b).expect("same name set");
        let backwards = distance(&b, &a).expect("same name set");

        assert_eq!(
            forwards, backwards,
            "the distance must not depend on argument order"
        );
        assert_ne!(
            forwards, 0.0,
            "this pair should not be at distance 0; the test would be vacuous"
        );
    }

    #[test]
    fn registry_order_does_not_change_the_distance() {
        // The mechanism: `Registry::for_pair` is itself argument-order independent, so
        // both directions of the test above are computed against the identical registry.
        let a = msa(&[("s1", "A-C--"), ("s2", "AC-G-"), ("s3", "-A-CG")]);
        let b = msa(&[("s1", "--A-C"), ("s2", "-ACG-"), ("s3", "AC--G")]);

        assert_eq!(
            Registry::for_pair(&a, &b).expect("same name set"),
            Registry::for_pair(&b, &a).expect("same name set")
        );
    }

    // -----------------------------------------------------------------------------
    // Pinned values
    // -----------------------------------------------------------------------------

    #[test]
    fn the_fixture_pair_distance_is_hand_verified() {
        // test/test.fasta      test/test2.fasta
        // >1 AA--              >1 A-A-
        // >2 A--A              >2 A--A
        // >3 AAA-              >3 AAA-
        //
        // Registry is sorted by name: "1" -> 0, "2" -> 1, "3" -> 2. Writing r(s,p) for a
        // residue and g(s,p) for a gap (a gap's position is the index of the residue
        // preceding it in its row):
        //
        //   seq pos | A(r)                  | B(r)                  | |AΔB| | |A|+|B|
        //   --------+-----------------------+-----------------------+-------+--------
        //   0   0   | {r(1,0), r(2,0)}      | {r(1,0), r(2,0)}      |   0   |   4
        //   0   1   | {g(1,0), r(2,1)}      | {g(1,0), r(2,2)}      |   2   |   4
        //   1   0   | {r(0,0), r(2,0)}      | {r(0,0), r(2,0)}      |   0   |   4
        //   1   1   | {g(0,1), g(2,2)}      | {g(0,1), g(2,2)}      |   0   |   4
        //   2   0   | {r(0,0), r(1,0)}      | {r(0,0), r(1,0)}      |   0   |   4
        //   2   1   | {r(0,1), g(1,0)}      | {g(0,0), g(1,0)}      |   2   |   4
        //   2   2   | {g(0,1), g(1,0)}      | {r(0,1), g(1,0)}      |   2   |   4
        //   --------+-----------------------+-----------------------+-------+--------
        //                                                       sum:   6        28
        //
        // 6 / 28 = 3/14 = 0.21428571428571427.
        //
        // This is the *same* number the pre-merge binary emitted under the old `2·|A|`
        // denominator, and that is not a coincidence: |A| = |B| = num_seqs - 1 = 2 in
        // every row of the table, so 2·|A| and |A| + |B| agree term for term. Numerator
        // 6 and denominator 28 are unchanged. See the module docs for why that holds in
        // general, and for the change that would break it.
        let distance = distance_from_files("test/test.fasta", "test/test2.fasta")
            .expect("fixtures should compare");
        assert_eq!(distance, 0.21428571428571427);
    }

    #[test]
    fn an_alignment_against_itself_is_exactly_zero() {
        let m = read_msa("test/test.fasta").expect("fixture should read");
        assert_eq!(distance(&m, &m).expect("same name set"), 0.0);
    }

    #[test]
    fn reordered_records_give_exactly_zero() {
        // `test/test_reordered.fasta` is `test/test.fasta` with the records written
        // 2, 1, 3. Same alignment, so distance 0. The old positional keying reported
        // 0.5714285714285714 here (`CODE_REVIEW.md` §0).
        assert_eq!(
            distance_from_files("test/test.fasta", "test/test_reordered.fasta")
                .expect("fixtures should compare"),
            0.0
        );
        // And in the other direction, since that is the whole point of this stage.
        assert_eq!(
            distance_from_files("test/test_reordered.fasta", "test/test.fasta")
                .expect("fixtures should compare"),
            0.0
        );
    }

    // -----------------------------------------------------------------------------
    // Errors rather than wrong numbers
    // -----------------------------------------------------------------------------

    #[test]
    fn a_mismatched_name_set_errs_end_to_end() {
        // Previously this produced a number: the iteration ran to A's sequence count and
        // the missing sequences were simply skipped, giving d(3-seq, 2-seq) = 0.25 and
        // d(2-seq, 3-seq) = 0.5 for the same pair. The failure has to surface as an
        // `Err`, and it has to surface from the same entry point `main.rs` calls.
        let a = msa(&[("x", "AA--"), ("y", "A--A"), ("z", "AAA-")]);
        let b = msa(&[("x", "A-A-"), ("y", "A--A")]);

        let err = match distance(&a, &b) {
            Ok(d) => {
                panic!("alignments over different name sets must not produce a distance, got {d}")
            }
            Err(e) => e.to_string(),
        };
        assert!(
            err.contains("'z'"),
            "the error must name the missing sequence, got: {err}"
        );

        // Both directions, so neither ordering can sneak a number out.
        assert!(distance(&b, &a).is_err());
    }

    #[test]
    fn an_undefined_ratio_errs_rather_than_returning_nan() {
        // One sequence per alignment: a residue's homology set is its column minus
        // itself, so every set is empty and the ratio is 0/0. The old code divided
        // anyway and returned NaN, which `main.rs` would have written to the CSV.
        let a = msa(&[("only", "AC-G")]);
        let b = msa(&[("only", "A-CG")]);

        assert!(
            distance(&a, &b).is_err(),
            "a 0/0 ratio must be an error, not NaN"
        );
    }

    #[test]
    fn all_gap_alignments_err_rather_than_returning_nan() {
        // No residues at all, so no slot is ever `Some` and the denominator is 0 for a
        // different reason than the test above.
        let a = msa(&[("s1", "----"), ("s2", "----")]);
        let b = msa(&[("s1", "...."), ("s2", "----")]);

        assert!(
            distance(&a, &b).is_err(),
            "a 0/0 ratio must be an error, not NaN"
        );
    }

    // -----------------------------------------------------------------------------
    // Shape handling
    // -----------------------------------------------------------------------------

    #[test]
    fn alignments_of_different_widths_compare_over_the_shared_residues() {
        // Same sequences, same residues, padded to different widths. Positions are
        // indexed among a row's *residues*, not its columns, so the extra trailing gap
        // columns add no residues and the two are at distance 0 — but only if the
        // iteration covers the wider view's slots without indexing past the narrower.
        let narrow = msa(&[("s1", "AC"), ("s2", "AC")]);
        let wide = msa(&[("s1", "AC---"), ("s2", "AC---")]);

        assert_eq!(distance(&narrow, &wide).expect("same name set"), 0.0);
        assert_eq!(distance(&wide, &narrow).expect("same name set"), 0.0);
    }

    #[test]
    fn every_homology_set_has_size_num_seqs_minus_one() {
        // The invariant the module docs lean on, pinned: it is what makes `2·|A|` and
        // `|A| + |B|` agree today, so if it ever stops holding, the fixture value above
        // is expected to move and this test says why.
        let m = msa(&[("a", "A-BC-"), ("b", "-AB-C"), ("c", "AB-C-")]);
        let registry = Registry::for_pair(&m, &m).expect("same name set");
        let view = homology_view(&m, &registry).expect("view should build");

        for sequence in 0..registry.len() {
            for position in 0..view.width() {
                let Some(column) = view.column_of(sequence, position) else {
                    continue;
                };
                assert_eq!(
                    column.len() - 1,
                    m.num_seqs() - 1,
                    "a residue's homology set is its column minus itself"
                );
            }
        }
    }

    /// The distance computed the slow, literal way: materialise every residue's
    /// homology set (its column *minus itself*) and work from those.
    ///
    /// This is the definition the metric is stated in, and the shape the code had
    /// before `CODE_REVIEW.md` §2 was fixed. `compute_symmetric_difference` now works
    /// from the shared column sets instead, which is only valid because the excluded
    /// element is the same value in both alignments and so cancels out of the symmetric
    /// difference. This function exists to check that claim rather than assume it.
    fn distance_from_materialised_sets(a: &Msa, b: &Msa) -> Option<f64> {
        let registry = Registry::for_pair(a, b).ok()?;
        let view_a = homology_view(a, &registry).ok()?;
        let view_b = homology_view(b, &registry).ok()?;

        let width = view_a.width().max(view_b.width());
        let (mut numerator, mut denominator) = (0usize, 0usize);

        for sequence in 0..registry.len() {
            for position in 0..width {
                let (Some(set_a), Some(set_b)) = (
                    view_a.homology_set_of(sequence, position),
                    view_b.homology_set_of(sequence, position),
                ) else {
                    continue;
                };
                numerator += set_a.symmetric_difference(&set_b).count();
                denominator += set_a.len() + set_b.len();
            }
        }

        if denominator == 0 {
            return None;
        }
        Some(numerator as f64 / denominator as f64)
    }

    /// Every alignment of `num_seqs` x `width` over `{A, -}`, exhaustively.
    fn all_alignments(num_seqs: usize, width: usize) -> Vec<Msa> {
        let names: Vec<String> = (0..num_seqs).map(|i| format!("s{i}")).collect();
        (0..(1u32 << (num_seqs * width)))
            .map(|bits| {
                let rows = (0..num_seqs)
                    .map(|s| {
                        (0..width)
                            .map(|c| {
                                if bits >> (s * width + c) & 1 == 1 {
                                    b'A'
                                } else {
                                    b'-'
                                }
                            })
                            .collect()
                    })
                    .collect();
                Msa::new(names.clone(), rows).expect("generated alignment should be valid")
            })
            .collect()
    }

    #[test]
    fn sharing_the_column_sets_computes_exactly_the_materialised_distance() {
        // *** The correctness argument for the `CODE_REVIEW.md` §2 memory fix. ***
        //
        // `compute_symmetric_difference` reads the shared column sets, which include
        // the residue itself, and subtracts one per side from the denominator only. The
        // claim is that this is not an approximation but exactly equal to working from
        // materialised homology sets, because the excluded element is
        // `{sequence, position, gap: false}` in *both* alignments and therefore lies in
        // both columns, so it cannot appear in the symmetric difference.
        //
        // Checked exhaustively rather than by sampling: every alignment of 2x3 and 3x2
        // over {A, -}, against every other, is 8192 comparisons.
        for (num_seqs, width) in [(2, 3), (3, 2)] {
            let alignments = all_alignments(num_seqs, width);
            for a in alignments.iter() {
                for b in alignments.iter() {
                    let shipped = distance(a, b).ok();
                    let materialised = distance_from_materialised_sets(a, b);

                    assert_eq!(
                        shipped,
                        materialised,
                        "shared-column and materialised-set distances disagree for\n  a: {:?}\n  b: {:?}",
                        a.rows()
                            .iter()
                            .map(|r| String::from_utf8_lossy(r).into_owned())
                            .collect::<Vec<_>>(),
                        b.rows()
                            .iter()
                            .map(|r| String::from_utf8_lossy(r).into_owned())
                            .collect::<Vec<_>>(),
                    );
                }
            }
        }
    }

    #[test]
    fn a_column_set_is_stored_once_no_matter_how_many_residues_share_it() {
        // The memory property itself, as opposed to the numeric one: storage is one set
        // per column, not one per residue. Before the fix this alignment held 6 sets of
        // 3 elements; it now holds 3 sets of 3.
        let m = msa(&[("a", "ABC"), ("b", "ABC"), ("c", "AB-")]);
        let registry = Registry::for_pair(&m, &m).expect("same name set");
        let view = homology_view(&m, &registry).expect("view should build");

        assert_eq!(view.column_sets().len(), m.width());
        for column in view.column_sets() {
            assert_eq!(column.len(), m.num_seqs(), "one element per sequence");
        }
    }

    #[test]
    fn case_only_differences_do_not_change_the_distance() {
        // Identity excludes the raw byte, so this is 0 without any normalisation pass.
        let upper = msa(&[("s1", "AC-G"), ("s2", "A-CG")]);
        let lower = msa(&[("s1", "ac-g"), ("s2", "a-cg")]);

        assert_eq!(distance(&upper, &lower).expect("same name set"), 0.0);
    }
}
