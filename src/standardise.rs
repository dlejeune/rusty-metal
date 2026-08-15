//! Canonicalisation of an alignment: sort the sequences by name, then permute the
//! columns into a canonical order, proving via a residue hash that no residue content
//! moved relative to any other residue in its own row.
//!
//! Originally ported from the standalone `standardise-msa` tool (`src/main.rs` @
//! `97978c9`), re-expressed over the shared [`Msa`] type: the original transposed the
//! alignment into a `Vec<Column>` and bubble-sorted that, where this works over column
//! *indices* and reads through into the rows, avoiding a second full copy of every
//! alignment (see `CODE_REVIEW.md` §2 — the crate already has a memory ceiling and
//! cannot afford another one).
//!
//! The ordering rule itself is **no longer the original's**. See
//! [`canonical_columns`]: the original's comparator was inverted relative to its own
//! stated intent, keyed on a single number that tied constantly, and — being a sort
//! over a non-transitive relation — did not determine a unique answer at all. It is
//! replaced by a canonical topological construction. This is a deliberate,
//! output-changing divergence, recorded in `INTEGRATION_NOTES.md` under Stage 6.
//!
//! # Why standardisation changes the metric
//!
//! The distance treats a gap's identity as "the index of the residue preceding it in
//! its row". Permuting columns moves residues within a row, so the same gap acquires a
//! different identity, hence different homology sets. Residues themselves are never
//! affected, because the comparator below never swaps two columns that both hold a
//! residue in the same row — which is precisely the invariant [`Msa::residue_hash`]
//! verifies. Standardising both inputs therefore strips gap-placement artefacts out of
//! the distance.
//!
//! # Intentional divergences from the original
//!
//! - The original hard-coded `const GAP: u8 = b'-'` and treated `.` as a residue. This
//!   uses the shared [`is_gap`], which counts both `-` and `.` (`CODE_REVIEW.md` §3):
//!   `.` is a real gap character in Pfam-derived FASTA, and the two tools previously
//!   disagreed about it.
//! - The column ordering rule is replaced outright; see [`canonical_columns`].
//! - All-gap columns are **dropped** rather than retained. The original pinned them in
//!   place — despite a comment claiming the opposite — so a column carrying no residue
//!   at all could still hold the rest of the alignment apart. A column of pure gaps
//!   states no homology relationship, so it is removed.
//!
//! Output is therefore **not** byte-compatible with `standardise-msa`, and distances
//! computed after standardisation differ from those of any earlier build of this crate.

use crate::msa::{is_gap, Msa};
use anyhow::{bail, Context, Result};
use std::cmp::Ordering;
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

/// Whether two columns are **pinned**: whether they are forbidden from moving past
/// each other.
///
/// This is a direct port of `Column::eq` in the original tool, where `eq` returning
/// `true` meant "pinned", not "equal" in any ordinary sense.
///
/// **This predicate is not an equivalence relation and is deliberately
/// non-transitive.** A may be movable past B, and B past C, while A is pinned against
/// C — a residue in a row that A and C share but B does not is enough. Consult
/// `column_order` before reaching for a standard sort.
///
/// Two columns are pinned when any row holds a residue in *both* of them: swapping
/// would reverse those two residues within that row and corrupt the sequence.
/// Whether two columns are **pinned**: whether they are forbidden from moving past
/// each other.
///
/// Two columns are pinned when any row holds a residue in *both* of them: swapping
/// would reverse those two residues within that row and corrupt the sequence. That is
/// the whole of the rule — see `canonical_columns` for what replaced the original
/// tool's four other branches.
///
/// **This predicate is not an equivalence relation and is deliberately
/// non-transitive.** A may be movable past B, and B past C, while A is pinned against
/// C — a residue in a row that A and C share but B does not is enough. It is a
/// precedence constraint, not an ordering, which is why `canonical_columns` treats it
/// as a DAG rather than handing it to a sort.
fn columns_are_pinned(msa: &Msa, a: usize, b: usize) -> bool {
    msa.rows()
        .iter()
        .any(|row| !is_gap(row[a]) && !is_gap(row[b]))
}

/// Whether a column holds no residue at all.
fn column_is_all_gaps(msa: &Msa, col: usize) -> bool {
    msa.rows().iter().all(|row| is_gap(row[col]))
}

/// Orders two columns by the canonical key: read top to bottom, **a residue sorts
/// before a gap**. The lexicographically smaller column is the one that fills the
/// higher-numbered rows further to the left.
///
/// Compared lazily row by row rather than by materialising a key per column, so this
/// allocates nothing and the whole pass stays within the O(width) extra memory the
/// crate can afford (`CODE_REVIEW.md` §2).
fn compare_columns(msa: &Msa, a: usize, b: usize) -> Ordering {
    for row in msa.rows().iter() {
        // `false` (residue) sorts before `true` (gap).
        match is_gap(row[a]).cmp(&is_gap(row[b])) {
            Ordering::Equal => continue,
            other => return other,
        }
    }
    Ordering::Equal
}

/// Computes the canonical column selection for `msa`, in `Msa::select_columns` form:
/// `result[new_index] = old_index`.
///
/// All-gap columns are **dropped**, so the result is generally shorter than
/// `msa.width()`. Everything else is ordered so that residues fill the higher rows as
/// far to the left as the pinning constraints permit.
///
/// Does not touch the alignment.
///
/// # Why this is a topological sort and not a sort
///
/// *** DO NOT REPLACE THIS WITH `sort_by`, `sort_unstable_by`, OR A BUBBLE SORT. ***
///
/// The original tool bubble-sorted columns under a comparator that returned `Equal`
/// for a pinned pair. That is not merely a violation of the total-order contract that
/// `sort_by` requires — it does not determine an answer at all. When a column is free
/// to sit in several places (movable past both of two columns that are pinned to each
/// other), a comparator has nothing to say about which place is canonical, and a sort
/// that only swaps adjacent pairs simply leaves the column wherever the input file
/// happened to put it.
///
/// The consequence was that standardisation **was not confluent**: two files holding
/// the same alignment, differing only by a legal column permutation, standardised to
/// different layouts and compared as different. Measured over random alignments, 78%
/// of legal column shuffles changed the standardised layout and 36% produced a nonzero
/// distance between an alignment and itself, up to the maximum of 1.0. That defeats
/// the entire purpose of standardising.
///
/// The fix is to *construct* the order rather than sort into it. The pinning relation
/// is a property of the columns themselves, and a legal permutation never reorders a
/// pinned pair, so the precedence DAG is identical no matter which legal permutation
/// of an alignment we are handed. Emitting the smallest available column under a total
/// key therefore depends only on that DAG — the result is confluent by construction.
/// `standardisation_is_confluent_over_legal_permutations` is the property test.
///
/// # Why the key never ties
///
/// All-gap columns are dropped first, so every remaining column holds a residue. If
/// two distinct columns had the same gap pattern, then either some row holds a residue
/// in both — in which case they are pinned, and one is an ancestor of the other in the
/// DAG, so they are never both available for selection at the same time — or every row
/// gaps in both, in which case both were all-gap and already dropped. So the minimum
/// among the available columns is always unique, and no tiebreak is needed.
/// `available_columns_never_tie` pins that argument.
///
/// Known cost: O(width^2 * num_seqs), the same as the bubble sort it replaces, in
/// O(width) extra memory.
pub fn canonical_columns(msa: &Msa) -> Vec<usize> {
    // All-gap columns carry no alignment information: no residue sits in one, so no
    // homology relationship depends on where it is. Keeping them would mean the
    // metric could be moved by a column that says nothing.
    let kept: Vec<usize> = (0..msa.width())
        .filter(|&col| !column_is_all_gaps(msa, col))
        .collect();
    let n = kept.len();

    // Precedence edges run from lower to higher *current* index between pinned pairs;
    // preserving them is exactly what makes the resulting order legal. Only the
    // in-degree is stored, and the edges are recomputed on release, so this is
    // O(width) memory rather than the O(width^2) an adjacency list would cost.
    let mut indegree = vec![0usize; n];
    for i in 0..n {
        for j in (i + 1)..n {
            if columns_are_pinned(msa, kept[i], kept[j]) {
                indegree[j] += 1;
            }
        }
    }

    let mut placed = vec![false; n];
    let mut order = Vec::with_capacity(n);
    for _ in 0..n {
        let next = (0..n)
            .filter(|&k| !placed[k] && indegree[k] == 0)
            .min_by(|&x, &y| compare_columns(msa, kept[x], kept[y]))
            .expect(
                "the pinning DAG is acyclic (every edge runs from a lower to a higher \
                 column index), so some column always has in-degree zero",
            );
        placed[next] = true;
        // Anything `next` pinned is a successor and so was still unplaced.
        for j in (next + 1)..n {
            if columns_are_pinned(msa, kept[next], kept[j]) {
                indegree[j] -= 1;
            }
        }
        order.push(kept[next]);
    }

    order
}

/// An order-independent digest of the alignment's residue content: one hash of
/// `(name, gap-filtered residues)` per sequence, combined with XOR.
///
/// [`Msa::residue_hash`] walks the sequences in row order, so it necessarily changes
/// when the sequences are reordered — which is exactly what sorting by name does. This
/// digest ignores sequence order, so it can bracket the *whole* of `standardise`,
/// including the sort. XOR is safe as a combiner here because names are unique
/// (enforced by `Msa::new`), so no two per-sequence hashes can cancel each other out.
fn residue_multiset_digest(msa: &Msa) -> u64 {
    let mut acc: u64 = 0;
    for (name, row) in msa.names().iter().zip(msa.rows().iter()) {
        let mut hasher = DefaultHasher::new();
        name.hash(&mut hasher);
        let residues: Vec<u8> = row.iter().copied().filter(|b| !is_gap(*b)).collect();
        residues.hash(&mut hasher);
        acc ^= hasher.finish();
    }
    acc
}

/// Standardises an alignment in place: sorts the sequences by name, then rewrites the
/// columns into canonical order, dropping any that hold no residue.
///
/// The result is a genuine canonical form: any two alignments that differ only by a
/// legal column permutation standardise to byte-identical output.
///
/// Returns `Err` if the operation altered any residue content. Two checks bracket the
/// work, for the reason described on [`residue_multiset_digest`]:
///
/// - `residue_hash` before and after the *column rewrite*. This is the original tool's
///   check, and the one that matters: it is what proves the ordering only ever made
///   legal moves. It covers column dropping too — `residue_hash` filters gaps out, so
///   discarding an all-gap column leaves it untouched, while discarding a column that
///   held a residue would change it and be caught.
/// - an order-independent digest before and after the *whole* operation, so the sort
///   by name is covered too. Sorting moves whole `(name, row)` pairs, so this can only
///   fire if `sort_sequences_by_name` is ever broken — which is the point of having it.
pub fn standardise(msa: &mut Msa) -> Result<()> {
    let digest_before = residue_multiset_digest(msa);

    msa.sort_sequences_by_name();

    // Taken *after* the sort, matching the original, which hashed once the sequences
    // were already in name order. `residue_hash` is order-dependent, so a hash taken
    // before the sort would differ from the one after for any input that was not
    // already sorted, and standardisation would reject its own legal work.
    let hash_before_columns = msa.residue_hash();

    let order = canonical_columns(msa);
    msa.select_columns(&order)
        .context("Standardisation computed an invalid column selection")?;

    let hash_after_columns = msa.residue_hash();
    if hash_before_columns != hash_after_columns {
        bail!(
            "Standardisation altered residue content: residue hash was {:x} before the column \
             rewrite and {:x} after. Columns were reordered in a way that moved residues \
             relative to each other within a sequence, or a column holding a residue was \
             dropped; the resulting alignment has been discarded.",
            hash_before_columns,
            hash_after_columns
        );
    }

    let digest_after = residue_multiset_digest(msa);
    if digest_before != digest_after {
        bail!(
            "Standardisation altered residue content across the sort by name: digest was {:x} \
             before and {:x} after.",
            digest_before,
            digest_after
        );
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::msa::read_msa;

    /// Builds an `Msa` from `(name, row)` pairs written the way an alignment reads.
    fn msa(records: &[(&str, &str)]) -> Msa {
        Msa::new(
            records.iter().map(|(n, _)| n.to_string()).collect(),
            records.iter().map(|(_, r)| r.as_bytes().to_vec()).collect(),
        )
        .expect("test fixture should be a valid Msa")
    }

    /// The alignment as `(name, row)` pairs, for readable assertions.
    fn rows_as_strings(m: &Msa) -> Vec<(String, String)> {
        m.names()
            .iter()
            .zip(m.rows().iter())
            .map(|(n, r)| (n.clone(), String::from_utf8_lossy(r).into_owned()))
            .collect()
    }

    fn expected(records: &[(&str, &str)]) -> Vec<(String, String)> {
        records
            .iter()
            .map(|(n, r)| (n.to_string(), r.to_string()))
            .collect()
    }

    #[test]
    fn standardise_is_idempotent() {
        for fixture in ["test/test.fasta", "test/test2.fasta", "test/unsorted_names.fasta"] {
            let mut once = read_msa(fixture).expect("fixture should read");
            standardise(&mut once).expect("standardisation should succeed");

            let mut twice = once.clone();
            standardise(&mut twice).expect("re-standardisation should succeed");

            assert_eq!(once, twice, "standardisation of {fixture} is not idempotent");
        }
    }

    #[test]
    fn residue_hash_is_unchanged_by_standardisation_on_real_fixtures() {
        for fixture in ["test/test.fasta", "test/test2.fasta"] {
            let original = read_msa(fixture).expect("fixture should read");

            // The fixtures are already in name order, so `residue_hash` is directly
            // comparable before and after; sorting is a no-op for them.
            let mut standardised = original.clone();
            standardise(&mut standardised).expect("standardisation should succeed");

            assert_eq!(
                original.residue_hash(),
                standardised.residue_hash(),
                "standardisation changed the residue hash of {fixture}"
            );
        }
    }

    #[test]
    fn sequences_are_sorted_by_name() {
        // `test/test.fasta` and `test/test2.fasta` are already in sorted order, so they
        // cannot distinguish a working sort from a no-op. This fixture is not.
        let mut m = read_msa("test/unsorted_names.fasta").expect("fixture should read");
        assert_eq!(
            m.names(),
            &["gamma".to_string(), "alpha".to_string(), "beta".to_string()],
            "the fixture must start out unsorted for this test to mean anything"
        );

        standardise(&mut m).expect("standardisation should succeed");

        // Rows travel with their names, and the columns are then canonicalised against
        // the *sorted* row order — the ordering key reads rows top to bottom, so the
        // sort necessarily happens first.
        //
        //   after sorting      columns          canonical order [0, 3, 1, 2]
        //   alpha  A--A        c0 = (A,A,A)     alpha  AA--
        //   beta   AAA-        c1 = (-,A,A)     beta   A-AA
        //   gamma  AA--        c2 = (-,A,-)     gamma  A-A-
        //                      c3 = (A,-,-)
        //
        // c0 is pinned before all three. c1 is pinned before c2. c3 is free of both c1
        // and c2, and holds a residue in row `alpha`, so it goes as far left as the
        // pins allow: immediately after c0.
        assert_eq!(
            rows_as_strings(&m),
            expected(&[("alpha", "AA--"), ("beta", "A-AA"), ("gamma", "A-A-")])
        );
    }

    #[test]
    fn all_gap_columns_are_dropped() {
        // *** DELIBERATE CHANGE FROM THE ORIGINAL TOOL (Stage 6). ***
        //
        // The original pinned all-gap columns, so they never moved — despite a comment
        // claiming "we can move it around". Verified against the built original at
        // `97978c9`, which emitted this input completely unchanged:
        //
        //     in:  >a A--   out:  >a A--
        //          >b --B         >b --B
        //
        // A column holding no residue states no homology relationship, so keeping it
        // let a column that says nothing hold the rest of the alignment apart. It is
        // now removed entirely and the width shrinks.
        let mut m = msa(&[("a", "A--"), ("b", "--B")]);
        standardise(&mut m).expect("standardisation should succeed");

        assert_eq!(
            rows_as_strings(&m),
            expected(&[("a", "A-"), ("b", "-B")]),
            "the middle column holds no residue and must be dropped"
        );
        assert_eq!(m.width(), 2, "dropping a column must shrink the width");
    }

    #[test]
    fn an_all_gap_alignment_standardises_to_zero_width() {
        // The limit case of dropping: nothing survives. Left representable rather than
        // rejected — `Msa::select_columns` documents why, and the distance stage
        // already reports an empty comparison as an error.
        let mut m = msa(&[("a", "---"), ("b", "---")]);
        standardise(&mut m).expect("standardisation should succeed");

        assert_eq!(m.width(), 0);
        assert_eq!(rows_as_strings(&m), expected(&[("a", ""), ("b", "")]));
    }

    #[test]
    fn residues_fill_the_higher_rows_leftward() {
        // The canonical rule, in its smallest form. Column 0 is `[A, -]`, column 1 is
        // `[-, B]`; no row holds a residue in both, so they are free to move.
        //
        // Row `a` is the higher row and it has its residue in column 0, so column 0
        // sorts first and the alignment is already canonical.
        //
        // *** The original tool emitted `-A` / `B-` here — the exact mirror image. ***
        // It keyed on the row index of the first *gap*, ascending, so a column with a
        // residue up top sorted last. That was inverted relative to its intent.
        let mut m = msa(&[("a", "A-"), ("b", "-B")]);
        standardise(&mut m).expect("standardisation should succeed");

        assert_eq!(
            rows_as_strings(&m),
            expected(&[("a", "A-"), ("b", "-B")])
        );

        // And the mirror input must be pulled into that same canonical form.
        let mut m = msa(&[("a", "-A"), ("b", "B-")]);
        standardise(&mut m).expect("standardisation should succeed");

        assert_eq!(
            rows_as_strings(&m),
            expected(&[("a", "A-"), ("b", "-B")]),
            "both column orders of one alignment must reach the same canonical form"
        );
    }

    #[test]
    fn pinned_columns_do_not_move() {
        // Chosen so that the pin is the *only* thing holding the columns in place:
        // column 0 is `[A, A, -]` (first gap in row 2) and column 1 is `[B, -, -]`
        // (first gap in row 1), so on first-gap order alone 2 > 1 and they would swap.
        // Row `a` holds a residue in both, so swapping would emit `BA` where the input
        // said `AB` — reversing two residues within a sequence. They are pinned.
        let mut m = msa(&[("a", "AB"), ("b", "A-"), ("c", "--")]);
        let before = m.clone();
        standardise(&mut m).expect("standardisation should succeed");

        assert_eq!(m, before, "columns sharing a residue row must never swap");
        assert_eq!(
            rows_as_strings(&m),
            expected(&[("a", "AB"), ("b", "A-"), ("c", "--")])
        );
    }

    #[test]
    fn the_test_fixture_standardises_to_its_canonical_form() {
        // *** CHANGED IN STAGE 6. *** `test/test.fasta` used to be a fixed point: the
        // original tool emitted it unchanged, and this test asserted a no-op.
        //
        //   input          columns          canonical order [0, 1, 3, 2]
        //   1  AA--        c0 = (A,A,A)     1  AA--
        //   2  A--A        c1 = (A,-,A)     2  A-A-
        //   3  AAA-        c2 = (-,-,A)     3  AA-A
        //                  c3 = (-,A,-)
        //
        // c0 is pinned before everything. c1 is pinned before c2 (both hold a residue
        // in row 3). c3 is free of c1 and c2. c1 holds a residue in row 1 where c3
        // gaps, so c1 goes first; then c3 beats c2 on row 2.
        let mut m = read_msa("test/test.fasta").expect("fixture should read");
        standardise(&mut m).expect("standardisation should succeed");

        assert_eq!(
            rows_as_strings(&m),
            expected(&[("1", "AA--"), ("2", "A-A-"), ("3", "AA-A")])
        );
    }

    #[test]
    fn the_test2_fixture_standardises_to_its_canonical_form() {
        // *** CHANGED IN STAGE 6. *** This was the port's regression gate against the
        // built `standardise-msa` binary at `97978c9`, which emitted:
        //
        //     >1 A--A   >2 A-A-   >3 AA-A
        //
        // The new rule gives a different answer, and deliberately so — see
        // `canonical_columns`. `test2` is `test` with two *pinned* columns exchanged,
        // so the two fixtures are genuinely different alignments, not a legal
        // permutation of one another, and they do not standardise to the same thing.
        //
        //   input          columns          canonical order [0, 3, 1, 2]
        //   1  A-A-        c0 = (A,A,A)     1  A--A
        //   2  A--A        c1 = (-,-,A)     2  AA--
        //   3  AAA-        c2 = (A,-,A)     3  A-AA
        //                  c3 = (-,A,-)
        let mut m = read_msa("test/test2.fasta").expect("fixture should read");
        standardise(&mut m).expect("standardisation should succeed");

        assert_eq!(
            rows_as_strings(&m),
            expected(&[("1", "A--A"), ("2", "AA--"), ("3", "A-AA")])
        );
    }

    #[test]
    fn dot_is_treated_as_a_gap_intentional_divergence() {
        // The original hard-coded `-` as the only gap, so it would have seen `.` as a
        // residue and pinned these two columns. The shared `is_gap` counts `.`, so the
        // columns are free — `CODE_REVIEW.md` §3.
        //
        // Being free, they are then ordered by the canonical rule: row `a` holds its
        // residue in column 0, so column 0 leads and the input is already canonical.
        let mut m = msa(&[("a", "A."), ("b", ".B")]);
        standardise(&mut m).expect("standardisation should succeed");

        assert_eq!(
            rows_as_strings(&m),
            expected(&[("a", "A."), ("b", ".B")])
        );

        // The mirror form must converge on the same layout. Note the gap characters
        // travel with their columns, so the `.` lands in row `b` here.
        let mut m = msa(&[("a", ".A"), ("b", "B.")]);
        standardise(&mut m).expect("standardisation should succeed");

        assert_eq!(
            rows_as_strings(&m),
            expected(&[("a", "A."), ("b", ".B")])
        );
    }

    #[test]
    fn canonical_columns_returns_new_to_old_indices() {
        let m = msa(&[("1", "A-A-"), ("2", "A--A"), ("3", "AAA-")]);
        assert_eq!(canonical_columns(&m), vec![0, 3, 1, 2]);
    }

    #[test]
    fn canonical_columns_omits_all_gap_columns() {
        // Column 1 holds no residue, so it is absent from the selection entirely
        // rather than appearing somewhere in it.
        let m = msa(&[("a", "A-C"), ("b", "-.C")]);
        let cols = canonical_columns(&m);

        assert!(!cols.contains(&1), "an all-gap column must not be selected, got {cols:?}");
        assert_eq!(cols.len(), 2);
    }

    #[test]
    fn select_columns_accepts_a_subset_and_shrinks_the_width() {
        let mut m = msa(&[("a", "ABCD")]);
        m.select_columns(&[3, 0]).expect("a subset is a legitimate selection");

        assert_eq!(rows_as_strings(&m), expected(&[("a", "DA")]));
        assert_eq!(m.width(), 2);
    }

    #[test]
    fn select_columns_rejects_more_columns_than_exist() {
        let mut m = msa(&[("a", "ABCD")]);
        assert!(
            m.select_columns(&[0, 1, 2, 3, 3]).is_err(),
            "more indices than there are columns cannot be satisfied without a repeat"
        );
        assert_eq!(m.rows()[0], b"ABCD".to_vec(), "a rejected selection must not mutate");
    }

    #[test]
    fn select_columns_rejects_duplicate_index() {
        let mut m = msa(&[("a", "ABCD")]);
        assert!(
            m.select_columns(&[0, 1, 1, 3]).is_err(),
            "a repeated column index would duplicate a residue"
        );
        assert_eq!(m.rows()[0], b"ABCD".to_vec());
    }

    #[test]
    fn select_columns_rejects_out_of_range_index() {
        let mut m = msa(&[("a", "ABCD")]);
        assert!(
            m.select_columns(&[0, 1, 2, 4]).is_err(),
            "an index at or past the width must be rejected, not panic"
        );
        assert_eq!(m.rows()[0], b"ABCD".to_vec());
    }

    // ---------------------------------------------------------------------------
    // Property tests for confluence.
    //
    // These are the reason the ordering rule was replaced in Stage 6, and they are
    // the tests that would have caught the original. Randomised but deterministically
    // seeded, so a failure is reproducible.
    // ---------------------------------------------------------------------------

    /// xorshift64. A generator is needed here and pulling in `rand` for four tests is
    /// not worth a dependency.
    struct Rng(u64);

    impl Rng {
        fn new(seed: u64) -> Rng {
            Rng(seed | 1)
        }
        fn next_u64(&mut self) -> u64 {
            let mut x = self.0;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            self.0 = x;
            x
        }
        fn below(&mut self, n: usize) -> usize {
            (self.next_u64() % (n as u64)) as usize
        }
    }

    /// A random alignment, weighted toward gaps so that free columns are common.
    fn random_msa(rng: &mut Rng) -> Msa {
        let num_seqs = 2 + rng.below(5);
        let width = 2 + rng.below(8);
        let alphabet: &[u8] = if rng.below(4) == 0 {
            b"----..ACGT"
        } else {
            b"------ACGT"
        };
        let names = (0..num_seqs).map(|i| format!("s{i}")).collect();
        let rows = (0..num_seqs)
            .map(|_| (0..width).map(|_| alphabet[rng.below(alphabet.len())]).collect())
            .collect();
        Msa::new(names, rows).expect("generated alignment should be valid")
    }

    /// Shuffles columns by repeatedly swapping adjacent *unpinned* pairs — every
    /// resulting alignment is the same alignment, written differently.
    fn legal_shuffle(rng: &mut Rng, msa: &Msa, steps: usize) -> Msa {
        let mut order: Vec<usize> = (0..msa.width()).collect();
        for _ in 0..steps {
            let i = rng.below(msa.width() - 1);
            if !columns_are_pinned(msa, order[i], order[i + 1]) {
                order.swap(i, i + 1);
            }
        }
        let mut shuffled = msa.clone();
        shuffled
            .select_columns(&order)
            .expect("a reordering is a valid selection");
        shuffled
    }

    #[test]
    fn standardisation_is_confluent_over_legal_permutations() {
        // *** THE HEADLINE PROPERTY. ***
        //
        // Two files holding the same alignment, differing only in where the gaps were
        // placed, must standardise to byte-identical output. The original tool failed
        // this on 78% of legal shuffles; this asserts it on every one.
        let mut rng = Rng::new(0xc0ffee);
        let mut shuffles_that_changed_the_layout = 0;

        for _ in 0..4000 {
            let original = random_msa(&mut rng);
            let shuffled = legal_shuffle(&mut rng, &original, 15);
            if shuffled != original {
                shuffles_that_changed_the_layout += 1;
            }

            let mut a = original.clone();
            let mut b = shuffled.clone();
            standardise(&mut a).expect("standardisation should succeed");
            standardise(&mut b).expect("standardisation should succeed");

            assert_eq!(
                a,
                b,
                "standardisation is not confluent.\n  original: {:?}\n  shuffled: {:?}",
                rows_as_strings(&original),
                rows_as_strings(&shuffled)
            );
        }

        assert!(
            shuffles_that_changed_the_layout > 1000,
            "only {shuffles_that_changed_the_layout} shuffles actually moved a column; \
             the generator is not exercising the property"
        );
    }

    #[test]
    fn the_original_tools_counterexample_now_converges() {
        // Shrunk from the Stage 6 search. Columns 0 and 2 both hold a residue in row
        // `s1`, so they are pinned; column 1 is free of both and so has three legal
        // resting places. The original left it wherever the input put it, which is
        // precisely the ambiguity that made the metric report 0.625 for an alignment
        // against itself.
        let a = msa(&[("s0", "A--"), ("s1", "A-A"), ("s2", "-A-")]);
        let b = msa(&[("s0", "A--"), ("s1", "AA-"), ("s2", "--A")]);

        let mut sa = a.clone();
        let mut sb = b.clone();
        standardise(&mut sa).expect("standardisation should succeed");
        standardise(&mut sb).expect("standardisation should succeed");

        assert_eq!(sa, sb);
        assert_eq!(
            rows_as_strings(&sa),
            expected(&[("s0", "A--"), ("s1", "AA-"), ("s2", "--A")])
        );
    }

    #[test]
    fn available_columns_never_tie() {
        // The argument that `canonical_columns` needs no tiebreak: once all-gap
        // columns are dropped, any two columns with the same gap pattern must share a
        // residue row, and so are pinned — meaning one is an ancestor of the other and
        // they are never both available for selection at once.
        let mut rng = Rng::new(0x7a1e5);

        for _ in 0..4000 {
            let m = random_msa(&mut rng);
            let kept: Vec<usize> = (0..m.width())
                .filter(|&c| !column_is_all_gaps(&m, c))
                .collect();

            for (i, &a) in kept.iter().enumerate() {
                for &b in kept.iter().skip(i + 1) {
                    if compare_columns(&m, a, b) == Ordering::Equal {
                        assert!(
                            columns_are_pinned(&m, a, b),
                            "columns {a} and {b} of {:?} tie on the ordering key but are \
                             not pinned, so the canonical order would be ambiguous",
                            rows_as_strings(&m)
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn standardisation_never_alters_residues_on_random_input() {
        // `standardise` checks this internally and returns `Err`; this asserts the
        // check never fires, which is the stronger claim that the ordering only ever
        // makes legal moves.
        let mut rng = Rng::new(0xd15ea5e);

        for _ in 0..4000 {
            let original = random_msa(&mut rng);
            let mut standardised = original.clone();
            standardise(&mut standardised).expect("standardisation must not reject its own work");

            let residues = |m: &Msa| -> Vec<Vec<u8>> {
                m.rows()
                    .iter()
                    .map(|r| r.iter().copied().filter(|b| !is_gap(*b)).collect())
                    .collect()
            };
            let mut sorted = original.clone();
            sorted.sort_sequences_by_name();

            assert_eq!(
                residues(&sorted),
                residues(&standardised),
                "residue content changed for {:?}",
                rows_as_strings(&original)
            );
        }
    }

    #[test]
    fn standardisation_is_idempotent_on_random_input() {
        let mut rng = Rng::new(0x1dea1);

        for _ in 0..2000 {
            let mut once = random_msa(&mut rng);
            standardise(&mut once).expect("standardisation should succeed");
            let mut twice = once.clone();
            standardise(&mut twice).expect("re-standardisation should succeed");

            assert_eq!(once, twice, "standardisation is not idempotent");
        }
    }
}
