//! The homology view: for every residue in an alignment, the set of elements that
//! share its column.
//!
//! This is `utils::create_hashsets` rebuilt as a *derived view* over the shared
//! [`Msa`] type rather than as the parse target. Two things change, and both are
//! deliberate behaviour fixes rather than refactoring noise:
//!
//! - **Sequences are identified by name, not by their ordinal position in the file.**
//!   The old reader passed a record counter as the sequence id (`utils.rs:74`), so two
//!   files holding the same alignment with the records written in a different order had
//!   disjoint element identities and reported a distance of `0.571` instead of `0`
//!   (`CODE_REVIEW.md` §0). Ids now come from a [`Registry`] shared across the pair of
//!   alignments being compared, which also gives us a place to reject two alignments
//!   that are not over the same sequences at all.
//! - **Element identity excludes the raw character.** See [`Element`].
//!
//! A residue is addressed as `[sequence][residue position]`, with `None` meaning "this
//! sequence has no residue at that position" (i.e. past its end). The old `MsaHashSets`
//! stored a separate set at every one of those slots; [`HomologyView`] stores one set
//! per *column* and indexes into it, which is the `CODE_REVIEW.md` §2 fix.

use crate::msa::{Msa, is_gap};
use anyhow::{Result, bail};
use std::collections::{HashMap, HashSet};

/// One position in one sequence of an alignment: either a residue or a gap.
///
/// # Identity
///
/// `Element` is keyed on `(seq, position, gap)` and *nothing else*. Two points about
/// that, both of which are easy to get wrong:
///
/// - **`gap` is a stored field, not something inferred.** The old `SequenceElement`
///   carried a `Base` enum whose `GAP` variant is what distinguished a gap from a
///   residue. That enum is being retired, and without an explicit flag a residue at
///   position 0 and a gap whose preceding residue is position 0 would be
///   indistinguishable — they would collide in the hash set and the two very different
///   things would compare equal.
/// - **The residue character is deliberately not part of this type.** The old
///   `SequenceElement::eq` compared `base`, but the reader uppercased every input byte
///   first. [`crate::msa::read_msa`] no longer does (it has to round-trip sequences to
///   disk faithfully), so keying on the character would make the same residue written
///   `a` in one file and `A` in another compare unequal and inflate the distance. This
///   loses no information: for a fixed sequence and a fixed position the residue is
///   already determined, so the character never contributed anything to identity in the
///   first place. Dropping it makes the metric case-insensitive with no normalisation
///   pass.
///
/// # Position
///
/// `position` follows the old `Sequence::from_characters`
/// (`datastructures.rs:206-233`) exactly:
///
/// - a **residue**'s position is its index among the non-gap characters of its row;
/// - a **gap**'s position is the index of the residue *preceding* it in its row, i.e.
///   `count.checked_sub(1)`, which is `None` for a gap before any residue.
///
/// A gap therefore has an identity tied to where it sits relative to its own row's
/// residues, which is why permuting columns changes gap identities (and hence the
/// distance) while leaving residue identities alone — the property the standardisation
/// pass exists to exploit.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Element {
    seq: u32,
    position: Option<u32>,
    gap: bool,
}

// `Element` deliberately has no accessors. Its three fields exist only to give the
// value an identity inside a `HashSet`, and nothing outside this module has ever
// needed to read one back — the metric only asks whether two elements are equal.
// Accessors were written in Stage 3 and removed in Stage 6 as dead code; if a caller
// ever needs them, they are three one-line functions.

/// A name-to-id mapping shared by the two alignments in a comparison.
///
/// The ids handed out here are what makes the metric independent of record order in the
/// input files. Names are unique on every [`Msa`] (enforced by `Msa::new` since Stage
/// 2), so no de-duplication is needed here.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Registry {
    /// Registry index -> name. Sorted, so `for_pair(a, b)` and `for_pair(b, a)` produce
    /// the same registry and the ids do not depend on which alignment was passed first.
    names: Vec<String>,
    /// Name -> registry index.
    index: HashMap<String, u32>,
}

impl Registry {
    /// Builds the shared registry for a pair of alignments, requiring that they contain
    /// exactly the same set of sequence names.
    ///
    /// This is the check that tells a user their two alignments are not comparable, so
    /// the error names the sequences involved on both sides rather than just reporting
    /// that the sets differ. Previously nothing validated this at all: mismatched
    /// alignments were silently compared over a truncated overlap, and `d(a, b)` and
    /// `d(b, a)` disagreed (`CODE_REVIEW.md` §0).
    pub fn for_pair(a: &Msa, b: &Msa) -> Result<Registry> {
        let names_a: HashSet<&str> = a.names().iter().map(|n| n.as_str()).collect();
        let names_b: HashSet<&str> = b.names().iter().map(|n| n.as_str()).collect();

        if names_a != names_b {
            let mut only_in_a: Vec<&str> = names_a.difference(&names_b).copied().collect();
            let mut only_in_b: Vec<&str> = names_b.difference(&names_a).copied().collect();
            only_in_a.sort_unstable();
            only_in_b.sort_unstable();

            bail!(
                "The two alignments are not over the same sequences, so their homology sets \
                 cannot be compared.\n  \
                 {} sequence(s) present only in the first alignment (missing from the second): {}\n  \
                 {} sequence(s) present only in the second alignment (missing from the first): {}",
                only_in_a.len(),
                format_names(&only_in_a),
                only_in_b.len(),
                format_names(&only_in_b),
            );
        }

        let mut names: Vec<String> = a.names().to_vec();
        names.sort();

        if u32::try_from(names.len()).is_err() {
            bail!(
                "The alignments contain {} sequences, which exceeds the {} that a sequence id \
                 can address",
                names.len(),
                u32::MAX
            );
        }

        let mut index: HashMap<String, u32> = HashMap::with_capacity(names.len());
        for (i, name) in names.iter().enumerate() {
            // `i` fits in u32: checked above.
            index.insert(name.clone(), i as u32);
        }

        Ok(Registry { names, index })
    }

    /// The number of sequences the pair has in common — which, because
    /// [`Registry::for_pair`] rejects anything else, is the number of sequences in each.
    pub fn len(&self) -> usize {
        self.names.len()
    }

    /// Registry index for a sequence name, or `None` if the name is not in the pair.
    pub fn index_of(&self, name: &str) -> Option<u32> {
        self.index.get(name).copied()
    }

    /// The name a registry index refers to.
    ///
    /// `#[cfg(test)]` because that is the honest scope: the tests below use it to
    /// check the id assignment, and no production caller does. `Registry::for_pair`
    /// builds its own error messages from the `Msa` names directly, before a registry
    /// exists to ask. Drop the attribute the moment something real needs it.
    #[cfg(test)]
    pub fn name_of(&self, index: u32) -> Option<&str> {
        self.names.get(index as usize).map(|n| n.as_str())
    }

    /// All names, in registry-index order. `#[cfg(test)]` for the same reason as
    /// [`Registry::name_of`].
    #[cfg(test)]
    pub fn names(&self) -> &[String] {
        &self.names
    }
}

fn format_names(names: &[&str]) -> String {
    if names.is_empty() {
        "(none)".to_string()
    } else {
        names
            .iter()
            .map(|n| format!("'{n}'"))
            .collect::<Vec<String>>()
            .join(", ")
    }
}

/// The homology sets of one alignment.
///
/// A residue's homology set is *its whole column minus itself*. This stores one set per
/// **column**, shared by every residue in it, plus an index from
/// `[registry sequence index][residue position]` to the column that residue sits in.
/// The "minus itself" is never materialised — see [`HomologyView::column_of`] for why
/// it does not need to be.
///
/// The sequence dimension is indexed by **registry index**, so two views built against
/// the same [`Registry`] are row-aligned regardless of record order in either file.
///
/// # Why this is not one set per residue
///
/// It used to be, and that was `CODE_REVIEW.md` §2: cloning the column set for every
/// residue costs O(num_seqs² × width) live memory per alignment, both alignments are
/// held at once, and that is multiplied again by the number of pairs running
/// concurrently. A 500 × 5000 alignment came to roughly 1.2×10⁹ set entries. It was a
/// hard ceiling on usable input size rather than a micro-optimisation.
///
/// Sharing the column set brings that to O(num_seqs × width) — the size of the
/// alignment itself.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HomologyView {
    /// One set per column of the alignment, in column order. Each holds exactly one
    /// element per sequence, so `len()` is the alignment's sequence count.
    columns: Vec<HashSet<Element>>,
    /// `[registry sequence index][residue position]` -> index into `columns`.
    ///
    /// `None` means that sequence has no residue at that position, either because it
    /// has fewer residues than the alignment is wide or because the slot was never
    /// filled. Gaps never get a slot: they appear only *inside* the columns of the
    /// residues they sit alongside.
    slots: Vec<Vec<Option<usize>>>,
}

impl HomologyView {
    /// The column containing the residue at `(sequence, position)`, or `None` if there
    /// is no residue there.
    ///
    /// **This is the residue's homology set *including* the residue itself**, which is
    /// deliberate and is what lets the sets be shared. Two facts make the difference
    /// cancel when two views are compared:
    ///
    /// - A residue's own element is `{sequence, position, gap: false}`, and that is the
    ///   same value in *both* alignments — a residue's identity does not depend on
    ///   where the gaps around it sit. Call it `x`.
    /// - `x` is therefore in both columns, so it is in neither symmetric difference:
    ///   `(Ca \ {x}) Δ (Cb \ {x}) == Ca Δ Cb`.
    ///
    /// So the numerator can be taken from the columns directly. Only the denominator
    /// has to remember the exclusion, by subtracting one from each column's size — see
    /// `distance::compute_symmetric_difference`.
    pub fn column_of(&self, sequence: usize, position: usize) -> Option<&HashSet<Element>> {
        let column = (*self.slots.get(sequence)?.get(position)?)?;
        self.columns.get(column)
    }

    /// The number of residue positions addressable per sequence — the width of the
    /// alignment this view was built from.
    pub fn width(&self) -> usize {
        self.slots.first().map_or(0, |positions| positions.len())
    }

    /// The homology set of the residue at `(sequence, position)`: its column **minus
    /// the residue itself**, which is the definition the metric is stated in terms of.
    ///
    /// Allocates, because it materialises the exclusion that [`HomologyView::column_of`]
    /// exists to avoid. `#[cfg(test)]` for that reason: the distance path must not call
    /// it, or `CODE_REVIEW.md` §2 comes straight back. It is here so the tests can
    /// assert homology sets as they are defined rather than as they are stored.
    #[cfg(test)]
    pub fn homology_set_of(&self, sequence: usize, position: usize) -> Option<HashSet<Element>> {
        let column = self.column_of(sequence, position)?;
        let mut set = column.clone();
        set.remove(&Element {
            seq: sequence as u32,
            position: Some(position as u32),
            gap: false,
        });
        Some(set)
    }

    /// Every column in the view, for tests that need to reason about sharing rather
    /// than about a particular residue.
    #[cfg(test)]
    pub fn columns(&self) -> &[HashSet<Element>] {
        &self.columns
    }
}

/// Builds the [`HomologyView`] for `msa`, resolving sequence ids through `registry`.
///
/// Returns `Err` if `msa` contains a name the registry does not know, or if the
/// alignment is wider than a `u32` position can address.
pub fn homology_view(msa: &Msa, registry: &Registry) -> Result<HomologyView> {
    let width = msa.width();
    if u32::try_from(width).is_err() {
        bail!(
            "The alignment is {} columns wide, which exceeds the {} that a residue position \
             can address",
            width,
            u32::MAX
        );
    }

    // Row index -> registry index. Resolved once up front so a name that is not in the
    // registry is reported before any work is done.
    let mut seq_ids: Vec<u32> = Vec::with_capacity(msa.num_seqs());
    for name in msa.names() {
        match registry.index_of(name) {
            Some(id) => seq_ids.push(id),
            None => bail!(
                "Sequence '{}' is not present in the shared name registry for this comparison; \
                 the registry holds {} sequence(s)",
                name,
                registry.len()
            ),
        }
    }

    // Elements laid out exactly as the alignment is: `elements[row][column]`.
    let elements: Vec<Vec<Element>> = msa
        .rows()
        .iter()
        .zip(seq_ids.iter())
        .map(|(row, &seq)| row_elements(seq, row))
        .collect();

    // Sized by the registry rather than by this alignment's row count so that both
    // views in a comparison have the same outer length. They are in practice equal —
    // `Registry::for_pair` guarantees the same name set — but indexing by registry
    // index means nothing here may assume row index == registry index.
    let mut slots: Vec<Vec<Option<usize>>> = vec![vec![None; width]; registry.len()];
    let mut columns: Vec<HashSet<Element>> = Vec::with_capacity(width);

    for col in 0..width {
        // One set per column, built once and pointed at by every residue in it. See
        // the `HomologyView` docs: this is the fix for `CODE_REVIEW.md` §2, which used
        // to clone this set per residue.
        columns.push(elements.iter().map(|row| row[col]).collect());

        for (row_idx, row) in elements.iter().enumerate() {
            let item = row[col];
            // Gaps get no slot of their own; they only ever appear *inside* the columns
            // of the residues they share a column with.
            if let (false, Some(position)) = (item.gap, item.position) {
                let seq_slot = seq_ids[row_idx] as usize;
                let position_slot = position as usize;
                match slots
                    .get_mut(seq_slot)
                    .and_then(|s| s.get_mut(position_slot))
                {
                    Some(slot) => *slot = Some(col),
                    // Unreachable: `seq_slot` came from the registry and `position` is
                    // an index among this row's residues, so it is < width. Handled
                    // rather than indexed so a future change cannot turn it into a
                    // panic inside a rayon worker.
                    None => bail!(
                        "Internal error: no slot at [sequence {}][position {}] for an alignment \
                         of {} sequence(s) and width {}",
                        seq_slot,
                        position_slot,
                        registry.len(),
                        width
                    ),
                }
            }
        }
    }

    Ok(HomologyView { columns, slots })
}

/// Turns one row of the alignment into elements, assigning positions exactly as the old
/// `Sequence::from_characters` did (`datastructures.rs:206-233`).
///
/// The caller must have checked that `row.len()` fits in a `u32`; `count` is bounded by
/// the number of residues in the row and so cannot overflow.
fn row_elements(seq: u32, row: &[u8]) -> Vec<Element> {
    let mut elements: Vec<Element> = Vec::with_capacity(row.len());
    let mut count: u32 = 0;

    for &byte in row {
        if is_gap(byte) {
            elements.push(Element {
                seq,
                // `None` for a gap that precedes every residue in this row.
                position: count.checked_sub(1),
                gap: true,
            });
        } else {
            elements.push(Element {
                seq,
                position: Some(count),
                gap: false,
            });
            count += 1;
        }
    }

    elements
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::distance::compute_symmetric_difference;
    use crate::msa::read_msa;

    fn msa(records: &[(&str, &str)]) -> Msa {
        Msa::new(
            records.iter().map(|(n, _)| n.to_string()).collect(),
            records.iter().map(|(_, r)| r.as_bytes().to_vec()).collect(),
        )
        .expect("test fixture should be a valid Msa")
    }

    fn residue(seq: u32, position: u32) -> Element {
        Element {
            seq,
            position: Some(position),
            gap: false,
        }
    }

    fn gap(seq: u32, position: Option<u32>) -> Element {
        Element {
            seq,
            position,
            gap: true,
        }
    }

    fn set(elements: &[Element]) -> HashSet<Element> {
        elements.iter().copied().collect()
    }

    /// Reads a pair of FASTA files and runs the whole current pipeline over them:
    /// shared registry, two views, and the Stage 4 distance. Deliberately skips
    /// `standardise`, which is not wired in until Stage 5.
    fn distance_via_new_path(path_a: &str, path_b: &str) -> f64 {
        let msa_a = read_msa(path_a).expect("fixture a should read");
        let msa_b = read_msa(path_b).expect("fixture b should read");
        let registry = Registry::for_pair(&msa_a, &msa_b).expect("fixtures share a name set");
        let view_a = homology_view(&msa_a, &registry).expect("view a should build");
        let view_b = homology_view(&msa_b, &registry).expect("view b should build");

        compute_symmetric_difference(&view_a, &view_b, &registry).expect("fixtures should compare")
    }

    // ---------------------------------------------------------------------------
    // The regression gate
    // ---------------------------------------------------------------------------

    #[test]
    fn new_path_reproduces_the_baseline_distance_exactly() {
        // `0.21428571428571427` is the number the pre-merge binary emits for this pair
        // (recorded in INTEGRATION_NOTES.md, "Baseline"). Reproducing it bit-for-bit
        // through the new reader, the new homology view and the new distance is what
        // proves Stages 1-3 were a rewrite and not a behaviour change.
        //
        // The two intentional divergences of the new path cannot fire here: these
        // fixtures contain no `.` (so the widened gap definition is inert) and no
        // lowercase (so dropping the uppercasing pass is inert). If this assertion ever
        // fails, something unintended moved — do not update the constant.
        //
        // **Stage 4 was expected to move this number and did not.** The denominator
        // changed from `2·|A|` to `|A| + |B|`, but a residue's homology set is its whole
        // column minus itself and every element in a column carries a distinct `seq`, so
        // `|A| = |B| = num_seqs - 1` for every residue and the two expressions agree term
        // for term. Old value 0.21428571428571427 = 6/28; new value 0.21428571428571427
        // = 6/28, same numerator and same denominator. The arithmetic is worked out
        // residue by residue in `distance::tests::the_fixture_pair_distance_is_hand_verified`.
        let distance = distance_via_new_path("test/test.fasta", "test/test2.fasta");
        assert_eq!(
            distance, 0.21428571428571427,
            "the new path must reproduce the pre-merge baseline exactly"
        );
    }

    // ---------------------------------------------------------------------------
    // The headline fix: record order no longer matters
    // ---------------------------------------------------------------------------

    #[test]
    fn reordered_records_give_distance_zero() {
        // *** This is the test that justifies the whole stage. ***
        //
        // `test/test_reordered.fasta` holds exactly the alignment in `test/test.fasta`
        // — same names, same rows, same columns — with the records written out in the
        // order 2, 1, 3 instead of 1, 2, 3. It is the same alignment, so its distance
        // from `test/test.fasta` must be 0.
        //
        // Under the old positional keying it was not. The old reader used each record's
        // ordinal index as its sequence id (`utils.rs:74`), so sequence "2" was element
        // id 0 in one file and id 1 in the other; elements that describe the same
        // residue of the same sequence failed to compare equal, and the tool
        // confidently reported 0.571 for an alignment against itself. That is the
        // number quoted in `CODE_REVIEW.md` §0, and this record order is the one that
        // produces it — the old answer depended on *which* permutation was used, giving
        // 0.5 for two of the other four and 0.167 for a full reversal. Arbitrary
        // sensitivity to record order is the bug; 0.571 is just one of its faces.
        //
        // This test used to assert the old path's 0.5714285714285714 alongside, so that
        // it demonstrated the fix rather than merely asserting the fixed state. That half
        // went with `utils::create_hashsets` in Stage 4; the old numbers survive in this
        // comment only.
        let new_distance = distance_via_new_path("test/test.fasta", "test/test_reordered.fasta");
        assert_eq!(
            new_distance, 0.0,
            "the same alignment with its records reordered must be at distance 0"
        );
    }

    #[test]
    fn registry_indices_are_independent_of_record_order() {
        // The mechanism behind the test above, isolated: the registry hands the same id
        // to the same name whichever file it came from, and whichever way round the
        // pair is passed.
        let a = msa(&[("one", "AA--"), ("two", "A--A"), ("three", "AAA-")]);
        let b = msa(&[("three", "AAA-"), ("one", "AA--"), ("two", "A--A")]);

        let registry = Registry::for_pair(&a, &b).expect("same name set");
        let reversed = Registry::for_pair(&b, &a).expect("same name set");
        assert_eq!(
            registry, reversed,
            "the registry must not depend on argument order"
        );

        let view_a = homology_view(&a, &registry).expect("view a");
        let view_b = homology_view(&b, &registry).expect("view b");
        assert_eq!(
            view_a, view_b,
            "the same alignment written in a different record order must give the same view"
        );
    }

    // ---------------------------------------------------------------------------
    // Case insensitivity
    // ---------------------------------------------------------------------------

    #[test]
    fn case_only_differences_do_not_affect_the_view_or_the_distance() {
        // Stage 1 stopped uppercasing input. Identity excludes the character, so the
        // metric is case-insensitive without any normalisation pass. If the raw byte
        // were ever added back to `Element`, this fails.
        let upper = msa(&[("s1", "AC-G"), ("s2", "A-CG")]);
        let lower = msa(&[("s1", "ac-g"), ("s2", "a-cg")]);

        let registry = Registry::for_pair(&upper, &lower).expect("same name set");
        let view_upper = homology_view(&upper, &registry).expect("view");
        let view_lower = homology_view(&lower, &registry).expect("view");

        assert_eq!(view_upper, view_lower);
        assert_eq!(
            compute_symmetric_difference(&view_upper, &view_lower, &registry).expect("comparable"),
            0.0
        );
    }

    // ---------------------------------------------------------------------------
    // Registry validation
    // ---------------------------------------------------------------------------

    #[test]
    fn for_pair_errs_when_the_second_alignment_is_missing_a_name() {
        let a = msa(&[("x", "AA"), ("y", "A-"), ("z", "-A")]);
        let b = msa(&[("x", "AA"), ("y", "A-")]);

        let err = match Registry::for_pair(&a, &b) {
            Ok(_) => panic!("alignments over different name sets must be rejected"),
            Err(e) => e.to_string(),
        };
        assert!(
            err.contains("'z'"),
            "the error must name the missing sequence, got: {err}"
        );
        assert!(
            err.contains("only in the first"),
            "the error must say which side 'z' is on, got: {err}"
        );
    }

    #[test]
    fn for_pair_errs_when_the_second_alignment_has_an_extra_name() {
        let a = msa(&[("x", "AA"), ("y", "A-")]);
        let b = msa(&[("x", "AA"), ("y", "A-"), ("z", "-A")]);

        let err = match Registry::for_pair(&a, &b) {
            Ok(_) => panic!("alignments over different name sets must be rejected"),
            Err(e) => e.to_string(),
        };
        assert!(
            err.contains("'z'"),
            "the error must name the extra sequence, got: {err}"
        );
        assert!(
            err.contains("only in the second"),
            "the error must say which side 'z' is on, got: {err}"
        );
    }

    #[test]
    fn for_pair_reports_both_directions_at_once() {
        // A name set difference is usually a mistake about *which* files were passed,
        // so both halves of the difference are worth showing in one message.
        let a = msa(&[("shared", "AA"), ("only-a", "A-")]);
        let b = msa(&[("shared", "AA"), ("only-b", "-A")]);

        let err = match Registry::for_pair(&a, &b) {
            Ok(_) => panic!("alignments over different name sets must be rejected"),
            Err(e) => e.to_string(),
        };
        assert!(err.contains("'only-a'"), "got: {err}");
        assert!(err.contains("'only-b'"), "got: {err}");
        assert!(
            !err.contains("'shared'"),
            "shared names are not the problem, got: {err}"
        );
    }

    #[test]
    fn for_pair_accepts_the_same_name_set_in_different_orders() {
        let a = msa(&[("x", "AA"), ("y", "A-")]);
        let b = msa(&[("y", "A-"), ("x", "AA")]);

        let registry = Registry::for_pair(&a, &b).expect("same name set, different order");
        assert_eq!(registry.len(), 2);
        assert_eq!(registry.names(), &["x".to_string(), "y".to_string()]);
        assert_eq!(registry.index_of("x"), Some(0));
        assert_eq!(registry.index_of("y"), Some(1));
        assert_eq!(registry.name_of(1), Some("y"));
        assert_eq!(registry.index_of("nope"), None);
    }

    #[test]
    fn homology_view_errs_on_a_name_outside_the_registry() {
        let a = msa(&[("x", "AA"), ("y", "A-")]);
        let registry = Registry::for_pair(&a, &a).expect("same name set");

        let stranger = msa(&[("x", "AA"), ("q", "A-")]);
        let err = match homology_view(&stranger, &registry) {
            Ok(_) => panic!("a sequence outside the registry must be rejected"),
            Err(e) => e.to_string(),
        };
        assert!(
            err.contains("'q'"),
            "the error must name the sequence, got: {err}"
        );
    }

    // ---------------------------------------------------------------------------
    // Element identity
    // ---------------------------------------------------------------------------

    #[test]
    fn a_gap_after_residue_zero_is_not_a_residue_at_position_zero() {
        // This is the collision the explicit `gap` field exists to prevent. With `Base`
        // retired and the flag inferred rather than stored, these two would be equal,
        // hash to the same bucket, and collapse into one entry in a homology set.
        let residue_zero = residue(7, 0);
        let gap_after_residue_zero = gap(7, Some(0));

        assert_ne!(residue_zero, gap_after_residue_zero);

        let mut both: HashSet<Element> = HashSet::new();
        both.insert(residue_zero);
        both.insert(gap_after_residue_zero);
        assert_eq!(
            both.len(),
            2,
            "a residue and a gap must not collide in a homology set"
        );
    }

    #[test]
    fn identity_distinguishes_sequence_and_position_but_not_the_character() {
        assert_ne!(residue(0, 0), residue(1, 0), "different sequences");
        assert_ne!(residue(0, 0), residue(0, 1), "different positions");
        assert_ne!(
            gap(0, None),
            gap(0, Some(0)),
            "leading gap vs gap after residue 0"
        );
        // Nothing in `Element` can express "an A" versus "a C" — that is the point.
        assert_eq!(residue(0, 0), residue(0, 0));
    }

    #[test]
    fn a_leading_gap_has_no_position() {
        // `count.checked_sub(1)` with `count == 0`: there is no preceding residue.
        let elements = row_elements(3, b"--A-");
        assert_eq!(elements[0], gap(3, None));
        assert_eq!(elements[1], gap(3, None));
        assert_eq!(elements[2], residue(3, 0));
        assert_eq!(elements[3], gap(3, Some(0)));
    }

    // ---------------------------------------------------------------------------
    // Positions and set contents, worked by hand
    // ---------------------------------------------------------------------------

    #[test]
    fn positions_follow_the_old_from_characters_rules() {
        //        col:  0 1 2 3 4
        //        a:    A - B C -
        //        b:    - A B - C
        //
        // Row a: residues at columns 0, 2, 3 take positions 0, 1, 2. The gap at column 1
        // follows residue 0, so position Some(0); the gap at column 4 follows residue 2,
        // so Some(2).
        // Row b: the gap at column 0 precedes every residue, so None. Residues at
        // columns 1, 2, 4 take positions 0, 1, 2. The gap at column 3 follows residue 1,
        // so Some(1).
        assert_eq!(
            row_elements(0, b"A-BC-"),
            vec![
                residue(0, 0),
                gap(0, Some(0)),
                residue(0, 1),
                residue(0, 2),
                gap(0, Some(2)),
            ]
        );
        assert_eq!(
            row_elements(1, b"-AB-C"),
            vec![
                gap(1, None),
                residue(1, 0),
                residue(1, 1),
                gap(1, Some(1)),
                residue(1, 2),
            ]
        );
    }

    #[test]
    fn homology_sets_match_a_hand_worked_example() {
        // Same alignment as the test above. Registry order is sorted by name, so
        // a -> 0 and b -> 1.
        let m = msa(&[("a", "A-BC-"), ("b", "-AB-C")]);
        let registry = Registry::for_pair(&m, &m).expect("same name set");
        let view = homology_view(&m, &registry).expect("view should build");

        assert_eq!(registry.index_of("a"), Some(0));
        assert_eq!(registry.index_of("b"), Some(1));

        // Sequence a, position 0 (column 0): shares its column with b's leading gap,
        // which has no position at all.
        // Sequence a, position 1 (column 2): shares its column with b's residue 1.
        // Sequence a, position 2 (column 3): shares its column with b's gap that
        // follows b's residue 1.
        // Positions 3 and 4 are past the end of a's three residues.
        let row_of = |sequence: usize| -> Vec<Option<HashSet<Element>>> {
            (0..5).map(|p| view.homology_set_of(sequence, p)).collect()
        };

        assert_eq!(
            row_of(0),
            vec![
                Some(set(&[gap(1, None)])),
                Some(set(&[residue(1, 1)])),
                Some(set(&[gap(1, Some(1))])),
                None,
                None,
            ]
        );

        // And the other row, for completeness.
        assert_eq!(
            row_of(1),
            vec![
                Some(set(&[gap(0, Some(0))])),
                Some(set(&[residue(0, 1)])),
                Some(set(&[gap(0, Some(2))])),
                None,
                None,
            ]
        );
    }

    #[test]
    fn a_residue_is_never_in_its_own_homology_set() {
        let m = msa(&[("a", "AC"), ("b", "AC"), ("c", "A-")]);
        let registry = Registry::for_pair(&m, &m).expect("same name set");
        let view = homology_view(&m, &registry).expect("view should build");

        // Column 0 holds three residues; each one's set is the other two.
        assert_eq!(
            view.homology_set_of(0, 0),
            Some(set(&[residue(1, 0), residue(2, 0)]))
        );
        assert_eq!(
            view.homology_set_of(1, 0),
            Some(set(&[residue(0, 0), residue(2, 0)]))
        );
        assert_eq!(
            view.homology_set_of(2, 0),
            Some(set(&[residue(0, 0), residue(1, 0)]))
        );

        // Column 1 holds two residues and one gap; the gap appears in the residues'
        // sets but has no set of its own.
        assert_eq!(
            view.homology_set_of(0, 1),
            Some(set(&[residue(1, 1), gap(2, Some(0))]))
        );
        assert_eq!(
            view.homology_set_of(2, 1),
            None,
            "sequence c has only one residue"
        );

        // The stored form shares one set per column rather than holding one per
        // residue — the whole point of the `CODE_REVIEW.md` §2 fix. The residue *is*
        // in the stored column; it is only absent from the derived homology set.
        assert_eq!(view.columns().len(), m.width());
        assert!(
            view.column_of(0, 0)
                .expect("a residue")
                .contains(&residue(0, 0))
        );
    }

    #[test]
    fn an_identical_pair_is_at_distance_zero() {
        let m = read_msa("test/test.fasta").expect("fixture should read");
        let registry = Registry::for_pair(&m, &m).expect("same name set");
        let view = homology_view(&m, &registry).expect("view should build");

        assert_eq!(
            compute_symmetric_difference(&view, &view, &registry).expect("comparable"),
            0.0
        );
    }

    #[test]
    fn dot_and_dash_gaps_are_interchangeable() {
        // The shared `is_gap` counts both, so an alignment written with `.` produces the
        // identical view to the same alignment written with `-`. Under the old `Base`
        // mapping the `.` was a residue and every position after it was shifted
        // (`CODE_REVIEW.md` §3). Not exercised by the fixtures, hence this test.
        let dashes = msa(&[("s1", "A--C"), ("s2", "-AC-")]);
        let dots = msa(&[("s1", "A..C"), ("s2", ".AC.")]);

        let registry = Registry::for_pair(&dashes, &dots).expect("same name set");
        assert_eq!(
            homology_view(&dashes, &registry).expect("view"),
            homology_view(&dots, &registry).expect("view")
        );
    }
}
