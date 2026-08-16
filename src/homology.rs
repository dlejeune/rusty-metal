//! The homology view: for every residue in an alignment, the set of elements that
//! share its column.
//!
//! A view is derived from an [`Msa`] rather than parsed from a file, so the same
//! alignment can be read once and viewed against several partners.
//!
//! Two properties of element identity govern everything here:
//!
//! - **Sequences are identified by name, not by their position in the file.** Ids come
//!   from a [`Registry`] shared across the pair of alignments being compared, so two
//!   files holding the same alignment with the records written in a different order
//!   produce the same element identities. The registry is also where a pair that is not
//!   over the same set of sequences is rejected.
//! - **Element identity excludes the raw character.** An element is
//!   `(sequence, position, is-gap)` and nothing else, so the same residue written `a` in
//!   one file and `A` in another compares equal, and the metric is case-insensitive
//!   without a normalisation pass over the input. See [`residue_code`] for how the
//!   identity is stored and the test-only `Element` type for the reasoning.
//!
//! A residue is addressed as `[sequence][residue position]`, with `None` meaning the
//! sequence has no residue at that position, i.e. it is past the end of that sequence's
//! residues.
//!
//! # Storage
//!
//! A column holds exactly one element per sequence, so a column is stored as a run of
//! `u32` codes rather than as a set: see [`residue_code`] and [`gap_code`] for the
//! encoding, and `distance::compute_symmetric_difference` for the identity that makes
//! comparing two runs give the same answer as comparing two homology sets. Every residue
//! in a column shares the one stored copy, reached through a slot table.

use crate::msa::{Msa, is_gap};
use anyhow::{Result, bail};
use std::collections::{HashMap, HashSet};

/// One position in one sequence of an alignment: either a residue or a gap.
///
/// **This is the decoded form, and it exists for the tests.** The stored form is a
/// `u32` code: within a column, an element's `seq` is its offset, so only
/// `(position, gap)` needs storing and both fit in one word. The tests work in
/// `Element`s because that is the vocabulary the metric is defined in, and `decode`
/// converts. Nothing on the distance path builds one.
///
/// # Identity
///
/// `Element` is keyed on `(seq, position, gap)` and nothing else. Two points about that:
///
/// - **`gap` is stored rather than inferred.** Without an explicit flag, a residue at
///   position 0 and a gap whose preceding residue is at position 0 would be
///   indistinguishable, and those two very different things would compare equal.
/// - **The residue character is not part of this type.** [`crate::msa::read_msa`]
///   preserves case, because the alignment has to round-trip to disk faithfully, so
///   keying on the character would make the same residue written `a` in one file and `A`
///   in another compare unequal and inflate the distance. Excluding it loses no
///   information: for a fixed sequence and a fixed position the residue is already
///   determined, so the character never contributed to identity.
///
/// # Position
///
/// - a **residue**'s position is its index among the non-gap characters of its row;
/// - a **gap**'s position is the index of the residue *preceding* it in its row, i.e.
///   `count.checked_sub(1)`, which is `None` for a gap before any residue.
///
/// A gap therefore has an identity tied to where it sits relative to its own row's
/// residues, which is why permuting columns changes gap identities (and hence the
/// distance) while leaving residue identities alone — the property the standardisation
/// pass exists to exploit.
#[cfg(test)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Element {
    seq: u32,
    position: Option<u32>,
    gap: bool,
}

// `Element` has no accessors. Its three fields exist to give the value an identity, and
// the metric only asks whether two elements are equal.

// ---------------------------------------------------------------------------------
// The stored form: one u32 per (sequence, column)
// ---------------------------------------------------------------------------------

/// A gap that precedes every residue in its row — the `position: None` case.
///
/// `u32::MAX` is odd, and every residue code is even, so this can never be mistaken for
/// a residue. It cannot collide with a real gap code either: the largest of those is
/// `2·width - 1`, and [`MAX_WIDTH`] keeps that below `u32::MAX`.
const LEADING_GAP: u32 = u32::MAX;

/// The widest alignment the encoding can address.
///
/// A gap after the last residue of a full-width row encodes as `2·(width-1) + 1`, which
/// must stay below [`LEADING_GAP`]. That caps width at `2³¹ - 1` rather than at the
/// `u32::MAX` a bare position would allow. The check is real but the bound is not: an
/// alignment that wide would need two billion columns.
const MAX_WIDTH: usize = (u32::MAX >> 1) as usize;

/// A residue at `position`, encoded.
///
/// Residues take the even codes and gaps the odd ones, so the low bit is the gap flag.
/// Keeping the two apart matters for the reason
/// `a_gap_after_residue_zero_is_not_a_residue_at_position_zero` gives: a residue at
/// position 0 and the gap that follows it are different elements, and conflating them
/// would silently shrink a homology set.
const fn residue_code(position: u32) -> u32 {
    position << 1
}

/// A gap, encoded. `None` means it precedes every residue in its row.
const fn gap_code(position: Option<u32>) -> u32 {
    match position {
        Some(position) => (position << 1) | 1,
        None => LEADING_GAP,
    }
}

/// Whether a code is a residue rather than a gap.
const fn is_residue(code: u32) -> bool {
    // `LEADING_GAP` is odd, so it fails this without a special case.
    code & 1 == 0
}

/// The residue position a code carries. Only meaningful when [`is_residue`] holds.
const fn code_position(code: u32) -> u32 {
    code >> 1
}

/// The decoded form of a code, for tests. `seq` is the element's offset in its column,
/// which the code does not carry and the caller therefore supplies.
#[cfg(test)]
fn decode(seq: u32, code: u32) -> Element {
    if code == LEADING_GAP {
        Element {
            seq,
            position: None,
            gap: true,
        }
    } else {
        Element {
            seq,
            position: Some(code_position(code)),
            gap: !is_residue(code),
        }
    }
}

// ---------------------------------------------------------------------------------
// Why a column of codes is a homology set
// ---------------------------------------------------------------------------------
//
// The elements sharing a column are one per *sequence*, and an element's identity is
// `(seq, position, gap)`. Inside a column the `seq` component is exactly the offset into
// the column, so storing it would be storing the index alongside the value. What is left
// — `(position, gap)` — is what the codes above encode.
//
// The set operations the metric asks for therefore reduce to positional comparison: see
// `distance::compute_symmetric_difference` for `|A Δ B| = 2·#{s : a[s] ≠ b[s]}`.
//
// The invariant this rests on is "one element per sequence per column". It is pinned by
// `distance::tests::every_homology_set_has_size_num_seqs_minus_one` and checked
// exhaustively against materialised `HashSet`s by
// `distance::tests::sharing_the_column_sets_computes_exactly_the_materialised_distance`.
// A move to residues-only or gap-filtered homology sets would break it — the same caveat
// that applies to the `|A| + |B|` denominator, and for the same reason.

/// A name-to-id mapping shared by the two alignments in a comparison.
///
/// The ids handed out here are what makes the metric independent of record order in the
/// input files. `Msa::new` enforces unique names on every [`Msa`], so nothing is
/// de-duplicated here.
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
    /// the error names the sequences involved on both sides rather than reporting only
    /// that the sets differ. Without it a mismatched pair would be compared over
    /// whichever sequences happened to overlap, and `d(a, b)` and `d(b, a)` could
    /// disagree.
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
/// A residue's homology set is *its whole column minus itself*. This stores each column
/// once, shared by every residue in it, plus an index from
/// `[registry sequence index][residue position]` to the column that residue sits in.
/// The "minus itself" is never materialised — see [`HomologyView::column_of`] for why
/// it does not need to be.
///
/// The sequence dimension is indexed by **registry index**, so two views built against
/// the same [`Registry`] are row-aligned regardless of record order in either file.
///
/// # Layout
///
/// Two flat `Vec<u32>`s, both `num_seqs × width`, so a view costs `8 · num_seqs · width`
/// bytes and allocates twice however large the alignment is. `codes` is **column-major**
/// because the hot path reads whole columns: `column_of` hands out a contiguous slice and
/// the comparison walks two of them in step.
///
/// # Why codes rather than sets
///
/// A homology set per residue would cost O(num_seqs² × width) per alignment, and both
/// alignments of a pair are held at once, multiplied again by the number of pairs running
/// concurrently. Sharing one set per column brings that to O(num_seqs × width) but keeps
/// a hash table per column, most of whose buckets are empty at these sizes. Since a
/// column holds one element per sequence and the sequence is the offset into the column,
/// the codes above carry the same information in a flat slice.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HomologyView {
    /// Column-major element codes: `codes[column * num_seqs + sequence]`.
    ///
    /// Every `(column, sequence)` cell is filled — a column holds one element per
    /// sequence, gap or residue alike — so there is no "absent" code and no `Option`.
    codes: Vec<u32>,
    /// `slots[sequence * width + position]` -> the column that residue sits in, or
    /// [`NO_RESIDUE`].
    ///
    /// [`NO_RESIDUE`] means that sequence has no residue at that position, either
    /// because it has fewer residues than the alignment is wide or because the slot was
    /// never filled. Gaps never get a slot: they appear only *inside* the columns of the
    /// residues they sit alongside.
    ///
    /// A column index is `< width ≤ MAX_WIDTH < u32::MAX`, so the sentinel is
    /// unambiguous, and a flat `Vec<u32>` needs neither an `Option` nor one allocation
    /// per sequence.
    slots: Vec<u32>,
    /// The sequence dimension, from the registry rather than from the alignment's own
    /// row count. Held explicitly because both `Vec`s are flat.
    num_seqs: usize,
    /// The alignment's width, i.e. the stride of `slots` and the number of columns.
    width: usize,
}

/// `slots` entry for "this sequence has no residue at this position".
///
/// Distinct from every real column index by the [`MAX_WIDTH`] bound.
const NO_RESIDUE: u32 = u32::MAX;

impl HomologyView {
    /// The column containing the residue at `(sequence, position)`, or `None` if there
    /// is no residue there.
    ///
    /// **This is the residue's homology set *including* the residue itself**, which is
    /// what lets a column be shared by the residues in it. Two facts make the difference
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
    ///
    /// The returned slice is one code per sequence, indexed by registry index, so two
    /// views built against the same [`Registry`] hand back slices that line up
    /// element-for-element. `len()` is the sequence count, which is what the denominator
    /// reads.
    pub fn column_of(&self, sequence: usize, position: usize) -> Option<&[u32]> {
        // Both bounds are checked before the flat index is formed. `slots` is strided by
        // `width`, so an out-of-range `position` would otherwise silently read the next
        // sequence's row rather than miss.
        if sequence >= self.num_seqs || position >= self.width {
            return None;
        }
        let column = self.slots[sequence * self.width + position];
        if column == NO_RESIDUE {
            return None;
        }
        let start = column as usize * self.num_seqs;
        Some(&self.codes[start..start + self.num_seqs])
    }

    /// The number of residue positions addressable per sequence — the width of the
    /// alignment this view was built from.
    pub fn width(&self) -> usize {
        self.width
    }

    /// The column at `(sequence, position)` as a set of [`Element`]s, i.e. the stored
    /// form decoded into the form the metric is defined in.
    ///
    /// `#[cfg(test)]` because materialising this is exactly what the layout exists to
    /// avoid; it is here so the tests can assert against the definition rather than
    /// against the representation.
    #[cfg(test)]
    pub fn column_set_of(&self, sequence: usize, position: usize) -> Option<HashSet<Element>> {
        Some(self.decode_column(self.column_of(sequence, position)?))
    }

    /// The homology set of the residue at `(sequence, position)`: its column **minus
    /// the residue itself**, which is the definition the metric is stated in terms of.
    ///
    /// Allocates, for the same reason as [`HomologyView::column_set_of`], and so is not
    /// for the distance path: materialising a set per residue is what the layout exists
    /// to avoid.
    #[cfg(test)]
    pub fn homology_set_of(&self, sequence: usize, position: usize) -> Option<HashSet<Element>> {
        let mut set = self.column_set_of(sequence, position)?;
        set.remove(&Element {
            seq: sequence as u32,
            position: Some(position as u32),
            gap: false,
        });
        Some(set)
    }

    /// Every column in the view, decoded, for tests that need to reason about sharing
    /// rather than about a particular residue.
    #[cfg(test)]
    pub fn column_sets(&self) -> Vec<HashSet<Element>> {
        self.codes
            .chunks(self.num_seqs)
            .map(|column| self.decode_column(column))
            .collect()
    }

    /// Decodes one column's codes into elements, restoring the `seq` component each
    /// code omits from its offset in the column.
    #[cfg(test)]
    fn decode_column(&self, column: &[u32]) -> HashSet<Element> {
        column
            .iter()
            .enumerate()
            .map(|(seq, &code)| decode(seq as u32, code))
            .collect()
    }
}

/// Builds the [`HomologyView`] for `msa`, resolving sequence ids through `registry`.
///
/// Returns `Err` if `msa` contains a name the registry does not know, or if the
/// alignment is wider than a `u32` position can address.
pub fn homology_view(msa: &Msa, registry: &Registry) -> Result<HomologyView> {
    let width = msa.width();
    if width > MAX_WIDTH {
        bail!(
            "The alignment is {} columns wide, which exceeds the {} that a residue position \
             can address",
            width,
            MAX_WIDTH
        );
    }

    // The sequence dimension comes from the registry, and every column must hold one
    // element per registry entry or the "column length == sequence count" invariant the
    // denominator and the symmetric-difference identity both rest on would not hold.
    // `Registry::for_pair` already guarantees this for any pair that reaches here; the
    // check is so that a future caller building a registry some other way gets an error
    // rather than a quietly wrong distance.
    if msa.num_seqs() != registry.len() {
        bail!(
            "The alignment has {} sequence(s) but the shared registry holds {}; a homology \
             view must have one row per registry entry",
            msa.num_seqs(),
            registry.len()
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

    // Both dimensions are sized by the registry rather than by this alignment's row
    // order, so nothing below may assume row index == registry index.
    let num_seqs = registry.len();
    let mut codes: Vec<u32> = vec![0; num_seqs * width];
    let mut slots: Vec<u32> = vec![NO_RESIDUE; num_seqs * width];

    // One pass per row, scattering into the column-major `codes`, rather than building a
    // transposed copy of the whole alignment and reading that column-wise. `scratch`
    // holds one row's codes and is reused across rows.
    let mut scratch: Vec<u32> = Vec::with_capacity(width);
    for (row, &seq) in msa.rows().iter().zip(seq_ids.iter()) {
        row_codes(row, &mut scratch);
        let seq = seq as usize;

        for (col, &code) in scratch.iter().enumerate() {
            // Column-major: the write stride is `num_seqs`, so this scatters. Reads on
            // the hot path are the contiguous direction, which is the one that matters —
            // this runs once per cell, `column_of` runs once per residue per comparison.
            codes[col * num_seqs + seq] = code;

            // Gaps get no slot of their own; they only ever appear *inside* the columns
            // of the residues they share a column with.
            if is_residue(code) {
                // `position < width`, because it counts residues in a row of `width`
                // cells, so the index is in range and the slot is written at most once.
                slots[seq * width + code_position(code) as usize] = col as u32;
            }
        }
    }

    Ok(HomologyView {
        codes,
        slots,
        num_seqs,
        width,
    })
}

/// Encodes one row of the alignment into `out`, one code per column: a residue's position
/// is its index among the row's non-gap characters, and a gap's is the index of the
/// residue preceding it in the row.
///
/// `out` is cleared first and reused across rows by the caller, so the build allocates
/// one row's worth of scratch for the whole alignment rather than one `Vec` per row.
///
/// The caller must have checked `row.len() <= MAX_WIDTH`; `count` is bounded by the
/// number of residues in the row and so cannot overflow.
fn row_codes(row: &[u8], out: &mut Vec<u32>) {
    out.clear();
    out.reserve(row.len());
    let mut count: u32 = 0;

    for &byte in row {
        if is_gap(byte) {
            // `checked_sub` is `None` for a gap that precedes every residue in this row.
            out.push(gap_code(count.checked_sub(1)));
        } else {
            out.push(residue_code(count));
            count += 1;
        }
    }
}

/// One row as decoded [`Element`]s, for the tests that state positions in those terms.
///
/// Goes through [`row_codes`] rather than duplicating the position rules, so the
/// hand-worked position tests still exercise the encoder the build actually uses.
#[cfg(test)]
fn row_elements(seq: u32, row: &[u8]) -> Vec<Element> {
    let mut codes = Vec::new();
    row_codes(row, &mut codes);
    codes.into_iter().map(|code| decode(seq, code)).collect()
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

    /// Reads a pair of FASTA files and runs the distance path over them: shared registry,
    /// two views, one distance. Skips `standardise`, so these tests pin the metric on the
    /// files as written.
    fn distance_between_files(path_a: &str, path_b: &str) -> f64 {
        let msa_a = read_msa(path_a).expect("fixture a should read");
        let msa_b = read_msa(path_b).expect("fixture b should read");
        let registry = Registry::for_pair(&msa_a, &msa_b).expect("fixtures share a name set");
        let view_a = homology_view(&msa_a, &registry).expect("view a should build");
        let view_b = homology_view(&msa_b, &registry).expect("view b should build");

        compute_symmetric_difference(&view_a, &view_b, &registry).expect("fixtures should compare")
    }

    // ---------------------------------------------------------------------------
    // Pinned values
    // ---------------------------------------------------------------------------

    #[test]
    fn the_fixture_pair_distance_is_pinned() {
        // 6/28 on the 3x4 fixtures, worked out residue by residue in
        // `distance::tests::the_fixture_pair_distance_is_hand_verified`.
        //
        // Neither of the metric's two case-dependent rules can fire on these fixtures:
        // they contain no `.`, so the gap definition is not exercised, and no lowercase,
        // so case insensitivity is not exercised either. A change in this number
        // therefore means something else moved.
        let distance = distance_between_files("test/test.fasta", "test/test2.fasta");
        assert_eq!(distance, 0.21428571428571427);
    }

    // ---------------------------------------------------------------------------
    // Record order does not matter
    // ---------------------------------------------------------------------------

    #[test]
    fn reordered_records_give_distance_zero() {
        // `test/test_reordered.fasta` holds exactly the alignment in `test/test.fasta` —
        // same names, same rows, same columns — with the records written out in the order
        // 2, 1, 3 instead of 1, 2, 3. It is the same alignment, so the distance is 0.
        //
        // This is what keying sequence ids on names rather than on a record's ordinal
        // index buys. Keyed positionally, sequence "2" would be id 0 in one file and id 1
        // in the other, elements describing the same residue of the same sequence would
        // not compare equal, and an alignment would have a nonzero distance from itself —
        // by an amount that depends on which permutation the file happens to use.
        let distance = distance_between_files("test/test.fasta", "test/test_reordered.fasta");
        assert_eq!(
            distance, 0.0,
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
        // Input is not uppercased on read, and element identity excludes the character,
        // so the metric is case-insensitive without a normalisation pass. Adding the raw
        // byte back to `Element` fails this.
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

    // ---------------------------------------------------------------------------
    // The code encoding
    // ---------------------------------------------------------------------------

    #[test]
    fn codes_round_trip_through_decode() {
        // The encoding is only allowed to drop `seq`, which the column offset supplies.
        // Everything else must survive, including the leading-gap `None`.
        for position in [0u32, 1, 2, 7, 1_000_000, (MAX_WIDTH as u32) - 1] {
            assert_eq!(
                decode(4, residue_code(position)),
                residue(4, position),
                "residue at {position}"
            );
            assert_eq!(
                decode(4, gap_code(Some(position))),
                gap(4, Some(position)),
                "gap after residue {position}"
            );
        }
        assert_eq!(decode(4, gap_code(None)), gap(4, None));
    }

    #[test]
    fn a_residue_code_is_never_a_gap_code() {
        // The encoding-level form of `a_gap_after_residue_zero_is_not_a_residue_at_
        // position_zero`: the low bit carries what the `Element::gap` field used to, so
        // a residue and the gap that follows it stay distinct values rather than
        // colliding into one.
        assert_ne!(residue_code(0), gap_code(Some(0)));
        assert!(is_residue(residue_code(0)));
        assert!(!is_residue(gap_code(Some(0))));

        // And the leading-gap sentinel is not mistaken for either. It is odd, so it
        // reads as a gap, and it cannot equal a real gap code because `MAX_WIDTH` caps
        // the largest of those at `2·MAX_WIDTH - 1 < u32::MAX`.
        assert!(!is_residue(LEADING_GAP));
        assert_ne!(LEADING_GAP, gap_code(Some((MAX_WIDTH as u32) - 1)));
        assert!(
            gap_code(Some((MAX_WIDTH as u32) - 1)) < LEADING_GAP,
            "the widest allowed alignment must still leave the sentinel free"
        );
    }

    #[test]
    fn every_code_in_a_row_is_distinct_from_every_other_of_the_same_kind() {
        // Within one row, no two residues share a position and no two gaps share a
        // preceding residue *except* the leading gaps, which are deliberately equal —
        // they are the same element identity, `(seq, None, gap)`.
        let mut codes = Vec::new();
        row_codes(b"--A-A--", &mut codes);
        assert_eq!(
            codes,
            vec![
                LEADING_GAP,
                LEADING_GAP,
                residue_code(0),
                gap_code(Some(0)),
                residue_code(1),
                gap_code(Some(1)),
                gap_code(Some(1)),
            ]
        );
    }

    #[test]
    fn homology_view_errs_when_the_row_count_does_not_match_the_registry() {
        // Unreachable through `Registry::for_pair`, which rejects mismatched name sets
        // outright. Checked anyway because the whole encoding rests on "a column holds
        // one element per registry entry": a short alignment would build short columns,
        // and the distance would come out quietly wrong rather than absent.
        let three = msa(&[("x", "AA"), ("y", "A-"), ("z", "-A")]);
        let registry = Registry::for_pair(&three, &three).expect("same name set");

        let two = msa(&[("x", "AA"), ("y", "A-")]);
        let err = match homology_view(&two, &registry) {
            Ok(_) => panic!("a row count below the registry's must be rejected"),
            Err(e) => e.to_string(),
        };
        assert!(err.contains("2 sequence(s)"), "got: {err}");
        assert!(err.contains("registry holds 3"), "got: {err}");
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

        // The stored form holds each column once rather than one copy per residue. The
        // residue *is* in the stored column; it is only absent from the derived homology
        // set.
        assert_eq!(view.column_sets().len(), m.width());
        assert!(
            view.column_set_of(0, 0)
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
        // identical view to the same alignment written with `-`. Treating `.` as a
        // residue instead would shift every position after it. No fixture contains a
        // `.`, hence this test.
        let dashes = msa(&[("s1", "A--C"), ("s2", "-AC-")]);
        let dots = msa(&[("s1", "A..C"), ("s2", ".AC.")]);

        let registry = Registry::for_pair(&dashes, &dots).expect("same name set");
        assert_eq!(
            homology_view(&dashes, &registry).expect("view"),
            homology_view(&dots, &registry).expect("view")
        );
    }
}
