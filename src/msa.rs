//! Shared multiple sequence alignment (MSA) representation.
//!
//! This is the row-major, byte-preserving representation that both the distance
//! pipeline and the (forthcoming) standardisation pass build on. Unlike the older
//! `datastructures::MultipleSequenceAlignment` / `Base` machinery in
//! `src/datastructures.rs`, this type:
//!
//! - keeps sequence names (the old reader discarded them entirely), and
//! - keeps raw bytes: no uppercasing, and no lossy mapping through an alphabet enum
//!   that collapses unrecognised characters (e.g. IUPAC ambiguity codes like `B`) to
//!   `X`.
//!
//! Losing either of those is fine for a pure distance metric, but not for a tool that
//! writes alignments back to disk, which is why this type exists.

use anyhow::{bail, Context, Result};
use seq_io::fasta::{Reader, Record};
use std::collections::hash_map::DefaultHasher;
use std::collections::HashSet;
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::io::{BufWriter, Write};
use std::path::Path;

/// A multiple sequence alignment: parallel names and rows, all rows the same width.
///
/// Bytes in `rows` are exactly what was read from the input FASTA (case preserved,
/// no character remapping). The only interpretation this module imposes is [`is_gap`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Msa {
    names: Vec<String>,
    rows: Vec<Vec<u8>>,
    width: usize,
}

impl Msa {
    /// Builds an `Msa` from names and rows, validating that:
    /// - `names` and `rows` have the same length,
    /// - there is at least one sequence,
    /// - every row has the same length (the alignment width), and
    /// - no sequence name occurs twice.
    ///
    /// Ragged input is reported as an `Err` naming the offending record and both
    /// lengths involved, rather than panicking later when the rows are indexed.
    ///
    /// Duplicate names are rejected here so that every later stage can assume names
    /// are unique. Two stages depend on that: [`Msa::sort_sequences_by_name`] would
    /// otherwise produce an order that depends on sort stability rather than on the
    /// data, and the forthcoming name-keyed homology view would silently collapse two
    /// records into one registry slot. Real FASTA does contain repeated ids, so this
    /// has to be a checked error rather than an assumption.
    pub fn new(names: Vec<String>, rows: Vec<Vec<u8>>) -> Result<Msa> {
        if names.len() != rows.len() {
            bail!(
                "Msa::new: names and rows have different lengths ({} names, {} rows)",
                names.len(),
                rows.len()
            );
        }
        if rows.is_empty() {
            bail!("Msa::new: an MSA must contain at least one sequence, got zero");
        }
        let width = rows[0].len();
        for (i, row) in rows.iter().enumerate() {
            if row.len() != width {
                bail!(
                    "Msa::new: record '{}' (index {}) has length {}, but the first record ('{}') has length {}; all rows in an MSA must be the same length",
                    names[i],
                    i,
                    row.len(),
                    names[0],
                    width
                );
            }
        }
        let mut seen: HashSet<&str> = HashSet::with_capacity(names.len());
        for name in names.iter() {
            if !seen.insert(name.as_str()) {
                bail!(
                    "Msa::new: sequence name '{}' occurs more than once; sequence names must be unique",
                    name
                );
            }
        }
        Ok(Msa { names, rows, width })
    }

    pub fn names(&self) -> &[String] {
        &self.names
    }

    pub fn rows(&self) -> &[Vec<u8>] {
        &self.rows
    }

    pub fn width(&self) -> usize {
        self.width
    }

    pub fn num_seqs(&self) -> usize {
        self.rows.len()
    }

    /// Rebuilds every row from `cols`, where `cols[new_index] = old_index`.
    ///
    /// This both **reorders and selects**: `cols` need not cover every column, so a
    /// column omitted from it is dropped from the alignment and the width shrinks to
    /// `cols.len()`. Standardisation uses that to discard all-gap columns, which carry
    /// no alignment information (see `standardise::canonical_columns`).
    ///
    /// Every index must be in range and no index may repeat. A repeated index would
    /// duplicate a residue — silently corrupting the alignment — so it is rejected with
    /// an `Err` rather than tolerated. Omission, by contrast, is legitimate and
    /// unchecked here: dropping a column that still holds a residue would be caught
    /// downstream by [`Msa::residue_hash`], which is what that check is for.
    ///
    /// Note that this is deliberately *unrestricted* as to ordering: it will happily
    /// apply a permutation that moves residues relative to each other within a row.
    /// Whether a particular ordering is legal for standardisation is the caller's
    /// business, and is again what [`Msa::residue_hash`] verifies afterwards.
    ///
    /// Passing an empty `cols` yields a zero-width alignment. That is reachable — an
    /// alignment of nothing but gaps standardises to no columns at all — and is left
    /// representable rather than rejected, so the caller decides what it means. The
    /// distance stage already reports an empty comparison as an error.
    pub fn select_columns(&mut self, cols: &[usize]) -> Result<()> {
        if cols.len() > self.width {
            bail!(
                "select_columns: {} columns requested, but the alignment is only {} columns wide; \
                 with no index repeated that cannot be satisfied",
                cols.len(),
                self.width
            );
        }
        let mut seen = vec![false; self.width];
        for (new_index, &old_index) in cols.iter().enumerate() {
            if old_index >= self.width {
                bail!(
                    "select_columns: cols[{}] = {} is out of range for an alignment {} columns wide",
                    new_index,
                    old_index,
                    self.width
                );
            }
            if seen[old_index] {
                bail!(
                    "select_columns: column {} appears more than once (again at cols[{}]); \
                     repeating a column would duplicate a residue",
                    old_index,
                    new_index
                );
            }
            seen[old_index] = true;
        }

        // One scratch row at a time rather than a transposed copy of the whole
        // alignment: peak extra memory is O(width), not O(width * num_seqs).
        let mut scratch: Vec<u8> = Vec::with_capacity(cols.len());
        for row in self.rows.iter_mut() {
            scratch.clear();
            scratch.extend(cols.iter().map(|&old_index| row[old_index]));
            row.clear();
            row.extend_from_slice(&scratch);
        }
        self.width = cols.len();

        Ok(())
    }

    /// Sorts the sequences by name, keeping `names` and `rows` parallel.
    ///
    /// Names are unique (enforced by [`Msa::new`]), so the resulting order is total
    /// and does not depend on the stability of the sort. Residue content is untouched:
    /// whole `(name, row)` pairs are moved together.
    pub fn sort_sequences_by_name(&mut self) {
        let names = std::mem::take(&mut self.names);
        let rows = std::mem::take(&mut self.rows);

        let mut pairs: Vec<(String, Vec<u8>)> = names.into_iter().zip(rows).collect();
        // Moves the row `Vec`s rather than copying their contents.
        pairs.sort_by(|a, b| a.0.cmp(&b.0));

        let (names, rows): (Vec<String>, Vec<Vec<u8>>) = pairs.into_iter().unzip();
        self.names = names;
        self.rows = rows;
    }

    /// Hashes this MSA's residue content: each sequence's name plus its gap-filtered
    /// bytes, in row order.
    ///
    /// This is independent of gap placement (gaps are filtered out before hashing),
    /// but does depend on residue content, residue order within each sequence, names,
    /// and the order of the sequences themselves. It is the invariant that later
    /// proves standardisation never altered the actual residues.
    ///
    /// It is *not* invariant under arbitrary column permutations — reversing the
    /// columns, say, reverses each row's residues and changes the hash. Invariance
    /// holds exactly for standardisation-legal permutations, those that never move two
    /// columns past each other when both hold a residue in the same row. That is the
    /// whole point: a hash invariant under any permutation would be worthless as a
    /// safety net.
    ///
    /// Because sequences are hashed in row order, this also changes when the sequences
    /// are reordered — see `standardise` for how that interacts with sorting by name.
    pub fn residue_hash(&self) -> u64 {
        let mut hasher = DefaultHasher::new();
        for (name, row) in self.names.iter().zip(self.rows.iter()) {
            name.hash(&mut hasher);
            let residues: Vec<u8> = row.iter().copied().filter(|b| !is_gap(*b)).collect();
            residues.hash(&mut hasher);
        }
        hasher.finish()
    }
}

/// The single shared gap definition for the whole crate: `-` and `.` are gaps,
/// everything else is a residue.
///
/// `.` is a real gap character in Pfam-derived FASTA; earlier code only recognised
/// `-`, which silently shifted residue positions for any sequence containing a `.`.
pub fn is_gap(b: u8) -> bool {
    b == b'-' || b == b'.'
}

/// Reads an MSA from a FASTA file, preserving record names, byte case, and any
/// unrecognised characters exactly as they appear in the input.
///
/// Returns `Err` (not a panic) if the file has no records, or if the records are not
/// all the same length.
pub fn read_msa<P: AsRef<Path>>(path: P) -> Result<Msa> {
    log::info!("Reading msa file: {}", path.as_ref().display());
    let mut reader = Reader::from_path(&path)
        .with_context(|| format!("Failed to open MSA file {}", path.as_ref().display()))?;

    let mut names: Vec<String> = Vec::new();
    let mut rows: Vec<Vec<u8>> = Vec::new();

    while let Some(record) = reader.next() {
        let record = record
            .with_context(|| format!("Failed to parse a record in {}", path.as_ref().display()))?;
        let name = String::from_utf8(record.id_bytes().to_owned()).with_context(|| {
            format!(
                "Record name in {} is not valid UTF-8",
                path.as_ref().display()
            )
        })?;
        names.push(name);
        rows.push(record.owned_seq());
    }

    if names.is_empty() {
        bail!(
            "MSA file {} contains no records (empty FASTA)",
            path.as_ref().display()
        );
    }

    Msa::new(names, rows)
        .with_context(|| format!("Invalid MSA in {}", path.as_ref().display()))
}

/// Writes an MSA back out as FASTA.
///
/// The writer is flushed explicitly and any flush error is propagated, rather than
/// relying on the implicit flush-on-drop (which discards its error).
pub fn write_msa<P: AsRef<Path>>(path: P, msa: &Msa) -> Result<()> {
    let file = File::create(&path)
        .with_context(|| format!("Failed to create {}", path.as_ref().display()))?;
    let mut writer = BufWriter::new(file);

    for (name, row) in msa.names.iter().zip(msa.rows.iter()) {
        seq_io::fasta::write_to(&mut writer, name.as_bytes(), row).with_context(|| {
            format!(
                "Failed to write record '{}' to {}",
                name,
                path.as_ref().display()
            )
        })?;
    }

    writer
        .flush()
        .with_context(|| format!("Failed to flush {}", path.as_ref().display()))?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicU64, Ordering};

    /// Returns a fresh path under the system temp dir for this test run. Cargo runs
    /// tests in parallel, so include a counter (in addition to the process id) to
    /// avoid collisions between tests in this file.
    fn temp_path(label: &str) -> std::path::PathBuf {
        static COUNTER: AtomicU64 = AtomicU64::new(0);
        let n = COUNTER.fetch_add(1, Ordering::Relaxed);
        std::env::temp_dir().join(format!(
            "rusty-metal-msa-test-{}-{}-{}.fasta",
            std::process::id(),
            n,
            label
        ))
    }

    #[test]
    fn round_trip_preserves_names_rows_and_width() {
        let original = read_msa("test/test.fasta").expect("test/test.fasta should read");

        let out_path = temp_path("roundtrip");
        write_msa(&out_path, &original).expect("write_msa should succeed");
        let read_back = read_msa(&out_path).expect("re-reading the written file should succeed");

        assert_eq!(original.names(), read_back.names());
        assert_eq!(original.rows(), read_back.rows());
        assert_eq!(original.width(), read_back.width());

        std::fs::remove_file(&out_path).expect("cleanup should succeed");
    }

    #[test]
    fn names_are_actually_read() {
        let msa = read_msa("test/test.fasta").expect("test/test.fasta should read");
        assert_eq!(
            msa.names(),
            &["1".to_string(), "2".to_string(), "3".to_string()]
        );
    }

    #[test]
    fn case_and_unusual_characters_survive_round_trip() {
        // Lowercase letters and a real IUPAC ambiguity code (B = Asx), with no
        // trailing newline on the file. The old `Base`-based reader would have
        // uppercased this and turned `B` into `X`.
        let fixture = "test/case_and_ambiguity.fasta";

        let msa = read_msa(fixture).expect("fixture should read");
        assert_eq!(msa.rows()[0], b"abBcXy-".to_vec());
        assert_eq!(msa.rows()[1], b"AbBcxY.".to_vec());

        let out_path = temp_path("case-roundtrip");
        write_msa(&out_path, &msa).expect("write_msa should succeed");
        let read_back = read_msa(&out_path).expect("re-reading should succeed");
        assert_eq!(msa.rows(), read_back.rows());

        std::fs::remove_file(&out_path).expect("cleanup should succeed");
    }

    #[test]
    fn ragged_fasta_returns_err_not_panic() {
        let result = read_msa("test/ragged.fasta");
        assert!(
            result.is_err(),
            "ragged FASTA must return Err, not panic or succeed"
        );
    }

    #[test]
    fn empty_fasta_returns_err_not_panic() {
        let result = read_msa("test/empty.fasta");
        assert!(
            result.is_err(),
            "empty FASTA must return Err, not underflow/panic"
        );
    }

    #[test]
    fn is_gap_covers_dash_and_dot_and_rejects_residues() {
        assert!(is_gap(b'-'));
        assert!(is_gap(b'.'));
        assert!(!is_gap(b'A'));
        assert!(!is_gap(b'a'));
        assert!(!is_gap(b'B'));
        assert!(!is_gap(b'X'));
    }

    #[test]
    fn residue_hash_ignores_gap_placement() {
        let a = Msa::new(
            vec!["s1".to_string()],
            vec![b"A--A".to_vec()],
        )
        .expect("valid Msa");
        let b = Msa::new(
            vec!["s1".to_string()],
            vec![b"AA--".to_vec()],
        )
        .expect("valid Msa");

        assert_eq!(a.residue_hash(), b.residue_hash());
    }

    #[test]
    fn residue_hash_ignores_column_order() {
        // `b`'s columns 1 and 2 are `a`'s columns 2 and 1 swapped. That swap is a
        // legal column-reordering move (in each row, at most one of the two columns
        // holds a residue, so no row's residue *order* changes) even though the
        // gap sits in a different column of each row afterwards.
        let a = Msa::new(
            vec!["s1".to_string(), "s2".to_string()],
            vec![b"A-CG".to_vec(), b"AC-G".to_vec()],
        )
        .expect("valid Msa");
        let b = Msa::new(
            vec!["s1".to_string(), "s2".to_string()],
            vec![b"AC-G".to_vec(), b"A-CG".to_vec()],
        )
        .expect("valid Msa");

        assert_eq!(a.residue_hash(), b.residue_hash());
    }

    #[test]
    fn residue_hash_differs_when_residue_changes() {
        let a = Msa::new(vec!["s1".to_string()], vec![b"A--A".to_vec()]).expect("valid Msa");
        let b = Msa::new(vec!["s1".to_string()], vec![b"A--C".to_vec()]).expect("valid Msa");

        assert_ne!(a.residue_hash(), b.residue_hash());
    }

    #[test]
    fn residue_hash_differs_when_name_changes() {
        let a = Msa::new(vec!["s1".to_string()], vec![b"A--A".to_vec()]).expect("valid Msa");
        let b = Msa::new(vec!["s2".to_string()], vec![b"A--A".to_vec()]).expect("valid Msa");

        assert_ne!(a.residue_hash(), b.residue_hash());
    }

    #[test]
    fn msa_new_rejects_mismatched_names_and_rows_lengths() {
        let result = Msa::new(vec!["only-one".to_string()], vec![]);
        assert!(result.is_err());
    }

    #[test]
    fn msa_new_rejects_empty_input() {
        let result = Msa::new(vec![], vec![]);
        assert!(result.is_err());
    }

    #[test]
    fn msa_new_rejects_duplicate_names() {
        // Sorting by name is ambiguous when names repeat, and the name-keyed homology
        // view in the next stage would collide two records into one slot. Uniqueness
        // is established here so nothing downstream has to re-check it.
        let result = Msa::new(
            vec!["dup".to_string(), "other".to_string(), "dup".to_string()],
            vec![b"AA".to_vec(), b"AC".to_vec(), b"CC".to_vec()],
        );
        let err = match result {
            Ok(_) => panic!("duplicate sequence names must be rejected"),
            Err(e) => e.to_string(),
        };
        assert!(
            err.contains("dup"),
            "the error should name the duplicate, got: {err}"
        );
    }

    #[test]
    fn select_columns_applies_new_index_to_old_index_mapping() {
        let mut msa = Msa::new(
            vec!["s1".to_string(), "s2".to_string()],
            vec![b"ABCD".to_vec(), b"abcd".to_vec()],
        )
        .expect("valid Msa");

        // cols[new] = old, so the new column 0 is the old column 3.
        msa.select_columns(&[3, 0, 2, 1]).expect("valid selection");

        assert_eq!(msa.rows()[0], b"DACB".to_vec());
        assert_eq!(msa.rows()[1], b"dacb".to_vec());
        assert_eq!(msa.width(), 4);
    }

    #[test]
    fn select_columns_can_drop_columns() {
        let mut msa = Msa::new(
            vec!["s1".to_string(), "s2".to_string()],
            vec![b"ABCD".to_vec(), b"abcd".to_vec()],
        )
        .expect("valid Msa");

        msa.select_columns(&[2, 0]).expect("a subset is a valid selection");

        assert_eq!(msa.rows()[0], b"CA".to_vec());
        assert_eq!(msa.rows()[1], b"ca".to_vec());
        assert_eq!(msa.width(), 2, "the width must follow the selection");
    }

    #[test]
    fn select_columns_can_empty_the_alignment() {
        // Reachable: an alignment of nothing but gaps standardises to no columns.
        let mut msa = Msa::new(vec!["s1".to_string()], vec![b"ABCD".to_vec()]).expect("valid Msa");

        msa.select_columns(&[]).expect("an empty selection is representable");

        assert_eq!(msa.width(), 0);
        assert_eq!(msa.rows()[0], Vec::<u8>::new());
        assert_eq!(msa.num_seqs(), 1, "dropping columns must not drop sequences");
    }

    #[test]
    fn sort_sequences_by_name_keeps_names_and_rows_parallel() {
        let mut msa = Msa::new(
            vec!["c".to_string(), "a".to_string(), "b".to_string()],
            vec![b"CC".to_vec(), b"AA".to_vec(), b"BB".to_vec()],
        )
        .expect("valid Msa");

        msa.sort_sequences_by_name();

        assert_eq!(msa.names(), &["a".to_string(), "b".to_string(), "c".to_string()]);
        assert_eq!(
            msa.rows(),
            &[b"AA".to_vec(), b"BB".to_vec(), b"CC".to_vec()]
        );
    }
}
