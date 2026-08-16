# Code review — rusty-metal

Reviewed at commit `6f73a81` (2026-08-15). 12 findings. All were verified against the
source, and the crash / wrong-result claims were checked by building the binary and
running it against constructed inputs.

Findings are grouped by what they cost you, not by file.

## Status: all 12 findings discharged as of `b33b7c2`

Re-verified at `b33b7c2` by running the release binary against the inputs each finding
names, not only by reading the code. `just check` (fmt, `clippy -D warnings`, 93 tests)
is clean.

**The file and line references in the body below are pre-merge and mostly dangle.**
`src/utils.rs` and `src/datastructures.rs` no longer exist; the tree is now `msa.rs`,
`homology.rs`, `standardise.rs`, `distance.rs`, `main.rs`. Where each finding ended up:

| § | Finding | Fixed in |
|---|---------|----------|
| 0 | Sequences matched by file position | `homology.rs:98` `Registry::for_pair` — names read at `msa.rs:241`, ids assigned per pair |
| 0 | Distance not symmetric | `distance.rs:76` — bounds come from the registry, denominator is `\|A\| + \|B\|` |
| 0 | Per-column `log::warn!` spam | Gone: a mismatched sequence set is one `Err` from `Registry::for_pair`, reported once per pair |
| 1 | Ragged FASTA panic | `msa.rs:52` `Msa::new` returns `Err` naming the record and both lengths |
| 1 | Empty FASTA underflow | `msa.rs:253` — `read_msa` reports "contains no records" |
| 1 | Non-UTF-8 path `unwrap` | `main.rs:443` `label()` — `to_string_lossy`, falling back to the whole path |
| 2 | Cloned `HashSet` per residue | `homology.rs:205` `HomologyView` — one set per column, indexed, O(num_seqs × width) |
| 3 | `.` not treated as a gap | `msa.rs:221` `is_gap` — the crate's single gap definition, `-` and `.` |
| 3 | Output silently truncated | `main.rs:542` `write_csv` — explicit flush, `write_all`, errors returned not `expect`ed |
| 4 | `Iterator for Sequence` never advances | Removed with `datastructures.rs` |
| 4 | `compute_jaccard_distance` panics | Removed (was dead code) |
| 5 | CSV fields unquoted | `main.rs:506` `csv_escape` — RFC 4180 |
| 5 | `Display for Base` prints quotes | Removed with `datastructures.rs` |
| — | 8 build warnings | Zero warnings; `clippy -D warnings` is enforced by `just check` |

Behavioural checks run at `b33b7c2`, each the scenario the finding describes:

```
empty.fasta            → Error: MSA file test/empty.fasta contains no records   (was: subtract overflow)
ragged.fasta           → Error: record '3' has length 3, first record has 4     (was: index out of bounds)
test vs test_reordered → 0                                                      (was: 0.571)
3-seq vs 2-seq         → error both ways, empty distance field, run continues   (was: 0.25 / 0.5)
`.` gaps vs `-` gaps   → 0                                                      (was: silently wrong)
path containing a comma→ "…/run,v2/aln.fasta",test/test2.fasta,0.5              (was: four fields)
d(a,b) vs d(b,a)       → 0.5 / 0.5
```

The rest of this document is the review as originally written, kept for the reasoning.

---

## 0. Deferred — handled by the incoming merge

These two are the most serious problems in the current tree, but they are already solved
in the other repository being folded in here. Recorded for completeness; **no action
needed in this pass.**

### Sequences are matched by file position, not by name
`src/utils.rs:74`

`read_msa` passes `counter` (the record's ordinal index) as the sequence id, and
`record.id_bytes()` is commented out at `src/utils.rs:81`, so names are never read.
`SequenceElement::eq` (`src/datastructures.rs:160`) compares `sequence_id`, so homology
sets only line up when records appear in identical order in both files.

Two files holding the *same alignment* with records reordered report `0.571` instead of
`0`.

### The distance is not symmetric
`src/distance.rs:44`

The iteration bound is `msa_a.num_seqs` only (passed at `src/main.rs:48`), and the
denominator at `src/distance.rs:87` is `homology_set_size * 2` = `2·|A|` rather than
`|A| + |B|`. Nothing validates that the two alignments contain the same sequence set, so
mismatched counts are silently computed over a truncated overlap:

```
d(3-seq, 2-seq) = 0.25
d(2-seq, 3-seq) = 0.5     # same file pair
```

Related, and *not* obviously covered by the merge: the `log::warn!` at
`src/distance.rs:71` is unreachable for the reason its message states —
`homology_set_a.get(x)` can never be `None` because `x < length`. It actually fires when
*B* lacks the sequence, and it fires once per column, so it spams the log with one line
per width unit. Worth fixing whichever way the merge lands.

---

## 1. Panics that abort the entire run

All three sit inside the `par_bridge` in `process`, so the panic propagates out of
`collect()` and kills the process. That bypasses the `Err(e) => log::warn!(…)` handler at
`src/main.rs:88` which is meant to log one failed comparison and continue — and it
happens before any CSV is written, so every other completed comparison in the run is lost.

### Ragged FASTA — index out of bounds
`src/utils.rs:102`

The transpose loop uses `width` derived from the **first** sequence only, then indexes
every sequence with it:

```rust
column.push(sequences[seq_idx].elements[col]);
```

Any record shorter than the first (unaligned or truncated FASTA) panics with
`index out of bounds: the len is 3 but the index is 3`.

**Fix:** validate that every sequence has the same length and return `Err` instead.

### Empty FASTA — usize underflow
`src/utils.rs:84`

```rust
counter -= 1;
```

With zero records `counter` is `0`, and the subtraction underflows: `attempt to subtract
with overflow`. This fires *before* the graceful `sequences.get(0).with_context(…)` error
on `src/utils.rs:86` gets a chance to run.

The variable is never read afterwards — `num_seqs` uses `sequences.len()` — so the line is
pure dead code.

**Fix:** delete the line. That alone restores the intended error path.

### Non-UTF-8 path — unwrap on the success path
`src/main.rs:79`

```rust
pair[0].file_stem().unwrap().to_str().unwrap(),
```

Panics on a path with no stem, or a non-UTF-8 path (reachable on Linux, including the
Docker image this repo builds). Note this is on the **success** branch, so it destroys a
run that otherwise computed fine. `display()` is already used safely two lines below in
the error arm.

---

## 2. Will OOM on realistic input

### A cloned HashSet per residue
`src/utils.rs:34`

```rust
let mut item_hashset: HashSet<&SequenceElement> = column_hashset.clone();
item_hashset.remove(&item);
```

This materialises one `HashSet` of `num_seqs` entries for *every residue*, i.e.
O(num_seqs² × width) live memory per MSA. Both MSAs are held simultaneously, and that is
multiplied again by the number of pairs running concurrently under the `par_bridge` at
`src/main.rs:68`.

A modest 500-sequence × 5000-column alignment is ~1.2×10⁹ set entries per MSA. This is a
hard ceiling on usable input size, not a micro-optimisation.

**Fix:** keep the single per-column set and a "remove self" view over it, or store column
indices instead of cloned sets.

---

## 3. Silently wrong results

### `.` is not treated as a gap
`src/datastructures.rs:68`

```rust
_ => Base::UNKNOWN,
```

The catch-all accepts any character as a residue. `.` — used as a gap character by several
aligners and by Pfam-derived FASTA — is not mapped to `Base::GAP`, so it falls into the
non-gap arm of `from_characters` (`src/datastructures.rs:220`), receives a real `position`,
and increments `count`. Every position after the first `.` in that sequence is shifted by
one, the homology sets are built against the wrong columns, and **the run completes
successfully and reports a confidently wrong distance.**

**Fix:** reject unrecognised characters, or map `.` to `Base::GAP`.

### Output can be silently truncated
`src/main.rs:131`

The `BufWriter` is never explicitly flushed. It flushes on drop at the end of `main`, where
the error is discarded — so a failure on the final flush (disk full, quota) yields a
truncated CSV **with exit code 0**.

Two smaller issues in the same block:
- `write` (lines 132, 135) instead of `write_all`, so a short write silently drops part of
  a row.
- The per-row `.expect("Failed to write result")` panics rather than returning the
  `Result` that `main` already supports.

---

## 4. Latent traps

### `Iterator for Sequence` never advances
`src/datastructures.rs:140`

```rust
fn next(&mut self) -> Option<Self::Item> {
    self.elements.iter().next().cloned()
}
```

A fresh iterator is built on every call, so this returns `elements[0]` forever and never
yields `None`. Any `for e in sequence` or `.filter(…).collect()` over a `Sequence` hangs and
grows memory without bound — which is exactly the pattern in the commented-out `load_msas`
at `src/utils.rs:52-56`. Because `Sequence` also implements `FromIterator`, this is easy to
fall back into when that code is revived.

**Fix:** consume the elements (delegate to `into_iter` on the `Vec`), or remove the impl.

### `compute_jaccard_distance` has three panics waiting
`src/distance.rs:14`

Currently dead code (the compiler confirms `never used`), so none of this fires today — but
all of it fires the moment the Jaccard metric is wired up:

- `homology_set_a[seq_idx][seq_element_idx]` is indexed without a bounds check inside the
  `while let`. For a sequence with no gaps every slot is `Some`, so `seq_element_idx`
  reaches `width` and panics instead of the loop terminating.
- Line 15's `.unwrap()` on the B side panics whenever the two MSAs' gap patterns differ.
- Line 24 divides by `union_sum`, producing `NaN` for an empty alignment.

---

## 5. Cosmetic

### CSV fields are unquoted
`src/main.rs:136`

An input path containing a comma (or a quote, or a newline) produces a malformed row.
`results/run,v2/aln.fasta` emits four fields where three are expected, and downstream
parsers misattribute the distance.

### `Display for Base` prints quotes
`src/datastructures.rs:186`

```rust
write!(f, "{:?}", c)   // c is a char
```

Emits `'A'` — with the single quotes — rather than `A`. Concatenating a sequence yields
`'A''C''G'`.

---

## Build warnings

The build is clean but emits 8 warnings, including unused imports (`log::log` and
`HashSet` in `distance.rs`, `PathBuf` in `utils.rs`) and the dead
`compute_jaccard_distance`. Worth clearing so that future real warnings stand out.
