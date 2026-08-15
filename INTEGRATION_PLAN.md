# Integration plan — folding `standardise-msa` into `rusty-metal`

Source: `github.com/dlejeune/standardise-msa` @ `97978c9` (v1.2.0).
Target: this repo @ `6f73a81`.

## Goal

Standardisation stops being a separate tool and becomes the first **stage** of the
distance pipeline. The standardised alignment is a first-class output that can be
emitted on its own, emitted alongside the distances, or not emitted at all.

```
rusty-metal -o dist.csv a.fa b.fa c.fa                       # standardise internally, distances only
rusty-metal -o dist.csv --emit-standardised out/ a.fa b.fa   # both
rusty-metal --standardise-only --emit-standardised out/ a.fa # standardise only, no distances
rusty-metal -o dist.csv --no-standardise a.fa b.fa           # escape hatch: pre-merge behaviour
```

Decisions taken (see conversation):

- Nothing invokes `standardise-msa` by name → no compat binary, no second `[[bin]]`.
- No residue-identity metric on the roadmap → the `Base` alphabet is removed rather
  than kept as a view.

`standardise-msa` introduces **no new dependencies** — `anyhow`, `clap`, `seq_io`,
`log`, `colored`, `simple_logger` are all already here.

## Why this is load-bearing, not cosmetic

Column standardisation genuinely changes the metric. A gap element is assigned
`position: count.checked_sub(1)` (`datastructures.rs:217`) — the index of the residue
preceding it in its row. Permuting columns moves residues within a row, so the same gap
acquires a different `position`, hence a different identity, hence different homology
sets. Residues themselves are unaffected, because the column comparator never swaps two
columns that both hold a residue in the same row — which is exactly the invariant the
residue hash check verifies.

So standardising both inputs strips gap-placement artefacts out of the distance. That is
the point of the merge.

---

## Stage 1 — one shared MSA type (`src/msa.rs`)

Replaces `datastructures.rs` and the reading half of `utils.rs`.

```rust
pub struct Msa {
    names: Vec<String>,   // parallel to rows
    rows:  Vec<Vec<u8>>,  // raw bytes: case- and character-preserving
    width: usize,
}
```

Row-major owns the data: FASTA is row-major in *and* out, standardise writes rows back to
disk, and sorting by name is a row operation.

**`Base` is deleted.** `SequenceElement::eq` compares `base`, `position`, `sequence_id` —
but for a fixed sequence and position the base is already determined, so it contributes
nothing to element identity. The only thing the metric needs from the alphabet is *is this
a gap*. Removing it also removes two things that are actively harmful now that standardise
round-trips sequences to disk:

- `to_ascii_uppercase()` in `read_msa` would destroy input case.
- the `_ => Base::UNKNOWN` catch-all writes `X` where the input said `B`.

In its place, one shared predicate:

```rust
pub fn is_gap(b: u8) -> bool { b == b'-' || b == b'.' }
```

This is also the fix for `CODE_REVIEW.md` §3 — the two tools currently disagree about `.`,
and this collapses them to a single definition both stages read.

Also in this stage, from `CODE_REVIEW.md` §1:

- Read the record name (`record.id_bytes()`, currently commented out at `utils.rs:81`).
- Validate every row has the same length → `Err`, not the index-out-of-bounds panic at
  `utils.rs:102`.
- Delete `counter -= 1` (`utils.rs:84`); it is dead and underflows on empty input,
  pre-empting the graceful error two lines below.
- Delete `datastructures::Sequence` and its `Iterator` impl, which never advances (§4).

Carried over from `standardise-msa`: `residue_hash(&self) -> u64`, hashing names plus
gap-filtered residues.

## Stage 2 — the standardise pass (`src/standardise.rs`)

Ported from `standardise-msa/src/main.rs`, operating on `Msa`.

1. `sort_sequences_by_name(&mut Msa)`
2. `column_order(&Msa) -> Vec<usize>` — bubble sort a permutation of **column indices**,
   with the comparator reading through into `rows`.
3. `apply_permutation(&mut Msa, &[usize])`
4. `standardise(&mut Msa) -> Result<()>` — hash before, hash after, `Err` on mismatch.

Sorting indices rather than materialising `Vec<Column>` avoids a second full copy of every
alignment, which matters given the memory ceiling in `CODE_REVIEW.md` §2.

> **Trap worth writing down in the code:** the column comparator is deliberately
> **non-transitive**. "Equal" means *these two columns may not move past each other*,
> which is not an equivalence relation — A can be movable past B, and B past C, while A is
> pinned against C. It therefore cannot be handed to `sort_by`, which assumes a total
> order and would produce an unspecified result. Bubble sort is the correct algorithm
> here precisely because it only ever swaps *adjacent* movable pairs. Keep it, and say why.

Known cost: O(width² × num_seqs) comparisons. Acceptable for now; optimising it is out of
scope for this pass.

## Stage 3 — homology view (`src/homology.rs`)

From `utils::create_hashsets`, but built as a *derived view* over an already-standardised
`Msa` rather than as the parse target.

```rust
struct Element { seq: u32, position: Option<u32>, gap: bool }
```

Two things to get right:

- **`gap` must be a real field**, not inferred. With `base` gone, a residue at position 0
  and a gap whose preceding residue is 0 would otherwise be indistinguishable and collide
  in the hash set.
- **Identity must exclude the raw byte.** Since `read_msa` no longer uppercases, including
  the character would make `A` and `a` compare unequal across two files. Keying on
  `(seq, position, gap)` makes the metric case-insensitive for free.

`seq` becomes an index into a **shared name registry** built across the pair, not the file
ordinal. Validate that both alignments carry the same name set and `Err` with the
difference if not. This is the fix for the first deferred finding in `CODE_REVIEW.md` §0 —
two files holding the same alignment with records reordered currently report `0.571`
instead of `0`. Standardisation sorts by name so the orders would mostly line up anyway,
but key explicitly rather than depend on the sort.

## Stage 4 — distance (`src/distance.rs`)

- Iterate over the shared registry length rather than `msa_a.num_seqs`
  (`main.rs:48`), and use `|A| + |B|` as the denominator rather than
  `homology_set_size * 2` — the second deferred finding in §0. **This changes numeric
  output**; it is the item to defer if you want the merge to land distance-neutral.
- Fix the `log::warn!` at `distance.rs:71`: it does not fire for the reason it states, and
  it fires once per column, spamming one line per width unit.
- Delete `compute_jaccard_distance` — dead, and holds three latent panics (§4).

## Stage 5 — two-phase pipeline and CLI (`src/main.rs`)

`compare_alignment_pair` currently re-reads both files for every pair. Standardising
inside it would standardise each file once *per pair it appears in*, and would have
several rayon threads writing the same `--emit-standardised` path concurrently. So:

**Phase 1** — read and standardise all inputs once, in parallel over files, hash-checking
each; emit here if asked (one writer per file, no race).
**Phase 2** — pairwise distances over the already-standardised set.

This also removes the redundant O(n²) re-reads that exist today. Raw `Msa`s are roughly
file-sized so holding all N is cheap; the homology views remain per-pair and are the real
memory consumers.

CLI changes:

- `--emit-standardised <DIR>`, `--standardise-only`, `--no-standardise`.
- `--standardise-only` without `--emit-standardised` is a hard error — the run would do
  nothing observable.
- Relax `num_args = 2..` to `1..`, applying the ≥2 check only when distances are on.
- Emitted filenames: `<stem>.standardised.fasta`. **Check for stem collisions** before
  writing — inputs from different directories can share a stem and would silently
  overwrite each other.

CSV output fixes (`CODE_REVIEW.md` §3, §5): `write_all` instead of `write`, an explicit
flush whose error is propagated, quoted fields, and returning `Result` rather than
`.expect("Failed to write result")` inside the loop. Also drop the two `.unwrap()`s at
`main.rs:79`, which panic on the *success* path for a non-UTF-8 path.

## Stage 6 — cleanup and verification

- Clear the 8 build warnings so real ones stand out.
- Tests, using the existing `test/` fixtures:
  - standardisation is idempotent (`standardise(standardise(x)) == standardise(x)`)
  - residue hash is invariant across standardisation
  - same alignment with records reordered → distance `0`
  - ragged FASTA and empty FASTA → `Err`, not panic
  - a known pair's distance, to pin the numeric change from Stage 4
- `justfile`: add a `test` recipe.
- Version bump, README covering the merged CLI.

**Verification before merging:** run the current binary and the new one over the `test/`
fixtures with `--no-standardise`, and confirm the numbers match pre-merge. Then run with
standardisation on and record the deltas — those are the gap-placement artefacts the merge
is meant to remove, and they should be explainable.

---

## Loose ends

- **History.** The code is being restructured rather than copied, so a subtree merge buys
  little. Suggest referencing source commit `97978c9` in the merge commit and archiving
  the GitHub repo. Say if you'd rather have the history reachable and I'll do
  `--allow-unrelated-histories` instead.
- **Package name.** `rusty-metAL` in `Cargo.toml` is unusual casing, and the justfile
  builds `dlejeune/pipeline-utils-rs`. Worth settling on one name while the crate is being
  reshaped anyway.
- **No Dockerfile in the repo**, though the justfile has `build-docker` recipes that
  reference one.
- **Deferred:** the per-residue `HashSet` clone (`CODE_REVIEW.md` §2, O(num_seqs² × width)
  live memory per MSA). Stage 3 rewrites this file, so it is the natural place to switch
  to a shared per-column set with a remove-self view — but it changes the
  symmetric-difference computation and is better as its own change.
