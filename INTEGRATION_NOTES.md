# Integration notes

Running log for the `standardise-msa` fold-in. Plan lives in `INTEGRATION_PLAN.md`.
Newest stage last. Each entry records what changed, what was decided in the moment, and
anything the next stage needs to know.

## Baseline — before any changes

Commit `6f73a81`. Builds clean apart from 8 warnings.

```
rusty-metAL -o baseline.csv test/test.fasta test/test2.fasta
→ test/test.fasta,test/test2.fasta,0.21428571428571427
```

Fixtures are tiny (3 sequences × 4 columns):

```
test.fasta      test2.fasta
>1 AA--         >1 A-A-
>2 A--A         >2 A--A
>3 AAA-         >3 AAA-
```

`0.21428571428571427` is the number `--no-standardise` must still reproduce at the end of
Stage 6.

*(Originally this paragraph warned that Stage 4's denominator fix would move the value.
It does not — `|A|` and `|B|` are always equal, so `2·|A|` and `|A| + |B|` are the same
number on every input that reaches the computation. The value survived Stages 1-4
untouched. See the Stage 4 entry.)*

`test.fasta` has no trailing newline. Worth keeping as a parser edge case.

---

## Decisions taken during execution

### Stage 1 adds rather than replaces (deviation from the plan)

`INTEGRATION_PLAN.md` has Stage 1 deleting `Base` and `datastructures::Sequence`. Doing
that in Stage 1 would break `distance.rs`, `utils.rs` and `main.rs` simultaneously and
leave the tree uncompilable across three stages.

Instead `src/msa.rs` lands *alongside* the existing types, and the old ones are deleted
when their last consumer migrates (Stage 3 for `create_hashsets`, Stage 4 for
`distance.rs`). Cost is some temporary dead-code warnings; benefit is that
`cargo build && cargo test` is green at the end of every stage, and the baseline number
stays reproducible while the replacement is built underneath it.

### The baseline check has to happen before Stage 4, not after

> **This entry's premise turned out to be wrong. Superseded by the Stage 4 entry —
> the denominator change moves nothing. Kept because the revised verification order
> below was still worth doing, and because the reasoning error is instructive.**

Noticed while writing the baseline entry above. The plan says Stage 6 should verify that
`--no-standardise` reproduces `0.21428571428571427`. It cannot: Stage 4 changes the
denominator from `2·|A|` to `|A| + |B|`, which changes that number regardless of whether
standardisation runs. `--no-standardise` disables the standardise stage, not the Stage 4
fix.

*(What this missed: `|A|` and `|B|` are always equal, so the two denominators are the same
number. A homology set is a column minus its own element, and every element in a column
carries a distinct sequence id, so no two collapse in the `HashSet` and the size is always
exactly `num_seqs - 1` on both sides. The argument above reasoned about the formula's
shape without checking the sets' cardinality. See Stage 4.)*

Revised verification order:

1. **End of Stage 3** — old and new parse paths coexist. Assert the new path reproduces
   `0.21428571428571427` exactly. This is the real regression gate, and it isolates
   "did the rewrite change anything?" from "did the deliberate fix change anything?".
2. **Stage 4** — record the new number and account for the delta by hand on the 3×4
   fixtures, which are small enough to work through on paper.
3. **Stage 6** — pin both numbers in a test so neither moves again silently.

Consequence: `--no-standardise` means "skip standardisation", *not* "reproduce pre-merge
output". It is an escape hatch for isolating the effect of standardisation on a real
dataset, not a bug-compatibility switch. Worth saying so in `--help` text so nobody
expects otherwise.

---

## Stage 1 — shared `Msa` type — DONE

`src/msa.rs` added; `mod msa;` wired into `main.rs`. Nothing else touched. 12 unit tests,
all passing. Old binary still emits `0.21428571428571427`.

API: `Msa::new` (validating), `read_msa`, `write_msa`, `is_gap`, `residue_hash`, plus
`names()` / `rows()` / `width()` / `num_seqs()` with private fields.

Fixes from `CODE_REVIEW.md` that landed here: names are now read (`id_bytes()`, §0),
ragged FASTA errors instead of panicking (§1), empty FASTA errors instead of underflowing
(§1), `.` counts as a gap (§3), and `write_msa` flushes explicitly with the error
propagated (§3).

New fixtures: `test/ragged.fasta`, `test/empty.fasta`, `test/case_and_ambiguity.fasta`.

### Correction to the plan's residue-hash spec

The plan said the residue hash should be "independent of column ordering". That is too
strong, and the agent hit it as a failing test before correcting the fixture. An
*arbitrary* column permutation (reversing every column, say) reorders residues within a
row and therefore must, and does, change the hash. Invariance holds only for
**standardisation-legal** permutations — those that never swap two columns both holding a
residue in the same row.

This is not a caveat, it is the whole point: the hash is precisely the check that
standardisation only ever made legal moves. If it were invariant under arbitrary
permutations it would be worthless as a safety net. The same invariant governs Stage 2's
comparator.

### Open gap: duplicate sequence names

Neither the plan nor the Stage 1 brief said anything about duplicate names, and nothing
currently rejects them. This matters in two places:

- Sorting by name (Stage 2) is ambiguous when names repeat — the result depends on sort
  stability rather than on the data.
- Name-keyed matching (Stage 3) silently collides: two records would map to one registry
  slot.

Real FASTA does contain duplicate ids. **Assigned to Stage 2**, to be rejected in
`Msa::new` so every downstream stage can assume uniqueness.

### Incidental findings from the agent

- `seq_io` strips trailing `\r` throughout, and the existing fixtures are CRLF with no
  final newline. Handled correctly by old and new code alike — not a new hazard, but the
  fixtures are not as plain as they look.
- `record.id_bytes()` is a trait method needing `use seq_io::fasta::Record;`. The old
  `utils.rs` imported only `Reader`, which is very likely why that line was commented out
  rather than fixed — a missing import, not a deliberate choice.
- `Msa` derives `PartialEq`/`Eq`, so Stage 2's idempotence test can compare values
  directly.

---

## Stage 2 — the standardise pass — DONE

`src/standardise.rs` added; `mod standardise;` wired into `main.rs` (the only change to
that file). Mutation primitives and duplicate-name rejection added to `src/msa.rs`.
28 unit tests, all passing, including all 12 from Stage 1 unchanged. Binary still emits
`0.21428571428571427`.

API added to `Msa`: `permute_columns(&[usize]) -> Result<()>` (`perm[new] = old`,
validated) and `sort_sequences_by_name()`. `Msa::new` now also rejects duplicate names,
closing the open gap recorded above.

API in `standardise`: `column_order(&Msa) -> Vec<usize>` and
`standardise(&mut Msa) -> Result<()>`.

### The port is bug-compatible, and that was checked rather than assumed

The original tool was built and run over 400 randomly generated alignments (1–5
sequences × 1–9 columns, alphabet weighted toward gaps, no `.`), alongside the port.
All 400 outputs were byte-identical. That corpus is throwaway; the interesting cases
from it are now fixtures or inline test alignments.

### All-gap columns never move — comment says otherwise

`standardise-msa/src/main.rs:53-58` returns `true` (= pinned) for an all-gap column
under a comment reading "If this column is only gaps, then we can move it around".
Confirmed empirically against the built original:

```
in:  >a A--        out:  >a A--
     >b --B              >b --B
```

Had the comment been the truth, column 0 would have sorted past both others and the
output would have been `--A` / `-B-`. It did not move. The port reproduces this;
`gap_only_column_does_not_move_disputed_behaviour` pins it and says it is disputed.

**Worth deciding deliberately at some point.** Gap-only columns are exactly the columns
whose placement is pure alignment artefact, so pinning them is the case where
standardisation does the least of what it is for. Changing it would change distances on
real data, so it is not a Stage 2 call.

### The residue hash cannot bracket the sequence sort

The brief asked for `residue_hash()` to be verified across the whole operation. It
cannot be: `residue_hash` walks sequences in row order, so sorting by name changes it
for any input not already sorted, and standardisation would reject its own legal work.
The original sidestepped this by hashing *after* the sort.

Resolved with two checks: `residue_hash` brackets the column permutation (the original's
check, and the one that carries the real invariant), and an order-independent XOR digest
of per-sequence `(name, residues)` brackets the whole operation including the sort. The
misleading "independent of column ordering" line in `residue_hash`'s doc comment — the
one flagged in the Stage 1 notes — was corrected at the same time.

### What Stage 3 can now assume

- Sequence names are unique, on every `Msa`, from construction onward. No re-check
  needed when building the name registry.
- After `standardise`, sequences are in byte-wise `String` order by name. Two
  standardised alignments over the same name set are therefore already row-aligned —
  but key on names anyway, as the plan says, because `--no-standardise` skips this.
- `Msa` fields stay private; reordering goes through `permute_columns` /
  `sort_sequences_by_name`.

### Carried forward

- Dead-code warnings are up from 8 to 19 (`standardise.rs` has no caller until Stage 5,
  same situation as `msa.rs` in Stage 1). Stage 6 clears them.
- New fixture: `test/unsorted_names.fasta` — the only fixture whose records are not
  already in name order, so the only one that can catch a no-op sort.
- The comparator is O(width² × num_seqs) and every comparison rescans both columns.
  Fine at fixture scale, and the plan already defers this.

---

## Stage 3 — homology view — DONE

`src/homology.rs` added; `mod homology;` wired into `main.rs` (the only change to that
file). Nothing else touched — `datastructures.rs`, `distance.rs`, `utils.rs` and
`standardise.rs` are all untouched, so the old path still runs and still emits
`0.21428571428571427`. 46 unit tests, all passing, including all 28 from Stages 1 and 2
unchanged.

API: `Element` (private fields, `seq()` / `position()` / `is_gap()` accessors),
`Registry::for_pair(&Msa, &Msa) -> Result<Registry>` with `len` / `is_empty` /
`index_of` / `name_of` / `names`, `type HomologyView = Vec<Vec<Option<HashSet<Element>>>>`
and `homology_view(&Msa, &Registry) -> Result<HomologyView>`.

### The regression gate reproduced the baseline on the first run

`new_path_reproduces_the_baseline_distance_exactly` asserts `0.21428571428571427` from
the new reader plus the new view, with the *old* distance arithmetic held fixed. It
passed first time, unmodified. `old_and_new_paths_agree_on_the_fixture_pair` makes the
stronger version of the same claim: it runs `utils::read_msa` + `create_hashsets` and
`msa::read_msa` + `homology_view` through one shared, generic copy of the old formula and
requires an identical `f64`. Nothing about the rewrite changed the number.

The old formula lives in the *test* module (`old_formula_distance`), not in the library,
because it encodes the two bugs Stage 4 fixes. It is generic over the element type,
which is the only reason it can be run over both `MsaHashSets` and `HomologyView`.

### The headline fix, with the old number pinned alongside

`reordered_records_give_distance_zero`: `test/test_reordered.fasta` is `test/test.fasta`
with the records written 2, 1, 3. New path → **0**. Old path → **0.5714285714285714**,
asserted in the same test so it demonstrates the fix rather than just asserting the fixed
state.

Worth knowing: `0.571` is not *the* wrong answer, it is one of several. All six record
orderings were measured through the old path — identity gives 0, two orderings give 0.5,
one gives 0.5714285714285714 (the figure in `CODE_REVIEW.md` §0, which is why the fixture
uses that ordering), and a full reversal gives 0.16666666666666666. The bug is arbitrary
sensitivity to record order; any single number understates it.

### Decisions taken

- **Registry ids are sorted by name**, so `for_pair(a, b) == for_pair(b, a)` and the ids
  do not depend on which file was passed first. Distances are unaffected by the choice —
  symmetric difference only cares that both views agree — but a canonical registry makes
  the view comparable across calls, which the tests rely on.
- **The view's outer dimension is sized by `registry.len()`, not by the alignment's row
  count**, and is indexed by registry index. Row index and registry index are equal only
  by accident, and nothing in `homology.rs` assumes they are.
- **The residue character is not stored on `Element` at all**, rather than stored and
  excluded from `PartialEq`/`Hash`. Keeping a field out of a derived `Hash` while it sits
  in the struct is exactly the kind of thing a later `#[derive]` change silently breaks.
  Case insensitivity therefore holds by construction.

### What Stage 4 needs to know

- `distance.rs` should iterate `0..registry.len()`, which removes the asymmetry at
  `CODE_REVIEW.md` §0 by construction: `Registry::for_pair` has already rejected any pair
  whose name sets differ, so `d(a, b)` and `d(b, a)` can no longer disagree, and the
  `log::warn!` at `distance.rs:71` becomes genuinely unreachable rather than
  misdescribed. It should be deleted, not moved.
- Three tests in `homology.rs` reach into the old types and must be deleted *with* them:
  `old_and_new_paths_agree_on_the_fixture_pair`, the second half of
  `reordered_records_give_distance_zero`, and `old_formula_distance` itself once no test
  uses it. They are marked in place with `*** Stage 4 ... ***`.
- `new_path_reproduces_the_baseline_distance_exactly` pins a number the new denominator
  will change. Replace it with the recomputed value; do not delete the pin.
- The per-residue `HashSet` clone is preserved verbatim, with a `TODO` pointing at
  `CODE_REVIEW.md` §2. It was kept identical on purpose so this stage could be shown to
  be behaviour-preserving; that argument no longer needs it after Stage 4.
- Dead-code warnings are up from 19 to 28 — the whole of `homology.rs` has no caller
  until Stage 4/5, same situation as the two stages before it.
- New fixture: `test/test_reordered.fasta`.

---

## Stage 4 — distance, and the deletion of the old path — DONE

`src/distance.rs` rewritten against `HomologyView` + `Registry`; `src/main.rs`'s
`compare_alignment_pair` repointed at it; `src/datastructures.rs` and `src/utils.rs`
deleted along with their `mod` lines. 57 unit tests, all passing. Build warnings 27 → 10,
test-build warnings 10 → 2; everything left is `standardise.rs` and the `Msa`/`Element`
accessors it needs, which Stage 5 wires up.

New signature:

```rust
pub fn compute_symmetric_difference(&HomologyView, &HomologyView, &Registry) -> Result<f64>
```

`compute_jaccard_distance` is gone (dead, three latent panics, `CODE_REVIEW.md` §4), as
is the `log::warn!` at the old `distance.rs:71` and the large commented-out block below
the return.

### The denominator fix does not move the number, and that is not luck

**This is the headline finding of the stage, and it contradicts what the plan and the
baseline note above both predicted.** The plan says the `2·|A|` → `|A| + |B|` change
"**changes numeric output**"; it does not, on any input that reaches the computation.

A residue's homology set is its whole column minus itself. Every element in a column
carries a distinct `seq`, so no two collapse in the `HashSet`, and therefore
`|A(r)| = |B(r)| = num_seqs - 1` for **every** residue `r`. `Registry::for_pair` has
already forced both alignments to the same sequence count. So `2·|A|` and `|A| + |B|`
agree term for term, not just in sum.

Worked by hand on the 3×4 fixtures before running anything (the full per-residue table is
in `distance::tests::the_fixture_pair_distance_is_hand_verified`): numerator 6,
denominator 28, `6/28 = 3/14 = 0.21428571428571427`. Old numerator 6, old denominator 28,
same value. The program agreed with the hand calculation on the first run, in both
argument orders.

Two consequences:

- The baseline note above ("Stage 4 changes the denominator … which changes that number")
  is **wrong**, and `--no-standardise` therefore *can* be expected to reproduce
  `0.21428571428571427` at the end of Stage 6 after all. Any delta seen there is
  standardisation, not this.
- `|A| + |B|` is still the correct thing to write, because the invariant it relies on is
  not one the metric should depend on. The deferred `CODE_REVIEW.md` §2 change (a shared
  per-column set with a remove-self view) — or any move to gap-filtered / residues-only
  sets — makes the two sides' sizes genuinely differ, and `2·|A|` would start producing
  asymmetric results the moment it lands.
  `distance::tests::every_homology_set_has_size_num_seqs_minus_one` pins the invariant so
  that whoever breaks it is told which pinned value is expected to move and why.

The *asymmetry* was real regardless: it came from `main.rs` passing `msa_a.num_seqs` as
the iteration bound over a truncated overlap, which is what produced the
`d(3-seq, 2-seq) = 0.25` / `d(2-seq, 3-seq) = 0.5` pair in `CODE_REVIEW.md` §0. That is
now an `Err` from `Registry::for_pair`, and `main.rs` logs it and continues to the next
pair (the CSV row gets an empty distance field), rather than reporting a number computed
over whichever sequences happened to overlap.

### 0/0 is now an error rather than `NaN`

Not asked for, but it fell out of touching the division. Two alignments with one sequence
each have an empty homology set for every residue, so the denominator is 0; likewise for
alignments made entirely of gaps. The old code divided anyway and returned `NaN`, which
`main.rs` writes into the CSV as `NaN` — indistinguishable from a real result to anything
downstream. It now returns `Err`, which the existing per-pair error handler already knows
how to report.

### Tests

Deleted: `old_and_new_paths_agree_on_the_fixture_pair` (the only whole test removed — it
ran `utils::read_msa` + `create_hashsets`, which no longer exist), the old-path half of
`reordered_records_give_distance_zero`, and the `old_formula_distance` helper. All three
were marked in place by Stage 3. The old numbers they asserted (`0.5714285714285714` for
the reordered pair) survive in the surviving test's comment.

`new_path_reproduces_the_baseline_distance_exactly` keeps its `0.21428571428571427` pin —
recomputed rather than carried over, and its comment now records both values and why they
are the same.

12 new tests in `distance.rs`, notably `the_distance_is_symmetric_on_the_fixture_pair` and
`the_distance_is_symmetric_for_differing_gap_patterns`, which assert exact `f64` equality
both ways round — the property that was broken.

### What Stage 5 needs to know

- `compare_alignment_pair` still re-reads both files per pair and does **not**
  standardise. Splitting it into the two phases is Stage 5's job and was deliberately not
  started here.
- Every `CODE_REVIEW.md` §1/§3/§5 item in `main.rs` is untouched: the two `.unwrap()`s at
  the file-stem log line, `write` vs `write_all`, the unflushed `BufWriter`, the
  `.expect("Failed to write result")`, and unquoted CSV fields.
- The `TODO(CODE_REVIEW.md §2)` per-residue `HashSet` clone in `homology.rs` is still
  there. Stage 3 kept it verbatim to argue behaviour preservation; that argument is now
  spent, so it is free to change — but see the note above about what it does to the
  denominator invariant.

---

## Stage 5 — two-phase pipeline and CLI — DONE

`src/main.rs` rewritten. Nothing else in `src/` was touched — `msa.rs`, `standardise.rs`,
`homology.rs` and `distance.rs` are byte-identical to their Stage 4 state, and no
dependency was added. 80 unit tests, all passing, including all 57 from Stages 1-4
unchanged. Build warnings 10 → 2, test-build warnings 2 → 2; both remaining ones are
`homology.rs` accessors (`Element::seq` / `position` / `is_gap`, `Registry::is_empty`)
that genuinely have no caller yet, which is Stage 6's list.

### The headline finding: the fixture distance moves, and it moves *up*

```
rusty-metAL -o d.csv --no-standardise test/test.fasta test/test2.fasta
→ test/test.fasta,test/test2.fasta,0.21428571428571427     (6/28, unchanged from baseline)

rusty-metAL -o d.csv test/test.fasta test/test2.fasta
→ test/test.fasta,test/test2.fasta,0.5                     (14/28)
```

Worked by hand on the 3×4 fixtures before running anything, and the program agreed on the
first run. `test/test.fasta` is already canonical — its records are in name order, and its
columns 2 and 3 have their first gap in the same row, so they compare `Equal` and never
swap — so standardisation is a no-op on it. `test/test2.fasta` is not: `column_order`
returns `[0, 1, 3, 2]`.

```
test.fasta      test2.fasta     standardise(test2)
>1 AA--         >1 A-A-         >1 A--A
>2 A--A         >2 A--A         >2 A-A-
>3 AAA-         >3 AAA-         >3 AA-A
```

Registry sorted by name (`"1"`→0, `"2"`→1, `"3"`→2); `r(s,p)` a residue, `g(s,p)` a gap
whose position is the index of the residue preceding it in its row:

```
seq pos | A(r)               | B(r) standardised   | |AΔB| | |A|+|B|
--------+--------------------+---------------------+-------+--------
0   0   | {r(1,0), r(2,0)}   | {r(1,0), r(2,0)}    |   0   |   4
0   1   | {g(1,0), r(2,1)}   | {g(1,1), r(2,2)}    |   4   |   4
1   0   | {r(0,0), r(2,0)}   | {r(0,0), r(2,0)}    |   0   |   4
1   1   | {g(0,1), g(2,2)}   | {g(0,0), g(2,1)}    |   4   |   4
2   0   | {r(0,0), r(1,0)}   | {r(0,0), r(1,0)}    |   0   |   4
2   1   | {r(0,1), g(1,0)}   | {g(0,0), g(1,0)}    |   2   |   4
2   2   | {g(0,1), g(1,0)}   | {r(0,1), g(1,1)}    |   4   |   4
--------+--------------------+---------------------+-------+--------
                                               sum:   14       28
```

Against the unstandardised table in
`distance::tests::the_fixture_pair_distance_is_hand_verified` (numerator 6), four of the
seven residue slots got worse and none got better. The denominator is untouched at 28, as
it must be: `|A(r)| = |B(r)| = num_seqs - 1` regardless of gap placement (Stage 4's
invariant), and standardisation changes neither the sequence count nor the residue count.

**Why up rather than down, since the merge is motivated as "strips gap-placement
artefacts out of the metric".** It does strip them, but stripping is not shrinking.
Standardisation moves each alignment to *its own* canonical column layout independently;
it does not move the two alignments toward each other. Here one input was already at its
canonical layout and the other was not, so the only effect was to rearrange test2's last
two columns *away* from test's. Every gap in test2's rows 0 and 1 acquired a new
preceding-residue index, so those gap identities stopped matching test's, and three of the
four differing slots went from partly shared to fully disjoint.

The case standardisation is actually *for* is a pair differing only by a legal column
permutation, and that case does go to 0. Two genuinely different alignments can move
either way.

> **The claim that a legal column permutation "does go to 0" was wrong, and was never
> tested — it was asserted from the shape of the algorithm. It failed on 78% of legal
> permutations. Superseded by the Stage 6 entry, which replaces the ordering rule and
> makes the claim true. Kept because the error is the reason Stage 6 exists.** `0.5` was not adjusted toward `0.214` and nothing was tuned to preserve the
old value; `end_to_end_with_standardisation_pins_the_new_fixture_distance` pins it with
the table above in a comment, and `end_to_end_without_standardisation_reproduces_the_baseline`
pins `0.21428571428571427` beside it so the two can never drift apart silently.

An incidental confirmation from a three-input run: `test.fasta` vs `test_reordered.fasta`
is `0` under standardisation, and `test2.fasta` vs `test_reordered.fasta` is `0.5` —
consistent, since `test_reordered` *is* `test` as an alignment.

### The `RunPlan` shape

Validation is a pure `RunPlan::from_args(&Args) -> Result<RunPlan>`, so the whole CLI
contract is unit-testable without `assert_cmd`. Illegal combinations are unrepresentable
rather than re-checked:

```rust
pub enum Standardisation { Skip, InMemory, Emit(PathBuf) }

pub enum Mode {
    StandardiseOnly { emit_dir: PathBuf },
    Distances { output_fp: PathBuf, standardisation: Standardisation },
}

pub struct RunPlan { inputs: Vec<PathBuf>, num_threads: usize, mode: Mode }
```

- `--standardise-only` without `--emit-standardised`: `emit_dir` is a `PathBuf`, not an
  `Option`, so the variant cannot be built without one.
- `--emit-standardised` with `--no-standardise`: `Skip` carries no directory, so there is
  nowhere to put one.
- `--standardise-only` with `--no-standardise`: disjoint variants.
- `-o` required only for distances: `Mode::Distances` holds it, `StandardiseOnly` has no
  field for it. `Args::output_fp` is `Option<PathBuf>` at the parser level and made
  mandatory in `from_args`.

`num_args` is relaxed to `1..` and the "at least 2" rule moved into `from_args`, applied
only in the distance branch. Accessors `emit_dir()` / `standardises()` / `output_fp()`
answer the questions the two phases actually ask, so no consumer matches on the enums.

`-o` alongside `--standardise-only` is accepted but ignored, with a `log::warn!` in `main`
saying so — coherent but pointless, and silently discarding an explicit output path would
be worse than mentioning it.

### Stem collisions are checked at plan time, before anything is opened

`check_stem_collisions` runs inside `from_args` whenever an emit directory is present, so
`a/aln.fasta b/aln.fasta` fails before the output directory is even created. The failure
mode it prevents is silent data loss — both map to `aln.standardised.fasta` and the second
write destroys the first — which is not something to discover halfway through a parallel
phase 1. The error lists every colliding input, grouped by stem in a `BTreeMap` so
multi-collision reports come out in a stable order. It deliberately does *not* fire when
nothing is emitted, where a shared stem is harmless.

Output paths are built by pushing `".standardised.fasta"` onto the stem's `OsString`
rather than going through `&str`, so a non-UTF-8 input path is named correctly instead of
being rejected.

### Rayon: a scoped pool per run, not `build_global`

`build_global` succeeds at most once per process, so the old `process()` — which called it
on every invocation — could only ever run the pipeline once, and
`the_pipeline_can_be_run_twice_in_one_process` would have failed on the second call for
reasons unrelated to what it tests.

Chosen: `ThreadPoolBuilder::new().num_threads(n).build()?` per run, with both phases inside
`pool.install(...)`. Rejected: tolerating the already-initialised error, and a `OnceLock` —
both work, but both silently ignore `--num-threads` on every run after the first, which is
a real behaviour change hiding inside a test-only fix. `install` makes the scoped pool the
one every nested `par_iter`/`par_bridge` uses, including the one inside
`distance::compute_symmetric_difference`, and it is torn down when `run` returns.
`num_threads(0)` still means "let rayon decide", matching what `--num-threads 0` documents.
The test asks for 2 threads on the first run and 3 on the second, so a `OnceLock`
regression would be caught rather than passing silently.

### Phase boundaries and the failure policy

Phase 1 (`load_inputs`) reads, standardises and emits, in parallel *across files*, once
each. `par_iter().map(...).collect::<Result<Vec<_>>>()` preserves input order, which
phase 2 relies on to pair an alignment back up with its path. Emission happens here, so
there is exactly one writer per output path — the race the plan warned about cannot arise.

The two phases fail differently on purpose, and both are tested:

- A **file** failure is hard. It invalidates every pair it appears in, so continuing would
  produce a result set silently missing an arbitrary subset of comparisons.
  `an_unreadable_input_is_a_hard_error` covers both a file that parses but is not a valid
  alignment (`test/ragged.fasta`) and one that does not exist, and asserts no CSV is
  written.
- A **pair** failure is soft, as before: logged, empty distance field, run continues.
  `a_single_failing_pair_is_logged_and_the_run_continues` uses
  `test/case_and_ambiguity.fasta`, whose names (`seq1`, `seq2`) do not match the fixtures'
  (`1`, `2`, `3`), so two of the three pairs are rejected by `Registry::for_pair` and the
  third still reaches the CSV.

`compare_alignment_pair` now takes `&Msa, &Msa` rather than two paths. The O(n²) re-reads
are gone: N files are read N times in total instead of N-1 times each.

### `CODE_REVIEW.md` §1/§3/§5 in `main.rs` — all cleared

- The two `.unwrap()`s on the file-stem log line are gone, replaced by a `label()` helper
  that falls back to the whole path. They were on the *success* branch inside a rayon
  worker, so a stemless or non-UTF-8 path destroyed an otherwise-good run and took every
  other completed comparison with it.
- `write` → `write_all` throughout.
- `write_results` returns `Result` instead of `.expect("Failed to write result")`.
- `write_csv` flushes the `BufWriter` explicitly and propagates the error, so a failing
  final flush is no longer a truncated CSV with exit code 0.
- Fields are escaped per RFC 4180 (`csv_escape`): quoted only when they contain `,`, `"`,
  `\r` or `\n`, with embedded quotes doubled. Ordinary paths are unchanged, so existing
  output is byte-identical.

`write_results` is generic over `W: Write` so
`a_row_with_an_awkward_path_stays_three_fields` can assert the exact bytes against a
`Vec<u8>` without touching the disk.

No `unwrap`/`expect` remains anywhere in non-test `main.rs`, for the reason §1 gives: a
panic inside `par_bridge` aborts the process and bypasses the per-pair handler entirely.

### Tests added (23)

Plan validation: `plan_accepts_distances_with_internal_standardisation`,
`plan_accepts_distances_and_emission_together`,
`plan_accepts_standardise_only_with_a_single_input_and_no_output`,
`plan_accepts_the_no_standardise_escape_hatch`,
`plan_rejects_standardise_only_without_emit_standardised`,
`plan_rejects_standardise_only_with_no_standardise`,
`plan_rejects_emit_standardised_with_no_standardise`,
`plan_rejects_missing_output_when_distances_are_computed`,
`plan_rejects_a_single_input_when_distances_are_computed`.

Stem collisions: `stem_collision_is_detected_across_directories`,
`distinct_stems_do_not_collide`, `a_stem_collision_is_caught_at_plan_time_when_emitting`,
`a_stem_collision_is_not_an_error_when_nothing_is_emitted`.

CSV: `csv_fields_are_escaped_per_rfc_4180`,
`a_row_with_an_awkward_path_stays_three_fields`.

End to end: `end_to_end_without_standardisation_reproduces_the_baseline`,
`end_to_end_with_standardisation_pins_the_new_fixture_distance`.

Phase 1: `phase_one_emits_one_standardised_file_per_input`,
`standardise_only_writes_no_csv`,
`standardised_output_path_names_the_file_from_the_stem`.

Failure policy: `an_unreadable_input_is_a_hard_error`,
`a_single_failing_pair_is_logged_and_the_run_continues`.

Rayon: `the_pipeline_can_be_run_twice_in_one_process`.

The plan tests go through `Args::try_parse_from` rather than constructing `Args` directly,
so they also pin the clap configuration itself — that `-o` is genuinely optional at the
parser level and that a single positional is genuinely accepted. Constructing `Args` by
hand would have let `required = true, num_args = 2..` survive untested.

### What Stage 6 still has

From `CODE_REVIEW.md`, still outstanding:

- **§2, the per-residue `HashSet` clone** in `homology.rs`, still marked with its `TODO`.
  This is the only remaining correctness-adjacent finding and it is the hard memory
  ceiling. Note that fixing it makes `|A|` and `|B|` genuinely differ, which is when
  Stage 4's `|A| + |B|` denominator starts mattering and when the two pinned fixture
  values above are expected to move —
  `distance::tests::every_homology_set_has_size_num_seqs_minus_one` is the tripwire.
- **The 2 remaining dead-code warnings**: `Element::seq` / `position` / `is_gap` and
  `Registry::is_empty` have no caller. `Registry::name_of` / `names`, `write_msa`,
  `standardise` and the rest of `standardise.rs` gained callers this stage and their
  warnings are gone.
- §4's `Iterator for Sequence` and `compute_jaccard_distance`, and §5's
  `Display for Base`, went with `datastructures.rs` / `distance.rs` in Stages 1-4 and
  need nothing further.

Also outstanding from the plan's Stage 6 list: the `justfile` `test` recipe, the version
bump, and a README covering the merged CLI. The `rusty-metAL` / `pipeline-utils-rs` naming
question in "Loose ends" is untouched.

Not a `CODE_REVIEW.md` item, but worth carrying: the disputed all-gap-column pinning from
Stage 2 now has a visible consequence, because standardisation is finally on the default
path. Gap-only columns are exactly the columns whose placement is pure artefact, so
pinning them is where standardisation does least of what it is for — and every user now
gets that behaviour by default rather than only those who ran the separate tool.

---

## Stage 6 â€” canonical ordering, and the confluence finding â€” DONE

`src/standardise.rs`'s ordering rule replaced outright; `Msa::permute_columns` generalised
to `Msa::select_columns`; three `homology.rs` accessors deleted and two scoped to tests;
`src/confluence_probe.rs` (the temporary Stage 6 Part A harness) removed. 91 unit tests,
all passing. **Build warnings 2 â†’ 0, test-build warnings 2 â†’ 0.**

Version bumped to `0.2.0` in `Cargo.toml` and the `justfile`; `just test` and `just check`
recipes added; `README.md` written, covering the merged CLI.

### The headline finding: standardisation was not confluent

Stage 5 asserted, in this file and in a test comment, that "a pair differing only by a
legal column permutation â€¦ does go to 0". **That was false**, and it was the whole
justification for folding `standardise-msa` in.

A property search over random alignments, restricted to column shuffles that are legal
(never exchanging two columns that both hold a residue in the same row â€” i.e. the same
alignment, written differently):

```
non-identity legal shuffles:            13405
standardised layouts disagreed:         10428   (78%)
  ...and produced a nonzero distance:    4837   (36%)
  ...with no all-gap column involved:    1703   (13%)
largest distance seen for an alignment against itself:  1.0
```

Shrunk to three sequences and three columns, nothing disputed involved:

```
   s0  A - -      c0 vs c2:  s1 holds a residue in both  -> pinned
   s1  A - A      c0 vs c1:  free
   s2  - A -      c1 vs c2:  free
       c0 c1 c2
```

`c1` is free to sit before `c0`, between, or after `c2`, and nothing chose. Two legal
orderings of this alignment standardised to two different layouts, distance **0.625**.

### It was three defects, not one, and only the third is interesting

Each was isolated by patching it out and re-running the search.

| Comparator | Counterexamples |
| --- | --- |
| as shipped | 10428 / 13405 |
| direction flipped | unchanged rate (22595 / 29446) |
| full-column key, ties broken | 8871 |
| ...and all-gap pinning removed | 2760 |
| topological sort (what landed) | **0 / 29446** |

1. **The direction was backwards.** The original keyed on the row index of the first
   *gap*, ascending, so a column with a residue at the top sorted *last* â€” gaps stacked
   up and to the left. On `s0 A- / s1 -A` the original emits `-A / A-`. That is the
   mirror image of the intent. Note this is the *second* inverted thing in that
   comparator; the all-gap branch already had a comment saying the opposite of what its
   code did (Stage 2). Flipping the sign fixes the direction and changes the
   counterexample rate not at all.
2. **The key was a single number**, so two columns tied constantly, and a tie meant "no
   swap", which meant "keep whatever order the file had".
3. **The real problem is structural.** Pinning is a precedence constraint, not an
   ordering, and a sort that only exchanges adjacent items has no unique fixed point
   under one. Even with a total key and no all-gap pinning, 2760 counterexamples
   survived. No comparator can fix this, because when a free column has three legal
   resting places there is no pairwise question whose answer picks one.

The irony worth recording: `column_order` carried a long, *correct* comment explaining
why `sort_by` must not be used here, because the relation is non-transitive. The
reasoning was sound and stopped one step short â€” the same non-transitivity means bubble
sort's answer is not unique either.

### The rule that replaced it

Decided by the repo owner, not derived here: **prefer the order that fills the higher
rows from the left.** Formally, order columns by reading top to bottom with a residue
sorting before a gap, and take the lexicographically smallest.

That is a total key, but it cannot simply be sorted on â€” pinning forbids the plain lex
order in **3894 of 5000** random alignments. So `canonical_columns` builds the order
instead: greedy topological sort over the pinning DAG, emitting the lex-smallest column
whose pinned predecessors are all placed.

**Why this is confluent, by construction rather than by testing:** pinning is a property
of the columns themselves, and a legal permutation never reorders a pinned pair, so the
precedence DAG is identical whichever legal ordering of an alignment you are handed. A
greedy topological sort that breaks every choice by a total key therefore depends only on
that DAG.

**Why no tiebreak is needed:** all-gap columns are dropped first, so every remaining
column holds a residue. Two distinct columns with the same gap pattern must share a
residue row, hence are pinned, hence one is an ancestor of the other and they are never
both available at once. `available_columns_never_tie` asserts this over random input.

Cost is O(widthÂ² Ã— num_seqs) in O(width) extra memory â€” the same as the bubble sort it
replaces. Only the in-degree vector is stored; edges are recomputed on release rather
than held in an O(widthÂ²) adjacency list.

### All-gap columns are now dropped

Also the owner's call, and it settles the disputed Stage 2 behaviour. The original
*pinned* them, so a column holding no residue could still hold the rest of the alignment
apart â€” while its comment claimed the opposite. A column of pure gaps states no homology
relationship, so it is removed rather than merely made movable.

An alignment of nothing but gaps therefore standardises to zero width. That is left
representable rather than rejected (`Msa::select_columns` documents why); the distance
stage already reports an empty comparison as an error rather than the old `NaN`.

### `permute_columns` â†’ `select_columns`

Dropping columns means the width changes, which the old primitive could not express. The
new one takes any duplicate-free, in-range list of column indices and rebuilds the rows
from it, so reordering and dropping are one operation. Repeats are still rejected â€” they
would duplicate a residue â€” but omission is now legitimate and unchecked, because
`residue_hash` catches a dropped column that held a residue.

### Numbers that moved

```
                                    before Stage 6        after
test.fasta standardised             AA-- / A--A / AAA-    AA-- / A-A- / AA-A
test2.fasta standardised            A--A / A-A- / AA-A    A--A / AA-- / A-AA
d(test, test2) standardised         0.5                   0.5
d(test, test2) --no-standardise     0.21428571428571427   0.21428571428571427
```

`test.fasta` used to be a fixed point and no longer is. The pair distance is unchanged at
`0.5`, but **by coincidence, not by construction** â€” the hand-worked homology table behind
it is entirely different, and has been redone in
`end_to_end_with_standardisation_pins_the_new_fixture_distance`. Both were recomputed by
hand before running anything, and the program agreed on the first run.

`--no-standardise` is untouched, as it must be: Stage 6 changed only the standardise stage.

The two fixtures are *not* a legal permutation of one another â€” `test2` is `test` with two
**pinned** columns exchanged, which is why nothing requires them to converge, and why
`0.5` is a legitimate answer rather than a residual bug. This was not obvious and cost
some time; it is the reason the new end-to-end test needed a new fixture.

### Tests

New fixture `test/test_legal_permutation.fasta`: `test/test.fasta` with column 3 moved
left past columns 2 and 1, both legal moves. It is the first fixture that is genuinely
the same alignment as another one, which is what the headline property needs.

- `standardisation_is_confluent_over_legal_permutations` â€” 4000 random alignments, each
  legally shuffled, standardised forms must be byte-identical. Asserts that at least 1000
  of the shuffles actually moved a column, so the property cannot pass vacuously.
- `a_legally_permuted_alignment_is_at_distance_zero` â€” the same claim end to end through
  `run` and the CSV, both argument orders, plus the `--no-standardise` value
  (`0.35714285714285715` = 10/28, worked by hand) to show the test is not trivially zero.
- `the_original_tools_counterexample_now_converges` â€” the shrunk 3Ã—3 case above.
- `available_columns_never_tie`, `standardisation_never_alters_residues_on_random_input`,
  `standardisation_is_idempotent_on_random_input` â€” 4000/4000/2000 random alignments.
- `all_gap_columns_are_dropped`, `an_all_gap_alignment_standardises_to_zero_width`,
  `residues_fill_the_higher_rows_leftward`, `canonical_columns_omits_all_gap_columns`,
  and four `select_columns_*` tests.

Renamed and re-derived: `gap_only_column_does_not_move_disputed_behaviour` â†’
`all_gap_columns_are_dropped`; `two_movable_columns_swap` â†’
`residues_fill_the_higher_rows_leftward`; `already_canonical_alignment_is_a_no_op` â†’
`the_test_fixture_standardises_to_its_canonical_form`;
`standardise_reproduces_the_original_tool_on_test2_fixture` â†’
`the_test2_fixture_standardises_to_its_canonical_form`. Each keeps the old expected value
in a comment, so what changed is visible rather than merely gone.

The randomised tests use a seeded xorshift rather than `rand`, so failures are
reproducible and no dependency was added. Stage 6 added no dependencies.

### Still outstanding

- **`CODE_REVIEW.md` Â§2, the per-residue `HashSet` clone** in `homology.rs`, still marked
  with its `TODO`. Unchanged by this stage and still the hard memory ceiling
  (O(num_seqsÂ² Ã— width) live per MSA). Fixing it makes `|A|` and `|B|` genuinely differ,
  which is when Stage 4's `|A| + |B|` denominator starts mattering;
  `every_homology_set_has_size_num_seqs_minus_one` is the tripwire.
- **The naming question** from the plan's Loose Ends: `rusty-metAL` in `Cargo.toml`,
  `dlejeune/pipeline-utils-rs` in the justfile, `rusty-metal` in the log line and the
  README. Untouched â€” it needs a decision, not a guess.
- **No Dockerfile**, though the justfile has `build-docker` recipes referencing one.
- The comparator is still O(widthÂ² Ã— num_seqs). Unchanged in kind by this stage, and
  still deferred.
- The merge is now ready to commit. `INTEGRATION_PLAN.md` Stage 6 also listed nothing
  further; the plan's "Verification before merging" step is satisfied by
  `end_to_end_without_standardisation_reproduces_the_baseline` (pre-merge number intact)
  and the Stage 6 table above (the deltas, explained).


### Postscript: `just check` does not include `cargo fmt`

`cargo clippy -- -D warnings` and `cargo test` are both clean as of `0.2.0` and are what
`just check` runs. `cargo fmt --check` is **not** in it: the tree predates rustfmt and has
never been run through it, so formatting it now would reformat far more than this merge
and bury the Stage 6 diff. There is a separate `just fmt` recipe; doing it deserves its
own commit.

---

## Stage 7 â€” the memory ceiling (`CODE_REVIEW.md` Â§2) â€” DONE

The last outstanding review finding. `HomologyView` changed from a type alias
(`Vec<Vec<Option<HashSet<Element>>>>`, one materialised set per residue) to a struct
holding **one set per column**, shared, plus an index from `[sequence][position]` to a
column. 93 unit tests, all passing; no build, clippy or rustdoc warnings.

Live memory per alignment: **O(num_seqsÂ² Ã— width) â†’ O(num_seqs Ã— width)** â€” the size of
the alignment itself. The review's example, a 500 Ã— 5000 alignment, went from roughly
1.2Ã—10â¹ set entries to 2.5Ã—10â¶.

### It is exact, and that is the whole argument

The stored column includes the residue itself; the homology set is the column minus that
one element. The two are reconciled differently in the two halves of the fraction, and
only one of them needs any reconciling at all:

- **Numerator needs none.** The excluded element is `x = {sequence, position, gap:
  false}`. A residue's identity does not depend on the gaps around it, so `x` is the
  *same value* in both alignments. It is therefore in both columns, and so in neither
  symmetric difference: `(Ca \ {x}) Î” (Cb \ {x}) == Ca Î” Cb`.
- **Denominator subtracts one per side**, `|A| = |Ca| - 1`.

That identity is what makes this a pure memory change rather than a numeric one. Every
pinned value in the suite â€” `0.21428571428571427`, `0.5`, `0`,
`0.35714285714285715` â€” was unchanged, and no test expectation was edited.

`sharing_the_column_sets_computes_exactly_the_materialised_distance` checks the identity
against a literal implementation that materialises every homology set, **exhaustively
over every 2Ã—3 and 3Ã—2 alignment over {A, -} compared against every other** â€” 8192
comparisons. Exhaustive rather than sampled because the alignment space at that size is
small enough to enumerate, which is a stronger claim than any number of random draws.

### Measured, not just reasoned about

Two random 400-sequence Ã— 2000-column alignments, release build:

```
--no-standardise   0.86 s
with standardisation   0.79 s
```

Under the old shape the same input needed roughly 3.6 GB of homology sets and would have
thrashed or died. Two incidental findings from that run:

- **Standardisation is nearly free on dense alignments, and a no-op on them.** With 400
  sequences at ~50% residue density, almost every pair of columns shares a residue row,
  so almost everything is pinned and nothing can move. The distance came back identical
  with and without standardisation. `columns_are_pinned` short-circuits on the first
  shared row, so the O(widthÂ²) pair scan is cheap in exactly the case that has the most
  pairs.
- **Sparse alignments are fast for the opposite reason.** 60 sequences Ã— 2000 columns at
  97% gaps ran in 0.16 s: most columns are all-gap, and dropping them shrinks the width
  before the quadratic work starts.

So the deferred O(widthÂ² Ã— num_seqs) comparator cost needs a specific shape to bite â€”
many columns, few sequences, and a gap density high enough that columns are mutually
free but low enough that they are not all-gap. Nothing pathological showed up on
realistic inputs. Still deferred, but now with a measurement behind the deferral rather
than an assumption.

### API notes

- `HomologyView::column_of(seq, pos)` is what the distance path uses. Its doc comment
  carries the cancellation argument, because the returned set is *not* the homology set
  and a caller assuming otherwise would be wrong by exactly one element.
- `HomologyView::homology_set_of(seq, pos)` materialises the exclusion and returns the
  homology set as defined. It is `#[cfg(test)]` on purpose: if the distance path ever
  calls it, Â§2 comes straight back. The existing hand-worked tests were repointed at it
  rather than weakened, so they still assert homology sets as the metric defines them.
- `every_homology_set_has_size_num_seqs_minus_one` still holds and still pins the
  invariant that keeps `2Â·|A|` and `|A| + |B|` numerically equal. Stage 4 predicted this
  fix would break that invariant and make the two denominators genuinely differ. **It
  does not** â€” sharing the column set does not change any set's cardinality. The
  prediction was about a different candidate fix (moving to residues-only or gap-filtered
  sets), which remains unimplemented.

### `CODE_REVIEW.md` is now fully discharged

Every finding in the review is addressed: Â§0 (both), Â§1 (all three panics), Â§2, Â§3
(both), Â§4 (both), Â§5 (both), and the build warnings. Nothing in that document is
outstanding.

Remaining known work, none of it from the review:

- The O(widthÂ² Ã— num_seqs) ordering cost, deferred with the measurements above.
- No Dockerfile, though the justfile has `build-docker` recipes referencing one.

