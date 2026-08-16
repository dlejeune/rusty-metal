# rusty-metAL

Pairwise distances between multiple sequence alignments, based on how much the
alignments disagree about which residues are homologous to which.

The crate and the binary are `rusty-metal`, all lower case; `rusty-metAL` is the
stylised form used in prose and help output.

As of `0.2.0` this also subsumes the standalone `standardise-msa` tool: column
standardisation is now the first stage of the pipeline rather than a separate program.

## What it computes

For every residue in an alignment, its **homology set** is the set of elements sitting
in the same column, excluding itself. An element is identified by `(sequence, position,
is-gap)`, where a gap's position is the index of the residue preceding it in its row.

The distance between two alignments is

```
sum over residues r of  |A(r) Δ B(r)|
--------------------------------------
sum over residues r of  |A(r)| + |B(r)|
```

— the total symmetric difference over the total set size, so `0` means the two
alignments agree everywhere and `1` means they agree nowhere.

Sequences are matched **by name**, not by their order in the file. Both alignments must
carry the same set of names; a pair that does not is reported as an error and skipped,
and the rest of the run continues. Residue characters are compared case-insensitively.
Both `-` and `.` count as gaps.

## Standardisation

Gap *placement* is largely an artefact of whichever aligner produced the file: the same
alignment can be written with its columns in many different orders, and those orders
change the gap identities and therefore the distance. Standardisation rewrites each
alignment into a canonical column layout so that this artefact does not reach the metric.

Two columns are **pinned** when some row holds a residue in both of them — swapping
those would reverse two residues within a sequence. Every other pair is free to move.
Standardisation:

1. sorts the sequences by name;
2. drops every column that holds no residue at all;
3. orders the rest so that residues fill the higher rows as far to the left as the
   pinning constraints allow.

The result is a genuine canonical form: **two files holding the same alignment, differing
only in where the gaps sit, standardise to byte-identical output and therefore compare at
distance 0.** Standardisation never changes residue content, and the result is verified
against a residue hash on every run.

Note that standardisation moves each alignment to *its own* canonical layout. It does not
move two genuinely different alignments toward each other, so a distance between two
different alignments can go either up or down relative to not standardising.

## Usage

```
rusty-metal [OPTIONS] --output-fp <FILE> <INPUT_FILES>...
```

Distances for every pair of inputs, standardising internally:

```
rusty-metal -o dist.csv a.fasta b.fasta c.fasta
```

Distances, and keep the standardised alignments too:

```
rusty-metal -o dist.csv --emit-standardised out/ a.fasta b.fasta
```

Standardise only, no distances — this is the replacement for `standardise-msa`:

```
rusty-metal --standardise-only --emit-standardised out/ a.fasta
```

Distances on the files exactly as they sit on disk:

```
rusty-metal -o dist.csv --no-standardise a.fasta b.fasta
```

### Options

| Option | Meaning |
| --- | --- |
| `-o`, `--output-fp <FILE>` | Where to write the distance CSV. Required unless `--standardise-only`. |
| `-n`, `--num-threads <N>` | Worker threads. `0` (the default) lets rayon choose. |
| `--emit-standardised <DIR>` | Write each standardised alignment to `<DIR>/<stem>.standardised.fasta`. Created if absent. |
| `--standardise-only` | Standardise and stop. Requires `--emit-standardised`. |
| `--no-standardise` | Skip standardisation and compare the files as they are. |

Illegal combinations (`--standardise-only` without `--emit-standardised`,
`--emit-standardised` with `--no-standardise`, and so on) are rejected before any file is
opened, as are two inputs whose filenames would produce the same output file.

`--no-standardise` is an escape hatch for isolating the effect of standardisation on a
real dataset. It is **not** a bug-compatibility switch: name-keyed matching, the
`|A|+|B|` denominator, `.` as a gap, and the ragged/empty-input errors all apply
regardless.

### Output

```csv
msa_a,msa_b,distance
a.fasta,b.fasta,0.21428571428571427
```

Fields are quoted per RFC 4180 when they contain a comma, quote or newline. A pair that
could not be compared gets an empty distance field and a logged warning.

## Errors

A **file-level** problem is fatal: a ragged or empty FASTA, or an unreadable path,
invalidates every pair it appears in, so the run stops before writing a CSV rather than
emitting a result set silently missing an arbitrary subset of comparisons.

A **pair-level** problem is not: mismatched sequence names, or a comparison with nothing
in it, is logged and the run continues.

## Performance

Every pair is compared independently, so the run parallelises across pairs (`-n`). Peak
memory, release build:

| Workload | Time | Peak memory |
| --- | --- | --- |
| 1 pair, 500 sequences × 5000 columns | 0.7 s | 49 MB |
| 6 pairs, 300 × 3000, `-n 1` | | 22 MB |
| 6 pairs, 300 × 3000, `-n 6` | | 92 MB |
| 1 pair, 60 × 2000 at 97% gaps | 0.16 s | |

The homology view of an alignment is about 4× the size of the FASTA that produced it, and
two views exist per comparison in flight — so the memory a run needs is set by the largest
pair times the thread count, not by the total input size. Standardisation is not a
significant cost at these sizes: the 500 × 5000 pair takes the same 49 MB with it on as
with `--no-standardise`.

Two things do scale in ways worth knowing about:

- **Every input alignment is held in memory for the whole run.** Files are read and
  standardised once, up front, rather than re-read for each pair they appear in. That is
  a large win in time — N reads rather than N-1 per file — but it means a run over a
  thousand 500 × 5000 alignments holds roughly 2.5 GB before the first comparison starts.
- **Choosing the canonical column order costs O(width² × num_seqs).** In practice it
  disappears at both extremes: dense alignments pin almost every column pair and the scan
  short-circuits immediately, while sparse ones lose most of their columns to the all-gap
  drop before the quadratic work begins. The shape that would bite is many columns, few
  sequences, and a gap density high enough that columns are mutually free but low enough
  that they are not all-gap.

## Development

```
just build      # cargo build
just test       # cargo test
just doc        # rustdoc, with warnings fatal
just doc-open   # ...and open it in a browser
just check      # fmt, clippy, docs and tests — what CI runs
just run -- -o dist.csv test/test.fasta test/test2.fasta
```

The internals are documented in some depth — each module opens with a header explaining
what it does and why the shape it has was chosen over the alternatives. CI publishes that
to **<https://dlejeune.github.io/rusty-metal/>** on every push to `main`, or run
`just doc-open` for the same thing locally.

`just check` is the gate. It includes `cargo doc` because a doc comment can rot in ways
nothing else catches: a link to an item that was later renamed still compiles and still
passes the tests.

### Docker

```
just build-docker              # dlejeune/rusty-metal:0.2.0
just run-docker-it             # interactive shell, cwd mounted at /data
just docker                    # build and push
```

The image has no entrypoint so that the interactive recipe can start a shell, which means
the binary is named explicitly when running it directly:

```
docker run --rm -v ./:/data dlejeune/rusty-metal:0.2.0 \
    rusty-metal -o /data/dist.csv /data/a.fasta /data/b.fasta
```

## Repository map

Alongside the source, three documents record how the tree got here. They are historical
records, not user documentation, and only this README describes current behaviour:

| File | What it is |
| --- | --- |
| `INTEGRATION_PLAN.md` | The plan for folding `standardise-msa` in, written before the work. Stage numbering throughout the other two refers to it. |
| `INTEGRATION_NOTES.md` | The running log of that work, one entry per stage, recording what was decided and what turned out to be wrong. Stages 6 and 8 are the substantial ones. |
| `CODE_REVIEW.md` | A review of the pre-merge tree at `6f73a81`. All 12 findings are discharged; the status header at the top maps each to where it was fixed. The body below it is the original review, and its file and line references are pre-merge and mostly dangle. |

## Compatibility

`0.2.0` is **not** output-compatible with `standardise-msa`, and distances computed with
standardisation differ from earlier builds of this crate. The ordering rule was replaced
because the original was not confluent — the same alignment written two different ways
standardised to two different layouts, which defeated the purpose of standardising at
all. `INTEGRATION_NOTES.md` (Stage 6) records the investigation and the numbers.
