```
                               __                                             __      ______  __
                             |  \                                           |  \    /      \|  \
  ______  __    __  _______ _| ▓▓_   __    __        ______ ____   ______  _| ▓▓_  |  ▓▓▓▓▓▓\ ▓▓
 /      \|  \  |  \/       \   ▓▓ \ |  \  |  \______|      \    \ /      \|   ▓▓ \ | ▓▓__| ▓▓ ▓▓
|  ▓▓▓▓▓▓\ ▓▓  | ▓▓  ▓▓▓▓▓▓▓\▓▓▓▓▓▓ | ▓▓  | ▓▓      \ ▓▓▓▓▓▓\▓▓▓▓\  ▓▓▓▓▓▓\\▓▓▓▓▓▓ | ▓▓    ▓▓ ▓▓
| ▓▓   \▓▓ ▓▓  | ▓▓\▓▓    \  | ▓▓ __| ▓▓  | ▓▓\▓▓▓▓▓▓ ▓▓ | ▓▓ | ▓▓ ▓▓    ▓▓ | ▓▓ __| ▓▓▓▓▓▓▓▓ ▓▓
| ▓▓     | ▓▓__/ ▓▓_\▓▓▓▓▓▓\ | ▓▓|  \ ▓▓__/ ▓▓      | ▓▓ | ▓▓ | ▓▓ ▓▓▓▓▓▓▓▓ | ▓▓|  \ ▓▓  | ▓▓ ▓▓_____
| ▓▓      \▓▓    ▓▓       ▓▓  \▓▓  ▓▓\▓▓    ▓▓      | ▓▓ | ▓▓ | ▓▓\▓▓     \  \▓▓  ▓▓ ▓▓  | ▓▓ ▓▓     \
 \▓▓       \▓▓▓▓▓▓ \▓▓▓▓▓▓▓    \▓▓▓▓ _\▓▓▓▓▓▓▓       \▓▓  \▓▓  \▓▓ \▓▓▓▓▓▓▓   \▓▓▓▓ \▓▓   \▓▓\▓▓▓▓▓▓▓▓
                                    |  \__| ▓▓
                                     \▓▓    ▓▓
                                      \▓▓▓▓▓▓
```

Pairwise distances between multiple sequence alignments, based on how much the
alignments disagree about which residues are homologous to which.

The metric is the one described in

> Blackburne, B. P. and Whelan, S. (2012). Measuring the distance between multiple
> sequence alignments. _Bioinformatics_ **28**(4), 495–502.
> <https://doi.org/10.1093/bioinformatics/btr701>

and the name follows the authors' own `metAL`. The crate and the binary are
`rusty-metal`, all lower case; `rusty-metAL` is the stylised form used in prose and help
output.

Column standardisation runs as the first stage of the pipeline, so the standalone
`standardise-msa` tool is no longer needed; `--standardise-only` does what it did.

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

Gap _placement_ is largely an artefact of whichever aligner produced the file: the same
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

The result is a canonical form: **two files holding the same alignment, differing only in
where the gaps sit, standardise to byte-identical output and therefore compare at distance 0.** Residue content is unchanged, which is verified against a residue hash on every run.

Note that standardisation moves each alignment to _its own_ canonical layout. It does not
move two genuinely different alignments toward each other, so a distance between two
different alignments can go either up or down relative to not standardising.

## Installation

Binaries and installer scripts can be found in the [Releases](https://github.com/dlejeune/rusty-metal/releases) tab on GitHub.
Alternatively, rusty-metAL can be built from scratch (see [Development](#development))

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

Standardise only, without computing distances:

```
rusty-metal --standardise-only --emit-standardised out/ a.fasta
```

Distances on the files exactly as they sit on disk:

```
rusty-metal -o dist.csv --no-standardise a.fasta b.fasta
```

### Options

| Option                      | Meaning                                                                                    |
| --------------------------- | ------------------------------------------------------------------------------------------ |
| `-o`, `--output-fp <FILE>`  | Where to write the distance CSV. Required unless `--standardise-only`.                     |
| `-n`, `--num-threads <N>`   | Worker threads. `0` (the default) lets rayon choose.                                       |
| `--emit-standardised <DIR>` | Write each standardised alignment to `<DIR>/<stem>.standardised.fasta`. Created if absent. |
| `--standardise-only`        | Standardise and stop. Requires `--emit-standardised`.                                      |
| `--no-standardise`          | Skip standardisation and compare the files as they are.                                    |

Illegal combinations (`--standardise-only` without `--emit-standardised`,
`--emit-standardised` with `--no-standardise`, and so on) are rejected before any file is
opened, as are two inputs whose filenames would produce the same output file.

`--no-standardise` disables the standardisation stage and nothing else, so it measures
what standardisation contributes to a distance on a given dataset. Every other rule
described here — name-keyed matching, `.` as a gap, the ragged and empty input errors —
applies either way.

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

Every pair is compared independently, so a run parallelises across pairs (`-n`). Measured
on a release build:

| Workload                             | Time   | Peak memory |
| ------------------------------------ | ------ | ----------- |
| 1 pair, 500 sequences × 5000 columns | 0.7 s  | 49 MB       |
| 1 pair, 60 × 2000 at 97% gaps        | 0.16 s | —           |
| 6 pairs, 300 × 3000, `-n 1`          | —      | 22 MB       |
| 6 pairs, 300 × 3000, `-n 6`          | —      | 92 MB       |

Peak memory scales with the largest pair in flight times the thread count. Standardisation
adds nothing measurable at these sizes: the 500 × 5000 pair peaks at the same 49 MB with
it on as under `--no-standardise`.

Each input is read and standardised once, up front, and held for the rest of the run. A
run over a thousand 500 × 5000 alignments therefore holds roughly 2.5 GB before the first
comparison starts.

## Development

```
just build      # cargo build
just test       # cargo test
just doc        # rustdoc, with warnings fatal
just doc-open   # ...and open it in a browser
just check      # fmt, clippy, docs and tests — what CI runs
just run -- -o dist.csv test/test.fasta test/test2.fasta
```

Each module opens with a header describing what it holds. CI publishes the rustdoc to
**<https://dlejeune.github.io/rusty-metal/>** on every push to `main`; `just doc-open`
builds the same pages locally.

`just check` is the gate CI runs. It includes `cargo doc` because a doc comment can rot in
ways nothing else catches: a link to an item that was later renamed still compiles and
still passes the tests.

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

## Compatibility

Standardised output is not byte-compatible with `standardise-msa`, and distances computed
with standardisation differ from those of `0.1.x`. The column ordering rule differs from
the one `standardise-msa` used, which did not produce a canonical form: an alignment
written two different ways could standardise to two different layouts, so a pair of files
holding the same alignment did not reliably compare at distance 0.

## License

This work is licensed under the MIT public license and can be found [here](LICENSE)
