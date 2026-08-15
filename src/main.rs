//! `rusty-metal` — pairwise symmetric-difference distances between multiple sequence
//! alignments, with standardisation as the first stage of the pipeline.
//!
//! # The two phases
//!
//! Stage 5 splits the run in two, and the split is load-bearing rather than tidy:
//!
//! - **Phase 1** reads and standardises every input file *once*, in parallel across
//!   files, and emits the standardised alignments here if asked.
//! - **Phase 2** computes the pairwise distances over the alignments phase 1 left in
//!   memory.
//!
//! Doing the standardisation inside the per-pair comparison instead would standardise
//! each file once per pair it appears in (O(n) redundant work per file), and — worse —
//! would have several rayon workers writing the *same* `--emit-standardised` path
//! concurrently for any file appearing in more than one pair. Phase 1 has exactly one
//! writer per output path, so that race cannot arise. It also removes the O(n²) re-reads
//! the previous shape performed: `compare_alignment_pair` used to open both files afresh
//! for every pair.
//!
//! Raw [`Msa`]s are roughly file-sized, so holding all N at once is affordable. The
//! memory-hungry structures are the homology views, and those stay per-pair.
//!
//! # Error handling policy
//!
//! The two phases fail differently, on purpose:
//!
//! - A **phase 1** failure — a file that will not read, or one whose standardisation
//!   residue-hash check fails — is a hard error that aborts the run. That file
//!   participates in every pair it appears in, so continuing would produce a result set
//!   silently missing an arbitrary subset of comparisons.
//! - A **phase 2** failure — one pair that cannot be compared, e.g. two alignments over
//!   different sequence sets — is logged and the run continues, with an empty distance
//!   field in that pair's CSV row. This is the pre-existing behaviour and is preserved.
//!
//! Note that a panic inside the rayon bridge aborts the whole process and bypasses the
//! per-pair handler entirely, so this file avoids `unwrap`/`expect` outside its tests.

use crate::distance::compute_symmetric_difference;
use crate::homology::{homology_view, Registry};
use crate::msa::{read_msa, write_msa, Msa};
use crate::standardise::standardise;
use anyhow::{bail, Context, Result};
use clap::builder::styling;
use clap::Parser;
use colored::Colorize;
use itertools::Itertools;
use rayon::prelude::*;
use std::collections::BTreeMap;
use std::ffi::OsString;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
mod distance;
mod homology;
mod msa;
mod standardise;

const VERSION: &str = env!("CARGO_PKG_VERSION");

/// The stylised form of the name, for output a human reads.
///
/// The binary itself is `rusty-metal`, all lower case — that is what gets typed, what
/// clap prints in its usage line, and what the Docker image is tagged with. This form
/// is for prose: the startup banner and the help text. Keep the two distinct rather
/// than making anything case-insensitive.
const DISPLAY_NAME: &str = "rusty-metAL";
const STYLES: styling::Styles = styling::Styles::styled()
    .header(styling::AnsiColor::Green.on_default().bold())
    .usage(styling::AnsiColor::Green.on_default().bold())
    .literal(styling::AnsiColor::Blue.on_default().bold())
    .placeholder(styling::AnsiColor::Cyan.on_default())
    .error(styling::AnsiColor::Red.on_default().bold());

/// The suffix appended to an input's file stem to name its standardised output.
const STANDARDISED_SUFFIX: &str = ".standardised.fasta";

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
#[command(styles = STYLES)]
pub struct Args {
    /// MSAs to process. Two or more are required when distances are computed; a single
    /// file is enough under --standardise-only.
    //
    // `num_args = 1..` rather than `2..`: the "at least two" rule belongs to the
    // distance mode, not to the parser, and is enforced by `RunPlan::from_args` so that
    // --standardise-only can accept one file.
    #[clap(required = true, num_args = 1..)]
    input_files: Vec<PathBuf>,

    /// Where to write the pairwise distance CSV. Required unless --standardise-only.
    //
    // Optional at the parser level for exactly that reason; `RunPlan::from_args` makes
    // it mandatory in the modes that produce distances.
    #[clap(short, long, value_name = "FILE")]
    output_fp: Option<PathBuf>,

    /// Worker threads to use. 0 lets rayon choose based on available parallelism.
    #[clap(short, long, default_value_t = 0)]
    num_threads: usize,

    /// Write each standardised alignment to this directory as
    /// `<input stem>.standardised.fasta`. The directory is created if it does not exist.
    #[clap(long, value_name = "DIR")]
    emit_standardised: Option<PathBuf>,

    /// Standardise the inputs and stop, without computing any distances. Requires
    /// --emit-standardised, since otherwise the run would produce nothing observable.
    #[clap(long)]
    standardise_only: bool,

    /// Skip the standardisation stage and compute distances on the alignments exactly as
    /// they appear on disk.
    ///
    /// This is an escape hatch for isolating the effect of standardisation on a real
    /// dataset. It is NOT a bug-compatibility switch and does not restore any pre-merge
    /// behaviour: name-keyed sequence matching, the |A|+|B| denominator, the widened gap
    /// definition (`.` counts as a gap) and the ragged/empty-input errors all apply
    /// regardless of this flag.
    #[clap(long)]
    no_standardise: bool,
}

// ---------------------------------------------------------------------------------
// The validated run plan
// ---------------------------------------------------------------------------------

/// What the standardisation stage should do, in the modes that compute distances.
///
/// The three states are separate variants rather than a `bool` plus an
/// `Option<PathBuf>` so that "emit the standardised alignments but do not standardise"
/// cannot be built: [`Standardisation::Skip`] carries no directory to emit into.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Standardisation {
    /// Do not standardise. `--no-standardise`.
    Skip,
    /// Standardise in memory; do not write the standardised alignments anywhere.
    InMemory,
    /// Standardise and write each result into this directory.
    Emit(PathBuf),
}

/// What the run produces.
///
/// [`Mode::StandardiseOnly`] holds a non-optional directory, so `--standardise-only`
/// without `--emit-standardised` is unrepresentable, and the variants are disjoint so
/// `--standardise-only --no-standardise` is unrepresentable too. Each illegal
/// combination is rejected once, in [`RunPlan::from_args`], rather than being re-checked
/// wherever it might matter.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Mode {
    /// Standardise every input and write it out. No distances, so no output CSV.
    StandardiseOnly { emit_dir: PathBuf },
    /// Compute pairwise distances into `output_fp`, standardising first per
    /// `standardisation`.
    Distances {
        output_fp: PathBuf,
        standardisation: Standardisation,
    },
}

/// A validated description of one run: the illegal flag combinations have already been
/// rejected, so nothing downstream re-checks them.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RunPlan {
    pub inputs: Vec<PathBuf>,
    pub num_threads: usize,
    pub mode: Mode,
}

impl RunPlan {
    /// Validates parsed [`Args`] into a plan, or explains what is wrong with them.
    ///
    /// This is deliberately a pure `&Args -> Result<RunPlan>` function with no I/O
    /// beyond inspecting the paths themselves, so the whole of the CLI's validation can
    /// be unit-tested without spawning the binary.
    pub fn from_args(args: &Args) -> Result<RunPlan> {
        if args.standardise_only && args.no_standardise {
            bail!(
                "--standardise-only and --no-standardise contradict each other: the first says \
                 standardisation is the only thing to do, the second says not to do it. Drop \
                 one of them."
            );
        }

        if args.emit_standardised.is_some() && args.no_standardise {
            bail!(
                "--emit-standardised has nothing to emit under --no-standardise, which skips \
                 the standardisation stage entirely. Drop --no-standardise to emit standardised \
                 alignments, or drop --emit-standardised to compute distances on the inputs as \
                 they are."
            );
        }

        let mode = if args.standardise_only {
            let emit_dir = match &args.emit_standardised {
                Some(dir) => dir.clone(),
                None => bail!(
                    "--standardise-only needs --emit-standardised <DIR>: with no distances to \
                     compute and nowhere to write the standardised alignments, the run would do \
                     nothing observable."
                ),
            };
            Mode::StandardiseOnly { emit_dir }
        } else {
            let output_fp = match &args.output_fp {
                Some(path) => path.clone(),
                None => bail!(
                    "-o/--output-fp <FILE> is required when distances are computed. Pass one, or \
                     pass --standardise-only --emit-standardised <DIR> to standardise without \
                     computing distances."
                ),
            };

            if args.input_files.len() < 2 {
                bail!(
                    "Computing distances needs at least 2 input alignments to form a pair, but \
                     {} was given. Pass another alignment, or use --standardise-only \
                     --emit-standardised <DIR> to standardise a single file.",
                    args.input_files.len()
                );
            }

            let standardisation = match (&args.emit_standardised, args.no_standardise) {
                (Some(dir), false) => Standardisation::Emit(dir.clone()),
                (None, false) => Standardisation::InMemory,
                // `(Some(_), true)` was rejected above.
                (_, true) => Standardisation::Skip,
            };

            Mode::Distances {
                output_fp,
                standardisation,
            }
        };

        let plan = RunPlan {
            inputs: args.input_files.clone(),
            num_threads: args.num_threads,
            mode,
        };

        // Checked here, before any file is opened or created, because the failure mode
        // is silent data loss: two inputs sharing a stem map to one output path and the
        // second write destroys the first.
        if plan.emit_dir().is_some() {
            check_stem_collisions(&plan.inputs)?;
        }

        Ok(plan)
    }

    /// The directory standardised alignments are written to, if any.
    pub fn emit_dir(&self) -> Option<&Path> {
        match &self.mode {
            Mode::StandardiseOnly { emit_dir } => Some(emit_dir.as_path()),
            Mode::Distances {
                standardisation: Standardisation::Emit(dir),
                ..
            } => Some(dir.as_path()),
            Mode::Distances { .. } => None,
        }
    }

    /// Whether phase 1 runs the standardisation pass.
    pub fn standardises(&self) -> bool {
        !matches!(
            self.mode,
            Mode::Distances {
                standardisation: Standardisation::Skip,
                ..
            }
        )
    }

    /// The CSV path, if this run computes distances at all.
    pub fn output_fp(&self) -> Option<&Path> {
        match &self.mode {
            Mode::Distances { output_fp, .. } => Some(output_fp.as_path()),
            Mode::StandardiseOnly { .. } => None,
        }
    }
}

/// The path a standardised alignment is written to: `<input stem>.standardised.fasta`
/// inside `dir`.
///
/// Built by pushing onto the stem's `OsString` rather than going through `&str`, so a
/// non-UTF-8 input path is named correctly instead of being rejected or mangled.
fn standardised_output_path(dir: &Path, input: &Path) -> Result<PathBuf> {
    let stem = input.file_stem().with_context(|| {
        format!(
            "Cannot name the standardised output for {}: the path has no file stem",
            input.display()
        )
    })?;
    let mut name = stem.to_os_string();
    name.push(STANDARDISED_SUFFIX);
    Ok(dir.join(name))
}

/// Rejects any two inputs whose file stems agree, because they would be written to the
/// same `<stem>.standardised.fasta` and the second would silently overwrite the first.
///
/// Different directories are no defence: `a/aln.fasta` and `b/aln.fasta` both emit
/// `aln.standardised.fasta`. The error lists every colliding input so the caller can see
/// which files to rename or split apart.
fn check_stem_collisions(inputs: &[PathBuf]) -> Result<()> {
    // `BTreeMap` so that a multi-collision report comes out in a stable, sorted order
    // rather than depending on hash iteration.
    let mut by_stem: BTreeMap<OsString, Vec<&PathBuf>> = BTreeMap::new();
    for input in inputs {
        let stem = input.file_stem().with_context(|| {
            format!(
                "Cannot name the standardised output for {}: the path has no file stem",
                input.display()
            )
        })?;
        by_stem.entry(stem.to_os_string()).or_default().push(input);
    }

    let collisions: Vec<(&OsString, &Vec<&PathBuf>)> = by_stem
        .iter()
        .filter(|(_, paths)| paths.len() > 1)
        .collect();

    if collisions.is_empty() {
        return Ok(());
    }

    let mut message = String::from(
        "Two or more inputs share a file stem, so their standardised outputs would collide \
         and overwrite each other:",
    );
    for (stem, paths) in collisions {
        message.push_str(&format!(
            "\n  {}{} would be written by: {}",
            Path::new(stem).display(),
            STANDARDISED_SUFFIX,
            paths
                .iter()
                .map(|p| p.display().to_string())
                .collect::<Vec<String>>()
                .join(", ")
        ));
    }
    message.push_str(
        "\nRename the inputs, or emit them into separate directories in separate runs.",
    );

    bail!(message)
}

// ---------------------------------------------------------------------------------
// Phase 1 — read, standardise and (optionally) emit, once per file
// ---------------------------------------------------------------------------------

/// Reads every input once, standardises it unless the plan says not to, and writes the
/// result out if the plan asks for it.
///
/// Runs in parallel across *files*, and the returned vector is in input order —
/// `par_iter().map(...).collect::<Result<Vec<_>>>()` preserves indexing, which phase 2
/// relies on to pair an alignment back up with its path.
///
/// Any failure here aborts the run: see the module docs for why a bad file is not
/// survivable the way a bad pair is.
fn load_inputs(plan: &RunPlan) -> Result<Vec<Msa>> {
    if let Some(dir) = plan.emit_dir() {
        std::fs::create_dir_all(dir).with_context(|| {
            format!(
                "Failed to create the --emit-standardised directory {}",
                dir.display()
            )
        })?;
    }

    plan.inputs
        .par_iter()
        .map(|path| {
            let mut msa = read_msa(path)?;

            if plan.standardises() {
                // `standardise` performs its own residue-hash check and returns `Err`
                // if the column permutation moved any residue, so there is nothing
                // extra to verify here.
                standardise(&mut msa).with_context(|| {
                    format!("Failed to standardise the alignment in {}", path.display())
                })?;
            }

            if let Some(dir) = plan.emit_dir() {
                let out_path = standardised_output_path(dir, path)?;
                log::info!(
                    "Writing standardised alignment to {}",
                    out_path.display().to_string().cyan()
                );
                write_msa(&out_path, &msa).with_context(|| {
                    format!(
                        "Failed to write the standardised form of {} to {}",
                        path.display(),
                        out_path.display()
                    )
                })?;
            }

            Ok(msa)
        })
        .collect()
}

// ---------------------------------------------------------------------------------
// Phase 2 — pairwise distances over the in-memory alignments
// ---------------------------------------------------------------------------------

/// The symmetric-difference distance between two already-loaded alignments.
///
/// The registry is what ties the two together: it rejects a pair that is not over the
/// same set of sequence names, and it assigns the sequence ids that both homology views
/// are keyed on, so the comparison does not depend on the order the records appear in
/// either file.
///
/// Neither reading nor standardisation happens here — both are phase 1's job, done once
/// per file rather than once per pair.
fn compare_alignment_pair(msa_a: &Msa, msa_b: &Msa) -> Result<f64> {
    let registry = Registry::for_pair(msa_a, msa_b)?;
    let view_a = homology_view(msa_a, &registry)?;
    let view_b = homology_view(msa_b, &registry)?;
    compute_symmetric_difference(&view_a, &view_b, &registry)
}

struct DistanceResult {
    msa_a: PathBuf,
    msa_b: PathBuf,
    distance: Option<f64>,
}

/// A short label for a path, for log lines only.
///
/// The file stem when there is one, falling back to the whole path. Replaces
/// `file_stem().unwrap().to_str().unwrap()` (`CODE_REVIEW.md` §1), which panicked on a
/// stemless or non-UTF-8 path — on the *success* branch, destroying a run that had
/// otherwise computed fine, and inside a rayon worker, so it took every other completed
/// comparison down with it.
fn label(path: &Path) -> String {
    match path.file_stem() {
        Some(stem) => stem.to_string_lossy().into_owned(),
        None => path.display().to_string(),
    }
}

/// Computes every pairwise distance over `msas`, which is parallel to `paths`.
///
/// A pair that cannot be compared is logged and yields a `None` distance; the remaining
/// pairs still run and still reach the CSV.
fn compare_all_pairs(paths: &[PathBuf], msas: &[Msa]) -> Vec<DistanceResult> {
    (0..msas.len())
        .combinations(2)
        .enumerate()
        .par_bridge()
        .map(|(curr_idx, pair)| {
            let (a, b) = (pair[0], pair[1]);
            DistanceResult {
                msa_a: paths[a].clone(),
                msa_b: paths[b].clone(),
                distance: match compare_alignment_pair(&msas[a], &msas[b]) {
                    Ok(distance) => {
                        log::info!("Comparison {curr_idx} complete.");
                        log::info!(
                            "{}",
                            format!(
                                "{} -> {}: {}",
                                label(&paths[a]),
                                label(&paths[b]),
                                distance
                            )
                            .green()
                            .bold()
                        );
                        Some(distance)
                    }
                    Err(e) => {
                        log::warn!(
                            "An error occurred while computing the distance between {} and {}. Continuing with remaining comparisons.",
                            paths[a].display(),
                            paths[b].display()
                        );
                        log::warn!("{}", e);
                        None
                    }
                },
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------------
// CSV output
// ---------------------------------------------------------------------------------

/// Escapes one CSV field per RFC 4180: fields containing a comma, a double quote, a
/// carriage return or a line feed are wrapped in double quotes, and any embedded double
/// quote is doubled. Everything else is emitted verbatim.
///
/// Paths really do contain commas — `results/run,v2/aln.fasta` used to emit four fields
/// where three were expected, and every downstream parser then misattributed the
/// distance (`CODE_REVIEW.md` §5).
fn csv_escape(field: &str) -> String {
    if field.contains([',', '"', '\n', '\r']) {
        format!("\"{}\"", field.replace('"', "\"\""))
    } else {
        field.to_string()
    }
}

/// Writes the CSV body to `writer`.
///
/// Every write is `write_all` (a short `write` used to drop part of a row silently) and
/// every error is returned rather than `expect`ed — `main` has always returned `Result`,
/// so a panic here was throwing away an error path that already existed
/// (`CODE_REVIEW.md` §3).
fn write_results<W: Write>(writer: &mut W, results: &[DistanceResult]) -> Result<()> {
    writer.write_all(b"msa_a,msa_b,distance\n")?;
    for result in results {
        let row = format!(
            "{},{},{}\n",
            csv_escape(&result.msa_a.display().to_string()),
            csv_escape(&result.msa_b.display().to_string()),
            match result.distance {
                Some(distance) => distance.to_string(),
                None => String::new(),
            }
        );
        writer.write_all(row.as_bytes())?;
    }
    Ok(())
}

/// Writes the results to `path`, flushing explicitly.
///
/// The explicit flush is the point: a `BufWriter` dropped at the end of `main` flushes
/// implicitly and *discards the error*, so a failure on the final flush (disk full,
/// quota) produced a truncated CSV with exit code 0 (`CODE_REVIEW.md` §3).
fn write_csv(path: &Path, results: &[DistanceResult]) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("Failed to create the output file {}", path.display()))?;
    let mut writer = BufWriter::new(file);

    write_results(&mut writer, results)
        .with_context(|| format!("Failed to write results to {}", path.display()))?;

    writer
        .flush()
        .with_context(|| format!("Failed to flush {}", path.display()))?;

    Ok(())
}

// ---------------------------------------------------------------------------------
// Driving the run
// ---------------------------------------------------------------------------------

/// Runs a validated plan inside its own rayon thread pool.
///
/// **A scoped pool, not `build_global`.** `rayon::ThreadPoolBuilder::build_global` may be
/// called at most once per process and returns `Err` on every subsequent call, so the
/// previous shape — building the global pool inside the pipeline function — meant the
/// pipeline could only ever be run once, and any test that ran it twice failed on the
/// second call for reasons that had nothing to do with what it was testing. Tolerating
/// the already-initialised error would work but silently ignores `--num-threads` on all
/// but the first run, and a `OnceLock` has the same problem. A pool built per run has
/// neither issue: `install` makes it the pool that every nested `par_iter` and
/// `par_bridge` in this crate uses, including the one inside
/// `distance::compute_symmetric_difference`, and it is torn down when `run` returns.
///
/// `num_threads(0)` means "let rayon decide from the available parallelism", which is
/// the same default the `--num-threads 0` flag documents.
pub fn run(plan: &RunPlan) -> Result<()> {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(plan.num_threads)
        .build()
        .context("Could not set up multithreading")?;

    pool.install(|| run_in_pool(plan))
}

fn run_in_pool(plan: &RunPlan) -> Result<()> {
    // Phase 1: read, standardise and emit — once per file.
    let msas = load_inputs(plan)?;

    // Phase 2: pairwise distances over what phase 1 produced.
    let Some(output_fp) = plan.output_fp() else {
        log::info!(
            "{}",
            "Standardisation complete; no distances requested."
                .green()
                .bold()
        );
        return Ok(());
    };

    let results = compare_all_pairs(&plan.inputs, &msas);

    log::info!("Writing output to {}", output_fp.display());
    write_csv(output_fp, &results)
}

fn main() -> Result<()> {
    simple_logger::SimpleLogger::new().env().init()?;
    let args = Args::parse();
    let plan = RunPlan::from_args(&args)?;

    log::info!(
        "This is {} version {}",
        DISPLAY_NAME.bold(),
        VERSION.bold().cyan()
    );
    if args.standardise_only && args.output_fp.is_some() {
        // Not an error — the flag combination is coherent, just pointless — but silently
        // ignoring an explicitly requested output path would be worse than saying so.
        log::warn!(
            "--standardise-only computes no distances, so -o/--output-fp is ignored and no CSV \
             will be written."
        );
    }
    match &plan.mode {
        Mode::StandardiseOnly { emit_dir } => log::info!(
            "Standardising the following MSAs into {} using {} threads. No distances will be computed.",
            emit_dir.display().to_string().bright_blue(),
            plan.num_threads.to_string().bright_red()
        ),
        Mode::Distances { output_fp, .. } => log::info!(
            "Computing the distance between the following MSAs, and writing the output to {} and using {} threads. Standardisation is {}.",
            output_fp.display().to_string().bright_blue(),
            plan.num_threads.to_string().bright_red(),
            if plan.standardises() { "on" } else { "off" }
        ),
    }
    for file in &plan.inputs {
        log::info!("> {}", file.display().to_string().cyan());
    }

    run(&plan)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicU64, Ordering};

    /// A fresh, empty directory under the system temp dir. Cargo runs tests in parallel,
    /// so the name carries a counter as well as the process id.
    fn temp_dir(label: &str) -> PathBuf {
        static COUNTER: AtomicU64 = AtomicU64::new(0);
        let n = COUNTER.fetch_add(1, Ordering::Relaxed);
        let dir = std::env::temp_dir().join(format!(
            "rusty-metal-main-test-{}-{}-{}",
            std::process::id(),
            n,
            label
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("temp dir should be creatable");
        dir
    }

    /// Parses an argv the way the binary would, so these tests exercise the clap
    /// configuration (which arguments are required, how many values they take) as well
    /// as `RunPlan::from_args`.
    fn parse(argv: &[&str]) -> Result<RunPlan> {
        let args = Args::try_parse_from(std::iter::once("rusty-metal").chain(argv.iter().copied()))
            .map_err(|e| anyhow::anyhow!("{e}"))?;
        RunPlan::from_args(&args)
    }

    fn err_of(result: Result<RunPlan>) -> String {
        match result {
            Ok(plan) => panic!("expected an error, got the plan {plan:?}"),
            Err(e) => e.to_string(),
        }
    }

    fn read_to_string(path: &Path) -> String {
        std::fs::read_to_string(path)
            .unwrap_or_else(|e| panic!("{} should be readable: {e}", path.display()))
    }

    // -----------------------------------------------------------------------------
    // The four valid invocations
    // -----------------------------------------------------------------------------

    #[test]
    fn plan_accepts_distances_with_internal_standardisation() {
        // rusty-metal -o dist.csv a.fa b.fa c.fa
        let plan = parse(&["-o", "dist.csv", "a.fa", "b.fa", "c.fa"]).expect("valid invocation");
        assert_eq!(
            plan.mode,
            Mode::Distances {
                output_fp: PathBuf::from("dist.csv"),
                standardisation: Standardisation::InMemory,
            }
        );
        assert_eq!(plan.inputs.len(), 3);
        assert!(plan.standardises());
        assert_eq!(plan.emit_dir(), None);
    }

    #[test]
    fn plan_accepts_distances_and_emission_together() {
        // rusty-metal -o dist.csv --emit-standardised out/ a.fa b.fa
        let plan = parse(&["-o", "dist.csv", "--emit-standardised", "out", "a.fa", "b.fa"])
            .expect("valid invocation");
        assert_eq!(
            plan.mode,
            Mode::Distances {
                output_fp: PathBuf::from("dist.csv"),
                standardisation: Standardisation::Emit(PathBuf::from("out")),
            }
        );
        assert!(plan.standardises());
        assert_eq!(plan.emit_dir(), Some(Path::new("out")));
        assert_eq!(plan.output_fp(), Some(Path::new("dist.csv")));
    }

    #[test]
    fn plan_accepts_standardise_only_with_a_single_input_and_no_output() {
        // rusty-metal --standardise-only --emit-standardised out/ a.fa
        //
        // Two things at once, both of which the pre-Stage-5 CLI made impossible: `-o` is
        // not supplied (it used to be unconditionally required) and there is only one
        // input (`num_args = 2..` used to reject it).
        let plan = parse(&["--standardise-only", "--emit-standardised", "out", "a.fa"])
            .expect("valid invocation");
        assert_eq!(
            plan.mode,
            Mode::StandardiseOnly {
                emit_dir: PathBuf::from("out"),
            }
        );
        assert_eq!(plan.inputs, vec![PathBuf::from("a.fa")]);
        assert!(plan.standardises());
        assert_eq!(plan.output_fp(), None);
    }

    #[test]
    fn plan_accepts_the_no_standardise_escape_hatch() {
        // rusty-metal -o dist.csv --no-standardise a.fa b.fa
        let plan =
            parse(&["-o", "dist.csv", "--no-standardise", "a.fa", "b.fa"]).expect("valid invocation");
        assert_eq!(
            plan.mode,
            Mode::Distances {
                output_fp: PathBuf::from("dist.csv"),
                standardisation: Standardisation::Skip,
            }
        );
        assert!(!plan.standardises());
        assert_eq!(plan.emit_dir(), None);
    }

    // -----------------------------------------------------------------------------
    // The four illegal combinations
    // -----------------------------------------------------------------------------

    #[test]
    fn plan_rejects_standardise_only_without_emit_standardised() {
        // Nothing would be written and no distances computed: the run does nothing
        // observable.
        let err = err_of(parse(&["--standardise-only", "a.fa"]));
        assert!(err.contains("--standardise-only"), "got: {err}");
        assert!(err.contains("--emit-standardised"), "got: {err}");
        assert!(
            err.contains("nothing observable"),
            "the error should say why, got: {err}"
        );
    }

    #[test]
    fn plan_rejects_standardise_only_with_no_standardise() {
        let err = err_of(parse(&[
            "--standardise-only",
            "--no-standardise",
            "--emit-standardised",
            "out",
            "a.fa",
        ]));
        assert!(err.contains("--standardise-only"), "got: {err}");
        assert!(err.contains("--no-standardise"), "got: {err}");
        assert!(
            err.contains("contradict"),
            "the error should say the flags contradict each other, got: {err}"
        );
    }

    #[test]
    fn plan_rejects_emit_standardised_with_no_standardise() {
        let err = err_of(parse(&[
            "-o",
            "dist.csv",
            "--emit-standardised",
            "out",
            "--no-standardise",
            "a.fa",
            "b.fa",
        ]));
        assert!(err.contains("--emit-standardised"), "got: {err}");
        assert!(err.contains("--no-standardise"), "got: {err}");
        assert!(
            err.contains("nothing to emit"),
            "the error should say why, got: {err}"
        );
    }

    #[test]
    fn plan_rejects_missing_output_when_distances_are_computed() {
        let err = err_of(parse(&["a.fa", "b.fa"]));
        assert!(err.contains("--output-fp"), "got: {err}");
        assert!(
            err.contains("required"),
            "the error should say it is required, got: {err}"
        );
        assert!(
            err.contains("--standardise-only"),
            "the error should point at the alternative, got: {err}"
        );
    }

    #[test]
    fn plan_rejects_a_single_input_when_distances_are_computed() {
        // The "at least 2" rule moved out of clap's `num_args` and into the plan, so it
        // has to still fire in the mode that needs it.
        let err = err_of(parse(&["-o", "dist.csv", "a.fa"]));
        assert!(err.contains("at least 2"), "got: {err}");
        assert!(
            err.contains("--standardise-only"),
            "the error should point at the mode that accepts one file, got: {err}"
        );
    }

    // -----------------------------------------------------------------------------
    // Stem collisions
    // -----------------------------------------------------------------------------

    #[test]
    fn stem_collision_is_detected_across_directories() {
        let inputs = vec![
            PathBuf::from("a/aln.fasta"),
            PathBuf::from("b/aln.fasta"),
            PathBuf::from("c/other.fasta"),
        ];
        let err = match check_stem_collisions(&inputs) {
            Ok(()) => panic!("two inputs sharing a stem must be rejected"),
            Err(e) => e.to_string(),
        };
        assert!(err.contains("aln.standardised.fasta"), "got: {err}");
        // Both colliding inputs must be named, since the point is to say which files to
        // rename.
        assert!(err.contains("a/aln.fasta") || err.contains("a\\aln.fasta"), "got: {err}");
        assert!(err.contains("b/aln.fasta") || err.contains("b\\aln.fasta"), "got: {err}");
        assert!(
            !err.contains("other.fasta"),
            "the non-colliding input is not the problem, got: {err}"
        );
    }

    #[test]
    fn distinct_stems_do_not_collide() {
        let inputs = vec![
            PathBuf::from("a/one.fasta"),
            PathBuf::from("b/two.fasta"),
            // Same stem *text* but a different extension is still a distinct stem only
            // if the stems differ — `one.fa` and `one.fasta` share the stem `one`, so
            // this list deliberately does not include such a pair.
            PathBuf::from("three.fa"),
        ];
        assert!(check_stem_collisions(&inputs).is_ok());
    }

    #[test]
    fn a_stem_collision_is_caught_at_plan_time_when_emitting() {
        let err = err_of(parse(&[
            "-o",
            "dist.csv",
            "--emit-standardised",
            "out",
            "a/aln.fasta",
            "b/aln.fasta",
        ]));
        assert!(err.contains("share a file stem"), "got: {err}");
    }

    #[test]
    fn a_stem_collision_is_not_an_error_when_nothing_is_emitted() {
        // Without --emit-standardised no file is written, so a shared stem is harmless.
        parse(&["-o", "dist.csv", "a/aln.fasta", "b/aln.fasta"])
            .expect("a stem collision only matters when emitting");
    }

    // -----------------------------------------------------------------------------
    // CSV escaping
    // -----------------------------------------------------------------------------

    #[test]
    fn csv_fields_are_escaped_per_rfc_4180() {
        // Plain: emitted verbatim, no quotes added.
        assert_eq!(csv_escape("results/aln.fasta"), "results/aln.fasta");
        // Comma: quoted, so it stays one field.
        assert_eq!(csv_escape("run,v2/aln.fasta"), "\"run,v2/aln.fasta\"");
        // Double quote: quoted, and the embedded quote is doubled.
        assert_eq!(csv_escape("say \"hi\".fasta"), "\"say \"\"hi\"\".fasta\"");
        // Newline: quoted, and the newline is preserved inside the quotes.
        assert_eq!(csv_escape("two\nlines.fasta"), "\"two\nlines.fasta\"");
        // Carriage return counts too, on its own.
        assert_eq!(csv_escape("cr\rhere"), "\"cr\rhere\"");
    }

    #[test]
    fn a_row_with_an_awkward_path_stays_three_fields() {
        let results = vec![
            DistanceResult {
                msa_a: PathBuf::from("run,v2/a.fasta"),
                msa_b: PathBuf::from("plain/b.fasta"),
                distance: Some(0.5),
            },
            DistanceResult {
                msa_a: PathBuf::from("q\"uote.fasta"),
                msa_b: PathBuf::from("plain/b.fasta"),
                distance: None,
            },
        ];

        let mut buffer: Vec<u8> = Vec::new();
        write_results(&mut buffer, &results).expect("writing to a Vec cannot fail");
        let csv = String::from_utf8(buffer).expect("output is UTF-8");

        assert_eq!(
            csv,
            "msa_a,msa_b,distance\n\
             \"run,v2/a.fasta\",plain/b.fasta,0.5\n\
             \"q\"\"uote.fasta\",plain/b.fasta,\n"
        );
    }

    // -----------------------------------------------------------------------------
    // End to end over the real fixtures
    // -----------------------------------------------------------------------------

    /// The distance the pipeline reports for one pair, read back out of the CSV it
    /// wrote. Exercises the whole of `run`, not just the arithmetic.
    fn distance_from_a_real_run(extra_flags: &[&str]) -> String {
        let dir = temp_dir("e2e");
        let csv = dir.join("dist.csv");
        let csv_arg = csv.display().to_string();

        let mut argv: Vec<&str> = vec!["-o", &csv_arg];
        argv.extend_from_slice(extra_flags);
        argv.extend_from_slice(&["test/test.fasta", "test/test2.fasta"]);

        let plan = parse(&argv).expect("valid invocation");
        run(&plan).expect("the run should succeed");

        let contents = read_to_string(&csv);
        let row = contents
            .lines()
            .nth(1)
            .unwrap_or_else(|| panic!("expected a data row, got:\n{contents}"))
            .to_string();
        let field = row
            .rsplit(',')
            .next()
            .expect("a row always has a last field")
            .to_string();

        std::fs::remove_dir_all(&dir).expect("cleanup should succeed");
        field
    }

    /// As `distance_from_a_real_run`, but for an arbitrary pair of fixtures.
    fn distance_between_fixtures(extra_flags: &[&str], a: &str, b: &str) -> String {
        let dir = temp_dir("pair");
        let csv = dir.join("dist.csv");
        let csv_arg = csv.display().to_string();

        let mut argv: Vec<&str> = vec!["-o", &csv_arg];
        argv.extend_from_slice(extra_flags);
        argv.extend_from_slice(&[a, b]);

        let plan = parse(&argv).expect("valid invocation");
        run(&plan).expect("the run should succeed");

        let contents = read_to_string(&csv);
        let field = contents
            .lines()
            .nth(1)
            .unwrap_or_else(|| panic!("expected a data row, got:\n{contents}"))
            .rsplit(',')
            .next()
            .expect("a row always has a last field")
            .to_string();

        std::fs::remove_dir_all(&dir).expect("cleanup should succeed");
        field
    }

    #[test]
    fn a_legally_permuted_alignment_is_at_distance_zero() {
        // *** THE POINT OF THE WHOLE MERGE, end to end. ***
        //
        // `test/test_legal_permutation.fasta` is `test/test.fasta` with column 3 moved
        // left past columns 2 and 1. Both moves are legal — no row holds a residue in
        // column 3 and either of the columns it crosses — so the two files hold the
        // *same alignment*, differing only in where the gaps sit:
        //
        //   test.fasta      test_legal_permutation.fasta   (columns [0, 3, 1, 2])
        //   >1 AA--         >1 A-A-
        //   >2 A--A         >2 AA--
        //   >3 AAA-         >3 A-AA
        //
        // Standardised, they must be identical, so the distance must be exactly 0.
        //
        // Under the pre-Stage-6 ordering rule this was NOT 0 — that rule left a free
        // column wherever the input file happened to put it, so the same alignment
        // written two ways standardised to two different layouts. See
        // `canonical_columns` for the full account.
        assert_eq!(
            distance_between_fixtures(&[], "test/test.fasta", "test/test_legal_permutation.fasta"),
            "0"
        );

        // Argument order must not matter either.
        assert_eq!(
            distance_between_fixtures(&[], "test/test_legal_permutation.fasta", "test/test.fasta"),
            "0"
        );

        // And without standardisation it is emphatically not 0 — which is what makes
        // this a test of standardisation rather than of the metric being trivial.
        //
        // 10/28, worked by hand: the four column-0 slots and the shared residues in
        // slots (0,1), (2,1) and (2,2) agree, and every gap identity that moved with
        // column 3 disagrees, contributing 2 + 4 + 2 + 2. This is the size of the
        // gap-placement artefact that standardisation removes on this pair.
        assert_eq!(
            distance_between_fixtures(
                &["--no-standardise"],
                "test/test.fasta",
                "test/test_legal_permutation.fasta"
            ),
            "0.35714285714285715"
        );
    }

    #[test]
    fn end_to_end_without_standardisation_reproduces_the_baseline() {
        // The pre-merge number, and the one Stages 1-4 preserved. `--no-standardise`
        // skips only the standardise stage, so this is the value that pins "nothing
        // except standardisation moved".
        assert_eq!(
            distance_from_a_real_run(&["--no-standardise"]),
            "0.21428571428571427"
        );
    }

    #[test]
    fn end_to_end_with_standardisation_pins_the_new_fixture_distance() {
        // *** This number is a RESULT of wiring standardisation in, not a target. ***
        //
        // Reworked by hand in Stage 6, when the canonical ordering rule changed. Both
        // sides now move; under the old rule `test.fasta` was a fixed point.
        //
        //   test/test.fasta   standardised   test/test2.fasta   standardised
        //   >1 AA--           >1 AA--        >1 A-A-            >1 A--A
        //   >2 A--A           >2 A-A-        >2 A--A            >2 AA--
        //   >3 AAA-           >3 AA-A        >3 AAA-            >3 A-AA
        //
        // Registry sorted by name ("1" -> 0, "2" -> 1, "3" -> 2); r(s,p) is a residue
        // and g(s,p) a gap whose position is the index of the residue preceding it in
        // its row. A residue's homology set is its whole column minus itself:
        //
        //   seq pos | A(r)               | B(r)               | |AΔB| | |A|+|B|
        //   --------+--------------------+--------------------+-------+--------
        //   0   0   | {r(1,0), r(2,0)}   | {r(1,0), r(2,0)}   |   0   |   4
        //   0   1   | {g(1,0), r(2,1)}   | {g(1,1), r(2,2)}   |   4   |   4
        //   1   0   | {r(0,0), r(2,0)}   | {r(0,0), r(2,0)}   |   0   |   4
        //   1   1   | {g(0,1), g(2,1)}   | {g(0,0), g(2,0)}   |   4   |   4
        //   2   0   | {r(0,0), r(1,0)}   | {r(0,0), r(1,0)}   |   0   |   4
        //   2   1   | {r(0,1), g(1,0)}   | {g(0,0), g(1,1)}   |   4   |   4
        //   2   2   | {g(0,1), g(1,1)}   | {r(0,1), g(1,1)}   |   2   |   4
        //   --------+--------------------+--------------------+-------+--------
        //                                                sum:   14       28
        //
        // 14 / 28 = 0.5, against 6/28 = 0.21428571428571427 unstandardised. The value
        // is unchanged from Stage 5 by coincidence, not by construction — the table
        // that produces it is entirely different. Every column-0 slot agrees (column 0
        // is all residues in both files and cannot move), and every slot that involves
        // a gap identity disagrees.
        //
        // The distance went *up* relative to not standardising, which is worth being
        // explicit about because the merge is motivated as "strips gap-placement
        // artefacts out of the metric". It does strip them — but stripping is not
        // shrinking. Standardisation moves each alignment to its *own* canonical layout;
        // it does not move two different alignments toward each other.
        //
        // These two fixtures are genuinely different alignments: `test2` is `test` with
        // two *pinned* columns exchanged, which is not a legal permutation, so nothing
        // requires them to converge. The case standardisation is for — a pair differing
        // only by a legal permutation — now provably goes to 0, which is what
        // `standardise::tests::standardisation_is_confluent_over_legal_permutations`
        // asserts and what the old rule failed to do.
        assert_eq!(distance_from_a_real_run(&[]), "0.5");
    }

    // -----------------------------------------------------------------------------
    // Phase 1
    // -----------------------------------------------------------------------------

    #[test]
    fn phase_one_emits_one_standardised_file_per_input() {
        let dir = temp_dir("emit");
        let out = dir.join("standardised");
        let out_arg = out.display().to_string();

        let plan = parse(&[
            "--standardise-only",
            "--emit-standardised",
            &out_arg,
            "test/test.fasta",
            "test/test2.fasta",
        ])
        .expect("valid invocation");

        // The directory does not exist yet: phase 1 must create it.
        assert!(!out.exists());
        run(&plan).expect("the run should succeed");

        let emitted_a = out.join("test.standardised.fasta");
        let emitted_b = out.join("test2.standardised.fasta");
        assert!(emitted_a.exists(), "expected {}", emitted_a.display());
        assert!(emitted_b.exists(), "expected {}", emitted_b.display());

        // Exactly two files, so nothing extra was written and nothing overwrote
        // anything.
        let count = std::fs::read_dir(&out)
            .expect("the emit dir should be readable")
            .count();
        assert_eq!(count, 2, "one output per input");

        // Both fixtures move under the Stage 6 canonical rule; the CRLF and missing
        // final newline of the fixtures are normalised on write. These two values are
        // derived in `standardise::tests::the_test_fixture_standardises_to_its_canonical_form`
        // and `..._the_test2_fixture_...`, which show the column working.
        assert_eq!(read_to_string(&emitted_a), ">1\nAA--\n>2\nA-A-\n>3\nAA-A\n");
        assert_eq!(read_to_string(&emitted_b), ">1\nA--A\n>2\nAA--\n>3\nA-AA\n");

        std::fs::remove_dir_all(&dir).expect("cleanup should succeed");
    }

    #[test]
    fn standardise_only_writes_no_csv() {
        let dir = temp_dir("nocsv");
        let out = dir.join("standardised");
        let out_arg = out.display().to_string();

        let plan = parse(&[
            "--standardise-only",
            "--emit-standardised",
            &out_arg,
            "test/test.fasta",
        ])
        .expect("valid invocation");
        run(&plan).expect("the run should succeed");

        assert!(out.join("test.standardised.fasta").exists());
        assert_eq!(plan.output_fp(), None);

        std::fs::remove_dir_all(&dir).expect("cleanup should succeed");
    }

    #[test]
    fn standardised_output_path_names_the_file_from_the_stem() {
        assert_eq!(
            standardised_output_path(Path::new("out"), Path::new("a/b/aln.fasta"))
                .expect("a stem exists"),
            PathBuf::from("out").join("aln.standardised.fasta")
        );
    }

    // -----------------------------------------------------------------------------
    // Failure policy: hard for a file, soft for a pair
    // -----------------------------------------------------------------------------

    #[test]
    fn an_unreadable_input_is_a_hard_error() {
        // A file that will not read invalidates every pair it appears in, so phase 1
        // aborts rather than dropping those pairs silently.
        let dir = temp_dir("hardfail");
        let csv = dir.join("dist.csv");
        let csv_arg = csv.display().to_string();

        // `test/ragged.fasta` parses as FASTA but is not a valid alignment.
        let plan = parse(&[
            "-o",
            &csv_arg,
            "test/test.fasta",
            "test/ragged.fasta",
            "test/test2.fasta",
        ])
        .expect("valid invocation");

        let err = match run(&plan) {
            Ok(()) => panic!("a file that fails to read must abort the run"),
            Err(e) => format!("{e:#}"),
        };
        assert!(err.contains("ragged.fasta"), "the error must name the file, got: {err}");
        assert!(
            !csv.exists(),
            "no CSV should be written when phase 1 fails"
        );

        // And a file that does not exist at all.
        let missing_plan = parse(&["-o", &csv_arg, "test/test.fasta", "test/nope.fasta"])
            .expect("valid invocation");
        assert!(run(&missing_plan).is_err(), "a missing file must abort the run");

        std::fs::remove_dir_all(&dir).expect("cleanup should succeed");
    }

    #[test]
    fn a_single_failing_pair_is_logged_and_the_run_continues() {
        // `test/case_and_ambiguity.fasta` reads fine, so phase 1 succeeds, but its
        // sequence names ("seq1", "seq2") do not match the fixtures' ("1", "2", "3"), so
        // `Registry::for_pair` rejects both pairs it takes part in. The third pair must
        // still be computed and must still reach the CSV.
        let dir = temp_dir("softfail");
        let csv = dir.join("dist.csv");
        let csv_arg = csv.display().to_string();

        let plan = parse(&[
            "-o",
            &csv_arg,
            "test/test.fasta",
            "test/test2.fasta",
            "test/case_and_ambiguity.fasta",
        ])
        .expect("valid invocation");
        run(&plan).expect("a failing pair must not fail the run");

        let contents = read_to_string(&csv);
        let rows: Vec<&str> = contents.lines().skip(1).collect();
        assert_eq!(rows.len(), 3, "every pair gets a row, got:\n{contents}");

        let with_a_distance: Vec<&&str> = rows.iter().filter(|r| !r.ends_with(',')).collect();
        assert_eq!(
            with_a_distance.len(),
            1,
            "only the one comparable pair should carry a distance, got:\n{contents}"
        );
        assert!(
            with_a_distance[0].ends_with(",0.5"),
            "the surviving pair is the standardised fixture pair, got: {}",
            with_a_distance[0]
        );

        std::fs::remove_dir_all(&dir).expect("cleanup should succeed");
    }

    // -----------------------------------------------------------------------------
    // The rayon global-pool trap
    // -----------------------------------------------------------------------------

    #[test]
    fn the_pipeline_can_be_run_twice_in_one_process() {
        // `rayon::ThreadPoolBuilder::build_global` can only succeed once per process,
        // and the previous `process()` called it on every invocation — so a second run
        // in the same process failed with "The global thread pool has already been
        // initialized". `run` builds a scoped pool instead, so this passes.
        //
        // Both runs also ask for a specific thread count, which is the part a
        // tolerate-the-error or `OnceLock` approach would silently drop on the second
        // call.
        let dir = temp_dir("twice");

        for (n, threads) in [(0usize, "2"), (1usize, "3")] {
            let csv = dir.join(format!("dist{n}.csv"));
            let csv_arg = csv.display().to_string();
            let plan = parse(&[
                "-o",
                &csv_arg,
                "-n",
                threads,
                "test/test.fasta",
                "test/test2.fasta",
            ])
            .expect("valid invocation");
            run(&plan).unwrap_or_else(|e| panic!("run {n} should succeed: {e:#}"));
            assert!(read_to_string(&csv).contains("0.5"));
        }

        std::fs::remove_dir_all(&dir).expect("cleanup should succeed");
    }
}
