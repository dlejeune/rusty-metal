use crate::datastructures::MsaHashSets;
use crate::distance::compute_symmetric_difference;
use crate::utils::{create_hashsets, read_msa};
use anyhow::{Context, Result};
use clap::builder::styling;
use clap::Parser;
use colored::Colorize;
use itertools::Itertools;
use rayon::prelude::*;
use std::cmp::max;
use std::fmt::format;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
mod datastructures;
mod distance;
mod utils;

const VERSION: &str = env!("CARGO_PKG_VERSION");
const STYLES: styling::Styles = styling::Styles::styled()
    .header(styling::AnsiColor::Green.on_default().bold())
    .usage(styling::AnsiColor::Green.on_default().bold())
    .literal(styling::AnsiColor::Blue.on_default().bold())
    .placeholder(styling::AnsiColor::Cyan.on_default())
    .error(styling::AnsiColor::Red.on_default().bold());
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
#[command(styles = STYLES)]
struct Args {
    // MSAs to get distances for
    #[clap(required = true, num_args=2..)]
    input_files: Vec<PathBuf>,

    #[clap(short, long)]
    output_fp: PathBuf,

    #[clap(short, long, default_value_t = 0)]
    num_threads: usize,
}

fn compare_alignment_pair(first_msa_fp: &PathBuf, second_msa_fp: &PathBuf) -> Result<f64> {
    let msa_a = read_msa(first_msa_fp)?;
    let msa_b = read_msa(second_msa_fp)?;
    let hom_set_a = create_hashsets(&msa_a)?;
    let hom_set_b = create_hashsets(&msa_b)?;
    let width = max(msa_a.width, msa_b.width);
    let symmetric_difference =
        compute_symmetric_difference(&hom_set_a, &hom_set_b, msa_a.num_seqs, width)?;
    Ok(symmetric_difference)
}

struct DistanceResult {
    msa_a: PathBuf,
    msa_b: PathBuf,
    distance: Option<f64>,
}
fn process(files: Vec<PathBuf>, n_threads: usize) -> Result<Vec<DistanceResult>> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(n_threads)
        .build_global()
        .context("Could not set up multithreading")?;

    let output: Vec<DistanceResult> = files
        .iter()
        .combinations(2)
        .into_iter()
        .enumerate()
        .par_bridge()
        .map(|(curr_idx, pair)| DistanceResult {
            msa_a: pair[0].clone(),
            msa_b: pair[1].clone(),
            distance: match compare_alignment_pair(pair[0], pair[1]) {
                Ok(distance) => {
                    log::info!("Comparison {curr_idx} complete.");
                    log::info!(
                        "{}",
                        format!(
                            "{} -> {}: {}",
                            pair[0].file_stem().unwrap().to_str().unwrap(),
                            pair[1].file_stem().unwrap().to_str().unwrap(),
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
                        pair[0].as_path().display(),
                        pair[1].as_path().display()
                    );
                    log::warn!("{}", e);
                    None
                }
            },
        })
        .collect();

    Ok(output)
}

fn main() -> Result<()> {
    simple_logger::SimpleLogger::new().env().init()?;
    let args = Args::parse();

    log::info!("This is rusty-metal version {}", VERSION.bold().cyan());
    log::info!(
        "Computing the distance between the following MSAs, and writing the output to {} and using {} threads.",
        args.output_fp.display().to_string().bright_blue(),
        args.num_threads.to_string().bright_red()
    );
    for file in &args.input_files {
        log::info!(
            "> {}",
            file.file_name()
                .with_context(|| format!("Expected the file {} to have a name", {
                    file.display()
                }))?
                .display()
                .to_string()
                .cyan()
        );
    }

    let results = process(args.input_files, args.num_threads)?;

    // let distance = compare_alignment_pair(&args.first_msa, &args.second_msa)?;
    log::info!("Writing output to {}", args.output_fp.display());
    let mut writer = BufWriter::new(File::create(args.output_fp)?);
    writer.write("msa_a,msa_b,distance\n".as_bytes())?;
    results.iter().for_each(|result| {
        writer
            .write(
                format!(
                    "{},{},{}\n",
                    result.msa_a.display(),
                    result.msa_b.display(),
                    match result.distance {
                        Some(distance) => distance.to_string(),
                        None => "".to_string(),
                    }
                )
                .as_bytes(),
            )
            .expect("Failed to write result");
    });

    Ok(())
}
