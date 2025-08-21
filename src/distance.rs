use crate::datastructures::{IsGap, MsaHashSets, MultipleSequenceAlignment, SequenceElement};
use anyhow::Result;
use colored::Colorize;
use itertools::Itertools;
use log::log;
use rayon::prelude::*;
use std::collections::HashSet;

fn compute_jaccard_distance(homology_set_a: &MsaHashSets, homology_set_b: &MsaHashSets) -> f32 {
    let mut union_sum: usize = 0;
    let mut intersection_sum: usize = 0;
    for seq_idx in 0..homology_set_a.len() {
        let mut seq_element_idx = 0;
        while let Some(element_homology_set_a) = &homology_set_a[seq_idx][seq_element_idx].clone() {
            let element_homology_set_b = &homology_set_b[seq_idx][seq_element_idx].clone().unwrap();
            union_sum += element_homology_set_a.union(element_homology_set_b).count();
            intersection_sum += element_homology_set_a
                .intersection(element_homology_set_b)
                .count();
            seq_element_idx += 1;
        }
    }
    println!("{}/{}", intersection_sum, union_sum);
    1f32 - (intersection_sum as f32) / (union_sum as f32)
}

pub fn compute_symmetric_difference(
    homology_set_a: &MsaHashSets,
    homology_set_b: &MsaHashSets,
    length: usize,
    width: usize,
) -> Result<f64> {
    log::info!("{}", "Computing Symmetric Difference".bold().purple());

    struct SymmetricDistanceResult {
        distance: usize,
        homology_set_size: usize,
    }
    // Note: the reason why we may sometimes fail to get a column or a hashset is because
    // the hashset vectors are pre-populated by None values before hand. getting a None element
    // hashet is indicative of reaching the end of the sequence.
    // Getting a none hashset on the column is because we iterate over the width of the longest MSA
    // which will always be longer than the sequences (or at least the same size).
    let symmetric_distances: Vec<Option<SymmetricDistanceResult>> = (0..length)
        .cartesian_product(0..width)
        .par_bridge()
        .map(|(x, y)| {
            if let (Some(sequence_homology_sets_a), Some(sequence_homology_sets_b)) =
                (homology_set_a.get(x), &homology_set_b.get(x))
            {
                if let (Some(column_homology_sets_a), Some(column_homology_sets_b)) = (
                    sequence_homology_sets_a.get(y),
                    sequence_homology_sets_b.get(y),
                ) {
                    if let (Some(element_homology_set_a), Some(element_homology_set_b)) =
                        (column_homology_sets_a, column_homology_sets_b)
                    {
                        Some(SymmetricDistanceResult {
                            distance: element_homology_set_a
                                .symmetric_difference(element_homology_set_b)
                                .count(),
                            homology_set_size: element_homology_set_a.len(),
                        })
                    } else {
                        None
                    }
                } else {
                    None
                }
            } else {
                log::warn!("The sequence with id {x} failed to be found in this hashset.");
                None
            }
        })
        .collect();

    let total_distance: usize = symmetric_distances
        .iter()
        .map(|d| match d {
            Some(d) => d.distance,
            None => 0,
        })
        .sum();
    let total_length: usize = symmetric_distances
        .iter()
        .map(|d| match d {
            Some(d) => d.homology_set_size * 2,
            None => 0,
        })
        .sum();
    // let total_chars = homology_set_a.iter()
    // let mut symmetric_difference_sum = 0;
    // let mut num_chars = 0;
    // let mut symmetric_diff_seqs: Vec<f64> = Vec::with_capacity(homology_set_a.len());
    //
    // for seq_idx in 0..homology_set_a.len() {
    //     let mut seq_dist = 0;
    //     let mut seq_dist_sum = 0;
    //     let mut seq_dist_total_chars = 0;
    //     let mut seq_element_idx = 0;
    //
    //     while let (Some(next_hom_set_a), Some(next_hom_set_b)) = (
    //         homology_set_a[seq_idx].get(seq_element_idx),
    //         homology_set_b[seq_idx].get(seq_element_idx),
    //     ) {
    //         if let (Some(element_homology_set_a), Some(element_homology_set_b)) =
    //             (next_hom_set_a, next_hom_set_b)
    //         {
    //             seq_dist = element_homology_set_a
    //                 .symmetric_difference(element_homology_set_b)
    //                 .count();
    //
    //             seq_dist_sum += seq_dist;
    //
    //             seq_dist_total_chars += 2 * element_homology_set_a.len();
    //             seq_element_idx += 1;
    //         } else {
    //             break;
    //         }
    //     }
    //     symmetric_difference_sum += seq_dist_sum;
    //     num_chars += seq_dist_total_chars;
    //     symmetric_diff_seqs.push((seq_dist as f64) / (seq_dist_total_chars as f64));
    //     // log::info!(
    //     //     "Seq {}: {}/{} = {}",
    //     //     seq_idx,
    //     //     seq_dist_sum,
    //     //     seq_dist_total_chars,
    //     //     seq_dist_sum as f64 / seq_dist_total_chars as f64,
    //     // )
    // }
    // let distance = (symmetric_difference_sum as f64) / (num_chars as f64);
    let distance: f64 = (total_distance as f64) / (total_length as f64);
    // log::info!(
    //     "{}",
    //     format!("{}/{} = {}", total_distance, total_length, distance)
    //         .bold()
    //         .green()
    // );
    Ok(distance)
}
