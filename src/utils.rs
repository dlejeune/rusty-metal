use crate::datastructures::{
    ElementFromString, IsGap, MsaHashSets, MultipleSequenceAlignment, Sequence, SequenceElement,
};
use anyhow::{Context, Result};
use seq_io::fasta::Reader;
use std::collections::HashSet;
use std::path::{Path, PathBuf};

// pub fn seq_from_iterator<T: IntoIterator<Item = SequenceElement>>(
//     iter: T,
//     seq_name: String,
// ) -> Sequence {
//     Sequence {
//         name: seq_name,
//         elements: iter.into_iter().collect(),
//     }
// }
pub fn create_hashsets<'a>(msa: &'a MultipleSequenceAlignment) -> Result<MsaHashSets<'a>> {
    let mut msa_homology_sets: MsaHashSets<'a> = Vec::with_capacity(msa.num_seqs);
    for seq_idx in 0..msa.num_seqs {
        msa_homology_sets.push(Vec::with_capacity(msa.width));
        for _col_idx in 0..msa.width {
            msa_homology_sets[seq_idx].push(None);
        }
    }

    for col_idx in 0..msa.width {
        let column = &msa.columns[col_idx];
        // Add this if doing ssp metric .filter(|x| !x.is_gap())
        let column_hashset: HashSet<&SequenceElement> = HashSet::from_iter(column.iter());
        for seq_idx in 0..column.len() {
            let item = column[seq_idx];
            if !item.is_gap() {
                let mut item_hashset: HashSet<&SequenceElement> = column_hashset.clone();
                item_hashset.remove(&item);
                msa_homology_sets[seq_idx][item.position.unwrap()] = Some(item_hashset);
            }
        }
    }

    Ok(msa_homology_sets)
}
//
// pub fn load_msas(files: Vec<PathBuf>) -> Result<MsaHashSets> {
//     // Create index of all sequence elements needed. We shouldn't need to re-assign after this.
//     let mut base_seq_elements: Vec<Sequence> = Vec::new();
//     let mut reader = Reader::from_path(&files[0])?;
//     let mut counter = 0;
//     let mut sequences: Vec<Sequence> = Vec::new();
//
//     while let Some(record) = reader.next() {
//         let record = record?;
//         base_seq_elements.push(
//             Sequence::from_characters(record.owned_seq().to_ascii_uppercase(), &counter)
//                 .filter(|c| !c.is_gap())
//                 .collect::<Sequence>(),
//         );
//
//         sequences.push(Sequence::from_characters(
//             record.owned_seq().to_ascii_uppercase(),
//             &counter,
//         ));
//         // ids.push(String::from_utf8(record.id_bytes().to_owned()).unwrap());
//         counter += 1
//     }
//     Ok()
// }
pub fn read_msa<P: AsRef<Path>>(path: P) -> Result<MultipleSequenceAlignment> {
    log::info!("Reading msa file: {}", path.as_ref().display());
    let mut reader = Reader::from_path(&path)?;
    let mut counter = 0;
    let mut sequences: Vec<Sequence> = Vec::new();

    while let Some(record) = reader.next() {
        let record = record?;

        sequences.push(Sequence::from_characters(
            record.owned_seq().to_ascii_uppercase(),
            &counter,
        ));
        // ids.push(String::from_utf8(record.id_bytes().to_owned()).unwrap());
        counter += 1;
    }
    counter -= 1;
    // We have to transpose the columns here
    let width = sequences
        .get(0)
        .with_context(|| {
            format!(
                "Tried to get the first sequence from the msa {} but failed.",
                path.as_ref().display()
            )
        })?
        .elements
        .len();
    let num_seqs = sequences.len();
    let mut columns: Vec<Vec<SequenceElement>> = Vec::with_capacity(width);

    for col in 0..width {
        let mut column: Vec<SequenceElement> = Vec::with_capacity(num_seqs);
        for seq_idx in 0..num_seqs {
            column.push(sequences[seq_idx].elements[col]);
        }
        // for seq in sequences {
        //     column.push(seq.elements.get(col).unwrap());
        // }
        columns.push(column);
    }

    Ok(MultipleSequenceAlignment {
        columns,
        num_seqs: sequences.len(),
        width,
    })
}
