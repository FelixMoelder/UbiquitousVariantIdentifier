mod cli;

use anyhow::Result;
use cli::ValidateArguments;
use rust_htslib::bcf;
use rust_htslib::bcf::record::{Buffer, BufferBacked};
use rust_htslib::bcf::Read;
use std::collections::{HashMap, HashSet};
use std::path::PathBuf;
use structopt::StructOpt;

fn main() -> anyhow::Result<()> {
    let opt = cli::UbiquitousVariantIdentifier::from_args();
    opt.validate()?;
    let mut variant_hashmap: HashMap<Vec<u8>, HashMap<Vec<u8>, usize>> = HashMap::new();
    for bcf_path in opt.bcf_paths.iter() {
        let bcf_variants_hashmap = gather_bcf_variants(bcf_path)?;
        update_global_variants(&mut variant_hashmap, bcf_variants_hashmap)?;
    }
    let bcf_count = opt.bcf_paths.len();
    print_ubiquitous_variants(variant_hashmap, bcf_count, opt.threshold)?;
    Ok(())
}

// Read all unique canonical variants from bcf file
fn gather_bcf_variants(bcf_path: &PathBuf) -> Result<HashMap<Vec<u8>, HashSet<Vec<u8>>>> {
    let mut bcf_file = bcf::Reader::from_path(bcf_path)?;
    let mut bcf_variants_hashmap = HashMap::new();
    for result in bcf_file.records() {
        let record = result?;
        let ann = record.info(b"ANN").string()?;
        let (gene, hgvsp_opt) = extract_canonical_transcript(ann)?;
        if let Some(hgvsp) = hgvsp_opt {
            let bcf_hgvsp_values: &mut HashSet<Vec<u8>> =
                bcf_variants_hashmap.entry(gene.unwrap()).or_default();
            bcf_hgvsp_values.insert(hgvsp);
        }
    }
    Ok(bcf_variants_hashmap)
}

// Extract canonical variants from bcf-record
fn extract_canonical_transcript(
    ann: Option<BufferBacked<'_, Vec<&[u8]>, Buffer>>,
) -> Result<(Option<Vec<u8>>, Option<Vec<u8>>)> {
    let mut gene = None;
    let mut hgvsp = None;
    if let Some(transcripts) = ann {
        for &transcript in transcripts.iter() {
            let entries: Vec<&[u8]> = transcript.split(|pos| *pos == b'|').collect();
            // Select gene and HGVSp from canonical transcript
            if entries[24] == b"YES" || !entries[26].is_empty() {
                gene = Some(entries[3].to_vec());
                if !entries[11].is_empty() {
                    hgvsp = entries[11].split(|x| *x == b':').last().map(|s| s.to_vec());
                };
                break;
            }
        }
    };
    Ok((gene, hgvsp))
}

// Update global variant hashmap with variants from current bcf record
fn update_global_variants(
    variant_hashmap: &mut HashMap<Vec<u8>, HashMap<Vec<u8>, usize>>,
    bcf_variant_hashmap: HashMap<Vec<u8>, HashSet<Vec<u8>>>,
) -> Result<()> {
    bcf_variant_hashmap
        .iter()
        .for_each(|(gene, bcf_hgvsp_hashset)| {
            bcf_hgvsp_hashset.iter().for_each(|hgvsp| {
                let global_hgvsp_hashmap = variant_hashmap.entry(gene.to_vec()).or_default();
                global_hgvsp_hashmap
                    .entry(hgvsp.to_vec())
                    .and_modify(|count| *count += 1)
                    .or_insert(1);
            });
        });
    Ok(())
}

// Output all ubiquitous variants that exceed threshold
fn print_ubiquitous_variants(
    variant_hashmap: HashMap<Vec<u8>, HashMap<Vec<u8>, usize>>,
    bcf_count: usize,
    threshold: f32,
) -> Result<()> {
    variant_hashmap.iter().for_each(|(gene, hgvsp_hashmap)| {
        hgvsp_hashmap.iter().for_each(|(hgvsp, variant_count)| {
            let ratio = *variant_count as f32 / bcf_count as f32;
            if ratio > threshold {
                println!(
                    "Gene: {:?}, HGVSP: {:?}, Variant Count: {}, Ratio: {}",
                    String::from_utf8_lossy(gene),
                    String::from_utf8_lossy(hgvsp),
                    variant_count,
                    ratio
                );
            }
        })
    });
    Ok(())
}
