use parse_args;
use std::io::{BufRead, BufReader};
use std::fs::File;
use std::collections::HashSet;

use filter::sv_locus_identifiers as sv_locus_identifiers;

const USAGE: &'static str = "
Usage:
  breakfast blacklist [options] <sv_files>...

Options:
  --freq-above=FREQ     Minimum samples(0.0 - 1.0) to consider [default: 0.0]
";

pub fn main() {

    let args = parse_args(USAGE);
    let min_frequency: f32 = args.get_str("--freq-above").parse()
        .expect("\tERROR: *** Invalid frequency! *** \n");
    let sv_files = args.get_vec("<sv_files>").to_vec();
    let mut sample_variants: Vec<Vec<String>> = vec![Vec::new(); sv_files.len()];

    for (s, sv_file) in sv_files.iter().enumerate() {

        let sv = BufReader::new(File::open(&sv_file).unwrap());
        let mut tmp: HashSet<String> = HashSet::new();

        for l in sv.lines() {
            let line = l.unwrap();
            if !line.starts_with("chr") { continue; }

            let tokens: Vec<&str> = line.split('\t').collect();
            let chrom = tokens[0];
            let pos: usize = tokens[2].parse().unwrap();
            let tmp1: Vec<String> = sv_locus_identifiers(chrom, pos, 5000).into_iter().collect();
            sample_variants[s].append(&mut tmp1.to_vec());

            let chrom = tokens[5];
            let pos: usize = tokens[7].parse().unwrap();

            let tmp2: Vec<String> = sv_locus_identifiers(chrom, pos, 5000).into_iter().collect();
            sample_variants[s].append(&mut tmp2.to_vec());

            //tmp = tmp1.union(&tmp2).cloned().collect();
            //println!("{}", tmp.len());
            //sample_variants[s] = tmp;
        }

    }


     let variants_cp = sample_variants.to_vec();
     let mut blacklist: HashSet<String> = HashSet::new();
     for loci in variants_cp {
        for l in loci {
         blacklist.insert(l);
        }
     }
     blacklist = natural_sorted(blacklist);
     eprintln!("blacklist loci before filtering {}", blacklist.len());

     let mut frequency: Vec<f32> = vec![0.0; blacklist.len()];
     for (k, bad_variant) in blacklist.iter().enumerate() {
        let mut bad_in_sample: Vec<usize> = vec![0; sample_variants.len()];
        for (n, loci) in sample_variants.iter().enumerate() {
            if loci.contains(bad_variant) {
                bad_in_sample[n] = 1;
                }
            }
        let sum: usize = bad_in_sample.iter().sum();
        let freq: f32  = sum as f32 / bad_in_sample.len() as f32;
        frequency[k]   = freq;
    }


    for (x, loci) in blacklist.iter().enumerate() {
        if frequency[x] >= min_frequency {
            println!("{}", loci);
        }
    }
}


pub fn natural_sorted(blacklist: HashSet<String>) -> HashSet<String> {
    let mut v: Vec<String> = blacklist.into_iter().collect();
    v.sort_by(|a, b| a.cmp(b));
    let out: HashSet<_> = v.into_iter().collect();
    out
}
