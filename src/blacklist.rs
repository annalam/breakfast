use std::io::{BufRead, BufReader};
use std::fs::File;
use std::collections::HashSet;

use filter::sv_locus_identifiers as sv_locus_identifiers;

pub fn natural_sorted(blacklist: HashSet<String>) -> HashSet<String> {
    let mut v: Vec<String> = blacklist.into_iter().collect();
    v.sort_by(|a, b| a.cmp(b)); 
    let out: HashSet<_> = v.into_iter().collect();
    out
}

pub fn generate_blacklist(sv_files: Vec<&str>, min_freq: usize) {
    println!("Total samples  {}", sv_files.len());

    let mut sample_variants: Vec<HashSet<String>> = Vec::with_capacity(sv_files.len());

    println!("Before {}", sample_variants.len());

    for (s, sv_file) in sv_files.iter().enumerate() {
        let sv = BufReader::new(File::open(&sv_file).unwrap());

        for l in sv.lines() {
            let line = l.unwrap();
            if !line.starts_with("chr") { continue; }

            let tokens: Vec<&str> = line.split('\t').collect();
            let chrom = tokens[0];
            let pos: usize = tokens[2].parse().unwrap();
            let tmp1: HashSet<_> = sv_locus_identifiers(chrom, pos, 5000).into_iter().collect();

            let chrom = tokens[5];
            let pos: usize = tokens[7].parse().unwrap();
            let tmp2: HashSet<_> = sv_locus_identifiers(chrom, pos, 5000).into_iter().collect();

            let tmp: HashSet<String> = tmp1.union(&tmp2).cloned().collect();
            sample_variants.push(tmp);
        }
    }

     let mut blacklist: HashSet<String> = HashSet::new();
     for loci in sample_variants {
        for l in loci {
         blacklist.insert(l);
        }
     }
   
     blacklist = natural_sorted(blacklist);
    
    let mut frequency: Vec<usize> = vec![0; blacklist.len()];
    for (k, bad_variant) in enumerate(blacklist) {
        let bad_in_sample: Vec<usize> = Vec::new();
        }
    }

    println!("{:?}", frequency);
    println!("{:?}", &blacklist.len());

}
