use std::io::{BufRead, BufReader};
use std::fs::File;
use std::collections::HashSet;

use filter::sv_locus_identifiers as sv_locus_identifiers;

pub fn natural_sorted() {
    unimplemented!();
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
            //println!("inner test {:?}", &tmp);
            sample_variants.push(tmp);

        }

        for loci in &sample_variants {
            println!("total variants in loci {:?}", loci);
        }


    }

}
