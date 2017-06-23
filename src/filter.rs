use std::io::{BufRead, BufReader};
use std::collections::HashSet;
use std::fs::File;

pub fn sv_locus_identifiers(chr: &str, pos: usize, resolution: i32) -> Vec<String> {
    let bin: i32 = pos as i32 / resolution;
	let mut bins: Vec<i32> = Vec::new();
    for i in bin-1..bin+2 as i32 {
    	bins.push(i * resolution);
    }

	let mut out: Vec<String> = Vec::new();
    for pos in bins {
    	out.push(format!("{}:{}", &chr, &pos.to_string()));
    }
	out
}

pub fn filter(sv_path: String, min_reads: usize, blacklist_path: String) {
	println!("{}\t{}", sv_path, 15);

	let mut sv = BufReader::new(File::open(&sv_path).unwrap());
	let mut bl =  BufReader::new(File::open(&blacklist_path).unwrap());
    let mut blacklist = HashSet::new();

	for b in bl.lines() {
       let line  = b.unwrap();
       blacklist.insert(line);
    }

	for l in sv.lines() {
		let line: String = l.unwrap();
		if !line.starts_with("chr") { continue; }

		let mut tokens: Vec<&str> = line.split('\t').collect();
		if tokens[10].parse::<i32>().unwrap() < min_reads as i32 { continue; }

        let mut chrom = tokens[0];
        let mut pos = tokens[2].parse::<usize>().unwrap();
        let mut loci_1: HashSet<_> = sv_locus_identifiers(chrom, pos, 5000).into_iter().collect();

                chrom = tokens[5];
                pos = tokens[7].parse::<usize>().unwrap();
        let mut loci_2: HashSet<_> = sv_locus_identifiers(chrom, pos, 5000).into_iter().collect();

        if loci_1.is_disjoint(&blacklist) || loci_2.is_disjoint(&blacklist) {
            println!("{}", line);
        }
	}
}
