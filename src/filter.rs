
use parse_args;
use ErrorHelper;
use std::io::{BufRead, BufReader};
use std::collections::HashSet;
use std::fs::File;

const USAGE: &'static str = "
Usage:
  breakfast filter [options] <sv_path>

Options:
  --min-reads=N     Minimum number of supporting reads [default: 0]
  --blacklist=PATH  Path to blacklist file
";

pub fn main() {
	let args = parse_args(USAGE);
	let sv_path = args.get_str("<sv_path>");
	let min_reads: usize = args.get_str("--min-reads").parse().unwrap();
	let blacklist_path = args.get_str("--blacklist");

	let mut blacklist = HashSet::new();
	if !blacklist_path.is_empty() {
		let bl = BufReader::new(File::open(&blacklist_path).on_error(
			&format!("Could not open blacklist file '{}'.", blacklist_path)));
		for line in bl.lines() { blacklist.insert(line.unwrap()); }
	}

	let mut sv_file = BufReader::new(File::open(&sv_path)
		.on_error("Could not open .sv file."));
	let mut header = String::new();
	sv_file.read_line(&mut header);
	print!("{}", header);

	for l in sv_file.lines() {
		let line: String = l.unwrap();
		let tokens: Vec<&str> = line.split('\t').collect();
		if tokens.len() < 9 { continue; }

		let num_reads = tokens[8].split(';').count();
		if num_reads < min_reads { continue; }

        let chrom = tokens[0];
        let pos: usize = tokens[2].parse().unwrap();
        let loci_1: HashSet<_> = sv_locus_identifiers(chrom, pos, 5000).into_iter().collect();

        let chrom = tokens[4];
        let pos = tokens[6].parse().unwrap();
        let loci_2: HashSet<_> = sv_locus_identifiers(chrom, pos, 5000).into_iter().collect();

        if loci_1.is_disjoint(&blacklist) || loci_2.is_disjoint(&blacklist) {
            println!("{}", line);
        }
	}
}


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
