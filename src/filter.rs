
use crate::common::{parse_args, FileReader};
use std::collections::HashSet;

const USAGE: &str = "
Usage:
  breakfast filter [options] <sv_path>

Options:
  --min-reads=N       Minimum number of supporting reads [default: 0]
  --min-distance=N    Minimum distance between breakpoints [default: 0]
  --blacklist=PATH    File containing blacklisted breakpoint signatures
  --merge-duplicates  Merge supporting reads with identical sequences
";

pub fn main() {
	let args = parse_args(USAGE);
	let sv_path = args.get_str("<sv_path>");
	let min_reads: usize = args.get_str("--min-reads").parse().unwrap();
	let min_distance: usize = args.get_str("--min-distance").parse().unwrap();
	let blacklist_path = args.get_str("--blacklist");
	let merge_duplicates = args.get_bool("--merge-duplicates");

	let mut line = String::new();

	let mut blacklist = HashSet::new();
	if !blacklist_path.is_empty() {
		let mut bl = FileReader::new(&blacklist_path);
		while bl.read_line(&mut line) {
			blacklist.insert(line.trim().to_string());
		}
	}

	let mut sv_file = FileReader::new(&sv_path);

	// Print the header
	sv_file.read_line(&mut line);
	print!("{}", line);

	while sv_file.read_line(&mut line) {
		let cols: Vec<&str> = line.split('\t').collect();
		assert!(cols.len() >= 10);

		let mut reads: Vec<&str> = cols[8].split(';').collect();
		if merge_duplicates {
			// Deduplicate based on last 20 bases of read
			reads.sort_by_key(|r| { let len = r.len(); &r[len-20..len] });
			reads.dedup_by_key(|r| { let len = r.len(); &r[len-20..len] });

			// Deduplicate based on first 20 bases of read
			reads.sort();
			reads.dedup_by_key(|r| &r[0..20]);
		}
		if reads.len() < min_reads { continue; }

		let signature = cols[9].trim();
		if blacklist.contains(signature) { continue; }

		// Filter rearrangements based on breakpoint distance
		if min_distance > 0 {
			let pos1: i64 = cols[2].parse().unwrap();
			let pos2: i64 = cols[6].parse().unwrap();
			if cols[0] == cols[4] && ((pos1 - pos2).abs() as usize) < min_distance {
				continue;
			}
		}

		print!("{}", cols[0]);
		for k in 1..8 { print!("\t{}", cols[k]); }
		print!("\t{}", reads.join(";"));
		for k in 9..cols.len() { print!("\t{}", cols[k]); }
		println!();

		// OLD CODE WITHOUT ALLOCATIONS
		//let num_reads = line.split('\t').nth(8).unwrap().split(';').count();
		//if num_reads < min_reads { continue; }

		//let signature = line.split('\t').nth(9).unwrap().trim();
		//if blacklist.contains(signature) { continue; }

		//println!("{}", line);
	}
}
