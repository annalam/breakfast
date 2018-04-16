
use parse_args;
use ErrorHelper;
use std::io::{BufRead, BufReader};
use std::collections::HashSet;
use std::fs::File;

const USAGE: &str = "
Usage:
  breakfast filter [options] <sv_path>

Options:
  --min-reads=N     Minimum number of supporting reads [default: 0]
  --blacklist=PATH  File containing blacklisted breakpoint signatures
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
		let line = l.unwrap();

		let num_reads = line.split('\t').nth(8).unwrap().split(';').count();
		if num_reads < min_reads { continue; }

		let signature = line.split('\t').nth(9).unwrap();
		if blacklist.contains(signature) { continue; }

		println!("{}", line);
	}
}
