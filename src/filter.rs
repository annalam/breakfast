
use common::{parse_args, FileReader};
use std::collections::HashSet;

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

	let mut line = String::new();

	let mut blacklist = HashSet::new();
	if !blacklist_path.is_empty() {
		let mut bl = FileReader::new(&blacklist_path);
		while bl.read_line(&mut line) { blacklist.insert(line.clone()); }
	}

	let mut sv_file = FileReader::new(&sv_path);

	// Print the header
	sv_file.read_line(&mut line);
	print!("{}", line);

	while sv_file.read_line(&mut line) {
		let num_reads = line.split('\t').nth(8).unwrap().split(';').count();
		if num_reads < min_reads { continue; }

		let signature = line.split('\t').nth(9).unwrap();
		if blacklist.contains(signature) { continue; }

		println!("{}", line);
	}
}
