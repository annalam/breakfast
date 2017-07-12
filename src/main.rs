
extern crate docopt;
extern crate bio;
extern crate rust_htslib;

use std::env;
use std::process::exit;
use docopt::{Docopt, ArgvMap};

mod detect;

const USAGE: &'static str = "
Breakfast is a software for detecting chromosomal rearrangements in DNA/RNA
sequencing data.

Usage:
  breakfast <subcommand>

Available subcommands:
  detect      Detect chromosomal rearrangements.
  filter      Filter rearrangements based on quality of evidence.
  blacklist   Construct a rearrangement blacklist based on various criteria.
";

fn main() {
	// TODO: Use match args.as_slice() { [_, "detect"] => ... } after slice
	// pattern matching stabilizes (see issue #23121).
	let args: Vec<String> = env::args().collect();

	if args.len() >= 2 && args[1] == "detect" { detect::main(); }
	//else if args.len() >= 2 && args[1] == "filter" { filter::main(); }
	//else if args.len() >= 2 && args[1] == "blacklist" { blacklist::main(); }
	else { println!("{}", USAGE); exit(-1); }
}

pub fn parse_args(usage: &str) -> ArgvMap {
	Docopt::new(usage).unwrap().parse()
		.unwrap_or_else(|_| {
			println!("Invalid arguments.\n{}", usage); exit(-1); })
}