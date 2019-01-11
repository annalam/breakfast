
extern crate docopt;
extern crate bio;
extern crate rust_htslib;
extern crate regex;
extern crate bitvec;

use std::env;

#[macro_use] mod common;
mod detect; mod filter; mod annotate; mod blacklist; mod matrix;

const USAGE: &str = "
Breakfast is a software for detecting chromosomal rearrangements in DNA/RNA
sequencing data.

Usage:
  breakfast <subcommand>

Available subcommands:
  detect      Detect chromosomal rearrangements.
  filter      Filter rearrangements based on quality of evidence.
  blacklist   Construct a rearrangement blacklist based on various criteria.
  annotate    Annotate genes adjacent to rearrangement breakpoints.
  matrix      Build a read count matrix for rearrangements.
";

fn main() {
	// TODO: Use match args.as_slice() { [_, "detect"] => ... } after slice
	// pattern matching stabilizes (see issue #23121).
	let args: Vec<String> = env::args().collect();

	if args.len() >= 2 && args[1] == "detect" { detect::main(); }
    else if args.len() >= 2 && args[1] == "filter" { filter::main(); }
	else if args.len() >= 2 && args[1] == "annotate" { annotate::main(); }
	else if args.len() >= 2 && args[1] == "blacklist" { blacklist::main(); }
	else if args.len() >= 2 && args[1] == "matrix" { matrix::main(); }
	else if args.len() == 1 { eprintln!("{}", USAGE); }
	else { error!("Invalid subcommand.\n\n{}", USAGE); }
}
