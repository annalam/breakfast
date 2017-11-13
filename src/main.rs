
extern crate docopt;
extern crate bio;
extern crate rust_htslib;
extern crate regex;

use std::env;
use std::process::exit;
use docopt::{Docopt, ArgvMap};

mod detect; mod filter; mod annotate; mod blacklist;

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
    else if args.len() >= 2 && args[1] == "filter" { filter::main(); }
	else if args.len() >= 2 && args[1] == "annotate" { annotate::main(); }
	else if args.len() >= 2 && args[1] == "blacklist" { blacklist::main(); }
	else { println!("{}", USAGE); exit(-1); }
}

pub fn parse_args(usage: &str) -> ArgvMap {
	Docopt::new(usage).unwrap().parse().
		on_error(&format!("Invalid arguments.\n{}", usage))
}

// Helper methods for error reporting
trait ErrorHelper<T> {
	fn on_error(self, msg: &str) -> T;
}

impl<T> ErrorHelper<T> for Option<T> {
	fn on_error(self, msg: &str) -> T {
		match self {
			Some(x) => x,
			None => { eprintln!("ERROR: {}\n", msg); exit(-1) }
		}
	}
}

impl<T, E> ErrorHelper<T> for Result<T, E> {
	fn on_error(self, msg: &str) -> T {
		match self {
			Ok(x) => x,
			Err(_) => { eprintln!("ERROR: {}\n", msg); exit(-1) }
		}
	}
}
