
use std::env;

#[macro_use] mod common;
mod detect; mod filter; mod annotate;

const USAGE: &str = "
Breakfast is a software for detecting chromosomal rearrangements in DNA/RNA
sequencing data.

Usage:
  breakfast <subcommand>

Available subcommands:
  detect      Detect chromosomal rearrangements.
  filter      Filter rearrangements based on quality of evidence.
  annotate    Annotate genes adjacent to rearrangement breakpoints.
";

fn main() {
	// TODO: Use match args.as_slice() { [_, "detect"] => ... } after slice
	// pattern matching stabilizes (see issue #23121).
	let args: Vec<String> = env::args().collect();

	if args.len() >= 2 && args[1] == "detect" { detect::main(); }
    else if args.len() >= 2 && args[1] == "filter" { filter::main(); }
	else if args.len() >= 2 && args[1] == "annotate" { annotate::main(); }
	//else if args.len() >= 2 && args[1] == "vaf" { vaf::main(); }
	else if args.len() == 1 { eprintln!("{}", USAGE); }
	else { error!("Invalid subcommand.\n\n{}", USAGE); }
}
