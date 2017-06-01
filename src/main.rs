// rust version of breakfast tool
// BreakFast is a toolkit for detecting chromosomal rearrangements
// based on RNA-seq data.

#[macro_use]
extern crate clap;
extern crate bio;

mod cli;
// struct for having breakfast options
struct BfOptions {
    anchor_len :    i32,
    max_frag_len:   i32,
    min_mapq:       i32,
    orientation:    &'static str,
    all_reads:      bool,
    discard_duplicates: &'static str,
    min_reads:      &'static str,
    freq_above:     i32,
}

fn main() {
  let defaults = BfOptions{ anchor_len : 0,
                            max_frag_len: 5000,
                            min_mapq: 15,
                            orientation:  "fr",
                            all_reads:  false,
                            discard_duplicates: "both-ends",
                            min_reads:    "0-0-0",
                            freq_above: 0};

  let mut anchor_len: i32   = defaults.anchor_len;
  let mut max_frag_len: i32 = defaults.max_frag_len;
  let mut min_mapq: i32     = defaults.min_mapq;
  let mut orientation       = String::from(defaults.orientation);
  let mut all_reads: bool   = defaults.all_reads;
  let mut discard_duplicates = String::from(defaults.discard_duplicates);
  let mut min_reads     = String::from(defaults.min_reads);
  let mut freq_above: i32 = defaults.freq_above;

  //arguments from src/cli.rs
  let matches = cli::build_cli().get_matches();
  if matches.is_present("anchor-len") {
      anchor_len = matches.value_of("anchor-len").unwrap().parse::<i32>().unwrap();
  }

  if matches.is_present("max-frag-len") {
      max_frag_len = matches.value_of("max-frag-len").unwrap().parse::<i32>().unwrap();
  }

  if matches.is_present("min-mapq") {
      min_mapq = matches.value_of("min-mapq").unwrap().parse::<i32>().unwrap();
  }

  if matches.is_present("orientation") {
      orientation = String::from(matches.value_of("orientation").unwrap());
  }

  if matches.is_present("all-reads") {
      all_reads = true;
  }

  if matches.is_present("discard-duplicates") {
      discard_duplicates = String::from(matches.value_of("discard-duplicates").unwrap());
  }

  if matches.is_present("min-reads") {
      min_reads = String::from(matches.value_of("min-reads").unwrap());
  }

  if matches.is_present("freq-above") {
      freq_above = matches.value_of("freq-above").unwrap().parse::<i32>().unwrap();
  }

  println!("Hello, world!");
}
