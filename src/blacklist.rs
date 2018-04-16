
use parse_args;
use ErrorHelper;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::collections::HashMap;

const USAGE: &str = "
Usage:
  breakfast blacklist [options] <sv_files>...

Options:
  --min-samples=N    Blacklist if present in N or more samples [default: 1]
";

pub fn main() {
    let args = parse_args(USAGE);
    let sv_paths = args.get_vec("<sv_files>").to_vec();
    let min_samples: usize = args.get_str("--min-samples").parse()
        .on_error("--min-samples must be numeric");

    // Mapping: rearrangement signature -> vector of evidence across samples
    let mut rearrangements: HashMap<String, Vec<bool>> = HashMap::new();

    for (s, sv_path) in sv_paths.iter().enumerate() {
        let mut sv_file = BufReader::new(File::open(&sv_path).unwrap());
        let mut header = String::new();
        sv_file.read_line(&mut header).unwrap();   // Skip the header
        for l in sv_file.lines() {
            let line = l.unwrap();
            let signature = line.split('\t').nth(9).unwrap().to_string();
            let entry = rearrangements.entry(signature)
            	.or_insert_with(|| vec![false; sv_paths.len()]);
            entry[s] = true;
        }
    }

    for (signature, evidence) in rearrangements.iter() {
    	if evidence.iter().filter(|x| **x == true).count() >= min_samples {
    		println!("{}", signature);
    	}
    }
}
