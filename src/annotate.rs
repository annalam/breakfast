
use parse_args;
use std::io::{BufRead, BufReader, stdin};
use std::fs::File;

const USAGE: &str = "
Usage:
  breakfast annotate [options] <sv_path> <bed_path>
";

fn distance(pos: u32, feature: &Feature) -> u32 {
	if pos < feature.start {
		feature.start - pos
	} else if pos > feature.end {
		pos - feature.end
	} else {
		0
	}
}

struct Feature {
	chr: String,
	start: u32,   // 1-based position of first base
	end: u32,     // 1-based position of last base
	name: String
}

pub fn main() {
	let args = parse_args(USAGE);
	let sv_path = args.get_str("<sv_path>");
	let bed_path = args.get_str("<bed_path>");

    let mut sv: Box<BufRead>;
    if sv_path == "-" {   // TODO: Make a function that handles this.
    	sv = Box::new(BufReader::new(stdin()));
    } else {
    	sv = Box::new(BufReader::new(File::open(&sv_path).unwrap()));
    }

	let mut header = String::new();
	sv.read_line(&mut header).unwrap();
	print!("{}", header);

    let bed = BufReader::new(File::open(&bed_path).unwrap());
    let mut features: Vec<Feature> = Vec::new();
    for l in bed.lines() {
        let line = l.unwrap();
        let mut cols = line.split('\t');
        features.push(Feature {
        	chr: cols.next().unwrap().to_string(),
        	start: cols.next().unwrap().parse::<u32>().unwrap() + 1,
        	end: cols.next().unwrap().parse().unwrap(),
        	name: cols.next().unwrap().to_string()
        });
    }

    for l in sv.lines() {
		let line = l.unwrap();
		if !line.starts_with("chr") { continue; }

		let cols: Vec<&str> = line.split('\t').collect();
        let chr_1 = cols[0];
        let pos_1: u32 = cols[2].parse().unwrap();
        let chr_2 = cols[4];
        let pos_2: u32 = cols[6].parse().unwrap();
        let reads = cols[8];

        let mut nearby_features_1: Vec<(u32, &Feature)> = Vec::new();
        for feature in &features {
        	if chr_1 != feature.chr { continue; }
        	let dist = distance(pos_1, feature);
        	if dist > 100_000 { continue; }
        	nearby_features_1.push((dist, feature));
        }
        nearby_features_1.sort_by_key(|x| x.0);

        let mut nearby_features_2: Vec<(u32, &Feature)> = Vec::new();
        for feature in &features {
        	if chr_2 != feature.chr { continue; }
        	let dist = distance(pos_2, feature);
        	if dist > 100_000 { continue; }
        	nearby_features_2.push((dist, feature));
        }
        nearby_features_2.sort_by_key(|x| x.0);

        print!("{}\t{}\t{}\t", chr_1, cols[1], pos_1);
        for (k, nf) in nearby_features_1.iter().enumerate() {
        	print!("{} ({})", nf.1.name, nf.0);
        	if k < nearby_features_1.len() - 1 { print!(", "); }
        }
        print!("\t");
        print!("{}\t{}\t{}\t", chr_2, cols[5], pos_2);
        for (k, nf) in nearby_features_2.iter().enumerate() {
        	print!("{} ({})", nf.1.name, nf.0);
        	if k < nearby_features_2.len() - 1 { print!(", "); }
        }
        println!("\t{}\t{}\t{}", reads, cols[9], cols[10]);
	}
}
