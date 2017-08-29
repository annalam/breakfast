
use parse_args;
use std::io::{BufRead, BufReader, stdin};
use std::fs::File;

const USAGE: &'static str = "
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
	sv.read_line(&mut header);
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

		let mut cols = line.split('\t');
        let chr_1 = cols.next().unwrap();
        let strand_1 = cols.next().unwrap();
        let pos_1: u32 = cols.next().unwrap().parse().unwrap();
        cols.next();
        let chr_2 = cols.next().unwrap();
        let strand_2 = cols.next().unwrap();
        let pos_2: u32 = cols.next().unwrap().parse().unwrap();
        cols.next();
        let reads = cols.next().unwrap();

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

        print!("{}\t{}\t{}\t", chr_1, strand_1, pos_1);
        for nf in nearby_features_1 { print!("{} ({}), ", nf.1.name, nf.0); }
        print!("\t");
        print!("{}\t{}\t{}\t", chr_2, strand_2, pos_2);
        for nf in nearby_features_2 { print!("{} ({}), ", nf.1.name, nf.0); }
        println!("\t{}", reads);
	}
}
