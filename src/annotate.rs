//src/annotate.rs

use std::io::{BufRead, BufReader};
use std::fs::File;
use regex::Regex;
use std::collections::HashMap;

pub fn distance_to_gene(sv_pos: usize, gene_pos: usize) ->usize {
    let mut dist: Vec<usize> = Vec::new();
    dist.push(0);
    dist.push(gene_pos - sv_pos);
    dist.push(sv_pos - gene_pos);
    dist.sort_by(|a, b| b.cmp(a));
    dist[1]
}


pub fn annotate(sv_path: String, bed_path: String) {
    let sv_file_header: String = "CHROM\tSTRAND\tPOSITION\tNEARBY_FEATURES\t\tCHROM\tSTRAND\tPOSITION\tNEARBY_FEATURES\t\tNUM_SPANNING_FRAGMENTS\tNUM_SPANNING_MATES\tSPANNING_MATE_SEQUENCES".to_uppercase();
    let mut sv  = BufReader::new(File::open(&sv_path).unwrap());
    let mut bed = BufReader::new(File::open(&bed_path).unwrap());

    let mut features: Vec<String> = Vec::new();

    for b in bed.lines() {
        let line = b.unwrap();
        let cols: Vec<&str> = line.split('\t').collect();
        let mut feat: Vec<&str> = Vec::new();
        feat.push(cols[0]); feat.push(cols[5]);
        feat.push(cols[1]);feat.push(cols[2]); feat.push(cols[3]);
        features.push(feat.join("\t"));
    }

    let re = Regex::new(r" \(ENSG.*?\)").unwrap();
    println!("{}", sv_file_header);
    for l in sv.lines() {
		let line: String = l.unwrap();
		if !line.starts_with("chr") { continue; }

		let mut tokens: Vec<&str> = line.split('\t').collect();
        let chr_1    = tokens[0];
        let strand_1 = tokens[1];
        let pos_1    = tokens[2].parse::<usize>().unwrap();

        let chr_2    = tokens[5];
        let strand_2 = tokens[6];
        let pos_2    = tokens[7].parse::<usize>().unwrap();

        let mut nearby_features_1 = HashMap::new();
        let mut nearby_features_2 = HashMap::new();

        for f in &features {
            let fe: Vec<&str> = f.split('\t').collect();
            if fe[0] == chr_1 {
                nearby_features_1.insert(re.replace(fe[4], "").to_string(), distance_to_gene(pos_1, fe[2].parse::<usize>().unwrap()));
            }
        }

		
        for f in &features {
            let fe: Vec<&str> = f.split('\t').collect();
            if fe[0] == chr_2 {
                nearby_features_2.insert(re.replace(fe[4], "").to_string(), distance_to_gene(pos_1, fe[2].parse::<usize>().unwrap()));
            }
        }

        nearby_features_1.retain(|ref x, &mut v| v < 100_000);
        nearby_features_2.retain(|ref x, &mut v| v < 100_000);
        
        let mut sorted_1: Vec<String> = nearby_features_1.keys().cloned().collect();
        		sorted_1.sort();
		let mut sorted_2: Vec<String> = nearby_features_2.keys().cloned().collect();
        		sorted_2.sort();
		
		let mut nb_feat1 = HashMap::new();
		let mut nb_feat2 = HashMap::new();
		
		for key in &sorted_1 {
			if nearby_features_1.contains_key(key) {
				nb_feat1.insert(key, nearby_features_1.get(key).unwrap());
			} else { continue; }
		}
		
		for key in sorted_2.iter() {
			if nearby_features_2.contains_key(key) {
				nb_feat2.insert(key, nearby_features_2.get(key).unwrap());
			} else { continue; }
		}
		
		
		let mut nb1 = String::new();
		let mut nb2 = String::new();
		
	}
}