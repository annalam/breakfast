
use common::{parse_args, FileReader};
use std::collections::HashMap;
use bitvec::*;
use rust_htslib::bam;
use rust_htslib::bam::Read;


const USAGE: &str = "
Usage:
  breakfast matrix [options] <sv_file> <bam_files>...

Options:
  --all-reads   Count all reads, not just unaligned ones.
";

// Each signature is 16 bp on both sides of breakpoint, for a total of
// 32 bp. With two bits per nucleotide, that amounts to 64 bits.
struct Rearrangement {
	signature: String,
	chromosome_left: String,
	position_left: usize,
	strand_left: char,
	chromosome_right: String,
	position_right: usize,
	strand_right: char,
	evidence: Vec<usize>
}

fn parse_strand(text: &str) -> char {
	match text {
		"+" => '+', "-" => '-',
		_ => error!("Invalid strand '{}' found.", text)
	}
}

fn ascii_nucleotide_to_0123(ascii: u8) -> usize {
	match ascii {
		b'A' => 0, b'C' => 1, b'G' => 2, b'T' => 3,
		_ => error!("Invalid base '{}'.", ascii)
	}
}


fn hash_bam_encoded_nucleotide(hash: u32, encoded: u8) -> u32 {
	((hash & 0xFFFF0FFF) << 2) + match encoded {
		1 => 0, 2 => 1, 4 => 2, 8 => 3,   // ACGT
		_ => 0x000F0000
	}
}

// Find the most frequent element in an unsorted vector. Operates in
// O(n log n) time.
fn most_frequent(elems: &Vec<String>) -> String {
	let mut sorted = elems.clone();
	sorted.sort_unstable();
	let mut most_frequent: usize = 0;
	let mut most_frequent_count: usize = 1;
	let mut curr_count: usize = 1;
	for k in 1..sorted.len() {
		if sorted[k] == sorted[k - 1] {
			curr_count += 1;
		} else {
			if curr_count > most_frequent_count {
				most_frequent = k - 1;
				most_frequent_count = curr_count;
			}
			curr_count = 1;
		}
	}
	if curr_count > most_frequent_count {
		most_frequent = sorted.len() - 1;
		most_frequent_count = curr_count;
	}
	sorted[most_frequent].clone()
}

pub fn main() {
	let args = parse_args(USAGE);
	let sv_path = args.get_str("<sv_file>");
	let bam_paths = args.get_vec("<bam_files>");

	let mut line = String::new();
	let mut rearrangements: Vec<Rearrangement> = Vec::new();

	// Read all rearrangement signatures into memory
	let mut sv_file = FileReader::new(&sv_path);
	while sv_file.read_line(&mut line) {
		if line.starts_with("CHROM\t") { continue; }
		let mut signatures: Vec<String> = Vec::new();
		let cols: Vec<&str> = line.split('\t').collect();
		let reads: Vec<&str> = cols[8].split(';').collect();
		for read in reads {
			let pipe = read.find('|').unwrap();
			if pipe < 20 { continue; }
			signatures.push(format!("{}{}",
				&read[pipe-20..pipe], &read[pipe+1..pipe+21]));
		}

		let mut signature = most_frequent(&signatures);
		signature.make_ascii_uppercase();
		if signature.contains('N') { continue; }

		rearrangements.push(Rearrangement {
			signature: signature,
			chromosome_left: cols[0].into(),
			strand_left: parse_strand(cols[1]),
			position_left: cols[2].parse().unwrap(),
			chromosome_right: cols[4].into(),
			strand_right: parse_strand(cols[5]),
			position_right: cols[6].parse().unwrap(),
			evidence: Vec::new()
		});
	}

	rearrangements.sort_unstable_by(|a, b| a.signature.cmp(&b.signature));
	/*for k in 1..rearrangements.len() {
		if rearrangements[k - 1].signature == rearrangements[k].signature {
			eprintln!("Found two rearrangements with identical signature.");
			break;
		}
	}*/
	rearrangements.dedup_by(|a, b| a.signature == b.signature);

	// Arrange junction signatures into a 65536-element table that is indexed
	// with the first 8 bp of the junction signature. This allows extremely
	// fast lookups.
	eprintln!("Building signature map...");
	let mut signature_exists = bitvec![0; 65536];
	let mut signature_map: Vec<Vec<u32>> =
		(0..65536).map(|x| Vec::new()).collect();
	for r in 0..rearrangements.len() {
		let mut hash: usize = 0;
		for base in rearrangements[r].signature[16..24].bytes() {
			hash *= 4; hash += ascii_nucleotide_to_0123(base);
		}
		signature_exists.set(hash, true);
		signature_map[hash].push(r as u32);
	}

	for bam_path in bam_paths {
		eprintln!("Analyzing {}...", bam_path);
		let mut bam = bam::Reader::from_path(&bam_path).unwrap_or_else(
			|_| error!("Could not open BAM file."));
		continue;

		// TODO: Make faster by reusing the same BAM record.
		for r in bam.records() {
			let read = r.unwrap();
			if read.is_unmapped() == false { continue; }
			let seq = read.seq();
			if seq.len() < 8 { continue; }
			let mut hash = 0u32;
			for k in 0..seq.len() {
				//bam_encoded_nucleotide_to_0123(seq.encoded_base(k))
			}
		}
	}

}

