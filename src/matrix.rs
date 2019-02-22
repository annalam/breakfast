
use common::{parse_args, FileReader};
use bitvec::*;
use rust_htslib::bam;
use rust_htslib::bam::Read;
use bio::alphabets::dna;

const USAGE: &str = "
Usage:
  breakfast matrix [options] <sv_file> <bam_files>...

Options:
  --all-reads   Count all reads, not just unaligned ones.
";

// Each signature is 20+20 bp, covering both sides of the breakpoint,
// for a total of 40 bp.
#[derive(Debug)]
struct Rearrangement {
	signature: String,
	signature_revcomp: String,
	first_8_cols: String,
	//chromosome_left: String,
	//position_left: usize,
	//strand_left: char,
	//chromosome_right: String,
	//position_right: usize,
	//strand_right: char,
	evidence: Vec<usize>
}

fn reverse_complement(seq: &str) -> String {
	String::from_utf8(dna::revcomp(seq.as_bytes())).unwrap()
}

//fn parse_strand(text: &str) -> char {
//	match text {
//		"+" => '+', "-" => '-',
//		_ => error!("Invalid strand '{}' found.", text)
//	}
//}

// This function takes the current hash as input, and modifies it according
// to the read sequence's next encoded base in the BAM file. The hash itself
// is the lower 16 bits of the u32. The upper 16 bits are used as error
// bits to indicate whether an ambiguous nucleotide was encountered.
// This way if the hash is larger than 65535, we know that the eight
// nucleotides used in the hash's calculation included some ambiguous ones.
fn hash_nucleotide(hash: u32, nuc: u8) -> u32 {
	((hash & 0b11111111_11111111_00111111_11111111u32) << 2) + match nuc {
		b'A' => 0u32, b'C' => 1u32, b'G' => 2u32, b'T' => 3u32,     // ACGT
		_ => 0b00000000_00000011_00000000_00000000u32   // Ambiguous nucleotide
	}
}

fn hash_8bp_sequence(seq: &str) -> u32 {
	assert!(seq.len() == 8);
	let mut hash: u32 = 0;
	for nuc in seq.bytes() { hash = hash_nucleotide(hash, nuc); }
	if hash & 0xFFFF0000u32 > 0 { error!("Invalid sequence '{}'.", seq); }
	hash
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

	// Convert BAM paths to sample names
	let mut samples: Vec<String> = Vec::new();
	for s in 0..bam_paths.len() {
		let start = if let Some(slash) =
			bam_paths[s].find('/') { slash + 1 } else { 0 };
		let end = if bam_paths[s].ends_with(".bam") {
			bam_paths[s].len() - 4 } else { bam_paths[s].len() };
		samples.push(bam_paths[s][start..end].into());
	}

	let mut line = String::new();
	let mut rearrangements: Vec<Rearrangement> = Vec::new();

	// Read all rearrangement signatures into memory
	let mut skipped_ambiguous = 0;
	let mut sv_file = FileReader::new(&sv_path);
	while sv_file.read_line(&mut line) {
		if line.starts_with("CHROM\t") { continue; }
		let mut signatures: Vec<String> = Vec::new();
		let cols: Vec<&str> = line.split('\t').collect();
		if cols.len() < 9 { continue; }
		let reads: Vec<&str> = cols[8].split(';').collect();
		for read in reads {
			let pipe = read.find('|').unwrap();
			if pipe < 20 { continue; }
			signatures.push(format!("{}{}",
				&read[pipe-20..pipe], &read[pipe+1..pipe+21]));
		}

		let mut signature = most_frequent(&signatures);
		signature.make_ascii_uppercase();
		if signature.chars().any(
			|b| b != 'A' && b != 'C' && b != 'G' && b != 'T') {
			skipped_ambiguous += 1;
			continue;
		}
		let signature_revcomp = reverse_complement(&signature);

		let first_8_cols = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
			cols[0], cols[1], cols[2], cols[3], cols[4], cols[5],
			cols[6], cols[7]);

		rearrangements.push(Rearrangement {
			signature, signature_revcomp, first_8_cols, 
			//chromosome_left: cols[0].into(),
			//strand_left: parse_strand(cols[1]),
			//position_left: cols[2].parse().unwrap(),
			//chromosome_right: cols[4].into(),
			//strand_right: parse_strand(cols[5]),
			//position_right: cols[6].parse().unwrap(),
			evidence: vec![0; bam_paths.len()]
		});
	}
	eprintln!("WARNING: Skipped {} rearrangements with signatures containing ambiguous nucleotides.", skipped_ambiguous);

	rearrangements.sort_unstable_by(|a, b| a.signature.cmp(&b.signature));
	for k in 1..rearrangements.len() {
		if rearrangements[k - 1].signature == rearrangements[k].signature &&
			rearrangements[k - 1].first_8_cols !=
			rearrangements[k].first_8_cols {
			eprintln!("WARNING: Found two rearrangements with signature {}:\n{}\n{}\n",
				rearrangements[k].signature,
				rearrangements[k - 1].first_8_cols,
				rearrangements[k].first_8_cols);
		}
	}
	rearrangements.dedup_by(|a, b| a.signature == b.signature);

	// Arrange junction signatures into a 65536-element table that is indexed
	// with the first 8 bp of the junction signature. This allows extremely
	// fast lookups.
	eprintln!("Building signature map for {} rearrangements...",
		rearrangements.len());
	let mut signature_exists = bitvec![0; 65536];
	let mut signature_map: Vec<Vec<u32>> =
		(0..65536).map(|_| Vec::new()).collect();
	for r in 0..rearrangements.len() {
		// Add the signatures for the two strands into the signature map.
		let hash = hash_8bp_sequence(&rearrangements[r].signature[16..24]);
		signature_exists.set(hash as usize, true);
		signature_map[hash as usize].push(r as u32);

		let hash = hash_8bp_sequence(
			&rearrangements[r].signature_revcomp[16..24]);
		signature_exists.set(hash as usize, true);
		signature_map[hash as usize].push(r as u32);
	}

	for s in 0..bam_paths.len() {
		eprintln!("Analyzing {}...", samples[s]);
		let mut bam = bam::Reader::from_path(&bam_paths[s]).unwrap_or_else(
			|_| error!("Could not open BAM file."));

		// TODO: Make faster by reusing the same BAM record.
		for r in bam.records() {
			let read = r.unwrap();
			if read.is_unmapped() == false { continue; }
			let seq = String::from_utf8(read.seq().as_bytes()).unwrap();

			// Start with some error bits set, so we only start checking
			// against the signature map once we have hashed at least eight
			// nucleotides.
			let mut hash = 0b00000000_00000011_00000000_00000000u32;
			'outer: for base in seq.bytes() {
				hash = hash_nucleotide(hash, base);
				if hash & 0xFFFF0000u32 > 0 { continue; }
				if signature_exists[hash as usize] == false { continue; }
				for ridx in &signature_map[hash as usize] {
					// This read contains the 4+4 bp junction signature.
					// Now check if the 20+20 bp junction is also found.
					let rearrangement = &mut rearrangements[*ridx as usize];
					if seq.contains(&rearrangement.signature) ||
						seq.contains(&rearrangement.signature_revcomp) {
						rearrangement.evidence[s] += 1;
						//eprintln!("Read {} supports rearrangement with signature {}.", seq, rearrangement.signature);
					}
					break 'outer;
				}
			}
		}

		//let mut supporting_reads = 0;
		//for r in &rearrangements { supporting_reads += r.evidence[s]; }
		//eprintln!("Found {} supporting reads for rearrangements.", supporting_reads);
	}

	print!("CHROM\tSTRAND\tPOSITION\tNEARBY FEATURES\t");
	print!("CHROM\tSTRAND\tPOSITION\tNEARBY FEATURES\t");
	print!("SUPPORTING READS\tSIGNATURE\tNOTES");
	for sample in &samples { print!("\t{}", sample); }
	println!();
	for r in rearrangements {
		print!("{}", r.first_8_cols);
		print!("\t\t{}|{}\t", &r.signature[0..20], &r.signature[21..]);
		for s in 0..samples.len() { print!("\t{}", r.evidence[s]); }
		println!();
	}
}

