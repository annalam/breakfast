
use crate::common::parse_args;
use std::mem::swap;
use std::{str, thread};
use std::process::{Command, Stdio};
use std::collections::HashMap;
use std::cmp::{min, max, Ordering};
use std::io::{BufReader, BufWriter, BufRead, Write};
use rust_htslib::bam;
use rust_htslib::bam::Read;
use bio::io::fasta;
use bio::alphabets::dna;

#[derive(Debug)]
struct Evidence {
	chr: String,
	pos: usize,               // Leftmost position of anchor #1 alignment
	strand: bool,
	mchr: String,
	mpos: usize,              // Leftmost position of anchor #2 alignment
	mstrand: bool,
	sequence: Vec<u8>,        // Full sequence of breakpoint overlapping read
	signature: Vec<u8>,       // Breakpoint signature (8 bp from both flanks)
	ref_signature_1: Vec<u8>, // Normal sequence around the 1st breakpoint
	ref_signature_2: Vec<u8>, // Normal sequence around the 2nd breakpoint
	frag_id: Vec<u8>          // Fragment QNAME from BAM file
}

const USAGE: &str = "
Usage:
  breakfast detect [options] <bam_file> <genome>

Options:
  --anchor-len=N       Anchor length for split read analysis [default: 30]
  --anchor-mm=N        Mismatches allowed in anchor alignments [default: 0]
  --max-frag-len=N     Maximum fragment length [default: 5000]
  --min-evidence=N     Minimum number of supporting DNA fragments [default: 2]
  --count-duplicates   Count also reads that have been flagged as duplicates
";

pub fn main() {
	let args = parse_args(USAGE);
	let sam_path = args.get_str("<bam_file>").to_string();
	let genome_path = args.get_str("<genome>");
	let anchor_len: usize = args.get_str("--anchor-len").parse().unwrap();
	//let anchor_mm: usize = args.get_str("--anchor-mm").parse().unwrap();
	let max_frag_len: usize = args.get_str("--max-frag-len").parse().unwrap();
	let min_evidence: usize = args.get_str("--min-evidence").parse().unwrap();
	let count_duplicates = args.get_bool("--count-duplicates");

	let fasta = fasta::Reader::from_file(format!("{}.fa", genome_path))
		.unwrap_or_else(|_| error!("Genome FASTA file {}.fa could not be read.", genome_path));
	eprintln!("Reading reference genome into memory...");

	let mut genome = HashMap::new();
	for entry in fasta.records() {
		let chr = entry.unwrap();
		genome.insert(chr.id().to_owned(), chr.seq().to_owned());
	}

	// TODO: Handle reads with multiple alignments...
	eprintln!("Splitting unaligned reads into {} bp anchors and aligning against the genome...", anchor_len);
	let bowtie = Command::new("bowtie")
		.args(&["-f", "-p1", "-v0", "-m1", "-B1", "--suppress", "5,6,7,8", &genome_path, "-"])
		.stdin(Stdio::piped()).stdout(Stdio::piped()).spawn().unwrap_or_else(
		|_| error!("Could not start Bowtie process."));

	let mut bowtie_in = BufWriter::new(bowtie.stdin.unwrap());
	let bowtie_out = BufReader::new(bowtie.stdout.unwrap());

	thread::spawn(move || {
		let mut bam = bam::Reader::from_path(&sam_path).unwrap_or_else(
			|_| error!("Could not open BAM file."));

		let mut num_reads_sent = 0;
		for r in bam.records() {
			let read = r.unwrap();
			if read.is_unmapped() == false { continue; }
			if read.is_duplicate() && count_duplicates == false { continue; }
			if read.seq().len() < anchor_len * 2 { continue; }

			// Any ':' characters in the fragment ID must be removed here
			// since ':' is used as a delimiter in our anchor descriptors.
			let mut frag_id = read.qname().to_vec();
			for k in 0..frag_id.len() {
				if frag_id[k] == b':' { frag_id[k] = b'_'; }
			}

			// Unaligned reads never need to be reverse-complemented
			let seq = read.seq().as_bytes();

			num_reads_sent += 1;

			// 5' anchor: >5p:READ#
			write!(bowtie_in, ">5p:{}\n", num_reads_sent).unwrap();
			bowtie_in.write_all(&seq[..anchor_len]).unwrap();

			// 3' anchor: >3p:READ#:FRAG_ID:FULL_SEQUENCE
			write!(bowtie_in, "\n>3p:{}:", num_reads_sent).unwrap();
			bowtie_in.write_all(&frag_id).unwrap();
			write!(bowtie_in, ":").unwrap();
			bowtie_in.write_all(&seq).unwrap();
			write!(bowtie_in, "\n").unwrap();
			bowtie_in.write_all(&seq[(seq.len() - anchor_len)..]).unwrap();
			writeln!(bowtie_in).unwrap();
		}
	});

	let mut evidence: Vec<Evidence> = Vec::new();
	let mut prev = String::new();
	let mut prev_read_num = 0;

	for l in bowtie_out.lines() {
		let line = l.unwrap();

		let read_num: usize = line.split(':').nth(1).unwrap().parse().unwrap();

		if line.starts_with("5p:") {
			prev = line;
			prev_read_num = read_num;
			continue;
		} else if line.starts_with("3p:") && read_num != prev_read_num {
			continue;
		}

		let mut anchor_info = line.split(':');
		let frag_id = anchor_info.nth(2).unwrap().as_bytes();
		let mut seq: Vec<u8> = anchor_info.next().unwrap().as_bytes().to_vec();
		let full_len: usize = seq.len();

		let mut cols = prev.split('\t');
		let mut mcols = line.split('\t');
		cols.next(); mcols.next();   // Skip first column

		let mut strand = cols.next().unwrap() == "+";
		let mut mstrand = mcols.next().unwrap() == "+";
		let mut chr = cols.next().unwrap();
		let mut mchr = mcols.next().unwrap();
		let mut pos: usize = cols.next().unwrap().parse().unwrap();
		let mut mpos: usize = mcols.next().unwrap().parse().unwrap();

		// Reorient the read so that anchor #1 has the lower coordinate.
		// This simplifies downstream analysis where we cluster the
		// rearrangement evidence by position.
		if chr > mchr || (chr == mchr && pos > mpos) {
			swap(&mut chr, &mut mchr);
			swap(&mut pos, &mut mpos);
			let tmp = strand; strand = !mstrand; mstrand = !tmp;
			seq = dna::revcomp(&seq);
		}

		// If the read is at the very edge of a chromosome, ignore it.
		if pos + full_len >= genome[chr].len() { continue; }
		if mpos + full_len >= genome[mchr].len() { continue; }

		let left_grch = if strand == true {
			genome[chr][pos-1..pos+full_len-1].to_vec()
		} else {
			dna::revcomp(&genome[chr][pos+anchor_len-full_len-1..pos+anchor_len-1].to_vec())
		};

		let right_grch = if mstrand == true {
			genome[mchr][mpos+anchor_len-full_len-1..mpos+anchor_len-1].to_vec()
		} else {
			dna::revcomp(&genome[mchr][mpos-1..mpos+full_len-1].to_vec())
		};

		// Identify the breakpoint location that minimizes the number of
		// nucleotide mismatches between the read and the breakpoint flanks.
		let mut mismatches: Vec<usize> = vec![0; full_len];
		for k in 0..anchor_len {
			mismatches[anchor_len] += (seq[k] != left_grch[k]) as usize;
		}
		for k in anchor_len..full_len {
			mismatches[anchor_len] += (seq[k] != right_grch[k]) as usize;
		}
		for k in anchor_len+1..full_len-anchor_len {
			mismatches[k] = mismatches[k-1] +
				(seq[k-1] != left_grch[k-1]) as usize -
				(seq[k-1] != right_grch[k-1]) as usize;
		}
		let mut bp = anchor_len;
		let mut least_mismatches = mismatches[anchor_len];
		for k in anchor_len+1..full_len-anchor_len {
			if mismatches[k] < least_mismatches {
				bp = k; least_mismatches = mismatches[k];
			}
		}

		let mut junction: Vec<u8> = vec![' ' as u8; full_len + 1];
		for k in 0..bp {
			junction[k] = if seq[k] == left_grch[k] { seq[k] } else { seq[k].to_ascii_lowercase() };
		}
		junction[bp] = '|' as u8;
		for k in bp..full_len {
			junction[k+1] = if seq[k] == right_grch[k] { seq[k] } else { seq[k].to_ascii_lowercase() };
		}

		// Construct a breakpoint signature, composed of 8 bp from both flanks
		let mut signature = vec![' ' as u8; 17];
		for k in 0..8 { signature[k] = left_grch[bp - 8 + k]; }
		signature[8] = '|' as u8;
		for k in 0..8 { signature[9 + k] = right_grch[bp + k]; }

		// Construct a reference signature for the first breakpoint.
		// TODO: Implement

		// Calculate the position of the first nucleotide immediately before
		// the breakpoint, on both sides of the junction.
		let left_bp_pos = if strand == true {
			pos + bp - 1
		} else {
			pos + anchor_len - bp
		};
		let right_bp_pos = if mstrand == true {
			mpos + anchor_len - (full_len - bp)
		} else {
			mpos + (full_len - bp)
		};

		evidence.push(Evidence {
			chr: chr.to_string(), pos: left_bp_pos, strand: strand,
			mchr: mchr.to_string(), mpos: right_bp_pos, mstrand: mstrand,
			sequence: junction, signature: signature,
			frag_id: frag_id.to_vec() });
	}

	eprintln!("Found {} rearrangement supporting reads.", evidence.len());

	eprintln!("Sorting rearrangement supporting reads by position...");
	evidence.sort_by(|a,b|
		if a.chr < b.chr { Ordering::Less }
		else if a.chr > b.chr { Ordering::Greater }
		else if a.pos < b.pos { Ordering::Less }
		else if a.pos > b.pos { Ordering::Greater }
		else { Ordering::Equal });

	eprintln!("Identifying rearrangements based on clusters of discordant reads...");
	println!("CHROM\tSTRAND\tPOSITION\tNEARBY FEATURES\tCHROM\tSTRAND\tPOSITION\tNEARBY FEATURES\tSUPPORTING READS\tSIGNATURE\tNOTES");
	let mut reported = vec![false; evidence.len()];
	for p in 0..evidence.len() {
		// We skip reads that were already incorporated into some cluster.
		if reported[p] { continue; }

		let read = &evidence[p];
		let mut cluster: Vec<&Evidence> = vec![read];
		for s in p+1..evidence.len() {
			// We try to add more reads into the cluster until we encounter
			// the first read that is so far that it cannot possibly belong
			// to the cluster. Then we terminate since we know that all
			// further reads are also too far away (since they are sorted).
			if evidence[s].chr != read.chr { break; }
			if evidence[s].pos - read.pos > max_frag_len { break; }

			// Before we add a read into the cluster, we check that both
			// anchors are consistent with other reads in the cluster.
			if evidence[s].mchr != read.mchr { continue; }
			if (evidence[s].mpos as i64 - read.mpos as i64).abs() > max_frag_len as i64 { continue; }
			if evidence[s].strand != read.strand { continue; }
			if evidence[s].mstrand != read.mstrand { continue; }
			if evidence[s].signature != read.signature { continue; }

			cluster.push(&evidence[s]);
			reported[s] = true;
		}

		/*if read.signature == "CAGAT|ACTTG".as_bytes() {
			for r in &cluster {
				println!("Template ID: {}\nFragment signature: {}\nSequence: {}\n", str::from_utf8(&r.frag_id).unwrap(), str::from_utf8(&r.frag_signature).unwrap(), str::from_utf8(&r.sequence).unwrap());
			}
		} else {
			continue;
		}*/

		if cluster.len() < min_evidence { continue; }
		cluster = remove_duplicates(cluster);
		if cluster.len() < min_evidence { continue; }

		print!("{}\t{}\t{}\t\t{}\t{}\t{}\t\t",
			read.chr, if read.strand { '+' } else { '-' }, read.pos,
			read.mchr, if read.mstrand { '+' } else { '-' }, read.mpos);
		for r in 0..cluster.len() {
			if r > 0 { print!(";"); }
			print!("{}", str::from_utf8(cluster[r].sequence.as_slice()).unwrap());
		}
		println!("\t{}\t", str::from_utf8(&read.signature).unwrap());
	}
}

fn remove_duplicates(evidence: Vec<&Evidence>) -> Vec<&Evidence> {
	let mut filtered: Vec<&Evidence> = Vec::new();
	let mut redundant_with = vec![-1i32; evidence.len()];
	for a in 0..evidence.len() {
		if redundant_with[a] >= 0 { continue; }
		redundant_with[a] = a as i32;
		let mut num_redundant = 1;
		for b in a+1..evidence.len() {
			if redundant_with[b] >= 0 { continue; }

			if evidence[a].frag_id == evidence[b].frag_id {
				redundant_with[b] = a as i32;
				num_redundant += 1;
			}
		}

		if num_redundant == 1 { filtered.push(evidence[a]); continue; }

		// Now we have a list of redundant reads. Since we can only keep
		// one of these reads as a source of evidence, we pick the "best"
		// one. Currently this means the one with the longest flanks.
		let mut best = 0;
		let mut best_score = 0;
		for k in 0..evidence.len() {
			if redundant_with[k] != a as i32 { continue; }
			let seq = &evidence[k].sequence;
			let left_len = seq.iter().position(|c| *c == '|' as u8).unwrap();
			let right_len = seq.len() - left_len - 1;

			// Score is primarily determined by the length of the shortest
			// flank, but if the shortest flank is equally long in both reads
			// then the length of the longest flank is considered also.
			let score = min(left_len, right_len) * 1000 + max(left_len, right_len);
			if score > best_score {
				best = k; best_score = score;
			}
		}
		filtered.push(evidence[best]);
	}
	filtered
}


