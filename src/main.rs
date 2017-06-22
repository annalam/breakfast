// rust version of breakfast tool
// BreakFast is a toolkit for detecting chromosomal rearrangements
// based on RNA-seq data.

use std::mem::swap;
use std::thread;
use std::str;
use std::process::{Command, Stdio};
use std::collections::HashMap;
use std::ascii::AsciiExt;
use std::cmp::{min, Ordering};
use std::io::{BufReader, BufWriter, BufRead, Read, Write, stderr};
use rust_htslib::bam;
use rust_htslib::bam::Read as HTSRead;

extern crate clap;
extern crate bio;
extern crate rust_htslib;

use bio::alphabets::dna;

mod cli;

struct Evidence {
	chr: String,
	pos: usize,            // Leftmost position of anchor #1 alignment
	strand: bool,
	mchr: String,
	mpos: usize,           // Leftmost position of anchor #2 alignment
	mstrand: bool,
	sequence: Vec<u8>,     // Full sequence of breakpoint overlapping read
	frag_id: String,       // Identifier of the DNA fragment in the BAM file
	signature: Vec<u8>,    // Breakpoint signature (5 bp from both flanks)
	frag_start_pos: usize  // Original DNA fragment boundary position
}

fn main() {
	let matches = cli::build_cli().get_matches();
	if let ("detect", Some(args)) = matches.subcommand() {
		detect_discordant_reads(
			args.value_of("bam_file").unwrap().to_string(),
			args.value_of("genome").unwrap().to_string(),
			args.value_of("anchor-len").unwrap().parse().unwrap(),
			args.value_of("max-frag-len").unwrap().parse().unwrap());
	} else if let ("filter", Some(args)) = matches.subcommand() {
		filter(
			args.value_of("sv_file").unwrap().to_string(),
			args.value_of("min-reads").unwrap().parse().unwrap());
	}
}


// DETECT //
//======= //

/*
fn detect_discordant_pairs(sam_path: String, out_prefix: String, max_frag_len: i32,
	min_mapq: i32, orientation: String) {
  println!("{:?}\t{:?}\t{:?}\t{:?}\t{:?}", sam_path, out_prefix, max_frag_len, min_mapq, orientation);

  let mut N = 0;
  let tmp = out_prefix.clone();
  let given_path = Path::new(&tmp);
  let mut sort_tmp_dir = Path::new("./");
  if given_path.is_absolute()  {
        sort_tmp_dir = given_path.parent().unwrap();
  }

  let extn = String::from(".discordant_pairs.tsv.gz");
  let out  = out_prefix + &extn;
  let mut file = File::create("test.txt").unwrap();

  lib::mkdir(&tmp);
  println!("Searching for discordant read pairs ...");

  let (code, stdout, stderr) = sh!("sam discordant pairs --min-mapq={} {} {} | sort -k1,1 -T {:?}", min_mapq, sam_path, max_frag_len, sort_tmp_dir);

  let mut prev : Vec<&str> = Vec::new();
  for line in stdout.lines() {
      let mut cols: Vec<&str> = line.split('\t').collect();
      if cols.len() < 9 { continue; }

      let r1 = Regex::new(r"N").unwrap();
      let r2 = Regex::new(r"S").unwrap();

      if r1.is_match(cols[5]) || r2.is_match(cols[5]) { continue; }
      else if cols[0].ends_with("/1") || cols[0].ends_with("/2") {
          let tmp : Vec<&str> = cols[0].split('/').collect();
          cols[0] = tmp[0];
      }
      if prev.len() <= 0 { prev = cols.clone(); }
      else if prev.len() > 0 && cols[0] != prev[0] { prev = cols.clone(); }

      let flags = &cols[1].parse::<i32>().unwrap();
      let mut chr   = String::from(cols[2]);
      let mut mchr  = String::from(prev[2]);
      let mut strand = '-';
      let mut mstrand = '-';
      let mut pos  = cols[3].parse::<i32>().unwrap();
      let mut mpos = prev[3].parse::<i32>().unwrap();
      let mut rlen = cols[9].clone().len() as i32;
      let mut mrlen= prev[9].clone().len() as i32;

      if chr.starts_with("chr") {  chr = chr; }
      else { let prefix = String::from("chr"); chr = prefix + &chr;}

      if mchr.starts_with("chr") {  chr = chr ;}
      else { let prefix = String::from("chr"); mchr = prefix + &mchr;}

      if chr.contains("chrM") || mchr.contains("chrM") { continue;}

      if orientation == "fr" {
        if chr > mchr || (chr == mchr && pos > mpos) {
          	swap(&mut chr, &mut mchr);
          	swap(&mut pos, &mut mpos);
          	swap(&mut rlen, &mut mrlen);
          	swap(&mut strand, &mut mstrand);
        }
        if mstrand == '+' { mstrand = '-';}
        else { mstrand = '+'; }
      }

      else if orientation == "rf" {
        if chr > mchr || (chr == mchr && pos > mpos) {
          	swap(&mut chr, &mut mchr);
          	swap(&mut pos, &mut mpos);
          	swap(&mut rlen, &mut mrlen);
          	swap(&mut strand, &mut mstrand);
        }
        if strand == '+' { strand = '-';}
        else { strand = '+'; }
      }

      else if orientation == "ff" {
        if chr > mchr || (chr == mchr && pos > mpos) {
          	swap(&mut chr, &mut mchr);
          	swap(&mut pos, &mut mpos);
          	swap(&mut rlen, &mut mrlen);

          if mstrand == '-' { mstrand = '+';} else { mstrand = '-'; }
          if strand  == '-' { strand  = '+';} else { strand  = '-'; }
          	swap(&mut strand, &mut mstrand);
        }
      } else { panic!("Unsupported read orientation detected."); }

      if strand  == '-' { pos  +=  rlen - 1;}
      if mstrand == '-' { mpos += mrlen - 1;}

      let mut out = String::new();
      out.push_str(&chr);out.push_str(&"\t");
      out.push_str(&strand.to_string());out.push_str(&"\t");
      out.push_str(&pos.to_string());out.push_str(&"\t");
      out.push_str(&mchr);out.push_str(&"\t");
      out.push_str(&mstrand.to_string());out.push_str(&"\t");
      out.push_str(&mpos.to_string());out.push_str(&"\t-\n");
      file.write_all(out.as_bytes());
      N += 1;
    }
    let (wcode, wstdout, wstderr) = sh!("gzip -c {:?}.txt > {}", given_path, out);
    let mut rm_tmp = tmp.to_string() + ".txt";
    fs::remove_file(rm_tmp);
    println!("Found {:?} discordant mate pairs.", N);
}
*/

fn remove_duplicate_evidence(evidence: Vec<&Evidence>) -> Vec<&Evidence> {
	let mut checked = vec![false; evidence.len()];
	let mut filtered: Vec<&Evidence> = Vec::new();

	for (e, read) in evidence.iter().enumerate() {
		// If the read was already found to be a redundant source of evidence,
		// we don't need to check for redundancy again.
		if checked[e] { continue; }

		// We consider two reads to be redundant sources of evidence for a
		// genomic rearrangement if they come from the same original DNA
		// fragment (e.g. two paired end mates from the same fragment).
		// Also, we consider them to be redundant if they share the same
		// start position (as in this case they likely arise from PCR or
		// optical duplicates).
		let mut todo_edges = vec![read.frag_start_pos];
		let mut todo_frag_ids = vec![read.frag_id.clone()];
		let mut redundant = vec![e];
		while !todo_edges.is_empty() || !todo_frag_ids.is_empty() {
			if let Some(edge) = todo_edges.pop() {
				for k in e+1..evidence.len() {
					if checked[k] { continue; }
					if evidence[k].frag_start_pos == edge {
						checked[k] = true;
						redundant.push(k);
						todo_frag_ids.push(evidence[k].frag_id.clone());
					}
				}
			}

			if let Some(frag_id) = todo_frag_ids.pop() {
				for k in e+1..evidence.len() {
					if checked[k] { continue; }
					if evidence[k].frag_id == frag_id {
						checked[k] = true;
						redundant.push(k);
						todo_edges.push(evidence[k].frag_start_pos);
					}
				}
			}
		}

		// Now we have a list of redundant reads. Since we can only keep
		// one of these reads as a source of evidence, we pick the "best"
		// one. Currently this means the one with the longest flanks.
		let mut best = 0;
		let mut best_score = 0;
		for k in redundant {
			let seq = &evidence[k].sequence;
			let left_len = seq.iter().position(|c| *c == '|' as u8).unwrap();
			let right_len = seq.len() - left_len - 1;
			let score = min(left_len, right_len);
			if score > best_score {
				best = k; best_score = score;
			}
		}
		filtered.push(evidence[best]);
	}

	/*if filtered.len() > 1 {
		for read in &evidence {
			println!("{}\t{}", read.frag_start_pos, read.frag_id);
		}
	}*/
	filtered
}

fn detect_discordant_reads(sam_path: String, genome_path: String, anchor_len: usize, max_frag_len: usize) {

	let fasta = bio::io::fasta::Reader::from_file(format!("{}.fa", genome_path)).unwrap();
	writeln!(stderr(), "Reading reference genome into memory...");
	// FIXME: Use eprintln!() once it stabilizes...

	let mut genome = HashMap::new();
	for entry in fasta.records() {
		let chr = entry.unwrap();
		genome.insert(chr.id().to_owned(), chr.seq().to_owned());
	}

    writeln!(stderr(), "Splitting unaligned reads into {} bp anchors and aligning against the genome...", anchor_len);
    let bowtie = Command::new("bowtie")
		.args(&["-f", "-p1", "-v0", "-m1", "-B1", "--suppress", "5,6,7,8", &genome_path, "-"])
        .stdin(Stdio::piped()).stdout(Stdio::piped())
        .spawn().unwrap();

	let mut bowtie_in = BufWriter::new(bowtie.stdin.unwrap());
	let mut bowtie_out = BufReader::new(bowtie.stdout.unwrap());

	thread::spawn(move || {
		let bam = bam::Reader::from_path(&sam_path).unwrap();
		for r in bam.records() {
			let read = r.unwrap();
			if read.is_unmapped() == false { continue; }
			if read.seq().len() < anchor_len * 2 { continue; }
			let read_id = str::from_utf8(read.qname()).unwrap();
			let seq = if read.is_reverse() {
				dna::revcomp(&read.seq().as_bytes())
			} else {
				read.seq().as_bytes()
			};
	 		let tail = seq.len() - anchor_len;
			write!(bowtie_in, ">{}#1\n{}\n>{}#{}\n{}", &read_id, str::from_utf8(&seq[..anchor_len]).unwrap(), &read_id, str::from_utf8(&seq).unwrap(), str::from_utf8(&seq[tail..]).unwrap());
	    }
    });

    let mut evidence: Vec<Evidence> = Vec::new();
    let mut prev = String::new();

    for l in bowtie_out.lines() {
        let line = l.unwrap();
        let anchor_info = line.split('\t').next().unwrap().to_owned();
        if anchor_info.ends_with("#1") {
        	prev = line;
        	continue;
        } else if prev.is_empty() {
        	continue;
        }

        let frag_id = anchor_info.split('#').next().unwrap();
        if !prev.starts_with(frag_id) { continue; }

        let mut cols = prev.split('\t');
        let mut mcols = line.split('\t');
        cols.next(); mcols.next();   // Skip first column

        let mut strand = cols.next().unwrap() == "+";
        let mut mstrand = mcols.next().unwrap() == "+";
        let mut chr = cols.next().unwrap();
        let mut mchr = mcols.next().unwrap();
        let mut pos: usize = cols.next().unwrap().parse().unwrap();
        let mut mpos: usize = mcols.next().unwrap().parse().unwrap();
        let mut seq: Vec<u8> = anchor_info.split('#').last().unwrap().as_bytes().to_vec();
        let full_len: usize = seq.len();

        // We store the fragment start position before we reorient the anchors.
        let frag_start_pos = pos;
            
       	// Do not report rearrangements involving mitochondrial DNA
	 	if chr.contains('M') || mchr.contains('M') { continue; }

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
		if pos < full_len { continue; }
		if mpos < full_len { continue; }

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
	 
	 	// Check that the read sequence is not too homologous on either side
		// of the breakpoint.
		/*let mut left_match: f32 = 0.0;
	  	for k in full_len-anchor_len+1..full_len { 
 			if seq[k] == left_grch[k]{ left_match += 1.0; }
   		}
  		left_match = left_match /anchor_len as f32;
  		
  		let mut right_match: f32 = 0.0;
	  	for k in 1..anchor_len { 
 			if seq[k] == right_grch[k]{ right_match += 1.0; }
   		}
  		right_match = right_match/anchor_len as f32;
		let max_homology = 0.7;
		
		if left_match >= max_homology || right_match >= max_homology { continue; }*/
				
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

		let mut signature = vec![' ' as u8; 7];   // 3 bp from both flanks
		signature[0] = left_grch[bp - 3];
		signature[1] = left_grch[bp - 2];
		signature[2] = left_grch[bp - 1];
		signature[3] = '|' as u8;
		signature[4] = right_grch[bp];
		signature[5] = right_grch[bp + 1];
		signature[6] = right_grch[bp + 2];

		// These were used for debugging junction printing.
		//println!(" Sequence: {}", str::from_utf8(&seq).unwrap());
		//println!(" Left ref: {}", str::from_utf8(&left_grch).unwrap());
		//println!("Right ref: {}", str::from_utf8(&right_grch).unwrap());
		//println!(" Junction: {}", str::from_utf8(&junction).unwrap());
		//println!("Mismatches: {:?}", mismatches);
		//println!("Least mismatches at: {}", bp);
		//println!("-----");

		evidence.push(Evidence {
			chr: chr.to_string(), pos: pos, strand: strand,
			mchr: mchr.to_string(), mpos: mpos, mstrand: mstrand,
			sequence: junction, frag_id: frag_id.to_string(),
			signature: signature, frag_start_pos: frag_start_pos });
    }

    writeln!(stderr(), "Found {} rearrangement supporting reads.", evidence.len());

    writeln!(stderr(), "Sorting rearrangement supporting reads by position...");
    evidence.sort_by(|a,b|
    	if a.chr < b.chr { Ordering::Less }
    	else if a.chr > b.chr { Ordering::Greater }
    	else if a.pos < b.pos { Ordering::Less }
    	else if a.pos > b.pos { Ordering::Greater }
    	else { Ordering::Equal });

    writeln!(stderr(), "Identifying rearrangements based on clusters of discordant reads...");
    println!("CHROM\tSTRAND\tPOSITION\tNEARBY FEATURES\tCHROM\tSTRAND\tPOSITION\tNEARBY FEATURES\tSUPPORTING READS");
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

    	if cluster.len() < 2 { continue; }
    	cluster = remove_duplicate_evidence(cluster);
    	if cluster.len() < 3 { continue; }

    	print!("{}\t{}\t{}\t\t{}\t{}\t{}\t\t",
    		read.chr, if read.strand { '+' } else { '-' }, read.pos,
    		read.mchr, if read.mstrand { '+' } else { '-' }, read.mpos);
    	for r in 0..cluster.len() {
    		if r > 0 { print!(";"); }
    		print!("{}", str::from_utf8(cluster[r].sequence.as_slice()).unwrap());
    	}
    	println!();
	}
}


fn filter(sv_path: String, min_reads: usize) {

}

