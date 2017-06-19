// rust version of breakfast tool
// BreakFast is a toolkit for detecting chromosomal rearrangements
// based on RNA-seq data.

use std::mem::swap;
use std::thread;
use std::path::Path;
use std::io::prelude::*;
use std::fs;
use std::str;
use std::fs::File;
use std::process::{Command, Stdio};
use std::collections::HashMap;
use std::io::{BufReader, BufWriter, BufRead, Read, Write};
use rust_htslib::bam;
use rust_htslib::bam::Read as HTSRead;

extern crate clap;
extern crate bio;
extern crate rust_htslib;

use bio::alphabets::dna;

mod cli;

// struct for having breakfast options
struct BfOptions {
    anchor_len :    usize,
    max_frag_len:   i32,
    min_mapq:       i32,
    orientation:    &'static str,
    all_reads:      bool,
    discard_duplicates: &'static str,
    min_reads:      &'static str,
    freq_above:     i32,
}

struct Evidence {
	chr: String,
	pos: u32,
	strand: bool,
	mchr: String,
	mpos: u32,
	mstrand: bool,
	sequence: String,     // Full sequence of breakpoint overlapping read
	frag_id: String,
	signature: String     // Breakpoint signature (5 bp from both flanks)
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

  let mut anchor_len: usize   = defaults.anchor_len;
  let mut max_frag_len: i32 = defaults.max_frag_len;
  let mut min_mapq: i32     = defaults.min_mapq;
  let mut orientation       = String::from(defaults.orientation);
  let mut all_reads: bool   = defaults.all_reads;
  let mut discard_duplicates = String::from(defaults.discard_duplicates);
  let mut min_reads     = String::from(defaults.min_reads);
  let mut freq_above: i32 = defaults.freq_above;

  //arguments from src/cli.rs
  let matches = cli::build_cli().get_matches();
  //arguments from detect subcommand
  if let ("detect", Some(detect)) = matches.subcommand() {
      let bam_file = String::from(detect.value_of("bam_file").unwrap());
      let genome   = String::from(detect.value_of("genome").unwrap());

      if detect.is_present("anchor-len") {
          anchor_len = detect.value_of("anchor-len").unwrap().parse().unwrap();
      }

      if detect.is_present("orientation") {
          orientation = String::from(detect.value_of("orientation").unwrap());
      }
      detect_discordant_reads(bam_file, genome, anchor_len);
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


fn detect_discordant_reads(sam_path: String, genome_path: String, anchor_len: usize) {

	let fastq = bio::io::fasta::Reader::from_file(format!("{}.fa", genome_path)).unwrap();
	println!("Reading reference genome into memory...");

	let mut genome = HashMap::new();
	for entry in fastq.records() {
		let chr = entry.unwrap();
		genome.insert(chr.id().unwrap().to_owned(), chr.seq().to_owned());
	}

    println!("Splitting unaligned reads into {} bp anchors and aligning against the genome...", anchor_len);
    let mut bowtie = Command::new("bowtie")
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
			let seq = read.seq().as_bytes();
			let tail = seq.len() - anchor_len;
			write!(bowtie_in, ">{}#1\n{}\n>{}#{}\n{}", &read_id, str::from_utf8(&seq[..anchor_len]).unwrap(), &read_id, str::from_utf8(&seq).unwrap(), str::from_utf8(&seq[tail..]).unwrap());
	    }
    });

    let mut evidence: Vec<Evidence> = Vec::new();
    let mut prev: String = String::new();

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
            
       	// Do not report rearrangements involving mitochondrial DNA
	 	if chr.contains('M') || mchr.contains('M') { continue; }

        if chr > mchr || (chr == mchr && pos > mpos) {
        	swap(&mut chr, &mut mchr);
        	swap(&mut pos, &mut mpos);
          	if mstrand == false { mstrand = true ;} else { mstrand = false ; }
          	if strand  == false { strand  = true ;} else { strand  = false ; }
          	swap(&mut strand, &mut mstrand);
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

	 	let right_grch = if strand == true {
	 		genome[chr][mpos+anchor_len-full_len..mpos+full_len-1].to_vec()
	 	} else {
	 		dna::revcomp(&genome[chr][mpos-1..mpos+full_len-1].to_vec())
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
		
	 	continue;
		
		// Identify the breakpoint location that minimizes the number of
		// nucleotide mismatches between the read and the breakpoint flanks.
		let mut mismatches: Vec<usize> = vec![0; full_len - anchor_len + 1];
//		println!("How mary can be there? {:?}", mismatches.capacity());	
				
		for k in 1..anchor_len+1 { if seq[k] != left_grch[k] {
		mismatches[anchor_len] += 1; }}

		for k in anchor_len..full_len+1 { if seq[k] != right_grch[k]{
		mismatches[anchor_len] += 1; }}
		

		for bp in anchor_len..full_len+1 - anchor_len {
			let mut lmatch:usize = 0 ; let mut rmatch:usize = 0; 
			if seq[bp] != left_grch[bp]{ lmatch += 1; }
			if seq[bp] !=right_grch[bp]{ rmatch += 1; }
			mismatches[bp] = mismatches[bp-1] + lmatch - rmatch;
		
		}
		let mut bp = 0;
		if !mismatches.is_empty() {
			//mismatches = mismatches.sort();
			bp = mismatches[0];
			}
	
		println!("{:?} breakpoint", bp);	
//		println!("How mary are there? {:?}", mismatches.len());	
		
    }
 	/* TEST VARAIBLES   
    let chr = "chr7";
    let seq  = "atcgtagtcgtacgtagctagatgctagatgctag";
 	let full_len = seq.len();
 	let pos = full_len / 4;
 	let left_grch = "atgctatcgtagtcgtacgtagctagatgctagatgctagacgtagctagat";
  	let mut left_match = 0;
  	for k in (full_len-anchor_len+1..full_len) { 
 		if seq.chars().nth(k) == left_grch.chars().nth(k){
  			left_match += 1;
  			println!("{:?}", seq.chars().nth(k).unwrap()); 
  		}
  	}
  	left_match = left_match/anchor_len;
 	println!("{:?}", left_match); */
}
