// rust version of breakfast tool
// BreakFast is a toolkit for detecting chromosomal rearrangements
// based on RNA-seq data.

use std::env;
use std::mem;
use std::path::Path;
use std::io::prelude::*;
use std::fs;
use std::fs::File;
use std::process::{Command, Stdio};
use std::io::{BufReader, BufRead, Read, Write};

#[macro_use]
extern crate clap;
#[macro_use]
extern crate shells;
#[macro_use]
extern crate bio;
#[macro_use]
extern crate regex;
#[macro_use]
extern crate subprocess;

use regex::Regex;
use subprocess::{Exec, Popen, PopenConfig, Redirection};

mod cli;
mod lib;

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
  //arguments from detect subcommand
  if let ("detect", Some(detect)) = matches.subcommand() {
      let bam_file = String::from(detect.value_of("bam_file").unwrap());
      let genome   = String::from(detect.value_of("genome").unwrap());
      let out_prefix  = String::from(detect.value_of("out_prefix").unwrap());

      if detect.is_present("anchor-len") {
          anchor_len = detect.value_of("anchor-len").unwrap().parse::<i32>().unwrap();
      }

      if detect.is_present("orientation") {
          orientation = String::from(detect.value_of("orientation").unwrap());
      }
      //detect_discordant_pairs(bam_file, out_prefix, max_frag_len, min_mapq, orientation);
      detect_discordant_reads(bam_file, genome, out_prefix, anchor_len);
  }

  //arguments from detectspecific subcommand
  if let ("detectspecific", Some(detectspecific)) = matches.subcommand() {
      let bam_file = String::from(detectspecific.value_of("bam_file").unwrap());
      let genome   = String::from(detectspecific.value_of("genome").unwrap());
      let out_prefix  = String::from(detectspecific.value_of("out_prefix").unwrap());
      let donors  = String::from(detectspecific.value_of("donors").unwrap());
      let acceptors  = String::from(detectspecific.value_of("acceptors").unwrap());
  }

  println!("Hello, world!");
}


// DETECT //
//======= //

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
          mem::swap(&mut chr, &mut mchr);
          mem::swap(&mut pos, &mut mpos);
          mem::swap(&mut rlen, &mut mrlen);
          mem::swap(&mut strand, &mut mstrand);
        }
        if mstrand == '+' { mstrand = '-';}
        else { mstrand = '+'; }
      }

      else if orientation == "rf" {
        if chr > mchr || (chr == mchr && pos > mpos) {
          mem::swap(&mut chr, &mut mchr);
          mem::swap(&mut pos, &mut mpos);
          mem::swap(&mut rlen, &mut mrlen);
          mem::swap(&mut strand, &mut mstrand);
        }
        if strand == '+' { strand = '-';}
        else { strand = '+'; }
      }

      else if orientation == "ff" {
        if chr > mchr || (chr == mchr && pos > mpos) {
          mem::swap(&mut chr, &mut mchr);
          mem::swap(&mut pos, &mut mpos);
          mem::swap(&mut rlen, &mut mrlen);

          if mstrand == '-' { mstrand = '+';} else { mstrand = '-'; }
          if strand  == '-' { strand  = '+';} else { strand  = '-'; }
          mem::swap(&mut strand, &mut mstrand);
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


fn detect_discordant_reads(sam_path: String, genome_path: String, out_prefix: String, anchor_len: i32) {

    //let mut N: i32 = 0;
    println!("Splitting unaligned reads into {} bp anchors and aligning against the genome...", anchor_len);

    println!("{:?}\t{:?}\t{:?}\t", &sam_path, &anchor_len, &genome_path);

    let mut sam = Command::new("samtools")
                    .args(&["fasta", "-f", "0x4", &sam_path])
                    .stdout(Stdio::piped())
                    .spawn()
                    .unwrap();

    let mut fas = Command::new("fasta")
                  .args(&["split", "interleaved", "-", &anchor_len.to_string()])
                  .stdin(Stdio::piped())
                  .stdout(Stdio::piped())
                  .spawn()
                  .unwrap();

    let mut bow = Command::new("bowtie")
                  .args(&["-f", "-p1", "-v0", "-m1", "-B1", "--suppress", "5,6,7,8", &genome_path, "-"])
                  .stdin(Stdio::piped())
                  .stdout(Stdio::piped())
                  .spawn()
                  .unwrap();
    /*
    if let Some(ref mut stdout)  = sam.stdout {
      if let Some(ref mut stdin) = fas.stdin {
          let mut buf: Vec<u8> = Vec::new();
          stdout.read_to_end(&mut buf).unwrap();
          stdin.write_all(&buf).unwrap();
          }
        }
    let res = fas.wait_with_output().unwrap().stdout;
    */
    let mut samres = String::new();
      match sam.stdout.unwrap().read_to_string(&mut samres) {
          Err(why) => panic!("samtools not running!"),
          Ok(_) => print!("samtools worked!"),
      }

    match fas.stdin.unwrap().write_all(samres.as_bytes()) {
        Err(why) => panic!("samtools results not reaching fasta"),
        Ok(_) => print!("fasta recvied input from samtools"),
    }

    let mut res = String::new();
      match fas.stdout.unwrap().read_to_string(&mut res) {
          Err(why) => panic!("fasta not running!"),
          Ok(_) => print!("fasta subcommand worked!"),
      }
      println!("{:?}", res.len());
}
