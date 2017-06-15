
use std::io::Read;
use std::path::Path;
use std::fs::File;
use std::process::Command;
use std::io::{BufReader, BufRead};
use bio::io::fastq;

extern crate bio;

pub fn mkdir(path: &String) {
  if Path::new(path).exists(){
     println!("{:?} already exists", path);
  } else {
    Command::new("mkdir")
            .arg(path)
            .spawn()
            .expect("failed to mkdir");
  }
}

pub fn zopen(path: &String, mode: String) {

}

pub fn read_whole_file<T: Read>(mut f: T) -> String {
        let mut content = String::new();
        f.read_to_string(&mut content).unwrap();
        content
}

pub fn fasta_split_interleaved() {
    let fastq = bio::io::fasta::Reader::from_file("test.fa").unwrap();
    for record in fastq.records() {
        println!("{:?}", record.unwrap().seq());
    }
}

pub fn shell_stdout() {
}
