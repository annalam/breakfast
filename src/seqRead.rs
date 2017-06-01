// src/seqRead.rs

extern crate bio;

use std::io;
use bio::io::fasta;

pub fn fasta_reader() -> Self {
        let reader = fasta::Reader::new(io::stdin());
    }
