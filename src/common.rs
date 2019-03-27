
use docopt::{Docopt, ArgvMap};
use std::process::{Command, Stdio};
use std::io::{stdin, BufRead, BufReader};
use std::fs::File;
use rust_htslib::bam;
use rust_htslib::bam::{Read, ReadError};

macro_rules! error {
	($($arg:tt)+) => ({
		use std::process::exit;
		eprint!("ERROR: "); eprintln!($($arg)+); exit(-1);
	})
}

pub fn parse_args(usage: &str) -> ArgvMap {
	Docopt::new(usage).unwrap().parse().unwrap_or_else(|_| {
		error!("Invalid arguments.\n{}", usage);
	})
}

pub struct FileReader {
	bufread: Box<BufRead>
}

impl FileReader {
	pub fn new(path: &str) -> FileReader {
		let bufread: Box<BufRead> = if path == "-" {
			Box::new(BufReader::new(stdin()))
		} else {
			let file = File::open(path).unwrap_or_else(
				|_| error!("Cannot open file {} for reading.", path));
			if path.ends_with(".gz") {
				Box::new(BufReader::new(Command::new("gunzip").arg("-c")
					.stdout(Stdio::piped()).stdin(file).spawn()
					.unwrap_or_else(|_| error!("Cannot start gunzip process."))
					.stdout.unwrap()))
			} else {
				Box::new(BufReader::new(file))
			}
		};
		FileReader { bufread }
	}

	pub fn read_line(&mut self, line: &mut String) -> bool {
		line.clear();
		match self.bufread.read_line(line) {
			Ok(len) => len > 0,
			_ => { error!("I/O error while reading from file."); }
		}
	}
}

// Function for reading BAM records, with proper user-friendly messages.
// Returns false after reading the last record, or if reading fails.
pub fn read_bam_record(bam: &mut bam::Reader, record: &mut bam::Record) -> bool {
	match bam.read(record) {
		Err(ReadError::NoMoreRecord) => false,
		Err(ReadError::Truncated) => error!("BAM file ended prematurely."),
		Err(ReadError::Invalid) => error!("Invalid BAM record."),
		Ok(_) => true
	}
}
