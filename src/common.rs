
use docopt::{Docopt, ArgvMap};
use std::process::{Command, Stdio};
use std::io::{stdin, BufRead, BufReader};
use std::fs::File;
use rust_htslib::bam::{self, Read};

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

pub fn create_file(path: &str) -> File {
	File::create(&path).unwrap_or_else(
		|_| error!("Cannot create file {}.", &path))
}

pub struct FileReader {
	bufread: Box<dyn BufRead>
}

impl FileReader {
	pub fn new(path: &str) -> FileReader {
		let bufread: Box<dyn BufRead> = if path == "-" {
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
		let Ok(bytes_read) = self.bufread.read_line(line) else {
			error!("I/O error while reading from file.");
		};
		if line.ends_with('\n') { line.pop(); }  // Remove trailing newline
		bytes_read > 0
	}
}

// Function for reading BAM records, with proper user-friendly messages.
// Returns false after reading the last record, or if reading fails.
pub fn read_bam_record(bam: &mut bam::Reader, record: &mut bam::Record) -> bool {
	match bam.read(record) {
		Some(Ok(())) => true,
		Some(Err(e)) => error!("Failed reading BAM record: {}", e),
		None => false
	}
}
