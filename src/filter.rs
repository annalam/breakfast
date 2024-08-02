
use crate::common::{parse_args, FileReader, create_file};
use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::fs;

const USAGE: &str = "
Usage:
  breakfast filter [options] <sample_sheet>

Options:
  --min-reads=N        Minimum required number of supporting reads [default: 4]
  --min-distance=N     Minimum distance between breakpoints [default: 10]
  --keep-duplicates    Do not merge supporting reads with identical sequences
  --output-dir=PATH    Directory for output files [default: annotated]
";

struct Evidence {
	sample: String,
	reads: usize
}

#[derive(Eq, Hash, PartialEq)]
struct Rearrangement {
	chr_1: String,
	pos_1: usize,
	strand_1: char,
	chr_2: String,
	pos_2: usize,
	strand_2: char
}

struct Sample {
	name: String,
	patient: String,
	germline: bool
}

pub fn main() {
	let args = parse_args(USAGE);
	let min_reads: usize = args.get_str("--min-reads").parse().unwrap_or_else(|_| error!("--min-reads must be a positive integer."));
	let min_distance: usize = args.get_str("--min-distance").parse().unwrap();
	let merge_duplicates = !args.get_bool("--keep-duplicates");
	let output_dir = args.get_str("--output-dir");

	eprintln!("Output files will be written in directory {}", &output_dir);

	eprintln!("Reading sample sheet into memory...");
	let cohort_path = args.get_str("<sample_sheet>");
	let mut sample_sheet = FileReader::new(&cohort_path);
	let mut line = String::new();
	sample_sheet.read_line(&mut line);
	let headers: Vec<&str> = line.split('\t').collect();
	let patient_col = headers.iter().position(|&x| x == "PATIENT").unwrap_or(99);
	let sample_col = headers.iter().position(|&x| x == "SAMPLE").unwrap_or(99);
	let type_col = headers.iter().position(|&x| x == "TYPE").unwrap_or(99);
	if patient_col == 99 || sample_col == 99 || type_col == 99 {
		error!("Required headers PATIENT, SAMPLE and TYPE were not all found on the first row of the sample sheet.");
	}

	let min_columns = *[patient_col, sample_col, type_col].iter().max().unwrap();
	let mut samples: Vec<Sample> = Vec::new();
	while sample_sheet.read_line(&mut line) {
		if line.trim().is_empty() { continue; }
		let cols: Vec<&str> = line.split('\t').collect();
		if cols.len() < min_columns {
			error!("Sample sheet line has too few columns:\n{}\n", &line);
		}
		let patient = cols[patient_col].trim();
		let sample = cols[sample_col].trim();
		let sample_type = cols[type_col].trim().to_lowercase();
		if patient == "" || sample == "" || sample_type == "" {
			error!("Missing PATIENT, SAMPLE or TYPE column on line:\n{}\n", &line);
		}

		let germline = sample_type == "wbc" || sample_type == "gdna" || sample_type == "germline";
		
		samples.push(Sample {
			name: sample.into(), patient: patient.into(), germline
		});
	}

	eprintln!("Reading rearrangements into memory...");
	let mut rearrangements: HashMap<Rearrangement, Vec<Evidence>> = HashMap::new();
	for sample in &samples {
		let mut sv_file = FileReader::new(&format!("{}.sv", sample.name));
		sv_file.read_line(&mut line);  // Discard the header line
		while sv_file.read_line(&mut line) {
			let cols: Vec<&str> = line.split('\t').collect();
			if cols.len() < 9 {
				error!("File {}.sv has a line with too few columns:\n{}\n",
					&sample.name, &line);
			}
			let r = read_rearrangement(&cols);
			let mut reads = cols[8].split(';').collect();
			if merge_duplicates { deduplicate_reads(&mut reads); }
			let evidence = rearrangements.entry(r).or_insert(Vec::new());
			evidence.push(Evidence { sample: sample.name.clone(), reads: reads.len() });
		}
	}


	eprintln!("Generating output files:");
	fs::create_dir_all(&output_dir); // Create the output directory if necessary
	for sample in &samples {
		let negative_controls: HashSet<String> = samples.iter()
			.filter(|s| s.patient != sample.patient || (sample.germline == false && s.germline)).map(|s| s.name.clone()).collect();
		if negative_controls.is_empty() {
			error!("No negative controls found for sample {}.", &sample.name);
		}

		let mut out = create_file(&format!("{}/{}.sv", &output_dir, &sample.name));

		let mut sv_file = FileReader::new(&format!("{}.sv", sample.name));
		sv_file.read_line(&mut line);
		write!(out, "{}\n", &line);
		while sv_file.read_line(&mut line) {
			let cols: Vec<&str> = line.split('\t').collect();
			let r = read_rearrangement(&cols);
			if r.chr_1 == r.chr_2 && r.strand_1 == r.strand_2 && r.pos_1.abs_diff(r.pos_2) < min_distance {
				continue;   // Breakpoints are too adjacent, skip
			}
			let mut reads = cols[8].split(';').collect();
			if merge_duplicates { deduplicate_reads(&mut reads); }
			if reads.len() < min_reads { continue; }
			let evidence = &rearrangements[&r];
			let max_control_reads = evidence.iter().filter(|e| negative_controls.contains(&e.sample)).map(|e| e.reads).max().unwrap_or(0);
			if reads.len() < 50 * max_control_reads { continue; }

			write!(out, "{}", cols[0]);
			for k in 1..8 { write!(out, "\t{}", cols[k]); }
			write!(out, "\t{}", reads.join(";"));
			for k in 9..cols.len() { write!(out, "\t{}", cols[k]); }
			write!(out, "\n");
		}
	}
}

fn read_rearrangement(cols: &Vec<&str>) -> Rearrangement {
	Rearrangement {
		chr_1: cols[0].into(), pos_1: cols[2].parse().unwrap(),
		strand_1: cols[1].chars().nth(0).unwrap(),
		chr_2: cols[4].into(), pos_2: cols[6].parse().unwrap(),
		strand_2: cols[5].chars().nth(0).unwrap()
	}
}

fn deduplicate_reads(reads: &mut Vec<&str>) {
	// Deduplicate based on last 20 bases of read
	reads.sort_by_key(|r| { let len = r.len(); &r[len-20..len] });
	reads.dedup_by_key(|r| { let len = r.len(); &r[len-20..len] });

	// Deduplicate based on first 20 bases of read
	reads.sort();
	reads.dedup_by_key(|r| &r[0..20]);
}
