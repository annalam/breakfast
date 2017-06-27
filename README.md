BreakFast is a toolkit for detecting chromosomal rearrangements
based on RNA-seq data.

Getting started
===============

Requirements
============
- [samtools](http://samtools.sourceforge.net/) 
- [bowtie](http://bowtie-bio.sourceforge.net/index.shtml)


Unix Binary
============
- Download release


Compile from source
===================
- install [Rust](https://www.rust-lang.org/en-US/)
```
	git clone https://github.com/mgvel/breakfast.rs.git
	
	cd breakfast.rs
	cargo build --release
```
Run with 
	
```
	./target/release/breakfast 
```
Or copy to your binaries location
	
```
	cp ./target/release/breakfast [DESTINATION-PATH]
```
	
	
Detailed description of algorithm
=================================

TODO
====
- Write "Getting started" documentation for the tool.
- Document and test how Breakfast can be run on CRAM files (use samtools view to convert to BAM, then pipe into Breakfast).
- Allow '-' to be used as BAM file name, to read from STDIN.
- Allow user to customize how many mismatches are allowed in a breakpoint-overlapping read (currently there is no limit).
- Change output format so that rearrangement coordinates represent the position of the last base immediately before the breakpoint.
- If the user specifies a Bowtie1 index as the reference genome and a FASTA file with the same prefix does not exist, read the chromosome sequences from the Bowtie1 index directly, using bowtie-inspect.
- Allow Breakfast to use Bowtie2 or BWA if Bowtie1 is not available.
- When aligning breakpoint-overlapping reads to putative rearrangement breakpoints, take indels into account in addition to mismatches.
- Consider changing the algorithm so that it uses DNA fragments (reconstructed based on paired end reads) as the source of evidence for rearrangements, and not individual reads. This is the only way to reliably infer which DNA fragments are PCR duplicates of one another.
- Include breakpoint-spanning read pairs (where neither mate directly overlaps the breakpoint) as source of evidence for rearrangements. Implement this in such a way that the input BAM file is only read through once (to allow reading from a pipe).
