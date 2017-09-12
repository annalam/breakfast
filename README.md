Breakfast
---------

Breakfast is a software for detecting genomic structural variants from DNA sequencing data. Its features include:
- Identifies structural variants based on breakpoint-overlapping reads
- Extremely fast, analyzes XX million reads per second
- Identifies PCR/optical duplicates and does not count them as independent sources of evidence
- Can be run on sorted or unsorted BAM files, or can read BAM input from a pipe
- Uses pre-existing Bowtie indexes to speed up alignment (does not require its own index)
- Provides tools for filtering out rearrangements that are present in control samples


Installation
------------

The easiest way to install Breakfast is to download one of the pre-built binary packages:
- Breakfast 0.1 (x86-64 Linux)

If a suitable binary package is not available, you can also build Breakfast directly from source code. Note that installing this way requires a Rust compiler and the Cargo build system to be available:
```
git clone https://github.com/mgvel/breakfast.rs.git
cd breakfast.rs
cargo install --release
```


Running Breakfast
-----------------

To run BreakFast, you need a BAM file containing sequenced reads (in this example, tumor.bam). You also need a Bowtie index and the Bowtie1 executable in your PATH. A Breakfast analysis begins with the "breakfast detect" command, which searches the BAM file for unaligned reads that support a genomic breakpoint:
```
breakfast detect tumor.bam bowtie_indexes/hg38 > tumor.sv
```




Detailed overview of the Breakfast algorithm
--------------------------------------------

Unaligned reads are split into two anchors of customizable size: one anchor from the 5' end of the read, and one anchor from the 3' end of the read. These anchors are then aligned against the reference genome. If both anchors align to the reference genome (but the read as a whole did not), the read is considered to support the existence of a genomic rearrangement.

Duplicate DNA fragments are identified based on "fragment signatures". When reporting candidate rearrangements, Breakfast shows the full sequence of all supporting reads.

TODO
----
- [ ] Write "Getting started" documentation for the tool.
- [ ] Document and test how Breakfast can be run on CRAM files (use samtools view to convert to BAM, then pipe into Breakfast).
- [ ] Allow '-' to be used as BAM file name, to read from STDIN.
- [ ] Allow user to customize how many mismatches are allowed in a breakpoint-overlapping read (currently there is no limit).
- [ ] Change output format so that rearrangement coordinates represent the position of the last base immediately before the breakpoint.
- [ ] If the user specifies a Bowtie1 index as the reference genome and a FASTA file with the same prefix does not exist, read the chromosome sequences from the Bowtie1 index directly, using bowtie-inspect.
- [ ] Allow Breakfast to use Bowtie2 or BWA if Bowtie1 is not available.
- [ ] When aligning breakpoint-overlapping reads to putative rearrangement breakpoints, take indels into account in addition to mismatches.
- [ ] Include breakpoint-spanning read pairs (where neither mate directly overlaps the breakpoint) as source of evidence for rearrangements. Implement this in such a way that the input BAM file is only read through once (to allow reading from a pipe).
