BreakFast is a toolkit for detecting chromosomal rearrangements
based on RNA-seq data.

TODO
====
- Allow '-' to be used as BAM file name, to read from STDIN.
- Allow user to customize how many mismatches are allowed in a breakpoint-overlapping read (currently there is no limit).
- When aligning breakpoint-overlapping reads to putative rearrangement breakpoints, take indels into account in addition to mismatches.
- Consider changing the algorithm so that it uses DNA fragments (reconstructed based on paired end reads) as the source of evidence for rearrangements, and not individual reads. This is the only way to reliably infer which DNA fragments are PCR duplicates of one another.
- Include breakpoint-spanning read pairs (where neither mate directly overlaps the breakpoint) as source of evidence for rearrangements. Implement this in such a way that the input BAM file is only read through once (to allow reading from a pipe).
