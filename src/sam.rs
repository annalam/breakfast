pub fn read_sam(sam_path: String, mode: Strig, min_quality: i32) {
	      let mut view_options = '';
	      let mut flag_on = 0x0;
	      let mut flag_off = 0x900;  //Ignore secondary and supplementary alignments

        match flag_on {
          A => 0x1,
          C => 0x3,
          u => 0x4,
          1 => 0x40,
          2 => 0x80,
          - => 0x10,
          _ => println!("No matching flag!");
        }

        match flag_off {
          a => 0x4,
          + => 0x10,
          D => 0x400,
          _ => println!("No matching flag!");
        }

	view_options += '-f 0x%x -F 0x%x ' % (flag_on, flag_off)

	if min_quality > 0: view_options += '-q%d ' % min_quality

	out = shell_stdout('samtools view %s %s' % (view_options, sam_path))
	for line in out:
		yield line.split('\t')
  }
