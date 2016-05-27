# clip-preprocess

This repository holds scripts for pre-processing of clip data.
This is not in a usable state for others!
This repository is not intended to be useful to others at the moment, as file paths are hardcoded, it is not generalized, ect.

First, the fastq files were split by barcode: 

	$ python src/split_by_barcode.py <fastq>

Several steps are done to get to uncollapsed bed files.
The script split_fastq_to_bed.py will move barcodes to read names, remove linkers, map and take the best map for each read, and finally create bams and beds from the sam output by STAR.
The sequence of folders is:

input dir -> temp_adapter_in_name/ -> temp_fastq/ -> temp_clipped/ -> sams/
 -> temp_sam_collapse_multi_map/ -> sams/ -> bed_uncollapsed/

bed_collapsed is the only output needed for everything else.

This is run by:

	$ python split_fastq_to_bed.py -i input_dir

That script does the following:

The barcodes were moved to the read header: 

	$ python move_barcode_to_name_in_fastq.py <fastq directory/>

The three prime linker was removed:

	$ python clip_adapter.py -t -i <in_dir> -o <out_dir>

The RT primer was removed:

	$ python clip_adapter.py -r -i <in_dir> -o <out_dir>

Reads were mapped to the genome:

	$ python map_with_star.py -i <in_dir>

WARNING: this version of STAR will output all loci for a multimapped read, and has no option to do otherwise, except to reject all multimapped reads.

Collapse multimapping reads to one:

	$ python collapse_multimapping_reads_to_one_sam.py <file or folder of sams>

This will output to the folder sam_collapse_multi_map/ by default.

TO DO: compare with the published STAR parameters from that comparison paper.

This outputs both bams to bams/ and beds to bed_uncollapsed/.

Commands to run the duplicates-removal script from the Zhang lab were printed with a call to:

	$ python write_commands_for_collapsing.py

The commands output by this script were then run, outputing to bed_collapsed/.

Bedgraphs were made and normalized with:

	$ python bed_to_wig.py -b bed_collapsed/ -o bedgraph_norm/ -u bedgraph_unnorm/

Bedgraphs can be normalized separately with:

	$ python normalize_bedgraph.py -b bed_folder -i bedgraph_folder -o output_bedgraph_folder

So to summarize the essential chain of commands:

	$ python split_fastq_to_bed.py -i input_dir
	$ python write_commands_for_collapsing.py
	$ python bed_to_wig.py -b bed_collapsed/ -o bedgraph_norm/ -u bedgraph_unnorm/

Peaks are called by importing peaks-by-permutations and running:

	$ generate_config_file.py .  # Edit output as needed.
	$ find_peaks_by_permutations.py -x -c auto.ini

Peak calling is covered elsewhere.
