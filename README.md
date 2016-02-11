# clip-preprocess

This repository will hold all of the scripts and output data for the sp/oo FBF CLIP-seq data.

First, the fastq files were split by barcode: 

	$ python src/split_by_barcode.py <fastq>

The barcodes were moved to the read header: 

	$ python move_barcode_to_name_in_fastq.py <fastq directory/>

The three prime linker was removed:

	$ python clip_adapter.py -t -i <in_dir> -o <out_dir>

The RT primer was removed:

	$ python clip_adapter.py -r -i <in_dir> -o <out_dir>

Reads were mapped to the genome:

	$ python map_with_star.py -i <in_dir>

TO DO: compare with the published STAR parameters from that comparison paper.

This outputs both bams to bams/ and beds to bed_uncollapsed/.

Bedgraphs were made with:

	$ python bed_to_wig.py <in_dir> <out_dir>

Bedgraphs were normalized with:

	$ python normalize_bedgraph.py <in_dir> <out_dir> 
