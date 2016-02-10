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


