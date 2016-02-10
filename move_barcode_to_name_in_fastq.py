__author__ = 'dp'
import HTSeq
import os
import sys
import re
import glob
import argparse


def move_barcode_to_name_in_fastq(filename, out_dir):
    if not os.path.exists(out_dir): os.system('mkdir ' + out_dir)
    outf = open(out_dir + '/%s' % os.path.basename(filename), 'w')
    fastq = HTSeq.FastqReader(filename)
    obs_let = set()
    # phred 33
    for read in fastq:
        if len(read.seq) < 14: continue
        if min(read.qualstr[:9]) < 53: continue
        n_read = HTSeq.SequenceWithQualities(
            read.seq[9:],
            read.name.partition(' ')[0] + '#' + read.seq[0:9],
            read.qualstr[9:])
        n_read.write_to_fastq_file(outf)
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''Move barcodes to read header.''')
    parser.add_argument('-i', '--input_dir',
                     help='Input directory.')
    args = parser.parse_args()
    for filename in glob.glob(args.input_dir + '/*.fastq'):
        move_barcode_to_name_in_fastq(filename, 'adapter_moved_to_name/')


