import glob
import argparse
import sys
import os
import re

def normalize_wig(bed_folder, input_bedgraphs, output_bedgraphs):
    if not os.path.exists(output_bedgraphs): os.system('mkdir ' + output_bedgraphs)
    for wig in glob.glob(input_bedgraphs + '/*.wig'):
        if re.search('all', os.path.basename(wig)):
            continue
        bed = '{a}/{b}.bed'.format(
            a=bed_folder, b=re.sub('_[\+-].wig', '', os.path.basename(wig)))
        print "{a} > {b}".format(a=wig, b=bed)
        out_bed = output_bedgraphs + '/' + os.path.basename(wig)
        print "Normalizing {a} / {b} to {c}".format(
            a=wig, b = bed, c=out_bed)

        if os.path.exists(out_bed):
            print "{A} file exists, skipping..".format(A=out_bed)
            continue
        total_reads = float(sum([1 for line in open(bed)]) - 1)
        wig_fh = open(wig, 'r')
        outli = next(wig_fh)  # Skip header.
        for line in wig_fh:
            s = line.rstrip('\n').split('\t')
            outli += '\t'.join(s[:3]) + '\t' 
            outli += str(1e6 * float(s[3])/total_reads) + '\n'
        open(out_bed, 'w').write(outli)        


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--input_bed',
        help='Folder of bed files for normalization.')
    parser.add_argument('-i', '--input_bedgraph',
        help='Folder of unnormalized bedgraph files (with .wig suffix) to normalize.')
    parser.add_argument('-o', '--output_bedgraph',
        help='Output folder for normalized bedgraphs.')
    args = parser.parse_args()
    normalize_wig(args.input_bed, args.input_bedgraph, args.output_bedgraph)


