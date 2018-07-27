"""
Creates bedgraph files of the read coverage, both normalized and raw.
"""

__author__ = 'dp'
import sys
import glob
import re
import HTSeq
import os
import normalize_bedgraph
import argparse

def add_to_ga(infile, global_ga):
    ga = HTSeq.GenomicArray('auto', stranded=True)

    above_len_cutoff = 0
    with open(infile, 'r') as f:
        i = 0
        for i, li in enumerate(f):
            s = li.rstrip('\n').split('\t')
            iv = HTSeq.GenomicInterval(
                s[0], int(s[1]), int(s[2]), s[5]
            )

            if (int(s[2]) - int(s[1])) > 101:
                #print('{}:'.format(i), li)
                #sys.stdin.readline()
                above_len_cutoff += 1
                continue
                
            ga[iv] += 1
            global_ga[iv] += 1
    print("\t...{i} reads in bed file {z}.".format(
        i=i, z=infile))
    return ga


def run(input_bed, output_bedgraph_unnorm, output_bedgraph_norm):
    if not os.path.exists(output_bedgraph_unnorm):
        os.system('mkdir ' + output_bedgraph_unnorm)
    if not os.path.exists(output_bedgraph_norm):
        os.system('mkdir ' + output_bedgraph_norm)
    ga_all_exp = HTSeq.GenomicArray('auto', stranded=True)
    ga_all_control = HTSeq.GenomicArray('auto', stranded=True)
    ga_other = HTSeq.GenomicArray('auto', stranded=True)
    for infile in glob.glob(input_bed + '/*.bed'):
        ga = HTSeq.GenomicArray('auto', stranded=True)
        if (re.match('.*fog.*', os.path.basename(infile)) is not None) or (
                re.match('.*exp.*', os.path.basename(infile)) is not None) or (
                re.match('.*fbf.*', os.path.basename(infile)) is not None):
            # if re.match('.*fbf1.*', os.path.basename(infile)) is not None:
            #     continue
            print(infile)
            ga = add_to_ga(infile, ga_all_exp)
        elif (re.match('.*control.*', os.path.basename(infile)) is not None) or (
            re.match('.*n2.*', os.path.basename(infile)) is not None):
            ga = add_to_ga(infile, ga_all_control)
        else:
            ga = add_to_ga(infile, ga_other)
        outname = "{d}/{b}".format(
            d=output_bedgraph_unnorm,
            b=os.path.basename(infile).partition('.bed')[0])
        print("Creating a bedgraph {c} from {a}...".format(
            c=outname, a=infile))
        outname_plus = outname + '_+.wig'
        ga.write_bedgraph_file(outname_plus, strand='+')
        outname_minus = outname + '_-.wig'
        ga.write_bedgraph_file(outname_minus, strand='-')
    ga_all_exp.write_bedgraph_file(
        output_bedgraph_unnorm + '/all_exp_+.wig', strand='+')
    ga_all_exp.write_bedgraph_file(
        output_bedgraph_unnorm + '/all_exp_-.wig', strand='-')
    ga_all_control.write_bedgraph_file(
        output_bedgraph_unnorm + '/all_control_+.wig', strand='+')
    ga_all_control.write_bedgraph_file(
        output_bedgraph_unnorm + '/all_control_-.wig', strand='-')
    normalize_bedgraph.normalize_wig(input_bed, output_bedgraph_unnorm,
                                     output_bedgraph_norm)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--input_bed',
        help='Folder of bed files.')
    parser.add_argument('-o', '--output_bedgraph_norm',
        help='Output folder for normalized bedgraphs.')
    parser.add_argument('-u', '--output_bedgraph_unnorm',
        help='Output folder for unnormalized bedgraphs.')
    args = parser.parse_args()
    run(args.input_bed, args.output_bedgraph_unnorm, args.output_bedgraph_norm)

