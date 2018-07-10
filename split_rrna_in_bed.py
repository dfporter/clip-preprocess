#!/usr/bin/python
import HTSeq
import sys
import os
import re
import argparse
import collections
import glob
import time
import collections



# Following:
# http://www.bioconductor.org/help/workflows/rnaseqGene/
#samfile = sys.argv[1]

# HTSeq segment from:
# http://www-huber.embl.de/users/anders/HTSeq/doc/counting.html#counting


def get_gtf(gtf_file=None):
    gtf = HTSeq.GFF_Reader(gtf_file)
    features = HTSeq.GenomicArrayOfSets('auto', stranded=True)
    for feature in gtf:
        if feature.type == 'exon':
            features[feature.iv] += feature.attr['gene_id']
    return features


def split_folder(bed_dir, features):
    for bed_filename in glob.glob(bed_dir + '/*.bed'):
        print "Assigning reads from %s..." % bed_filename
        map_bed(bed_filename, features)


def map_bed(bed_filename, features):
    rRNA = set([
        "WBGene00004512",
    "WBGene00004513",
    "WBGene00077468",
    "WBGene00077469",
    "WBGene00077466",
    "WBGene00077467",
    "WBGene00077465",
    "WBGene00004622",
    "WBGene00004567",
    "WBGene00014454",
    "WBGene00014472",
    "WBGene00077475",
    "WBGene00077474",
    "WBGene00077477",
    "WBGene00077476",
    "WBGene00077471",
    "WBGene00077470",
    "WBGene00077473",
    "WBGene00077472",
    "WBGene00235197",
    "WBGene00189966",
    "WBGene00014621"])

    start_time = time.clock()
    counts = collections.defaultdict(int)
    rrna, non_rrna = '', ''
    for line_num, li, in enumerate(open(bed_filename).readlines()):
        s = li.rstrip('\n').split('\t')
        read_iv = HTSeq.GenomicInterval(s[0], int(s[1]), int(s[2]), s[5])
        if line_num and (not line_num % 100000):
            print line_num
            print counts
            elapsed = time.clock() - start_time
            per_ten_million = 1e6 * elapsed/max([1, float(line_num)])
            print "Time elapsed: %f. Seconds per million reads %s." % (
                elapsed, per_ten_million)
        gene_ids = set()
        for iv, val in features[read_iv].steps():
            gene_ids |= val
        if len(rRNA & gene_ids) > 0:
            rrna += li
            counts['rrna'] += 1
        else:
            non_rrna += li
            counts['non_rrna'] += 1
    rrna_bed = 'bed_uncollapsed/rrna/' + os.path.basename(bed_filename)
    non_rrna_bed = 'bed_uncollapsed/no_rrna/' + os.path.basename(bed_filename)
    if not os.path.exists('bed_uncollapsed/'):
        os.system('mkdir bed_uncollapsed/')
    if not os.path.exists('bed_uncollapsed/rrna/'):
        os.system('mkdir bed_uncollapsed/rrna/')
    if not os.path.exists('bed_uncollapsed/no_rrna/'):
        os.system('mkdir bed_uncollapsed/no_rrna/')
    open(rrna_bed, 'w').write(rrna)
    open(non_rrna_bed, 'w').write(non_rrna)


def run(gtf_file=None, bed_folder=None):
    print "Loading gtf..."
    features = get_gtf(gtf_file=gtf_file)#'./genomes/Saccharomyces_cerevisiae.EF4.70.gtf')
    split_folder(bed_folder, features)


if __name__ == '__main__':
    """Create a combined_counts.txt file of the reads/gene for a 
    folder of .bed files. Requires only one value from lib, lib['gtf_raw'].
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bed',
                        help="""Folder of bed files.""")
    parser.add_argument('-g', '--gtf',
        default='/opt/lib/Caenorhabditis_elegans.WBcel235.78.noheader.gtf',
        help="""GTF file.""")
    args = parser.parse_args()
    run(bed_folder=args.bed, gtf_file=args.gtf)
