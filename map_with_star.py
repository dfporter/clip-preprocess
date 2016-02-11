"""
Calls the STAR mapper and outputs to sams/.
Then calls samtools view -q 20 to export filtered bam files
to bams/.
The unfilterefd sams are also converted to bedgraph in bedgraphs/.
"""

import glob
import sys
import os
import re
import collections
import argparse

def scribe(cmd):
    print cmd
    os.system(cmd)


def count_gaps_in_sam(samfile):
    hist = collections.defaultdict(int)
    gap_pat = re.compile('(\d+)N')
    gap = gap_pat.finditer
    header_lines = 0
    with open(samfile, 'r') as f:
        in_header = True
        while in_header:
            for li in f:
                if re.match('\A@.*', li) is not None:
                    header_lines += 1
                else:
                    in_header = 0
    with open(samfile, 'r') as f:
        for n in range(header_lines):
            next(f)
        for li in f:
            s = li.rstrip('\n').split('\t')
            # cigar = s[5]
            m = gap(s[5])
            if m is not None:
                num = 0
                for ins in m:
                    hist[min([1000,int(ins.group(1))])] += 1
                    num += 1
                    # print s[5]
                    # print ins.group(1)
                if num == 0:
                    hist[0] += 1
            else:
                hist[0] += 1
    keys = sorted(hist.keys(), key=lambda x: hist[x])
    for k in keys:
        print 'gap\t{k}:\t{n}'.format(k=k, n=hist[k])

def convert_sam_to_bed(samname):
    cmd = 'sam2bed {i} {d}/{b}'.format(
        i=samname, d='bed_uncollapsed/',
        b=os.path.basename(samname).partition('.sam')[0] + '.bed'
    )
    scribe(cmd)


def init():
    parser = argparse.ArgumentParser(
        description="Always writes to a ./sam/ directory.")
    parser.add_argument('-i', '--input_dir')
    args = parser.parse_args()
    os.system('mkdir bams/')
    os.system('mkdir sams/')
    os.system('mkdir bedgraphs/')
    os.system('mkdir bed_uncollapsed/')
    return args

def run():
    args = init()
    call_star(args)
    sort_and_index_sam(in_dir='sams/')


def call_star(args):
    in_dir = args.input_dir
    for fastq_filename in glob.glob(in_dir + '/*.fastq'):
        print fastq_filename
        bname = os.path.basename(fastq_filename).partition('.fastq')[0]
        cmd = '/groups/Kimble/Aman\ Prasad/clip/STAR-STAR_2.4.2a/bin/Linux_x86_64/STAR --runThreadN 5 '
        cmd += '--genomeDir /groups/Kimble/Aman\ Prasad/clip/STAR-STAR_2.4.2a/bin/Linux_x86_64/indexes/'
        cmd += ' --alignIntronMax 1000 --alignSJoverhangMin 100 --alignSJDBoverhangMin 15'
        cmd += ' --sjdbGTFfile /groups/Kimble/Common/fog_iCLIP/calls/lib/Caenorhabditis_elegans.WBcel235.78.noheader.gtf'
        cmd += ' --outFilterIntronMotifs RemoveNoncanonicalUnannotated'
        cmd += ' --outFilterScoreMin 20'
        cmd += ' --readMapNumber -1'
        cmd += ' --outReadsUnmapped Fastx'
        cmd += ' --outSAMstrandField intronMotif'
        cmd += ' --readFilesIn {rin} --outFileNamePrefix {prefx}'.format(
            rin=fastq_filename,
            prefx=bname + '_')
        outname = bname + '_Aligned.out.sam'
        fixed_name = outname.partition('_Aligned.out.sam')[0] + '.sam'
        print cmd
        os.system(cmd)
        os.system('mv {i} {o}'.format(i=outname, o=fixed_name))
        os.system('mv {o} sams/{o}'.format(o=fixed_name))
        convert_sam_to_bed('sams/{o}'.format(o=fixed_name))


def sort_and_index_sam(in_dir='sams/'):
    for sam in glob.glob('{a}/*.sam'.format(a=in_dir)):
        # bed_to_wig.run('sams/', 'bedgraphs/')
        bname = os.path.basename(sam).partition('.sam')[0]
        scribe('samtools view  -bS {i} > {o}'.format(
            i=sam, o=bname + '.bam'
        ))
        scribe('samtools sort {o} {s}'.format(
            o=bname + '.bam',
            s=bname
        ))
        scribe('samtools index {o}'.format(
            o=bname + '.bam'
        ))
        scribe('mv {i} bams/{o}'.format(
            i=bname + '.bam',
            o=bname + '.bam'
        ))

if __name__ == '__main__':
    run()
    
'''
for bam in glob.glob('bams/*.bam'):
    outname = re.sub('_Aligned\.out', '', bam)
    cmd = 'mv {i} {o}'.format(i=bam, o=outname)
    print cmd
    os.system(cmd)

for bai in glob.glob('*.bai'):
    outname = re.sub('_Aligned\.out', '', bai)
    cmd = 'mv {i} {o}'.format(i=bai, o=outname)
    print cmd
    os.system(cmd)

for samn in glob.glob('sams/*.sam'):
    outname = re.sub('_Aligned\.out', '', samn)
    cmd = 'mv {i} {o}'.format(i=bai, o=outname)
    print cmd
    os.system(cmd)


if not os.path.exists('bed_uncollapsed_star_q20/'):
    os.mkdir('bed_uncollapsed_star_q20')
for fname in glob.glob('bams/*.bam'):
#    os.system('bamToBed -split -i {i} > bed_uncollapsed_star_q10/{o}'.format(
#        i=fname, o=os.path.basename(fname).partition('.bam')[0] + '.bed'
#	))

    k = os.path.basename(fname).partition('.bam')[0] + '.bed'
    cmd = 'python src/bam_to_bed.py bed_uncollapsed_star_q10/{bedn} sams/{samn}'.format(
        bedn=k, samn=os.path.basename(fname).partition('.bam')[0] + '.sam'
        )
    print cmd
    os.system(cmd)
'''
