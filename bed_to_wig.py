"""
python bed_to_wig.py input_dir output_dir

Creates bedgraph files of the read coverage.
"""

__author__ = 'dp'
import sys
import glob
import re
import HTSeq
import os


def add_to_ga(infile, global_ga):
    ga = HTSeq.GenomicArray('auto', stranded=True)
    with open(infile, 'r') as f:
        for li in f:
            s = li.rstrip('\n').split('\t')
            iv = HTSeq.GenomicInterval(
                s[0], int(s[1]), int(s[2]), s[5]
            )
            if (int(s[2]) - int(s[1])) > 101:
                print li
                sys.stdin.readline()
            ga[iv] += 1
            global_ga[iv] += 1
    return ga


def run(indir, outdir):
    if not os.path.exists(outdir):
        os.system('mkdir ' + outdir)
    ga_all_fog = HTSeq.GenomicArray('auto', stranded=True)
    ga_all_control = HTSeq.GenomicArray('auto', stranded=True)
    for infile in glob.glob(indir +'*.bed'):
        print infile
        ga = HTSeq.GenomicArray('auto', stranded=True)
        if (re.match('.*fog.*', os.path.basename(infile)) is not None) or (
                re.match('.*exp.*', os.path.basename(infile)) is not None) or (
                re.match('.*fbf.*', os.path.basename(infile)) is not None):
            # if re.match('.*fbf1.*', os.path.basename(infile)) is not None:
            #     continue
            print infile
            ga = add_to_ga(infile, ga_all_fog)
        if (re.match('.*control.*', os.path.basename(infile)) is not None) or (
            re.match('.*n2.*', os.path.basename(infile)) is not None
        ):
            ga = add_to_ga(infile, ga_all_control)
        outname = "{d}/{b}".format(
            d=outdir,
            b=os.path.basename(infile).partition('.bed')[0])
        outname_plus = outname + '_+.wig'
        ga.write_bedgraph_file(outname_plus, strand='+')
        outname_minus = outname + '_-.wig'
        ga.write_bedgraph_file(outname_minus, strand='-')
    ga_all_fog.write_bedgraph_file(outdir + '/all_exp_+.wig', strand='+')
    ga_all_fog.write_bedgraph_file(outdir + '/all_exp_-.wig', strand='-')
    ga_all_control.write_bedgraph_file(outdir + '/all_control_+.wig', strand='+')
    ga_all_control.write_bedgraph_file(outdir + '/all_control_-.wig', strand='-')

'''
nohup python /groups/Kimble/Common/fog_iCLIP/pre_calls/src/collapse_duplicate_tags_for_cims.py bed_uncollapsed/fbf1_oo_lane2_rt6.bed bed_collapsed/ > fbf1_oo_lane2_rt6.bed.log 2>&1 < /dev/null &
nohup python /groups/Kimble/Common/fog_iCLIP/pre_calls/src/collapse_duplicate_tags_for_cims.py bed_uncollapsed/fbf1_oo_lane2_rt9.bed bed_collapsed/ > fbf1_oo_lane2_rt9.bed.log 2>&1 < /dev/null &
nohup python /groups/Kimble/Common/fog_iCLIP/pre_calls/src/collapse_duplicate_tags_for_cims.py bed_uncollapsed/fbf2_oo_lane2_rt2.bed bed_collapsed/ > fbf2_oo_lane2_rt2.bed.log 2>&1 < /dev/null &
nohup python /groups/Kimble/Common/fog_iCLIP/pre_calls/src/collapse_duplicate_tags_for_cims.py bed_uncollapsed/fbf2_sp_lane1_rt13.bed bed_collapsed/ > fbf2_sp_lane1_rt13.bed.log 2>&1 < /dev/null &
nohup python /groups/Kimble/Common/fog_iCLIP/pre_calls/src/collapse_duplicate_tags_for_cims.py bed_uncollapsed/fbf2_oo_lane2_rt13.bed bed_collapsed/ > fbf2_oo_lane2_rt13.bed.log 2>&1 < /dev/null &

nohup python /groups/Kimble/Common/fog_iCLIP/pre_calls/src/collapse_duplicate_tags_for_cims.py bed_uncollapsed/fbf2_sp_lane1_rt14.bed bed_collapsed/ > fbf2_sp_lane1_rt14.bed.log 2>&1 < /dev/null &
nohup python /groups/Kimble/Common/fog_iCLIP/pre_calls/src/collapse_duplicate_tags_for_cims.py bed_uncollapsed/n2_oo_lane1_rt16.bed bed_collapsed/ > n2_oo_lane1_rt16.bed.log 2>&1 < /dev/null &
nohup python /groups/Kimble/Common/fog_iCLIP/pre_calls/src/collapse_duplicate_tags_for_cims.py bed_uncollapsed/fbf1_oo_lane2_rt1.bed bed_collapsed/ > fbf1_oo_lane2_rt1.bed.log 2>&1 < /dev/null &
nohup python /groups/Kimble/Common/fog_iCLIP/pre_calls/src/collapse_duplicate_tags_for_cims.py bed_uncollapsed/fbf1_sp_lane1_rt1.bed bed_collapsed/ > fbf1_sp_lane1_rt1.bed.log 2>&1 < /dev/null &
nohup python /groups/Kimble/Common/fog_iCLIP/pre_calls/src/collapse_duplicate_tags_for_cims.py bed_uncollapsed/fbf2_oo_lane2_rt11.bed bed_collapsed/ > fbf2_oo_lane2_rt11.bed.log 2>&1 < /dev/null &
nohup python /groups/Kimble/Common/fog_iCLIP/pre_calls/src/collapse_duplicate_tags_for_cims.py bed_uncollapsed/fbf1_sp_lane1_rt6.bed bed_collapsed/ > fbf1_sp_lane1_rt6.bed.log 2>&1 < /dev/null &

nohup python /groups/Kimble/Common/fog_iCLIP/pre_calls/src/collapse_duplicate_tags_for_cims.py bed_uncollapsed/fbf1_sp_lane1_rt9.bed bed_collapsed/ > fbf1_sp_lane1_rt9.bed.log 2>&1 < /dev/null &
nohup python /groups/Kimble/Common/fog_iCLIP/pre_calls/src/collapse_duplicate_tags_for_cims.py bed_uncollapsed/fbf2_sp_lane1_rt2.bed bed_collapsed/ > fbf2_sp_lane1_rt2.bed.log 2>&1 < /dev/null &
nohup python /groups/Kimble/Common/fog_iCLIP/pre_calls/src/collapse_duplicate_tags_for_cims.py bed_uncollapsed/n2_oo_lane1_rt15.bed bed_collapsed/ > n2_oo_lane1_rt15.bed.log 2>&1 < /dev/null &

'''
if __name__ == '__main__':
    if len(sys.argv) != 3:
        print "Run as " \
              "python bed_to_wig.py <input_dir> <output_dir>"
        sys.exit()
    indir = sys.argv[1]
    outdir = sys.argv[2]
    run(indir, outdir)

