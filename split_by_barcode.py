import re
import sys
import HTSeq

"""
FOG-3:

    for bc in ['GGTT', 'TGGC', 'CGGA', 'GGCA']:
        outfiles[bc] = open('exp_%s.fastq' % bc, 'w')
    for bc in ['TTGT', 'CCGG', 'GGCA', 'TTAA']:
        outfiles[bc] = open('control_%s.fastq' % bc, 'w')
"""

def split_by_barcode(initial_filename):
    outfiles = {}
    outfiles['CAAT'] = open('n2_sp_lane3_rt3.fastq', 'w')
    outfiles['AATA'] = open('n2_sp_lane3_rt15.fastq', 'w')
    outfiles['TTAA'] = open('n2_sp_lane3_rt16.fastq', 'w')
    missingf = open('no_recognized_barcode_lane3.fastq', 'w')
    skip = """
    # Lane 1:
    outfiles['GGTT'] = open('fbf1_sp_lane1_rt1.fastq', 'w')
    outfiles['TTGT'] = open('fbf2_sp_lane1_rt2.fastq', 'w')
    outfiles['CAAT'] = open('n2_oo_lane1_rt3.fastq', 'w')
    outfiles['CCGG'] = open('fbf1_sp_lane1_rt6.fastq', 'w')
    outfiles['TGGC'] = open('fbf1_sp_lane1_rt9.fastq', 'w')
    outfiles['CGGA'] = open('fbf2_sp_lane1_rt13.fastq', 'w')
    outfiles['GGCA'] = open('fbf2_sp_lane1_rt14.fastq', 'w')
    outfiles['AATA'] = open('n2_oo_lane1_rt15.fastq', 'w')
    outfiles['TTAA'] = open('n2_oo_lane1_rt16.fastq', 'w')
    missingf = open('no_recognized_barcode_lane1.fastq', 'w')
    # Lane 2:
    outfiles['GGTT'] = open('fbf1_oo_lane2_rt1.fastq', 'w')
    outfiles['TTGT'] = open('fbf2_oo_lane2_rt2.fastq', 'w')
    outfiles['CCGG'] = open('fbf1_oo_lane2_rt6.fastq', 'w')
    outfiles['TGGC'] = open('fbf1_oo_lane2_rt9.fastq', 'w')
    outfiles['CGGA'] = open('fbf2_oo_lane2_rt13.fastq', 'w')
    outfiles['GGCA'] = open('fbf2_oo_lane2_rt11.fastq', 'w')
    missingf = open('no_recognized_barcode_lane2.fastq', 'w')
    outfiles['TGGC'] = open('exp_fbf1_TGGC.fastq', 'w')
    outfiles['CGGA'] = open('exp_fbf1_CGGA.fastq', 'w')
    # There is an irregularity here. The GEO dataset indicates a GCCA/TCCG barcode
    # (the rc of the above) and a GGTT barcode (so we expect AACC).
    # What we use in the sp/oo is GGTT, though.
    outfiles['AACC'] = open('exp_fbf1_AACC.fastq', 'w')
    outfiles['GGTT'] = open('exp_fbf1_GGTT.fastq', 'w')

    outfiles['CCGG'] = open('fbf1_n2_CCGG.fastq', 'w')
    outfiles['TTGT'] = open('fbf1_n2_TTGT.fastq', 'w')
    outfiles['GGCA'] = open('fbf1_n2_GGCA.fastq', 'w')
    missingf = open('no_recognized_barcode_fbf1.fastq', 'w')

"""
    skip = """
    outfiles['TGGC'] = open('exp_fbf2_TGGC.fastq', 'w')
    outfiles['CGGA'] = open('exp_fbf2_CGGA.fastq', 'w')
    outfiles['AACC'] = open('exp_fbf2_AACC.fastq', 'w')
    outfiles['GGTT'] = open('exp_fbf2_GGTT.fastq', 'w')

    outfiles['CCGG'] = open('fbf2_n2_CCGG.fastq', 'w')
    outfiles['TTGT'] = open('fbf2_n2_TTGT.fastq', 'w')
    outfiles['GGCA'] = open('fbf2_n2_GGCA.fastq', 'w')
"""
    fastq_file = HTSeq.FastqReader(initial_filename)
    total_reads = 0
    for read in fastq_file:
        total_reads += 1
        if(not (total_reads % 100000)):
            print "Read: %i " % (total_reads)
        found = False
        for bc in outfiles.keys():
            if read.seq[3:7] == bc:
            # if(re.match('\w{3}' + bc + '.*', read.seq)):
                read.write_to_fastq_file(outfiles[bc])
                found = True
        if not found:
            read.write_to_fastq_file(missingf)
    for bc in outfiles:
        outfiles[bc].close()


import collections
def count_barcodes(fname):
    fastq_file = HTSeq.FastqReader(initial_filename)
    barcode_counts = collections.defaultdict(int)
    barcode_reg_counts = collections.defaultdict(int)
    for read in fastq_file:
        barcode_counts[read.seq[3:7]] += 1
        barcode_reg_counts[read.seq[:8]] += 1
    with open('barcode_counts.txt', 'w') as f:
        for barcode in barcode_counts:
            f.write("{b}\t{c}\n".format(b=barcode, c=str(barcode_counts[barcode])))
    with open('barcode_region_counts.txt', 'w') as f:
        for barcode in barcode_reg_counts:
            f.write("{b}\t{c}\n".format(b=barcode, c=str(barcode_reg_counts[barcode])))


def count_diversity(bcs):
    pass


if __name__ == '__main__':
    initial_filename = sys.argv[1]
    split_by_barcode(initial_filename)
    count_barcodes(initial_filename)
