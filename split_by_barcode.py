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
    # Lane 1:
    skip = """
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
"""
    # Lane 2:
    outfiles['GGTT'] = open('fbf1_oo_lane2_rt1.fastq', 'w')
    outfiles['TTGT'] = open('fbf2_oo_lane2_rt2.fastq', 'w')
    outfiles['CCGG'] = open('fbf1_oo_lane2_rt6.fastq', 'w')
    outfiles['TGGC'] = open('fbf1_oo_lane2_rt9.fastq', 'w')
    outfiles['CGGA'] = open('fbf2_oo_lane2_rt13.fastq', 'w')
    outfiles['GGCA'] = open('fbf2_oo_lane2_rt11.fastq', 'w')
    missingf = open('no_recognized_barcode_lane2.fastq', 'w')
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

def count_barcodes(fname):
    fastq_file = HTSeq.FastqReader(initial_filename)
    barcode_counts = {}
    barcode_reg_counts = {}
    for read in fastq_file:
        barcode_reg = read.seq[:8]
        barcode = read.seq[3:7]
        barcode_counts.setdefault(barcode, 0)
        barcode_counts[barcode] += 1
        barcode_reg_counts.setdefault(barcode_reg, 0)
        barcode_reg_counts[barcode_reg] += 1
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
