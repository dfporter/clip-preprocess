import glob
import sys
import os


def normalize_wig(_folder, out_dir):
    if not os.path.exists(out_dir): os.system('mkdir ' + out_dir)
    for bed in glob.glob(_folder + '/*.wig'):
        out_bed = out_dir + '/' + os.path.basename(bed)
        if os.path.exists(out_bed):
            print "{A} file exists, skipping..".format(A=out_bed)
            continue
        total_reads = float(sum([1 for line in open(bed)]) - 1)
        bed_fh = open(bed, 'r')
        outli = next(bed_fh)  # Skip header.
        for line in bed_fh:
            s = line.rstrip('\n').split('\t')
            outli += '\t'.join(s[:3]) + str(1e6 * float(s[3])/total_reads)
            outli += '\n'
        open(out_bed, 'w').write(outli)        

if __name__ == "__main__":
    print "Usage: python normalize_wig.py <in_dir> <out_dir>"
    normalize_wig(sys.argv[1], sys.argv[2])


