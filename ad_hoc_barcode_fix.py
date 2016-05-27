import HTSeq
import glob
import sys
import os
def process_file(filename, output_filename):
    outf = open(output_filename, 'w')
    fastq = HTSeq.FastqReader(filename)
    for read in fastq:
        barcode = read.name[-9:]
        new_barcode = barcode[:3] + barcode[7:]
        n_read = HTSeq.SequenceWithQualities(
            read.seq,
            read.name[:-9] + new_barcode,
            read.qualstr)
        n_read.write_to_fastq_file(outf)
    outf.close()


def process_sam(filename, output_filename):
    outf = open(output_filename, 'w')
    sam = open(filename)
    outli = ''
    for read in sam:
        if read[0] == '#': outli += read
        s = read.split('/t')
        barcode = s[0][-9:]
        new_barcode = barcode[:3] + barcode[7:]
        s[0] = s[0][:-9] + new_barcode
        outli += "\t".join(s)
    outf.write(outli)
    outf.close()



input_dir = sys.argv[1]
try:
    output_dir = sys.argv[2]
except:
    output_dir = 'fixed_barcodes'
if not os.path.exists(output_dir):
    os.system('mkdir ' + output_dir)

for filename in glob.glob(input_dir + '/*.fastq'):
    process_file(filename, output_dir + '/' + os.path.basename(filename))
for filename in glob.glob(input_dir + '/*.sam'):
    process_sam(filename, output_dir + '/' + os.path.basename(filename))
