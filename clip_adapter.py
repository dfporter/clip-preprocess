import os
import glob
import sys


def scribe(cmd):
    print cmd
    os.system(cmd)


print "Usage: python clip_adapter.py input_directory/ ouput_directory/"
if len(sys.argv) != 3:
    print "Expected two command line arguments (input dir and output dir)."
    raise IOError
if sys.argv[1] == '.' or sys.argv[2] == '.':
    print "Can't pass the current directory (lazy)."
    raise IOError
if not os.path.exists(sys.argv[2]):
    cmd = 'mkdir ' + sys.argv[2]
    scribe(cmd)
for infile in glob.glob(sys.argv[1] + '/*.fastq'):
    print infile
    outfile = '{c}/{b}'.format(
        b=os.path.basename(infile),
        c=sys.argv[2])
    if os.path.exists(outfile):
        print "Output file {a} exists... Aborting...".format(a=outfile)
        raise IOError
    cmd = "fastx_clipper -a AGATCGGAAGAGCGGTTCAGCAGGAATGCC "
    cmd += "-i {a} -o {b} -Q 33".format(
        a=infile, b=outfile)
    scribe(cmd)
