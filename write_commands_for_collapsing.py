#!/usr/bin/python
import glob
import os
import sys

"""
export PERL5LIB=/groups/Kimble/Aman\ Prasad/clip/plib/

"""
usage = "USAGE: <input dir> <output dir"
print usage
in_dir = sys.argv[1]
if len(sys.argv) < 3:
    out_dir = 'bed_collapsed/rrna'
    print "Deafults to output dir {0}. Is this OK?\n".format(
        out_dir)
    if raw_input('Y/N:').upper()[0] != 'Y':
        sys.exit()
else:
    print "Input dir: {0}. Output dir: {1}".format(
        sys.argv[1], sys.argv[2])
    out_dir = sys.argv[2]

print "export PERL5LIB=/groups/Kimble/Aman\ Prasad/clip/plib/"
if not os.path.exists('bed_collapsed'): os.system('mkdir bed_collapsed')
if not os.path.exists(out_dir):
    os.system('mkdir ' + out_dir)
for fname in glob.glob(in_dir + '/*.bed'):
    outname = out_dir + '/' + os.path.basename(fname)
    logname = os.path.basename(fname) + '.log'
    n = 2
    while os.path.exists(logname):
         logname = logname.partition('.log')[0] + '.log{0}'.format(n)
         if n > 100:
              print "Too many log files for {0} already exist in this directory...".format(
                  fname)
              sys.exit()
    if os.path.exists(outname):
        continue
    cmd = 'nohup python {0}/collapse_duplicate_tags_for_cims.py'.format(
        os.path.dirname(os.path.realpath(__file__)))
    cmd += ' ' + fname
    cmd += ' {od}/ > {v} 2>&1 < /dev/null &'.format(
        od=out_dir, v=logname)
    print cmd

