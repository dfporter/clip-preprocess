#!/usr/bin/python
import glob
import sys
import os
import re

'''
export PERL5LIB=/groups/Kimble/Aman\ Prasad/clip/plib/
'''

def run(input_filename, out_dir):
    """Need to export the plib/ directory to PERL5LIB before running.
    export PERL5LIB=/groups/Kimble/Aman\ Prasad/clip/plib/
    export PERL5LIB=/Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/FBF_gendered_gl/plib/
    """ 

    output_filename = out_dir + '/'
    output_filename += os.path.basename(input_filename)

    cims_dir = '/groups/Kimble/Aman\ Prasad/clip/CIMS/'
    cims_dir = '//Users/dfporter/Desktop/macbook_air_Desktop/shared/sp_oo/cims/CIMS/'
    # if os.path.exists(output_filename): return
    cmd = 'perl {c}/tag2collapse.pl --verbose --random-barcode --seq-error-model fix=0.02 -EM 50 --keep-cache {i} {o}'.format(
        c=cims_dir, i=input_filename, o=output_filename)
    #cmd = 'perl ../clip/CIMS/tag2collapse.pl --random-barcode --seq-error-model fix=0.02 -EM 50 --keep-cache {i} {o}'.format(
    #    i=input_filename, o=output_filename)

    print(cmd)
    os.system(cmd)


if __name__ == '__main__':

    input_filename = sys.argv[1]
    out_dir = sys.argv[2]

    if input_filename[-1] == '/':

        for fname in glob.glob(input_filename + '/*.bed'):
            run(fname, out_dir)

    else:
        run(input_filename, out_dir)
