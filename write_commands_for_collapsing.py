
import glob
import os


"""
export PERL5LIB=/groups/Kimble/Aman\ Prasad/clip/plib/

"""
if not os.path.exists('bed_collapsed'): os.system('mkdir bed_collapsed/')
for fname in glob.glob('bed_uncollapsed/*.bed'):
    outname = 'bed_collapsed/' + os.path.basename(fname)
    if os.path.exists(outname):
        continue
    cmd = 'nohup python /groups/Kimble/Common/fog_iCLIP/pre_calls/src/collapse_duplicate_tags_for_cims.py'
    cmd += ' ' + fname
    cmd += ' bed_collapsed/ > {v} 2>&1 < /dev/null &'.format(
        v=os.path.basename(fname) + '.log')
    print cmd


