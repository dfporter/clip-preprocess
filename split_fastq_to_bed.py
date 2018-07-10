import os
import sys
import argparse
import glob
import HTSeq
import clip_adapter
import map_with_star
import move_barcode_to_name_in_fastq


def scribe(cmd):
    print(cmd)
    os.system(cmd)


def run(a_dir, lib=None):
    num_fastq = len(glob.glob(a_dir + '/*.fastq'))
    paths = {
        'star': '/opt/STAR-STAR_2.4.2a/bin/Linux_x86_64/STAR',
        'indexes': '/opt/STAR-STAR_2.4.2a/bin/Linux_x86_64/indexes/',
        'gtf': '/opt/lib/Caenorhabditis_elegans.WBcel235.78.noheader.gtf'}
    if lib is None:
        print('Using default paths:')
    else:
        paths = lib
    all_ok = True
    for path in list(paths.values()):
        if not os.path.exists(path):
            print("Missing %s" % path)
            all_ok = False
    if not all_ok:
        if input("Will not be able to run STAR. Continue anyway [y/n]?"
                     )[0].upper() != "Y":
            return False
    if num_fastq > 0:
        out_dir = 'temp_adapter_in_name/'
        move_barcode_to_name_if_not_present(a_dir, out_dir=out_dir)
        clip_adapters_if_not_already_clipped(
            out_dir, 'temp_clipped/', args)
        args.input_dir = 'temp_clipped/'
        map_with_star.run(args, paths=paths)
        scribe('python write_commands_for_collapsing.py')

# a_dir -> temp_adapter_in_name/ -> temp_fastq/ -> temp_clipped/ -> sams/
# -> temp_sam_collapse_multi_map/ -> sams/

def clear_directory_of_output_files(_dir):
    for subdir in [_dir + '/' + x for x in [
        'temp_adapter_in_name', 'temp_fastq',
        'temp_clipped', 'sams/', 'temp_sam_collapse_multi_map/']]:
       scribe('rm -r' + subdir)


def move_barcode_to_name_if_not_present(a_dir, out_dir='adapter_in_name'):
    for fname in glob.glob(a_dir + '*.fastq'):
        first_line = open(fname).readline()
        if first_line[0] != '@':
            print("Error: fastq file does not start with @...")
            continue
        if '#' in first_line: continue
        move_barcode_to_name_in_fastq.move_barcode_to_name_in_fastq(
            fname, out_dir)


def clip_adapters_if_not_already_clipped(in_dir, out_dir, args):
    n_lines_to_check = 4e3
    need_to_clip = True
    for _file in glob.glob(in_dir + '/*.fastq'):
        fastq_reader = HTSeq.FastqReader(_file)
        read_lens = set()
        for i, read in enumerate(fastq_reader):
            read_lens.add(len(read.seq))
            if i > n_lines_to_check: break
        if len([x for x in read_lens if x>20]) > 2:
            need_to_clip = False
            break
    if not need_to_clip:
        print("Adapters in {0} are apparently already clipped...".format(in_dir))
        return
    print("Clipping adapters in {0}...".format(in_dir))
    if in_dir == out_dir:
        print("Input and output dir can't be the same.")
        sys.exit()
    args.input_dir = in_dir
    args.output_dir = 'temp_fastq/'
    if os.path.exists('temp_fastq/'): os.system('rm -r temp_fastq/')
    args.adapter = ''
    args.three_prime_linker = True
    args.rt_primer = False
    clip_adapter.run(args)
    for k, v in list({'three_prime_linker': False, 'rt_primer': True,
    'input_dir': 'temp_fastq/', 'output_dir': out_dir, 'adapter': ''}.items()):
        setattr(args, k, v)
    clip_adapter.run(args)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input',
                        help='Input folder of fastq files,\
sam files or bed files')
    args = parser.parse_args()
    run(args.input)
