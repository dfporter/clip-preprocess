import os
import glob
import sys
import argparse

def scribe(cmd):
    print(cmd)
    os.system(cmd)

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_dir')
    parser.add_argument('-o', '--output_dir')
    parser.add_argument('-a', '--adapter', default='')
    parser.add_argument('-t', '--three_prime_linker', default=False,
                        action='store_true')
    parser.add_argument('-r', '--rt_primer', default=False,
                        action='store_true')
    args = parser.parse_args()
    return args

def init(args):
    if (args.adapter != '') and (args.three_prime_linker or args.rt_primer):
        print("Only one of the arguments -a, -t or -r can be used.")
        sys.exit()
    if args.three_prime_linker and args.rt_primer:
        print("Only one of the arguments -a, -t or -r can be used.")
        sys.exit()
    if not(args.adapter or args.three_prime_linker or args.rt_primer):
        print("One of the arguments -a, -t or -r must be used.")
        sys.exit()
    if args.three_prime_linker:
        args.adapter = 'TGAGATCGGAAGAGCGGTTCAG'
#        args.adapter = 'AGATCGGAAGAGCGGTTCAG'
    if args.rt_primer:
        args.adapter = 'AGATCGGAAGAGCGGTTCAGCAGGAATGCC'
    if args.input_dir == '.' or args.output_dir == '.':
        print("Can't pass the current directory (lazy).")
        raise IOError
    args.adapter = str(args.adapter)
    return args

def run(args):
    args = init(args)
    if not os.path.exists(args.output_dir):
        cmd = 'mkdir ' + args.output_dir
        scribe(cmd)
    for infile in glob.glob(args.input_dir + '/*.fastq'):
        print(infile)
        outfile = '{0}/{1}'.format(args.output_dir, os.path.basename(infile))
        if os.path.exists(outfile):
            print("Output file {a} exists... Skipping...".format(a=outfile))
            continue
        cmd = "fastx_clipper -n -a {seq} -i {a} -o {b} -Q 33".format(
            seq=args.adapter,
            a=infile, b=outfile)
        scribe(cmd)
    #    cmd = "fastx_clipper -a AGATCGGAAGAGCGGTTCAGCAGGAATGCC"
    #    cmd += " -i {a} -o {b} -Q 33".format(
    #        a=infile, b=outfile)
    #    scribe(cmd)


if __name__ == '__main__':
    args = get_args()
    run(args)
