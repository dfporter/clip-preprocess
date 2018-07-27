import re
import collections
import os
import sys
import argparse
import glob
import HTSeq


def get_auc(bedgraph_list):
	auc = collections.defaultdict(float)

	for fname in bedgraph_list:

		with open(fname) as fh:

			next(fh) # Skip header
			for li in fh:
				s = li.rstrip('\n').split('\t')
				auc[fname] += float(s[3]) * (int(s[2]) - int(s[1]))

		guess_negative = re.sub('_+', '_-', fname)

		if os.path.exists(guess_negative):
			print("Reading {} as a negative strand.".format(guess_negative))

			with open(guess_negative) as fh:
				next(fh) # Skip header
				for li in fh:
					s = li.rstrip('\n').split('\t')
					auc[fname] += float(s[3]) * (int(s[2]) - int(s[1]))

		print(fname, auc)

	return auc


parser = argparse.ArgumentParser()

parser.add_argument(
    '-i', '--input',
    help='Input folder of bedgraph files.')

args = parser.parse_args()

auc = get_auc(glob.glob('{}/*_+.wig'.format(args.input)))

for k, v in auc.items():
	print(k, v)

