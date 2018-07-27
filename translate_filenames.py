	import re, argparse

def translate_fnames(in_folder, out_folder, suffix):

	table = {
	'SRR3715251': 'exp_fbf1_sp_1',
	'SRR3715252': 'exp_fbf1_sp_2',
	'SRR3715253': 'exp_fbf1_sp_3',
	'SRR3715254': 'exp_fbf2_sp_1',
	'SRR3715255': 'exp_fbf2_sp_2',
	'SRR3715256': 'exp_fbf2_sp_3',
	'SRR3715257': 'exp_fbf1_oo_1',
	'SRR3715258': 'exp_fbf1_oo_2',
	'SRR3715259': 'exp_fbf1_oo_3',
	'SRR3715260': 'exp_fbf2_oo_1',
	'SRR3715261': 'exp_fbf2_oo_2',
	'SRR3715262': 'exp_fbf2_oo_3',
	'SRR3715263': 'control_oo_1',
	'SRR3715264': 'control_oo_2', 
	'SRR3715265': 'control_oo_3',
	'SRR3715266': 'control_sp_1',
	'SRR3715267': 'control_sp_2',
	'SRR3715268': 'control_sp_3',
	}


	for fname in glob.glob(in_folder + '/*' + suffix):
		for srr_string in table:
			if srr_string in fname:
				out_fname = os.path.basename(
					re.sub(srr_string, table[srr_string], fname))
				cmd = 'rsync {} {}/{}'.format(fname, out_folder, out_fname)
				print(cmd)
				os.system(cmd)


parser = argparse.ArgumentParser()

parser.add_argument(
    '-i', '--input',
    help='Input folder.')
parser.add_argument(
    '-o', '--output', default='',
    help='Output folder.')
parser.add_argument(
    '-s', '--suffix', default='',
    help='Suffix (optional).')

args = parser.parse_args()

if args.output == '':
	args.output = os.path.dirname(args.input) + '_name_change/'
	print("Using default output directory {}.".format(args.output))

if not os.path.exists(args.output):
	os.system('mkdir {}'.format(args.output))
translate_fnames(args.input, args.output, args.suffix)