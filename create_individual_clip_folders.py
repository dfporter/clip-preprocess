import glob, os, sys, collections, re, argparse, HTSeq

print(__file__)
sys.path.append(os.path.dirname(__file__))
import bed_to_wig

def _mk(dirname, **kwargs):

	if not os.path.exists(dirname):

		if kwargs['run_commands']:
			print('mkdir {}'.format(dirname))
			os.system('mkdir {}'.format(dirname))

		else:
			print('mkdir {}'.format(dirname))

	else:
		print("Output directory exists: {}".format(dirname))


def scan_filenames(top_dir):
	sp_fnames = [x for x in glob.glob(top_dir + '/*') if re.search('sp', os.path.basename(x).lower())]
	oo_fnames = [x for x in glob.glob(top_dir + '/*') if re.search('oo', os.path.basename(x).lower())]
#	print("SP filenames: {}".format(sp_fnames))
#	print("OO filenames: {}".format(oo_fnames))
	fbf1_fnames = [x for x in glob.glob(top_dir + '/*') if re.search('fbf1', os.path.basename(x).lower())]
	fbf2_fnames = [x for x in glob.glob(top_dir + '/*') if re.search('fbf2', os.path.basename(x).lower())]

	control_fnames = [x for x in glob.glob(top_dir + '/*') if re.search('control', os.path.basename(x).lower())]

	print("Looking for files in {}...".format(top_dir))
	individual_control_exps = [x for x in control_fnames if re.search('_\d.\w\w\w$', x)]

	sp_fnames = set(sp_fnames)
	oo_fnames = set(oo_fnames)
	fbf1_fnames = set(fbf1_fnames)
	fbf2_fnames = set(fbf2_fnames)
	control_fnames = set(control_fnames)
	individual_control_exps = set(individual_control_exps)

	print("Found {}".format(individual_control_exps))

	return sp_fnames, oo_fnames, fbf1_fnames, fbf2_fnames, control_fnames, individual_control_exps 


def create_folders(top_dir, out_top_dir, **kwargs):

	print("Creating folders...")
	print(' * ~ ' * 10)	

	for datatype_dirname in ['bed_collapsed/', 'bedgraph_unnorm/', 'bedgraph_norm/']:
		dir_in = '{}/{}'.format(top_dir, datatype_dirname)

		sp_fnames, oo_fnames, fbf1_fnames, fbf2_fnames, control_fnames, individual_control_exps = scan_filenames(dir_in)

		if datatype_dirname == 'bed_collapsed/':
			if concatenate_controls(sp_fnames, oo_fnames, individual_control_exps,
			 dir_in, '{}/'.format(top_dir), **kwargs):
				sp_fnames, oo_fnames, fbf1_fnames, fbf2_fnames, control_fnames, individual_control_exps = scan_filenames(dir_in)

		sp_fbf1 = (sp_fnames & fbf1_fnames) | (control_fnames & sp_fnames)
		sp_fbf2 = (sp_fnames & fbf2_fnames) | (control_fnames & sp_fnames)
		oo_fbf1 = (oo_fnames & fbf1_fnames) | (control_fnames & oo_fnames)
		oo_fbf2 = (oo_fnames & fbf2_fnames) | (control_fnames & oo_fnames)

		_mk(out_top_dir, **kwargs)
		
		for sample_dirname in ['sp_fbf1/', 'sp_fbf2/', 'oo_fbf1/', 'oo_fbf2/']:
			_mk('{}/{}'.format(out_top_dir, sample_dirname), **kwargs)
			_mk('{}/{}/{}'.format(out_top_dir, sample_dirname, datatype_dirname), **kwargs)

		for fnames, sample_dirname in zip(
			[sp_fbf1, sp_fbf2, oo_fbf1, oo_fbf2],
			['/sp_fbf1', '/sp_fbf2', '/oo_fbf1', '/oo_fbf2']):

			out_subdir = '{}/{}/{}'.format(out_top_dir, sample_dirname, datatype_dirname)

			copy_over(fnames, out_subdir, **kwargs)

#	_mk(top_dir + '/sp_fbf2/')
#	copy_over(sp_fbf2 )


def copy_over(fname_list, target_dir, **kwargs):

	print('*' * 14)
	print("Copying from {} to {}".format(os.path.dirname(list(fname_list)[0]), target_dir))

	for fname in fname_list:
		cmd = 'rsync {} {}/{}'.format(
			fname, target_dir, os.path.basename(fname))

		print(cmd)

		if kwargs['run_commands']:
			os.system(cmd)


def concatenate_controls(sp_fnames, oo_fnames, individual_control_exps, top_dir, bedgraphs_unnorm_dir, **kwargs):
	individual_control_exps_sp = sp_fnames & individual_control_exps
	combined_sp = '{}/control_sp.bed'.format(top_dir)

	files_already_existing = 0

	if os.path.exists(combined_sp):
		files_already_existing += 1
	else:
		cmd = 'cat ' + ' '.join([x for x in individual_control_exps_sp])
		cmd += ' > {}'.format(combined_sp)
		print(cmd)

		if kwargs['run_commands']:
			os.system(cmd)

	individual_control_exps_oo = oo_fnames & individual_control_exps
	combined_oo = '{}/control_oo.bed'.format(top_dir)

	if os.path.exists(combined_oo):
		files_already_existing += 1
	else:
		cmd = 'cat ' + ' '.join([x for x in individual_control_exps_oo])
		cmd += ' > {}'.format(combined_oo)
		print(cmd)

		if kwargs['run_commands']:
			os.system(cmd)

	make_bedgraphs_for_concatenated_controls(combined_oo, bedgraphs_unnorm_dir)

	if files_already_existing == 2:
		return False

	return True

def make_bedgraphs_for_concatenated_controls(infile, bedgraph_unnorm_dir):
	global_ga = HTSeq.GenomicArray('auto', stranded=True)

	ga = bed_to_wig.add_to_ga(infile, global_ga)

	for subdir in ['/sp_fbf1', '/sp_fbf2', '/oo_fbf1', '/oo_fbf2']:
		if not os.path.exists('{}/{}'.format(bedgraph_unnorm_dir, subdir)):
			os.system('mkdir {}/{}'.format(bedgraph_unnorm_dir, subdir))
		if not os.path.exists('{}/{}/bedgraph_unnorm'.format(bedgraph_unnorm_dir, subdir)):
			os.system('mkdir {}/{}/bedgraph_unnorm'.format(bedgraph_unnorm_dir, subdir))
			
		ga.write_bedgraph_file(
	        bedgraph_unnorm_dir + '{}/bedgraph_unnorm/{}_+.wig'.format(subdir, os.path.basename(infile).split('.')[0]), strand='+')
		ga.write_bedgraph_file(
	        bedgraph_unnorm_dir + '{}/bedgraph_unnorm/{}_-.wig'.format(subdir, os.path.basename(infile).split('.')[0]), strand='-')

parser = argparse.ArgumentParser()

parser.add_argument(
    '-i', '--input',
    help='Input folder.')
parser.add_argument(
    '-o', '--output', default='',
    help='Output folder.')
parser.add_argument(
    '-r', '--run_commands', action='store_true', default=False,
    help='Run commands, rather than just print them.')

args = parser.parse_args()

create_folders(args.input, args.output, run_commands=args.run_commands)