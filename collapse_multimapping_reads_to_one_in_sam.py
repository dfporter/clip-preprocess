import collections
import sys
import os
import re
import glob


def get_as_val(line):
    s = line.split('\t')
    as_val = re.search('AS:i:(\d+)\s', line.rstrip('\n'))
    return as_val.group(1), int(s[4])  # STAR score, MAPQ


def get_max_map_val(map_list):
    name_to_as = dict([(x, get_as_val(x)) for x in map_list])
    ranked_AS_vals = sorted(name_to_as.values(), key=lambda x: x[0])
    ranked_mapq_vals = sorted(name_to_as.values(), key=lambda x: x[1])
    if ranked_AS_vals[0][0] != ranked_AS_vals[-1][0]:
        best_map = sorted(map_list, key=lambda x: name_to_as[x][0])[-1]
    elif ranked_mapq_vals[0][1] != ranked_mapq_vals[-1][1]:
        best_map = sorted(map_list, key=lambda x: name_to_as[x][1])[-1]
    else:
        best_map = map_list[0]
    return best_map


def process(samfile, outfile):
    header = ''
    by_name = collections.defaultdict(list)
    with open(samfile) as f:
        for li in f:
            if li[0] == '#' or li[0] == '@':
                header += li
                continue
            s = li.split('\t')
            by_name[s[0]].append(li)
    uniques = collections.defaultdict(str)
    for name in by_name:
        if len(by_name[name]) == 1:
            uniques[name] = by_name[name][0]
        else:
            uniques[name] = get_max_map_val(by_name[name])
    open(outfile, 'w').write(
        header + "".join(uniques.values())
    )


def run(samfile_or_folder, out_dir='sam_collapse_multi_map/'):
    if not os.path.exists(out_dir): os.system('mkdir ' + out_dir)
    if os.path.isdir(samfile_or_folder):
        for samfile in glob.glob(samfile_or_folder +'/*.sam'):
            outsam = out_dir + '/' + os.path.basename(samfile)
            process(samfile, outsam)
    elif os.path.isfile(samfile_or_folder):
        outsam = out_dir + '/' + os.path.basename(samfile_or_folder)
        process(samfile_or_folder, outsam)
    else:
        print "Input was not a sam file or folder."


if __name__ == '__main__':
    run(sys.argv[1])
