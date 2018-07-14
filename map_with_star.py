#!/usr/bin/python
"""
Calls the STAR mapper and outputs to sams/.
Then calls samtools view -q 20 to export filtered bam files
to bams/.
"""

import glob
import sys
import os
import re
import collections
import argparse
import collapse_multimapping_reads_to_one_in_sam

def scribe(cmd):
    print(cmd)
    os.system(cmd)


def cleanup_star():
    os.system('mkdir idk')
    os.system('mv *out idk/')
    os.system('mv *tab idk/')
    os.system('mv *mate1 idk/')
    os.system('rm -r idk/*STARgenome')
    os.system('rm -r idk/*STARtmp')
    os.system('mv *STARgenome idk/')
    os.system('mv *STARtmp idk/')
    os.system('mv *.log idk/')


def count_gaps_in_sam(samfile):
    hist = collections.defaultdict(int)
    gap_pat = re.compile('(\d+)N')
    gap = gap_pat.finditer
    header_lines = 0
    with open(samfile, 'r') as f:
        in_header = True
        while in_header:
            for li in f:
                if re.match('\A@.*', li) is not None:
                    header_lines += 1
                else:
                    in_header = 0
    with open(samfile, 'r') as f:
        for n in range(header_lines):
            next(f)
        for li in f:
            s = li.rstrip('\n').split('\t')
            # cigar = s[5]
            m = gap(s[5])
            if m is not None:
                num = 0
                for ins in m:
                    hist[min([1000,int(ins.group(1))])] += 1
                    num += 1
                    # print s[5]
                    # print ins.group(1)
                if num == 0:
                    hist[0] += 1
            else:
                hist[0] += 1
    keys = sorted(list(hist.keys()), key=lambda x: hist[x])
    for k in keys:
        print('gap\t{k}:\t{n}'.format(k=k, n=hist[k]))


def convert_sam_to_bed(samname):
    cmd = 'sam2bed < {i} > {d}/{b}'.format(
        i=samname, d='bed_uncollapsed/',
        b=os.path.basename(samname).partition('.sam')[0] + '.bed'
    )
    scribe(cmd)


def mk_dirs():
    for _dir in ['bams/', 'sams/', 'bedgraphs/', 'bed_uncollapsed/']:
        if not os.path.exists(_dir):
            os.system('mkdir ' + _dir)
    

def run(args, paths=None):
    #        'star': '/groups/Kimble/Aman\ Prasad/clip/STAR-STAR_2.4.2a/bin/Linux_x86_64/STAR',
    if paths is None:
        paths={
#        'indexes': '/opt/STAR-STAR_2.4.2a/bin/Linux_x86_64/indexes/',
        'gtf': '/opt/lib/Caenorhabditis_elegans.WBcel235.78.noheader.gtf',
#        'star': '/Users/dfporter/Desktop/macbook_air_desktop/shared//star/STAR_2.4.2a/bin/MacOSX_x86_64/STAR',
        'star': '/Volumes/mond//STAR-STAR_2.4.2a/bin/MacOSX_x86_64/STAR',
#        'indexes': '/Volumes/Seagate/STAR-STAR_2.4.2a/bin/Linux_x86_64/indexes/',
        'indexes': '/Volumes/mond/STAR-STAR_2.4.2a/bin/Linux_x86_64/indexes/',

        'sjdb': '/opt/lib/sjdb.txt'}
#    paths['sjdb'] = '/opt/lib/sjdb.txt'
    create_splice_junctions_file(args, paths['gtf'], paths['sjdb'])

    paths['gtf'] = paths['sjdb']

    mk_dirs()

    call_star(args, paths=paths)

    split_multimapping_reads(
            in_dir='sams/', unique_dir='uniquely_mapping/',
            multi_dir='map_to_two/', unmapped_dir='unmapped/')

    os.system('mkdir uniquely_mapping_20_AS')

    for sam in glob.glob('uniquely_mapping/*.sam'):
        header = ''
        outli = ''

        for li in open(sam).readlines():

            if li[0] == '@':
                header += li
                continue

            s = li.rstrip('\n').split('\t')
            as_flag = re.match('.*AS:i:(\d+).*', str(s[11:]))

            if as_flag is not None:
                as_value = int(as_flag.group(1))
            else:
                print("No AS value?")

            if as_value < 20:
                continue

            outli += li

        open('uniquely_mapping_20_AS/{}'.format(os.path.basename(sam)), 'w').write(header + outli)

    if not os.path.exists('unfiltered_star_sams_output/'):
        os.system('mkdir unfiltered_star_sams_output/')

    os.system('rsync -r -v sams/ unfiltered_star_sams_output/')
    os.system('rsync -r -v uniquely_mapping_20_AS/ sams/')
    #collapse_multimapping_reads_to_one_in_sam.run('sams/', out_dir='temp_sam_collapse_multi_map/')
    #os.system('mv temp_sam_collapse_multi_map/* sams/')
    #os.system('rm -r temp_sam_collapse_multi_map/')
    sort_and_index_sam(in_dir='sams/')
    cleanup_star()


def create_splice_junctions_file(args, gtf_file, sjdb_file):
    outname = sjdb_file
    if os.path.exists(outname):
        print("Using splice junctions file " + outname)
        return
    else:
        print("Writing to splice junctions file " + outname)
    outli = ''
    with open(gtf_file) as f:
        for li in [x for x in f if x.split('\t')[2] == 'exon']:
            outli += li
    open(outname, 'w').write(outli)
    print("\t... Created splice juntions file in GFF format " + outname)


def call_star(args, paths=None):
    """
[dfporter@ad.wisc.edu@kimble-1 map_params_taken_from_sam_files_in_old_mapping_data_for_sp_oo]$ samtools view -H /mnt/Groups/KimbleGrp/General/Common/fbf_celltype/old_mapping_data_for_sp_oo/old/sams/n2_oo_lane1_rt15.sam 
@HD VN:1.4
@SQ SN:I    LN:15072434
@SQ SN:II   LN:15279421
@SQ SN:III  LN:13783801
@SQ SN:IV   LN:17493829
@SQ SN:V    LN:20924180
@SQ SN:X    LN:17718942
@SQ SN:MtDNA    LN:13794
@PG ID:STAR PN:STAR VN:STAR_2.4.2a  CL:/groups/Kimble/Aman Prasad/clip/STAR-STAR_2.4.2a/bin/Linux_x86_64/STAR   --runThreadN 5  
 --genomeDir "/groups/Kimble/Aman Prasad/clip/STAR-STAR_2.4.2a/bin/Linux_x86_64/indexes/"   --readFilesIn adapter_moved_to_name/n2_oo_lane1_rt15.fastq  
     --readMapNumber 18446744073709551615   --outFileNamePrefix n2_oo_lane1_rt15_   --outReadsUnmapped Fastx   --outSAMstrandField intronMotif  
      --outFilterScoreMin 20   --outFilterIntronMotifs RemoveNoncanonicalUnannotated   --alignIntronMax 1000   --alignSJoverhangMin 100  
       --alignSJDBoverhangMin 15   --sjdbGTFfile /groups/Kimble/Common/fog_iCLIP/calls/lib/Caenorhabditis_elegans.WBcel235.78.noheader.gtf
@CO user command line: /groups/Kimble/Aman Prasad/clip/STAR-STAR_2.4.2a/bin/Linux_x86_64/STAR --runThreadN 5 
--genomeDir "/groups/Kimble/Aman Prasad/clip/STAR-STAR_2.4.2a/bin/Linux_x86_64/indexes/" --alignIntronMax 1000 
--alignSJoverhangMin 100 --alignSJDBoverhangMin 15 --sjdbGTFfile /groups/Kimble/Common/fog_iCLIP/calls/lib/Caenorhabditis_elegans.WBcel235.78.noheader.gtf
 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outFilterScoreMin 20 --readMapNumber -1 --outReadsUnmapped Fastx 
 --outSAMstrandField intronMotif --readFilesIn adapter_moved_to_name/n2_oo_lane1_rt15.fastq --outFileNamePrefix n2_oo_lane1_rt15_


# Original Parameters
        cmd = '''STAR --runThreadN 8 \
        --genomeDir {indexes}\
         --alignIntronMax 1000 --alignSJoverhangMin 100 --alignSJDBoverhangMin 15\
         --sjdbGTFfile {gtf}\
         --outFilterIntronMotifs RemoveNoncanonicalUnannotated\
         --outFilterScoreMin 20\
         --outReadsUnmapped Fastx\
         --outSAMstrandField intronMotif\
         --readFilesIn {rin} --outFileNamePrefix {prefx}'''.format(
            rin=fastq_filename,
            prefx=bname + '_', indexes=paths['indexes'],
            gtf=paths['gtf'])''')

# CSEQ with --alignEndsType and --outFileterMultimapScoreRange changed
{star} --alignIntronMax 1 \
--sjdbGTFfile {sjdb} \
--genomeDir {indexes} \
--readFilesIn {rin} \
--outSAMunmapped "Within" \
--outFilterMultimapNmax 3 --outFilterMismatchNmax 2 \
--seedSearchStartLmax 6 --winAnchorMultimapNmax 10000 --alignEndsType Local \
--sjdbGTFtagExonParentTranscript transcript_id \
--outFilterMultimapScoreRange 0  --runThreadN 8 \
--outFileNamePrefix {prefix} \
--alignTranscriptsPerReadNmax 100000

"""
    in_dir = args.input_dir

    for fastq_filename in glob.glob(in_dir + '/*.fastq'):
        print(fastq_filename)
        bname = os.path.basename(fastq_filename).partition('.fastq')[0]
        # CSEQ parameters. http://psb.stanford.edu/psb-online/proceedings/psb16/kassuhn.pdf
        # Defaults changed:
        # --outFilterMultimapScoreRange is default 1
        # --alignEndsType EndToEnd for CSEQ.
        # We've added the --runThreadN for multithreading (faster speed, same output).
        #
        # Normally --runThreadN 8
CL:/groups/Kimble/Aman Prasad/clip/STAR-STAR_2.4.2a/bin/Linux_x86_64/STAR  
 --runThreadN 5   --genomeDir "/groups/Kimble/Aman Prasad/clip/STAR-STAR_2.4.2a/bin/Linux_x86_64/indexes/"  
  --readFilesIn adapter_moved_to_name/n2_oo_lane1_rt3.fastq    
    --readMapNumber 18446744073709551615   --outFileNamePrefix n2_oo_lane1_rt3_   
    --outReadsUnmapped Fastx   --outSAMstrandField intronMotif   --outFilterScoreMin 20 
     --outFilterIntronMotifs RemoveNoncanonicalUnannotated   --alignIntronMax 1000   
     --alignSJoverhangMin 100   --alignSJDBoverhangMin 15   
     --sjdbGTFfile /groups/Kimble/Common/fog_iCLIP/calls/lib/Caenorhabditis_elegans.WBcel235.78.noheader.gtf

        cmd = '''
{star} --alignIntronMax 1000
--sjdbGTFfile {sjdb} \
--genomeDir {indexes} \
--readFilesIn {rin} \
--alignSJoverhangMin 100 --alignSJDBoverhangMin 15 \
--outFilterIntronMotifs RemoveNoncanonicalUnannotated \
--outFilterScoreMin 20 --readMapNumber -1 \
--outReadsUnmapped Fastx --outSAMstrandField intronMotif \
--outFileNamePrefix {prefix} \
--runThreadN 16
'''.format(
        star=paths['star'], sjdb=paths['sjdb'],
        indexes=paths['indexes'], rin=fastq_filename, prefix=bname + '_')

    alt_cmd = '''
{star} --alignIntronMax 1 \
--sjdbGTFfile {sjdb} \
--genomeDir {indexes} \
--readFilesIn {rin} \
--outSAMunmapped "Within" \
--outFilterMultimapNmax 3 --outFilterMismatchNmax 2 \
--seedSearchStartLmax 6 --winAnchorMultimapNmax 10000 --alignEndsType Local \
--sjdbGTFtagExonParentTranscript transcript_id \
--outFilterMultimapScoreRange 0  --runThreadN 8 \
--outFileNamePrefix {prefix} \
--alignTranscriptsPerReadNmax 100000
'''.format(
        star=paths['star'], sjdb=paths['sjdb'],
        indexes=paths['indexes'], rin=fastq_filename, prefix=bname + '_')

        # Use --readMapNumber 10000 to only map the first 10,000 reads.
        outname = bname + '_Aligned.out.sam'
        fixed_name = outname.partition('_Aligned.out.sam')[0] + '.sam'
        if os.path.exists(outname) or os.path.exists('sams/{0}'.format(fixed_name)):
            print("{0} already exists, skipping...".format(outname))
            continue
        print(cmd)
        os.system(cmd)
        os.system('mv {i} {o}'.format(i=outname, o=fixed_name))
        os.system('mv {o} sams/{o}'.format(o=fixed_name))


"""
# CSEQ, except with EndToEnd changed to Local (allow soft-clipping,
# so we can get stuff with a polyA tail.)
STAR --alignIntronMax 1 --sjdbGTFfile annotation.gtf --outSAMunmapped "Within" \
--outFilterMultimapNmax 3 --outFilterMismatchNmax 2 --seedSearchStartLmax 6 \
--winAnchorMultimapNmax 10000 --alignEndsType Local

# FAST iCLIP parameters
        cmd = "STAR --genomeDir {indexes} --runThreadN 8 \
    --genomeLoad LoadAndKeep \
    --readFilesIn {reads} --outFileNamePrefix {prefix} --alignEndsType EndToEnd \
    --outFilterMismatchNoverLmax 0.08".format(
        star=paths['star'],
        indexes=paths['indexes'],
        #gtf=paths['gtf'],
        reads=fastq_filename,
        prefix=bname + '_',
        )
        scribe(cmd)

# Original Parameters
        cmd = '''STAR --runThreadN 8 \
        --genomeDir {indexes}\
         --alignIntronMax 1000 --alignSJoverhangMin 100 --alignSJDBoverhangMin 15\
         --sjdbGTFfile {gtf}\
         --outFilterIntronMotifs RemoveNoncanonicalUnannotated\
         --outFilterScoreMin 20\
         --outReadsUnmapped Fastx\
         --outSAMstrandField intronMotif\
         --readFilesIn {rin} --outFileNamePrefix {prefx}'''.format(
            rin=fastq_filename,
            prefx=bname + '_', indexes=paths['indexes'],
            gtf=paths['gtf'])''')"""


def sort_and_index_sam(in_dir='sams/'):
    for sam in glob.glob('{a}/*.sam'.format(a=in_dir)):
        print(sam)
        convert_sam_to_bed(sam)
        # bed_to_wig.run('sams/', 'bedgraphs/')
        bname = os.path.basename(sam).partition('.sam')[0]
        scribe('samtools view -bS {i} > {o}'.format(
            i=sam, o=bname + '.bam'
        ))
        scribe('samtools sort {o} -T temp -o {s}.bam'.format(
            o=bname + '.bam',
            s=bname
        ))
        scribe('samtools index {o}'.format(
            o=bname + '.bam'
        ))
        scribe('mv {i} bams/{o}'.format(
            i=bname + '.bam',
            o=bname + '.bam'
        ))
        scribe('mv {i} bams/{o}'.format(
            i=bname + '.bam.bai',
            o=bname + '.bam.bai'
        ))


def split_multimapping_reads(
    in_dir='sams',
    unique_dir='uniquely_mapping/', multi_dir='multimapping/',
    unmapped_dir='unmapped/'):

    for folder in [unique_dir, multi_dir, unmapped_dir]:
        if not os.path.exists(folder):
            os.system('mkdir {}'.format(folder))

    for sam in glob.glob('{a}/*.sam'.format(a=in_dir)):
        # bed_to_wig.run('sams/', 'bedgraphs/')
        bname = os.path.basename(sam).partition('.sam')[0]
        header = ''
        chosen_from_multi = ''

        # Write alignments with MAPQ=255 as uniquely mapping to a separate directory.
        with open('{}/{}'.format(unique_dir, os.path.basename(sam)), 'w') as uniquef:
            unique = get_lines_with_mapq_val(samn=sam, mapq='255')
            uniquef.write(unique)

            # Also write these to the multimappers directory (includes all reads).
            with open('{d}/{a}'.format(d=multi_dir, a=os.path.basename(sam)), 'w') as multif:
                multif.write(unique)
        
        # Write multimapping reads (MAPQ=3) to the multimappers directory.
        with open('{d}/{a}'.format(d=multi_dir, a=os.path.basename(sam)), 'a') as multif:
            multi = get_lines_with_mapq_val(samn=sam, mapq='3', include_header=False)
            multif.write(multi)

        # Write unmapped reads (MAPQ=0) to a directory of unmapped reads.
        with open('{d}/{a}'.format(d=unmapped_dir, a=os.path.basename(sam)), 'w') as unmapf:
            unmapped = get_lines_with_mapq_val(samn=sam, mapq='0')
            unmapf.write(unmapped)


def get_lines_with_mapq_val(samn='', mapq='-1', include_header=True):
    
    if type(mapq) != type(''):
        mapq = str(mapq)
    
    outli = ''
    with open(samn, 'r') as f:
        for li in f:
    
            if li[0] == '@':
                if include_header:
                    outli += li
                continue

            s = li.split('\t')
            
            if len(s) < 5:
                print("Syntax error in sam file: {0}".format(li))
                sys.exit()
            
            if s[4] == mapq:
                outli += li
    
    return outli


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Always writes to a ./sam/ directory.")
    parser.add_argument('-i', '--input_dir')
    args = parser.parse_args()
    run(args)
    
