#!/usr/bin/python3
# NGMSAT Next Generation Microsatellite Analysis Tool
# Purpose: Output statistics about number of repeats found across an Illumina Read Dataset
# Usage: ngmsat

import pprint
import argparse
import sys
import os
import gzip
import re
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import matplotlib.pyplot as plt
import numpy as np


def parse_arguments():
    """ Parse command line arguments """
    # Parse arguments
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("-i", dest="input_table", type=str, required=True, help="input manifest file, a csv with sample_name,for_primer,rev_primer,optional_repeat")
    arg_parser.add_argument("-r", dest="reads_dir", type=str, required=True, help="fastq directory contianing forward amd reverse reads")
    arg_parser.add_argument("-o", dest="outdir", type=str, required=False, default="ngmsat/", help="output directory")
    arg_parser.add_argument("-l", dest="outlog", type=str, required=False, default="ngmsat_log.txt", help="output log")
    arg_parser.add_argument("-f", dest="force", action="store_true", help="force rewrite of output dircetory")
    return arg_parser.parse_args()


def CreateHistoGram(data, save_path, header, x_label, y_label):
    """
    :param Data: list of values, must be a complete set [11,1,1,1,1,1,1,2,2,2,2,2,2,28,8,8,8,8,8,8,8,8,8,10,101,101101,10]
    :param save_path: exact name for file to be created
    :return: nothing
    """
    # An "interface" to matplotlib.axes.Axes.hist() method
    n, bins, patches = plt.hist(x=data, bins=np.arange(len(set(data))) - 0.5, color='#0504aa', alpha=0.7, rwidth=0.85)
    plt.grid(axis='y', alpha=0.75)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(header)
    #plt.text(23, 45, r'$\mu=15, b=3$')
    plt.xticks(range(len(set(data))))
    maxfreq = n.max()
    # Set a clean upper y-axis limit.
    #plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)
    plt.savefig(save_path, bbox_inches='tight')
    plt.clf()


def parse_manifest(manifest_file):
    """
    :param manifest_file is a csv file with sample names (identical to read samples), forward primer, reverse primer and repeat
    :return: table_lookup = {} # sample: [forward_primer, reverse_primer, repeat_sequence]
    """
    table_lookup = {}
    for line in open(manifest_file):
        sample_name, forward_primer, reverse_primer, repeat_seq = line.rstrip().split(',')
        table_lookup[sample_name] = [forward_primer, reverse_primer, repeat_seq]
    return table_lookup


def parse_read_dir(read_directory):
    """
    :param read_directory: a directory full of gzipped reads, paired end.
        reads should have this format SampleName_AACCGTTC-TGAGCTAG_L001_R1_001.fastq.gz
    :return: fastq lookup, based on sample names
    """
    fastq_lookup = {}  # sample: [forward, reverse]
    for fastq in os.listdir(read_directory):
        if '_R1_' in fastq and fastq.endswith('.fastq.gz'):
            sample = fastq.split('/')[-1].split('_')[0] # this is based on read format
            fastq_lookup.setdefault(sample, [read_directory + fastq, read_directory + fastq.replace('_R1_', '_R2_')])
    return fastq_lookup


def IdentifyPrimerRepeat(fastqgz, primer, repeat):
    """
    :param fastqgz: path to fastq file
    :param primer: primer string
    :param repeat: repeat string
    :return: read_dict contains all info about repeat and primer in read,
    :return: repeat dict contains counts for each number of reapeats. i.e. the repeat occurred seven times in X number of reads
    """
    # Initiat returned files
    read_dict = {}  # header: [repeat_count, primer_found, repeat_found]
    repeat_dict = {}  # repeat_count: total

    # initiate variables for various statistics
    repeat_length = len(repeat)
    num_unfound = 0  # number of sequences that contain 0 reads
    num_unfound_primers = 0 # sequences without the primer
    data_set = []  # list containing the number of repeats for each fastq

    # get sample name
    basename = os.path.splitext(os.path.basename(fastqgz))[0].split('_')[0]

    # determine read direction
    direction = 'reverse'
    if '_R1_' in fastqgz:
        direction = 'forward'

    # #### Setup saving specific fastq files
    # heter_info = sample_hetero_dict[basename]
    # output_most_abundant = open(basename + '-' + direction + '-' + str(heter_info[0]) + '.fasta', 'w')
    # output_second_abundant = open(basename + '-' + direction + '-' + str(heter_info[1]) + '.fasta', 'w')
    # output_one_abundant = open(basename + '-' + direction + '-' + str(1) + '.fasta', 'w')
    # most_count = 0
    # second_count = 0
    # one_count = 0
    # ######

    with gzip.open(fastqgz, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            header = str(record.id)
            sequence = str(record.seq)

            # Determine if primer is in sequence, must be exact match
            primer_match = re.findall(r"((" + primer + ")+)", sequence)
            if not primer_match:
                num_unfound_primers += 1

            # recover the longest repeat sequence
            match = ''
            matches = re.findall(r"((" + repeat + ")+)", sequence)
            for tup in matches: # loop through matches and get the largest hit
                for t in tup:
                    if len(t) > len(match):
                        match = t
            num_repeats = int(len(match)/repeat_length)
            if num_repeats == 0:
                num_unfound += 1

            # write read info into dictionaries
            data_set.append(num_repeats)
            read_dict[header] = [num_repeats, primer_match, match]
            if primer_match: # only output data if the primer was found
                if str(num_repeats) in repeat_dict.keys():
                    repeat_dict[str(num_repeats)] += 1
                else:
                    repeat_dict[str(num_repeats)] = 1

            # # Write specific reads to files
            # if num_repeats in heter_info or num_repeats == 1:
            #     if num_repeats == heter_info[0] and most_count < 5:
            #         output_most_abundant.writelines('>' + header + '\n' + sequence + '\n')
            #         most_count += 1
            #     if num_repeats == heter_info[1] and second_count < 5:
            #         output_second_abundant.writelines('>' + header + '\n' + sequence + '\n')
            #         second_count += 1
            #     if num_repeats == 1 and one_count < 5:
            #         output_one_abundant.writelines('>' + header + '\n' + sequence + '\n')
            #         one_count += 1

    # Create per sample histogram
    CreateHistoGram(data_set, basename+'-'+direction+'.png', basename + '-'+direction+'  Repeat Size Distribution', 'Number of Repeats', 'Frequency')
    # print statistics
    print ('#Repeat=', repeat, 'Reads=', len(read_dict.keys()), '#noprimer=', num_unfound_primers, "#norepeat=", num_unfound)
    return repeat_dict, read_dict


# parse arguments
args = parse_arguments()

# create output directory and log
if args.outdir[-1] != '/':  # make sure its consistent
    args.outdir = args.outdir + '/'
if os.path.exists(args.outdir):
    if not args.force:
        print ("Output directory '{}' already exists, please use a new name or use -f".format(args.outdir))
        sys.exit()
    else:
        os.rmdir(args.outdir)
        os.mkdir(args.outdir)
else:
    os.mkdir(args.outdir)

# fill input lookup tables
table_lookup = parse_manifest(args.input_table)
fastq_lookup = parse_read_dir(args.reads_dir)

# write input data info to log files
pp = pprint.PrettyPrinter(indent=4)
pp.pprint(table_lookup)
pp.pprint(fastq_lookup)

print("Sample_Name,Repeat,Repeat_count,#Reads_Forward, #Reads_Reverse")
for sample in table_lookup.keys():
    print('#working on sample,', sample)
    forward_reads, reverse_reads = fastq_lookup[sample]
    f_primer, r_primer, repeat = table_lookup[sample]

    print("\n",sample, f_primer, repeat,"\n")
    if repeat != 'complex': # skip complex repeats
        all_num_repeats = []
        f_dict, f_read_dict = IdentifyPrimerRepeat(forward_reads, f_primer, repeat) # parse forward read
        for k in f_dict.keys():
            all_num_repeats.append(int(k))
        my_dna = Seq(repeat, generic_dna)
        revcomp_repeat = str(my_dna.reverse_complement())
        r_dict, r_read_dict = IdentifyPrimerRepeat(reverse_reads, r_primer, revcomp_repeat) # parse reverse reads
        for k in f_dict.keys():
            if int(k) not in all_num_repeats:
                all_num_repeats.append(int(k))
        all_num_repeats = sorted(all_num_repeats)
        for key in all_num_repeats:
            key = str(key)
            try:
                forward_count = f_dict[key]
            except:
                forward_count = 0
            try:
                reverse_count = r_dict[key]
            except:
                reverse_count = 0
            print (sample, repeat, key, forward_count, reverse_count, sep=',')





# # Test Data Location
# ## INPUT DATA
# fastq_dir = "/home/genome/vanessa/Project_STRB/raw_reads" # directory containing forward and reverse reads, two files per sample
# # THO samples
# #fastq_dir = '/home/genome/vanessa/Project_VanTHOtest/raw_reads'
# table = "/home/genome/vanessa/Project_STRB/raw_reads/repeat_table.tsv" # table that includes repeat sequence per sample

# data for heterozygotes
sample_hetero_dict = {"CSF1POB":[12,11], "D13S317B":[12,13], "D16S39B":[9,13], 'D5S818B':[11,12], "D7S820B":[9,10], "TPOB":[8,11]}
