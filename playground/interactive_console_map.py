from Bio import SeqIO
from Bio.Seq import Seq
from PySAIS import sais
from bisect import bisect_right
from time import time
from smudgeplot.map import mapper
import pysam
import numpy as np
import gzip
import logging
import csv
import os
import time
import argparse
# argument 1
# args.genomefile
kmer_genome_file = 'tests/data/fake_genome_h1.fa'
output_pattern = 'tests/data/toy_middle_kmers'

parser = argparse.ArgumentParser(description='whatever')
args = parser.parse_args()
args.genomefile = kmer_genome_file
args.o = output_pattern
args.s = 0

kmer_map = mapper(args)
kmer_map.loadGenome()

##########################
###### SUFFIX ARRAY ######
##########################

kmer_map.constructSuffixArray()

kmer_file_name_s1 = output_pattern + '_in_smudge_1.txt'
with open(kmer_file_name_s1, 'r') as s1_kmer_file:
    kmers = [kmer.rstrip() for kmer in s1_kmer_file]


start = time.time()
mapping_list = [kmer_map.searchKmer(kmer) for kmer in kmers]
end = time.time()
print('Done in ' + str(round(end - start, 1)) + ' s')

########### SAVE BAM
header = { 'HD': {'VN': '1.3', 'SO': 'coordinate'},
           'SQ': [{'LN': kmer_map.scf_sizes[i], 'SN': scf} for i,scf in enumerate(kmer_map.scf_names)] }

name2index = dict()
for i,scf in enumerate(kmer_map.scf_names):
    name2index[scf] = i

hlaf_of_kmer = str(int((len(kmers[0]) - 1) / 2))
cigar = hlaf_of_kmer + '=1X' + hlaf_of_kmer + '='
bamfile_name = kmer_map.output_pattern + "_" + str(proc_smudge) + "_mapped.bam"

bamfile = pysam.AlignmentFile(bamfile_name, "wb", header = header)
for i, mapped_kmer in enumerate(mapping_list):
    a = pysam.AlignedSegment()
    a.query_name = 'k' + str(i + 1)
    a.query_sequence = kmers[i]
    # a.mapping_quality = 255
    # a.next_reference_id = 0
    # a.next_reference_start = 199
    # a.template_length = 167
    # a.query_qualities = pysam.qualitystring_to_array("<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<")
    # a.tags = (("NM", 1),
    #           ("RG", "L1"))
    for entry in mapped_kmer:
        a.reference_id = name2index[entry[0]]
        a.reference_start = entry[1]
        # a.cigar = ((0,10), (2,1), (0,25))
        # cigar
        if entry[2] == '+':
            a.flag = 0
        else:
            a.flag = 16
        bamfile.write(a)
    if not mapped_kmer:
        a.flag = 4
        bamfile.write(a)







######
# hit_numbers = [len(mapped_kmer) for mapped_kmer in mapping_list]
# hist = Counter(hit_numbers)
# fractions = [100 * round(hist[i] / len(mapping_list), 4) for i in range(0,10)]

import logging

logging.info("Missing in the assembly: " + str(fractions[0]) + "% of duplicated kmers.")
logging.info("Collapsed in the assembly: " + str(fractions[1]) + "% of duplicated kmers.")
logging.info("Correctly assembled: " + str(fractions[2]) + "% of duplicated kmers.")
logging.info("Assembled three times (once more than it should): " + str(fractions[3]) + "% of duplicated kmers.")
logging.info("Assembled more than three times: " + str(1 - sum(fractions[0:3])) + "% of duplicated kmers.")

import matplotlib.pyplot as plt
# plt.hist(, bins='auto')
bins = range(0, 10)
plt.hist(hit_numbers, bins=bins)
plt.show()

# # # # # # # # # # # # # #
# # scaffold processing # #
# # # # # # # # # # # # # #

# this will be somehow integrate in the class, but for now...

# how I can easily get scf length from scaffold name
scf2size = dict()
for i, seq in enumerate(kmer_map.sequences):
    scf2size[kmer_map.scf_names[i]] = len(seq)

### correctly assembled duplications
from collections import defaultdict
scaffold_dupl = defaultdict(int)
within_scaffold_dupl = defaultdict(int)
for mapped_kmer in mapping_list:
    if len(mapped_kmer) == 2:
        if mapped_kmer[0][0] == mapped_kmer[1][0]:
            within_scaffold_dupl[mapped_kmer[0][0]] += 1
        else:
            scaffold_dupl[mapped_kmer[0][0]] += 1
            scaffold_dupl[mapped_kmer[1][0]] += 1

### Collapsed assembled duplications
collapsed_dupl = defaultdict(int)
for mapped_kmer in mapping_list:
    if len(mapped_kmer) == 1:
        collapsed_dupl[mapped_kmer[0][0]] += 1


with open(kmer_map.output_pattern + "_scaffold_duplicates_list.txt", 'w') as dupl_file:
    for scf in kmer_map.scf_names:
        dupl_file.write(scf + "\t" + str(scaffold_dupl[scf]) + "\t" + str(scaffold_dupl[scf] / scf2size[scf]) + "\n")

with open(kmer_map.output_pattern + "_scaffold_collapsed_duplicates_list.txt", 'w') as dupl_file:
    for scf in kmer_map.scf_names:
        dupl_file.write(scf + "\t" + str(collapsed_dupl[scf]) + "\t" + str(collapsed_dupl[scf] / scf2size[scf]) + "\n")

### distributions

# ????

### Networks

dupl_network = defaultdict(int)
for mapped_kmer in mapping_list:
    if len(mapped_kmer) == 2:
        scf1 = mapped_kmer[0][0]
        scf2 = mapped_kmer[1][0]
        if scf1 < scf2:
            dupl_network[(scf1, scf2)] += 1
        else :
            dupl_network[(scf2, scf1)] += 1

## save it???

# # # # # # # # #
# # full kmer # #
# # # # # # # # #

l = 0
r = len(sa) - 1
while l <= r:
    m = (l + r) // 2
    eval_pos = sa[m]
    genome_kmer = genome[eval_pos:(eval_pos + 21)]
    if kmer == genome_kmer:
        print('jackpot at index: ', m, ' genomic position: ', eval_pos)
        break
    elif genome_kmer < kmer:
        l = m + 1
    else:
        r = m - 1
# This works
# jackpot at index:  421217008  genomic position:  480