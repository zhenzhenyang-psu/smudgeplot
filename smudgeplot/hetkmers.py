from collections import defaultdict
from itertools import combinations
import logging

###################
###  DEFINE     ###
###   FUNCTIONS ###
###################

def get_one_away_pairs(kmer_index_family, k):
  """kmer_index_family is a list of (kmer, index) pairs currently under consideration. k is the kmer length. get_one_away_pairs returns a list of pairs of indices where each pair of indices corresponds to a pair of kmers different in exactly one base."""

  #This is the base case for the recursion. Return every pair of indices where the kmers corresponding to those indices differ at exactly one base.
  if k == 1:
    return [(i,j) for ((kmer1,i),(kmer2,j)) in combinations(kmer_index_family, 2) if kmer1 != kmer2]

  #Initialize one_away_pairs, which will be returned by get_one_away_pairs.
  one_away_pairs = []

  #Initialize dictionaries in which the key is a kmer_half (kmer_L or kmer_R) and the value is a list of (other_kmer_half, index) pairs.
  kmer_L_to_index_family = defaultdict(list)
  kmer_R_to_index_family = defaultdict(list)

  #Get the locations for the two halves of the kmer.
  k_L = k // 2
  k_R = k-k_L
  i_L_L = 0
  i_L_R = k_L - 1
  i_R_L = k_L
  i_R_R = k-1

  #For each kmer and index calculate the corresponding left half and right half, then add the necessary (kmer_half, index) pair to the corresponding entries of the dictionary
  for kmer, i in kmer_index_family:
    kmer_L = kmer[i_L_L:i_L_R+1]
    kmer_R = kmer[i_R_L:i_R_R+1]
    kmer_L_to_index_family[kmer_L].append((kmer_R, i))
    kmer_R_to_index_family[kmer_R].append((kmer_L, i))

  #For each left half in which there are multiple kmers with that left half, find the list of pairs in which the right half differs by 1. (aka, if left half matches, recurse on right half).
  for kmer_L_index_family in kmer_L_to_index_family.values(): #same in left half
    if len(kmer_L_index_family) > 1:
      one_away_pairs.extend(get_one_away_pairs(kmer_L_index_family, k_R)) #differ by 1 in right half

  del kmer_L_to_index_family

  #For each right half in which there are multiple kmers with that same right half, find the list of pairs in which the left half differs by 1. (aka, if right half matches, recurse on left half).
  for kmer_R_index_family in kmer_R_to_index_family.values(): #same in right half
    if len(kmer_R_index_family) > 1:
      one_away_pairs.extend(get_one_away_pairs(kmer_R_index_family, k_L)) #differ by 1 in left half

  del kmer_R_to_index_family

  return(one_away_pairs)

#########################################
# WRAPPING FUNCTIONS OF THE TWO MODULES #
#########################################

def middle_one_away(args):
  logging.info('Extracting kmer pairs that differ in the middle nt')

  # file_one_away_pairs = open(args.o + '_one_away_pairs.tsv', 'w')
  file_coverages = open(args.o + '_coverages.tsv', 'w')
  file_kmers = open(args.o + '_sequences.tsv', 'w')

  duplicated = set()
  filtered = set()

  #Initialize a dictionary in which the key is the right kmer_half (not including the middle nucleotide), and the value is a list of (index, coverage) tuples corresponding to kmers that have that particular right kmer_half.
  kmer_R_to_index_family = defaultdict(list)

  # read the first line to get the length of the kmer
  with open(args.infile.name) as dump_file:
    kmer, coverage = dump_file.readline().split()
    k = len(kmer)

  #Get the locations for the two halves of the kmer.
  k_middle = k // 2
  i_L_L = 0
  i_L_R = k_middle - 1
  i_R_L = k_middle + 1
  i_R_R = k-1

  logging.info('Saving ' + args.o + '_coverages.tsv and ' + args.o + '_sequences.tsv files.')
  # Read each line of the input file in order to load the kmers and coverages and process the kmer halves.
  current_kmer_L = ""
  for i1, line in enumerate(args.infile):
    kmer, coverage1 = line.split()
    coverage1 = int(coverage1)

    new_kmer_L = kmer[i_L_L:i_L_R+1]
    kmer_R = kmer[i_R_L:i_R_R+1]
    if new_kmer_L == current_kmer_L:
      if kmer_R in kmer_R_to_index_family:
        if kmer_R in duplicated:
          filtered.discard(kmer_R)
        else:
          duplicated.add(kmer_R)
          filtered.add(kmer_R)
    else:
      for kmer_R in filtered:
        (i1, coverage1), (i2, coverage2) = kmer_R_to_index_family[kmer_R]
        if coverage2 < coverage1:
          # file_one_away_pairs.write(str(i2) + '\t' + str(i1) + '\n')
          file_coverages.write(str(coverage2) + '\t' + str(coverage1) + '\n')
        else:
          # file_one_away_pairs.write(str(i1) + '\t' + str(i2) + '\n')
          file_coverages.write(str(coverage1) + '\t' + str(coverage2) + '\n')
        file_kmers.write(current_kmer_L + 'N' + kmer_R + '\n')
      duplicated = set()
      filtered = set()
      kmer_R_to_index_family = defaultdict(list)
      current_kmer_L = new_kmer_L
    kmer_R_to_index_family[kmer_R].append((i1,coverage1))

  # file_one_away_pairs.close()
  file_coverages.close()
  file_kmers.close()

def all_one_away(args):
  #Initiate kmer and coverages lists.
  kmers = []
  coverages = []

  # Read each line of the input file in order to
  # load the kmers and coverages and process the kmer halves.
  for i, line in enumerate(args.infile):
    kmer, coverage = line.split()
    coverage = int(coverage)
    coverages.append(coverage)
    kmers.append(kmer)

  logging.info('Kmers and coverages loaded.')

  k = len(kmer) # all the kmers in the dump file have the same length, so I can just calc the number of nts in the last one
  # get_one_away_pairs is a recursive function that gatheres indices of all kmer 1 SNP from each other
  one_away_pairs = get_one_away_pairs([(kmer,i) for i,kmer in enumerate(kmers)], k)

  logging.info('Kmer pairs identified.')

  repeated = {}
  for (i1, i2) in one_away_pairs:
    repeated[i1] = i1 in repeated
    repeated[i2] = i2 in repeated

  logging.info('Kmers in unique kmer pairs identified.')

  with open(args.o + '_sequences.tsv', 'w') as file_seqs, open(args.o + '_coverages.tsv', 'w') as file_coverages:
    for (i1, i2) in one_away_pairs:
      if not repeated[i1] and not repeated[i2]:
        cov1 = coverages[i1]
        cov2 = coverages[i2]
        if cov1 < cov2:
          file_coverages.write(str(cov1) + '\t' + str(cov2) + '\n')
          file_seqs.write(kmers[i1] + '\t' + kmers[i2] + '\n')
        else:
          file_coverages.write(str(cov2) + '\t' + str(cov1) + '\n')
          file_seqs.write(kmers[i2] + '\t' + kmers[i1] + '\n')

  logging.info(args.o + '_families.tsv and ' + args.o + '_coverages.tsv files saved.')

