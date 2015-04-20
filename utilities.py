#!/usr/bin/env python

__author__ = "Jonathan Leff"
__version__ = "0.0.1"
__email__ = "leff.jonathan@gmail.com"

"""Utilities"""

from Bio import SeqIO
import random

 
def randomize_seq_order(input_fp,output_fp,datasets=None):
	""" Take a fasta file and randomize the order of the the sequences. This is
	used when doing clustering when there isn't a logical way to order the sequences
	for UCLUST. This is an effort to remove bias for clusters being created around one
	dataset or another.
	"""

	seqs = []
	datasetsTracker = []
	for seq_record in SeqIO.parse(open(input_fp,'rU'), "fasta"):
		seqs.append(seq_record)
	randIds = random.sample(range(0,len(seqs)),len(seqs))
	seqsRand = []
	for i in randIds:
		seqsRand.append(seqs[i])
		datasetsTracker.append(datasets[i])
	out = open(output_fp,'w')
	SeqIO.write(seqsRand, out, "fasta")
	out.close()
	return datasetsTracker