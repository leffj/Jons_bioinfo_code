#!/usr/bin/env python

__author__ = "Jonathan Leff"
__version__ = "0.0.1"
__email__ = "leff.jonathan@gmail.com"

"""Calculate GC content for every sequence in a fasta file"""

from Bio import SeqIO
from Bio.SeqUtils import GC
import argparse
import sys
from os.path import split, splitext


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--input_fp', required=True,
		help='The fastx file input filepath')
	parser.add_argument('-o', '--output_fp',
		help='The output fp')
	parser.add_argument('-t','--file_type',\
		help='Filetype detected from file extension, but can be overridden here.'\
		'The input sequences file type (either "fastq" or "fasta")')

	args = parser.parse_args()

	in_fp = args.input_fp

	# determine filetype
	if args.file_type:
		fileType = args.file_type
	else:
		fileType = splitext(in_fp)[1].split('.')[1]

	if fileType == 'fasta' or fileType == 'fa':
		ftype = 'fasta'
	elif fileType == 'fastq' or fileType == 'fq':
		ftype = 'fastq'

	# if the output fp isn't specified, create one
	output_fp = args.output_fp
	if not output_fp:
		input_file_basename, input_file_ext = \
			splitext(split(in_fp)[1])
		output_fp = '%s_gc.txt' % (input_file_basename)

	samples,seqIDs,GCs = calc_gc(in_fp,ftype)
	
	out = open(output_fp,'w')
	for i,j,k in zip(samples,seqIDs,GCs):
		out.write('%s\t%s\t%s\n' %(i,j,k))


def calc_gc(in_fp, ftype):
	# calculate gc content for each sequence in a given fastx file and return sampleID,
	# seqID, and GC content arrays
	in_f = open(in_fp, "U")
	samples = []
	seqIDs = []
	GCs = []
	for seq_rec in SeqIO.parse(in_f,ftype):
		nameSpl = seq_rec.name.split('_')
		samples.append(nameSpl[0])
		seqIDs.append(nameSpl[1])
		GCs.append(GC(seq_rec.seq))
	return samples,seqIDs,GCs




if __name__ == "__main__":
		main()


