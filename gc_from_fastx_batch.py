#!/usr/bin/env python

__author__ = "Jonathan Leff"
__version__ = "0.0.1"
__email__ = "leff.jonathan@gmail.com"

"""Calculate GC content for every sequence in a fasta file"""

from Bio import SeqIO
from Bio.SeqUtils import GC
import argparse
import sys
import os
from os.path import split, splitext
from gc_from_fastx import calc_gc


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--input_dir', required=True,
		help='The directory path contraining the fastx files')
	parser.add_argument('-o', '--output_fp',
		help='The output fp')
	parser.add_argument('-t','--file_type',\
		help='Filetype detected from file extension, but can be overridden here.'\
		'The input sequences file type (either "fastq" or "fasta")')

	args = parser.parse_args()

	in_dir = args.input_dir

	# if the output fp isn't specified, create one
	output_fp = args.output_fp
	if not output_fp:
		dirFP = in_dir.rstrip('/')
		output_fp = '%s_gc.txt' % (dirFP)
		out = open(output_fp,'w')

	for datafile in os.listdir(in_dir):
		# determine filetype
		if args.file_type:
			fileType = args.file_type
		else:
			fileType = splitext(datafile)[1].split('.')[1]

		if fileType == 'fasta' or fileType == 'fa':
			ftype = 'fasta'
		elif fileType == 'fastq' or fileType == 'fq':
			ftype = 'fastq'
		else:
			print "%s not valid filetype - ignoring." % datafile
			continue

		samples,seqIDs,GCs = calc_gc(in_dir + '/' + datafile,ftype)

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


