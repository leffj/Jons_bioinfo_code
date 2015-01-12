#!/usr/bin/env python

__author__ = "Jonathan Leff"
__version__ = "0.0.1"
__email__ = "jleff@symbiota.ag"

"""Reformat fasta/q in tab delimited"""

from Bio import SeqIO
import random
import argparse
import sys
from os.path import split, splitext


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--input_fp', required=True,
		help='The input filepath')
	parser.add_argument('-o', '--output_fp',
		help='The output filepath')
	parser.add_argument('-t','--file_type',\
		help='Filetype detected from file extension, but can be overridden here.'\
		'The input sequences file type (either "fastq" or "fasta")')

	args = parser.parse_args()

	input_fp = args.input_fp
	input_seqs = open(input_fp, "U")


	fileSize = get_file_size(input_fp)

	# determine filetype
	if args.file_type:
		fileType = args.file_type
	else:
		fileType = splitext(input_fp)[1].split('.')[1]

	if fileType == 'fasta' or fileType == 'fa':
		ftype = 'fasta'
	elif fileType == 'fastq' or fileType == 'fq':
		ftype = 'fastq'

	# if the output fp isn't specified, create one
	output_fp = args.output_fp
	if not output_fp:
		input_file_basename, input_file_ext = \
			splitext(split(input_fp)[1])
		output_fp = '%s.txt' % (input_file_basename)

	output = open(output_fp, "w")
	printcounter = 0
	seq_iterator = SeqIO.parse(input_seqs,ftype)
	for seq_record in seq_iterator:
		output.write("%s\t%s\n" %(seq_record.name,seq_record.seq))
		if printcounter == 1000:
			pos = input_seqs.tell()
			display_progress(pos, fileSize)
			printcounter = 0
		printcounter += 1


def get_file_size(FileName):
	File = open(FileName,'U')
	Pos = File.tell()
	File.seek(0, 2)
	FileSize = File.tell()
	File.seek(Pos)
	File.close()
	return FileSize

def display_progress(position,fileSize):
	pct = (100.0*position)/fileSize
	sys.stdout.write('\r%5.1f%%' % (pct))
	sys.stdout.flush()

 
if __name__ == "__main__":
		main()


