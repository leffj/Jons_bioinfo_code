#!/usr/bin/env python

__author__ = "Jonathan Leff"
__version__ = "0.0.1"
__email__ = "leff.jonathan@gmail.com"

"""Reformat fasta/q to tab delimited or 2 column tab delimited to fasta"""

from Bio import SeqIO
import random
import argparse
import sys
from os.path import split, splitext


def main():
	parser = argparse.ArgumentParser(description='Reformat fasta/q to tab delimited or 2 column'\
		'tab delimited to fasta. For tab-to-fasta, first column converts to headers.')
	parser.add_argument('-i', '--input_fp', required=True,
		help='The input filepath')
	parser.add_argument('-o', '--output_fp',
		help='The output filepath')
	parser.add_argument('-t','--file_type',\
		help='Filetype detected from file extension, but can be overridden here.'\
		'The input sequences file type ("fastq", "fasta", or "txt")')

	args = parser.parse_args()

	input_fp = args.input_fp
	input_seqs = open(input_fp, "U")


	fileSize = get_file_size(input_fp)

	ftype = get_file_type(input_fp, args.file_type)


	# if the output fp isn't specified, create one
	output_fp = args.output_fp
	if not output_fp:
		input_file_basename, input_file_ext = \
			splitext(split(input_fp)[1])
		if ftype == 'fasta' or ftype == 'fastq':
			output_fp = '%s.txt' % (input_file_basename)
		elif ftype == 'txt':
			output_fp = '%s.fasta' % (input_file_basename)

	output = open(output_fp, "w")

	if ftype == 'fasta' or ftype == 'fastq':
		fastx_to_tab(input_seqs, ftype, output)
	elif ftype == 'txt':
		tab_to_fasta(input_seqs, output)
	else:
		raise RuntimeError("Unknown error.")
	pos = input_seqs.tell()
	display_progress(pos, fileSize)



def fastx_to_tab(input_seqs, file_type, output_file):
	printcounter = 0
	seq_iterator = SeqIO.parse(input_seqs, file_type)
	for seq_record in seq_iterator:
		output_file.write('%s\t%s\n' %(seq_record.name, seq_record.seq))
		if printcounter == 1000:
			pos = input_seqs.tell()
			display_progress(pos, fileSize)
			printcounter = 0
		printcounter += 1


def tab_to_fasta(input_seqs, output_file):
	printcounter = 0
	for seq_record in input_seqs:
		if seq_record.startswith('#'):
			continue
		split_record = seq_record.strip().split('\t')
		header = split_record[0]
		seq = split_record[1]
		output_file.write('>%s\n%s\n' %(header, seq))
		if printcounter == 1000:
			pos = input_seqs.tell()
			display_progress(pos, fileSize)
			printcounter = 0
		printcounter += 1


def get_file_type(fp, specified_file_type):
	# determine filetype
	if specified_file_type:
		fileType = specified_file_type
	else:
		fileType = splitext(fp)[1].split('.')[1]
	if fileType == 'fasta' or fileType == 'fa':
		ftype = 'fasta'
	elif fileType == 'fastq' or fileType == 'fq':
		ftype = 'fastq'
	elif fileType == 'txt' or fileType == 'tab':
		ftype = 'txt'
	else:
		raise RuntimeError("Input file type is not recognized. Should be 'fasta', 'fastq', or 'txt'.")
	return ftype


def get_file_size(FileName):
	File = open(FileName,'U')
	Pos = File.tell()
	File.seek(0, 2)
	FileSize = File.tell()
	File.seek(Pos)
	File.close()
	return FileSize

def display_progress(position, fileSize):
	pct = (100.0*position)/fileSize
	sys.stdout.write('\r%5.1f%%\n' % (pct))
	sys.stdout.flush()

 
if __name__ == "__main__":
		main()


