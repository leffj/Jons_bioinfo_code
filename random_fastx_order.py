#!/usr/bin/env python

__author__ = "Jonathan Leff"
__version__ = "0.0.1"
__email__ = "jleff@symbiota.ag"

"""Filter fasta/q file to keep sequences from a list of samples"""

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
		output_fp = '%s_rand.%s' % (input_file_basename,ftype)


	randomize_seq_order(in_fp,output_fp,ftype)

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

 
def randomize_seq_order(input_fp,output_fp,file_type):
	""" Take a fasta file and randomize the order of the the sequences. This is
	used when doing clustering when there isn't a logical way to order the sequences
	for UCLUST. This is an effort to remove bias for clusters being created around one
	dataset or another.
	"""

	seqs = []
	for seq_record in SeqIO.parse(open(input_fp,'rU'), file_type):
		seqs.append(seq_record)
	randIds = random.sample(range(0,len(seqs)),len(seqs))
	seqsRand = []
	for i in randIds:
		seqsRand.append(seqs[i])
	out = open(output_fp,'w')
	SeqIO.write(seqsRand, out, "fasta")
	out.close()


if __name__ == "__main__":
		main()


