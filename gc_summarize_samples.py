#!/usr/bin/env python

__author__ = "Jonathan Leff"
__version__ = "0.0.1"
__email__ = "leff.jonathan@gmail.com"

"""Calculate GC content for every sequence in a fasta file"""

import pandas as pd
import argparse
import sys
from os.path import split, splitext


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--input_fp', required=True,
		help='The filepath of the output from gc_from_fastx[_batch].py')
	parser.add_argument('-o', '--output_fp',
		help='The output fp')

	args = parser.parse_args()

	in_fp = args.input_fp

	# if the output fp isn't specified, create one
	output_fp = args.output_fp
	if not output_fp:
		input_file_basename, input_file_ext = \
			splitext(split(in_fp)[1])
		output_fp = '%s_summary.csv' % (input_file_basename)

	out = open(output_fp,'w')

	df = pd.read_table(in_fp,sep='\t',header=None)
	grouped = df.groupby([0])
	means = grouped[2].mean().reset_index()
	stds = grouped[2].std().reset_index()
	smry = pd.merge(means,stds,on=0)
	smry.columns = ['sampleID','mean','std_dev']
	smry.to_csv(output_fp,index=False)



if __name__ == "__main__":
		main()


