#!/usr/bin/env python

__author__ = "Jonathan Leff"
__version__ = "0.0.1"
__email__ = "leff.jonathan@gmail.com"

""" Take multiple datasets and cluster them using UCLUST. The output will be a fasta file with
	centroids for each cluster and a text file listing each input sequence, the file origin, and
	which cluster (ie OTU) it was clustered into.
"""

import argparse
import os
import utilities
from subprocess import call



def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--input_fps', nargs='+', required=True,
		help='The input fasta filepaths. Can add multiple with spaces')
	parser.add_argument('-o', '--output_dir', type=str,
		help='The output directory')
	parser.add_argument('-s', '--similarity', type=float, default=0.97,
		help='The sequence similarity within clusters (Default = 0.97)')

	args = parser.parse_args()

	# create output directory
	if not args.output_dir:
		outdir = 'clustered_data'
	else:
		outdir = args.output_dir
	if check_output_dir(outdir) == 'abort':
		return


	# concatenate the datasets
	catOutfp = outdir + '/all_seqs.fa'
	DStracker = [] # keep track of which sequences come from which dataset
	with open(catOutfp, 'w') as outfile:
		for idx,fname in enumerate(args.input_fps):
			with open(fname, 'U') as infile:
				for line in infile:
					outfile.write(line)
					if line[0] == '>':
						DStracker.append(idx)

	# randomize the sequences
	randOutfp = outdir + '/all_seqs_random.fa'
	DStrackerRando = utilities.randomize_seq_order(catOutfp,randOutfp,DStracker)


	# run uclust. Assumes that usearch executable is called usearch7.
	print "\n\nClustering Sequences...\n\n"
	cmd = 'usearch7 -cluster_smallmem %s -id %f -centroids %s/centroids.fa -uc %s/clust_data.txt \
			-usersort -maxaccepts 0 -maxrejects 0' %(randOutfp,args.similarity,outdir,outdir)
	call(cmd,shell=True)
	print "\n"


	# parse the output
	print "Writing output...\n"
	clustsOut = open(outdir+'/seqs_and_cluster_ids.txt','w')
	clustsOut.write('input_seq\tdataset\tcluster\n')
	for idx,line in enumerate(open(outdir+'/clust_data.txt')):
		fields = line.split('\t')
		if fields[0] == 'C':
			continue
		clustsOut.write('%s\t%s\t%s\n' %(fields[8],DStrackerRando[idx],fields[1])) #seqID,dataset,clusterID, 
	clustsOut.close()


	print "Done."


def check_output_dir(output_dir):
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)
	else:
		print("Overwrite existing output? (Y/n)")
		response = raw_input()
		if response != 'Y':
			print("Specify new output directory.")
			return 'abort'



if __name__ == "__main__":
	main()