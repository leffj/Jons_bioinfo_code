#!/usr/bin/env python
from __future__ import division

__author__ = "Jonathan Leff"
__version__ = "0.0.1"

from os.path import split, splitext
import argparse
import re
import sys
from Bio import SeqIO
import random

# script uses the reservoir sampling approach


def main():
    parser = argparse.ArgumentParser(description="Get a random subsample of sequences from a "\
        "fastx file")
    parser.add_argument('-i','--input_seqs_fp',required=True,\
        help='The input sequences in fastq or fasta format')
    parser.add_argument('-n','--number_subsample',required=True,\
        help='Specify the number of sequences to subsample')
    parser.add_argument('-o', '--output_fp',
        help='The output fp')
    parser.add_argument('-t','--file_type',\
        help='Filetype detected from file extension, but can be overridden here.'\
        'The input sequences file type (either "fastq" or "fasta")')

    args = parser.parse_args()

    in_fp = args.input_seqs_fp

    nSub = int(args.number_subsample)
    if nSub <= 0:
        raise ValueError,('number_subsample must be >0')

    # determine filetype
    if args.file_type:
        fileType = args.file_type
    else:
        fileType = splitext(in_fp)[1].split('.')[1]

    if fileType == 'fasta' or fileType == 'fa':
        ftype = 'fasta'
    elif fileType == 'fastq' or fileType == 'fq':
        ftype = 'fastq'
    
    fileSize = get_file_size(in_fp)

    # if the output fp isn't specified, create one
    output_fp = args.output_fp
    if not output_fp:
        input_file_basename, input_file_ext = \
         splitext(split(in_fp)[1])
        output_fp = '%s_sub_%s.%s' % (input_file_basename,nSub,ftype)

    in_f = open(in_fp, "U")
    out_seqs = []
    printcounter = 0
    for i, seq_rec in enumerate(SeqIO.parse(in_f,ftype)):
        # fill reservoir array with first nSub sequences in input file
        if i < nSub:
            out_seqs.append(seq_rec)
        # starting with the nSub+1 sequence, use probability to determine if it will replace an
        # existing sequence
        else:
            j = random.randint(0,i)
            if j <= nSub:
                out_seqs[random.randint(0,nSub-1)] = seq_rec
        # for progress status
        if printcounter == 1000:
            pos = in_f.tell()
            display_progress(pos, fileSize)
            printcounter = 0
        printcounter += 1
    sys.stdout.write('\r%5.1f%%' % (100))
    sys.stdout.write('\n')

    out = open(output_fp, "w")
    SeqIO.write(out_seqs,out,ftype)
    out.close()


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
