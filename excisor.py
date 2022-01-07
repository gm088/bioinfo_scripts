
import numpy as np
import re
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq
import gzip
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--fname', type=str)
parser.add_argument('--direction', type=str)
parser.add_argument('--outname', type=str)
parser.add_argument('--outdir', type=str)
parser.add_argument('--minlen', type=int)
args = parser.parse_args()

tr = 0
untr = 0

if(args.direction=="read1"):
    #sequences = []
    with gzip.open(args.fname, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            rtwoumi = str(Seq(record.name.split("_")[1][-10:]).reverse_complement())
            z = re.search("(.*)({0})".format(rtwoumi[:5]),str(record.seq))
            try:
                if len(z.group(1))>args.minlen:
                    phred = record.letter_annotations['phred_quality']
                    record.letter_annotations = {}
                    record.seq = Seq(z.group(1))
                    record.letter_annotations = {'phred_quality': phred[:len(z.group(1))]}
                    tr += 1
            except AttributeError:
                if(len(record.seq)>args.minlen):
                	untr += 1
                	with open("{0}/{1}".format(args.outdir, args.outname), "a") as f:
                        SeqIO.write(record, f, "fastq")
                else:
                    continue
            with open("{0}/{1}".format(args.outdir, args.outname), "a") as f:
            	SeqIO.write(record, f, "fastq")

    #SeqIO.write(sequences, "{0}/{1}".format(args.outdir, args.outname), "fastq")

if(args.direction=="read2"):
    #sequences = []
    with gzip.open(args.fname, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            roneumi = str(Seq(record.name.split("_")[1][:12]).reverse_complement())
            z = re.search("(.*)({0})".format(roneumi[:5]),str(record.seq))
            try:
                if len(z.group(1))>args.minlen:
                    phred = record.letter_annotations['phred_quality']
                    record.letter_annotations = {}
                    record.seq = Seq(z.group(1))
                    record.letter_annotations = {'phred_quality': phred[:len(z.group(1))]}
                    #sequences.append(record)
            except AttributeError:
                if(len(record.seq)>args.minlen):
                	with open("{0}/{1}".format(args.outdir, args.outname), "a") as f:
                        SeqIO.write(record, f, "fastq")
                else:
                	continue
            with open("{0}/{1}".format(args.outdir, args.outname), "a") as f:
            	SeqIO.write(record, f, "fastq")

    #SeqIO.write(sequences, "{0}/{1}".format(args.outdir, args.outname), "fastq")

with open("{0}/{1}.txt".format(args.outdir, args.outname), "w") as ff:
	ff.write("Of seqs > {0} bp, {1} had some umi, {2} did not\n".format(args.minlen, tr, untr))




