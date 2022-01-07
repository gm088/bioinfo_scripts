#!/Library/Frameworks/Python.framework/Versions/3.9/bin/python3

import sys
sys.path.append('/Users/IEO5559/Desktop/misc/splice/script')
import spl
import numpy as np
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('bedfile', type=str)
parser.add_argument('direction', type=str)
parser.add_argument('outname', type=str)
parser.add_argument('--usename', type=str, default="y")
args = parser.parse_args()

if not args:
        print("give bedfile S/AS outname")
        sys.exit()

ends = pd.read_table('{0}'.format(args.bedfile), header=None)

if(args.direction=="S"):
    seqs = [i.seq for i in spl.getfastas(ends[ends[5]=="+"])] + [i.antisense for i in spl.getfastas(ends[ends[5]=="-"])]
    if(args.usename=="y"):
        namelist = ['{0}_(+)'.format(i) for i in ends[ends[5]=="+"].iloc[:,3]] + ['{0}_(-)'.format(i) for i in ends[ends[5]=="-"].iloc[:,3]]
elif(args.direction=="AS"):
    seqs = [i.seq for i in spl.getfastas(ends[ends[5]=="-"])] + [i.antisense for i in spl.getfastas(ends[ends[5]=="+"])]
    if(args.usename=="y"):
        namelist = ['{0}_(-)'.format(i) for i in ends[ends[5]=="-"].iloc[:,3]] + ['{0}_(+)'.format(i) for i in ends[ends[5]=="+"].iloc[:,3]]
else:
	print("bad direction")
if(args.usename!="y"):
    namelist = list(range(len(seqs)))

with open('{0}.fa'.format(args.outname),'w') as f:
    for i in range(len(seqs)):
        f.write(">{0}".format(namelist[i]))
        f.write("\n")
        f.write(seqs[i])
        f.write("\n")

