#!/usr/bin/python3

import pandas as pd
import numpy as np
import seaborn as sns
import pyfastx as pasta
import matplotlib.pyplot as plt
import sys
import os
import tensorflow as tf
from keras.models import load_model
from pkg_resources import resource_filename
from spliceai.utils import one_hot_encode
from pybedtools import BedTool
from maxentpy import maxent
from maxentpy import maxent_fast
from maxentpy.maxent import load_matrix5, load_matrix3


###read and parse data from junctions.bed file
###store processed stuff as attributes
### maybe add some plotting utility functions
###have genome and genome index file in directory, or point to it if its elsewhere

def allstop1(profiles):
	#takes i x j x k array

	berr = []
	for i in profiles:
		err = []
		for j in i:
			thold = 0.05*np.max(j)
			for k in range(50,len(j)):
				if(np.mean(j[k-10:k+10]) < thold):
					err.append(k)
					break
		berr.append(err)
	return berr

def ploot(list,labels=None):

	ctl = [list[i] for i in [0,2,4,6]]
	siz = [list[i] for i in [1,3,5,7]]
	labels = ["DMSO-ctl_for","DMSO-siZC3H4_for","DRB+10-ctl_for","DRB+10-siZC3H4_for","DRB+20-ctl_for","DRB+20-siZC3H4_for","DRB-ctl_for","DRB-siZC3H4_for"]
	clabels = [labels[i] for i in [0,2,4,6]]
	dlabels = [labels[i] for i in [1,3,5,7]]
	fig, axs = plt.subplots(1, 2, sharex=True, gridspec_kw={'hspace': 0}, figsize=[10,5], sharey=True)
	fig.suptitle("Coding Transcripts")
	for i in range(len(ctl)):
	    axs[0].plot(ctl[i], label = clabels[i])
	axs[0].legend()
	for i in range(len(siz)):
	    axs[1].plot(siz[i], label = dlabels[i])
	
	#axs[0].set_ylim([0,110])
	axs[1].legend()
	axs[0].set_title("siCtl")
	axs[1].set_title("siZC3H4")
	axs[0].set_ylabel("Relative Read Density")
	fig.text(0.5, 0.04, '(Kb from TSS) x 10', ha='center')
	fig.patch.set_facecolor('white')
	plt.show()

def allstop(profiles, bps, tholdcoeff=0.05):
	#takes i x j x k array

	indiv_stop = []
	for i in profiles:
		err = []
		for j in i:
			added=False
			#start=10
			thold = tholdcoeff*np.max(j)
			for h in range(5,len(j)):
				if(np.mean(j[h-5:h+5]) > thold):
					start = h
					break
			if(start>10):
				for k in range(start,len(j)):
					if(np.mean(j[k-10:k+10]) < thold):
						err.append(k-start)
						added=True
						break
			else:
				for k in range(10,len(j)):
					if(np.mean(j[k-10:k+10]) < thold):
						err.append(k)
						added=True
						break
			if not added:
				err.append(bps/2)
		indiv_stop.append(err)

	return indiv_stop

def spd(indices):
	result = []
	times = [np.nan,np.nan,10,10,20,20,1,1]
	for i in range(len(times)):
		result.append((indices[i])/(times[i]*10))
	return result

def ezspeed(arr):	#pass an array of averages(collapsed)
	a2 = [i - np.mean(arr[7]) for i in arr]
	a2 = [i - np.min(i[50:]) for i in a2]
	thresh = 0.05*np.max(np.array(a2).flatten(order='C'))

	CI = np.full((8,),np.nan)
	for i in range(len(a2)):
		for j in range(20,len(a2[i])):
			if(np.mean(a2[i][j-10:j+10]) < thresh):
				CI[i] = j
				break

	#print(CI)

	return CI, spd(CI)

def bedtool2df(bt):
	return(pd.read_table(bt.fn, header=None))

def df2bedtool(df):
	return BedTool(df.to_string(index=False,header=False),from_string=True)

def getfastas(df):
    obj = BedTool(df.to_string(index=False,header=False),from_string=True).sequence(fi='/home/gm089/misc/splice/genome.fa')
    return pasta.Fasta(obj.seqfn)

def writefasta(seq, fname):
    with open(fname,'w') as f:
        for i in range(len(seq)):
            f.write(">{0}".format(i))
            f.write("\n")
            f.write(seq[i].seq)
            f.write("\n")
    pass

def randomly_sample_bed(fname, num):
	bed = pd.read_table(fname,sep="\t|,",header=None)
	for i in range(num):
		ind = int(np.random.random_sample()*len(bed))
		try:
			sampledbed = sampledbed.append(bed[ind:ind+1])
		except NameError:
			sampledbed = pd.DataFrame(bed[ind:ind+1])
	sampledbed = sampledbed.sort_values(by=[0,1])
	#procbed = BedTool(sampledbed.to_string(index=False,header=False),from_string=True)
	#del sampledbed

	return sampledbed

#this function is taken from spliceai.utils
def one_hot_encode(seq):

    map = np.asarray([[1, 0, 0, 0],
                      [0, 1, 0, 0],
                      [0, 0, 1, 0],
                      [0, 0, 0, 1]])

    seq = seq.replace('A', '\x00').replace('C', '\x01')
    seq = seq.replace('G', '\x02').replace('T', '\x03')

    return map[np.fromstring(seq, np.int8) % 4]

def infomat(seqlist):

	enc = np.array([one_hot_encode(i) for i in seqlist])
	baseprobs = [i/len(seqlist) for i in np.sum(enc, axis=0)]
	u = [np.sum(-i*np.log2(i)) for i in baseprobs]
	for i in range(len(u)):
		if str(u[i])=='nan':
			u[i] = 0
	ic = 2 - np.array(u)

	return [baseprobs[i]*ic[i] for i in range(len(baseprobs))]

def maxentcalc(donors, acceptors):
	m5 = load_matrix5()
	m3 = load_matrix3()
	donor_scores = [maxent_fast.score5(i, matrix=m5) for i in donors]
	acceptor_scores = [maxent_fast.score3(i, matrix=m3) for i in acceptors]
	donor_probs = list(map(lambda x: 2**x/(2**x+1), donor_scores))
	acceptor_probs = list(map(lambda x: 2**x/(2**x+1), acceptor_scores))

	return donor_scores, acceptor_scores, donor_probs, acceptor_probs

def config_SAI():

	gpus = tf.config.list_physical_devices('GPU')
	if gpus:
		try:
		# Currently, memory growth needs to be the same across GPUs
			for gpu in gpus:
				tf.config.experimental.set_memory_growth(gpu, True)
		except RuntimeError as e:
		# Memory growth must be set before GPUs have been initialized
			print(e)

	context = 10000
	paths = ('models/spliceai{}.h5'.format(x) for x in range(1, 6))
	models = [load_model(resource_filename('spliceai', x)) for x in paths]

	return context, models

def SAI_seq(sequences, context, models):			#input_sequence has to be raw sequence

	dp = []
	ap = []

	for i in sequences:

		input_sequence = i
		x = one_hot_encode('N'*(context//2) + input_sequence + 'N'*(context//2))[None, :]
		#y = np.mean([models[m].predict(x) for m in range(5)], axis=0)
		x = tf.convert_to_tensor(x)
		y = models[1].predict(x)

		acceptor_prob = y[0, :, 1]     #mind these indices
		donor_prob = y[0, :, 2]         #mind these indices
		dp.append(donor_prob)
		ap.append(acceptor_prob)

	return dp, ap

def maxentbox(scores1, scores2):   #scores1 is donor
	fig,axs = plt.subplots(1,2, figsize = (5,5), sharey=True)
	fig.patch.set_facecolor('white')
	fig.suptitle('Splice Sites: MaxEnt')
	axs[0].boxplot(scores1, notch=True, showfliers=False)
	#axs[0].set_ylim([0,15])
	plt.sca(axs[0])
	plt.xticks(range(1,2),['First'])
	axs[1].boxplot(scores2, notch=True, showfliers=False)
	plt.sca(axs[1])
	plt.xticks(range(1,2),['Second'])
	plt.show()

class junctions:

	def __init__(self, fname, fromfile=True):		#fromfile false if you're passing an in-memory dataframe

		if fromfile:
			self.data = pd.read_table(fname, sep="\t|,", header=None)
		else:
			self.data = fname
		
		self.rawjuncs = None
		self.overhang_removed=False
		self.win_don = None
		self.win_acc = None
		#self.DS, self.AS, self.DP, self.AP = None


	def remove_overhang(self):

		if(self.overhang_removed==False):

			self.data[1] = self.data[1]+self.data[12]
			self.data[2] = self.data[2]-self.data[13]
			self.overhang_removed=True
		else:
			print("overhang already removed")


	def window_about_junction(self, d1=3, d2=9, a1=23, a2=3):    
		# d1 to the left of don, a1 to the left of acc (sense dirn)
	    
	    if(self.overhang_removed==False):
	    	print("call remove_overhang first")
	    	return 
	    ##first remove overhang
	    donors = self.data.copy(deep=True)              #for maxent, pass 3,9,23,3
	    #donors[1] = donors[1]+donors[12]                #for spliceport, 80,162,162,80
	    #donors[2] = donors[2]-donors[13]
	    donors_plus = donors[donors[5]=="+"]
	    donors_minus = donors[donors[5]=="-"]
	    donors_plus[1] = donors_plus[1]-d1
	    donors_plus[2] = donors_plus[1]+d2
	    donors_minus[2] = donors_minus[2]+d1
	    donors_minus[1] = donors_minus[2]-d2
	    
	    acceptors = self.data.copy(deep=True)
	    #acceptors[1] = acceptors[1]+acceptors[12]
	    #acceptors[2] = acceptors[2]-acceptors[13]
	    acceptors_plus = acceptors[acceptors[5]=="+"]
	    acceptors_minus = acceptors[acceptors[5]=="-"]
	    acceptors_plus[2] = acceptors_plus[2]+a2
	    acceptors_plus[1] = acceptors_plus[2]-a1
	    acceptors_minus[1] = acceptors_minus[1]-a2
	    acceptors_minus[2] = acceptors_minus[1]+a1

	    self.win_don = [i.seq for i in getfastas(donors_plus)] + [i.antisense for i in getfastas(donors_minus)]
	    self.win_acc = [i.seq for i in getfastas(acceptors_plus)] + [i.antisense for i in getfastas(acceptors_minus)]
	    
	def raw_junctions(self):			# maybe add some functionality here

		juncs = self.data.copy(deep=True)
		jplus = juncs[juncs[5]=="+"]
		jminus = juncs[juncs[5]=="-"]

		self.rawjuncs = [i.seq for i in getfastas(jplus)] + [i.antisense for i in getfastas(jminus)]

		#return [i.seq for i in getfastas(juncplus)] + [i.antisense for i in getfastas(juncminus)]

