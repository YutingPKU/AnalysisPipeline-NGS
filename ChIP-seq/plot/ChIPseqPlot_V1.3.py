#!/usr/bin/env python
# -*- encoding = utf8 -*-
# Date created : 2016.5.15
# Python version : 2.7.9
# Date Update  : 2017.3.15
#

import sys
import argparse
import pysam
import os
from threading import Thread


__auther__ = "Czh3"
__version__ = "1.2"
__email__ = "zhangchao3@hotmail.com"


def get_TSS_TTS(bed12):
	'''
	geneFeature[
		[strand, chr, start, end, geneName],
		[strand, chr, start, end, geneName],
		...
	]
	'''
	geneFeature = []

	for i in open(bed12):
		if i.startswith("#"):continue
		j = i.split()

		geneFeature.append([j[3], j[2], int(j[4]), int(j[5]), j[1]])

	return(geneFeature)


def getDepth(geneFeature, bam, boundary):
	'''
	len(geneFeatureDepth) = 201
	'''
	global geneFeatureDepth1, geneFeatureDepth2
	geneFeatureDepth1 = []
	geneFeatureDepth2 = []
	depth1 = depth2 = []
	boundary = int(boundary)
	
	if boundary % 100 != 0:
		sys.exit("Wrong boundary length!")

	step = boundary / 100

	DP = {}
	bamfile = pysam.AlignmentFile(bam, "rb")
	for i in geneFeature:
		if (i[3] - i[2] < 100) or i[2] < boundary :
			continue
		for j in bamfile.pileup(region = "%s:%d-%d" % (i[1], i[2]-boundary, i[3]+boundary+1)):
			DP[j.pos] = j.n

		depth1 = []
		for s in range(i[2]-boundary, i[2]+boundary+1, step):
			d = 0
			for j in range(s-10, s+10):
				if DP.has_key(j):
					d += DP[j]
			depth1.append(d/20.0)

		depth2 = []
		for s in range(i[3]-boundary, i[3]+boundary+1, step):
			d = 0
			for j in range(s-10, s+10):
				if DP.has_key(j):
					d += DP[j]
			depth2.append(d/20.0)

		if(i[0] != "+"):
			# reverse if gene in - strand
			tmp = depth1[::-1]
			depth1 = depth2[::-1]
			depth2 = tmp[::-1]

		geneFeatureDepth1.append([i[4]] + depth1)
		geneFeatureDepth2.append([i[4]] + depth2)

	return(geneFeatureDepth1, geneFeatureDepth2)


def splitInteger(m, n):
	assert n > 0
	quotient = m / n
	remainder = m % n
	arr = []
	if remainder > 0:
		arr = [quotient] * (n - remainder) + [quotient + 1] * remainder
	elif remainder < 0:
		arr = [quotient] * (n + remainder) + [quotient - 1] * -remainder
	else:
		arr =  [quotient] * n
	res = [0]
	for i in range(len(arr)):
		res.append(sum(arr[0:i+1]))
	return res

def body_split(start, end, boundary):
	res = []
	step = boundary / 50
	res = range(start-boundary,start,step) + [start+i for i in splitInteger(end-start, 99)] + range(end+step, end+boundary+step, step)
	return res

def getDepth_geneBody(geneFeature, bam, boundary):
	'''
	len(geneBodyDepth) = 200
	'''
	global geneBodyDepth
	geneBodyDepth = []
	bamfile = pysam.AlignmentFile(bam, "rb")
	depth = []
	DP = {}

	for i in geneFeature:
		if (i[3] - i[2] < 100) or i[2] < boundary :
			continue
		depth = []
		for j in bamfile.pileup(region = "%s:%d-%d" % (i[1], i[2]-boundary,i[3]+boundary+1)):
			DP[j.pos] = j.n

		body_step = body_split(i[2], i[3] , boundary)
		for s in body_step:
			d = 0
			for j in range(s-10, s+10):
				if DP.has_key(j):
					d += DP[j]
			depth.append(d/20.0)

		# reverse if gene in - strand
		if i[0] == '-':
			depth = depth[::-1]
		geneBodyDepth.append( [i[4]] + depth)

	return(geneBodyDepth)


def getBamReads(bam):
	'''
	get total reads from a bam file
	'''
	return(reduce(lambda x, y: x + y, [ int(l.split('\t')[2]) for l in pysam.idxstats(bam).split('\n')[0:-1]]))

######## main ###########

parser = argparse.ArgumentParser(description='chip seq plot. <zhangchao3@hotmail.com>',
                                formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-f', '--bed12',help="bed12 file", required=True)
parser.add_argument('-b', '--bam',help="sort index bam file", required=True)
parser.add_argument('-u', '--boundary',help="uptream, downstream boundary [2000]", type=int, default=2000)
parser.add_argument('-o', '--output',help="output dir [ChIP_plot]", default="ChIP_plot")

args = parser.parse_args()


# init 
geneFeatureDepth1 = geneFeatureDepth2 = geneBodyDepth = []

features = get_TSS_TTS(args.bed12)

# mutliprocessing
P1 = Thread(name = getDepth, target = getDepth, args = (features, args.bam, args.boundary,))
P1.start()
P2 = Thread(name = getDepth_geneBody, target = getDepth_geneBody, args = (features, args.bam, args.boundary,))
P2.start()

P1.join()
#P2.join()

#single processing
#(geneFeatureDepth1, geneFeatureDepth2) = getDepth(features, args.bam, args.boundary)
#geneBodyDepth = getDepth_geneBody(features, args.bam, args.boundary)

totolReads = getBamReads(args.bam)
totolReads = float(totolReads / 1000000)

#make output dir
if not os.path.isdir(args.output):
	os.mkdir(args.output)

tssTable = os.path.join(args.output, "tss.datatable")
tss = open(tssTable, 'w')
for i in geneFeatureDepth1:
	tss.write(i[0] + '\t' + '\t'.join([str(j / totolReads) for j in i[1:]]) + "\n")


ttsTable = os.path.join(args.output, "tts.datatable")
tts = open(ttsTable, 'w')
for i in geneFeatureDepth2:
	tts.write(i[0] + '\t' + '\t'.join([str(j / totolReads) for j in i[1:]]) + "\n")

geneBodyTable = os.path.join(args.output, "genebody.datatable")
geneBody = open(geneBodyTable, 'w')

for i in geneBodyDepth:
	geneBody.write(i[0] + '\t' + '\t'.join([str(j / totolReads) for j in i[1:]]) + "\n")





