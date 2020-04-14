#!/usr/local/bin/python

#---------------------------------------
# Python default modules
#---------------------------------------
import sys, getopt, os, re
from pygr import cnestedlist
from string import join
from datetime import datetime


#---------------------------------------
# GFP python modules
#---------------------------------------
from GFP.ExonMapper import *
from GFP.DictBuilder import *
from GFP.HomologyChecker import *
from GFP.DistEstimator import *
from GFP.Filter import *


class Fusion:
	def __init__(self, name, ctext):
		self.name= name
		self.ctext= ctext
		self.nPairs= 0
		self.nSpans= 0
		self.evidences= []
		self.minDist= ''


def main():

	als, als_chrDic, strDict= [], {}, {}

	#----------------------------------------------------------
	# Required parameters
	#----------------------------------------------------------
	infile= ''
	indir= ''
	outprefix= ''
	bl2seqPATH= ''

	#----------------------------------------------------------
	# Optional parameters
	#----------------------------------------------------------
	min_pair= 1			# Minimum # of discordant read-pairs
	min_span= 2			# Minimum # of fusion spanning reads
	min_cov= 10			# Minimum # of base-pairs for both genes.
	min_shift= 1		# Minimum shifting pattern(bp) around fusion point.

	if len(sys.argv) == 1:
		print
		print "GFP --- A tool to detect fusion genes using RNA-Seq"
		print "\nRequired parameters"
		print "\t-i <string>             GSNAP result file."
		print "\t-d <string>             Pre-built exon index directory."
		print "\t-o <string>             Output prefix."
		print "\t--bl2seq <string>       bl2seq excutable path."
		print
		print "Optional parameters"
		print "\t--mpair <integer>       Minimum # of discordant read-pairs, DEFAULT: %d." %min_pair
		print "\t--mspan <integer>       Minumum # of fusion spanning reads, DEFAULT: %d." %min_span
		print "\t--mcov <integer>        Minimum # of base-pairs for both genes, DEFAULT: %d." %min_cov
		print "\t--mshift <integer>      Minimum # of shifting pattern(bp), DEFAULT: %d." %min_shift
		print
		sys.exit(1)
	opts, args= getopt.getopt(sys.argv[1:], "i:d:o:", ["bl2seq=", "mpair=", "mspan=", "mcov=", "mshift="])
	for opt, arg in opts:
		if opt == "-i": infile= arg
		elif opt == "-d": indir= arg
		elif opt == "-o": outprefix= arg
		elif opt == "--bl2seq": bl2seqPATH= arg
		elif opt == "--mpair": min_pair= int(arg)
		elif opt == "--mspan": min_span= int(arg)
		elif opt == "--mcov": min_cov= int(arg)
		elif opt == "--mshift": min_shift= int(arg)

	# Generate strand dictionary
	for line in open(os.path.join(indir, "transcript.bed")).readlines():
		fields= line.rstrip().split("\t")
		strDict[fields[3]]= fields[-1]

	# Preprocessing for pygr
	print str(datetime.now())+"\tPreprocessing pygr requirements..."
	bedfiles= []
	for file in os.listdir(indir):
		if file.endswith(".bed") and file != "transcript.bed": bedfiles.append(file)
	for i in range(len(bedfiles)):
		als.append(cnestedlist.NLMSA(os.path.join(indir, bedfiles[i].split('.')[0]), 'r', pairwiseMode= True))
		als_chrDic[bedfiles[i].split('.')[0]]= i

	# Read GSNAP result
	print str(datetime.now())+"\tReading GSNAP result & extracting fusion evidence..."
	read1Exons, read2Exons= '', ''
	aligns1, aligns2= [], []
	poss1, poss2= '', ''
	strand1, strand2= '', ''
	fi= open(infile, 'r')
	fo= open(outprefix+"_raw.txt", 'w')
	while 1:
		line= fi.readline()
		if not line: break

		if line.startswith('>'):	# Read1
			aligns1= []
			read1Exons= ''
			if line.split()[1] != '1': continue	# Skip multiply-mapped read
			while 1:
				line= fi.readline()
				if line == "\n": break
				aligns1.append(line)
			poss1, read1Exons, strand1= exonMapper(aligns1, fo, indir, als, als_chrDic, strDict)

		if line.startswith('<'):	# Read2
			aligns2= []
			read2Exons= ''
			if line.split()[1] != '1': continue
			while 1:
				line= fi.readline()
				if line == "\n": break
				aligns2.append(line)
			poss2, read2Exons, strand2= exonMapper(aligns2, fo, indir, als, als_chrDic, strDict)

			if read1Exons != '' and read2Exons != '':
				sameFlg= False
				sep= re.compile("[;|]+")
				for read1exon in re.split(sep, read1Exons):
					for read2exon in re.split(sep, read2Exons):
						if read1exon.split('.')[1] == read2exon.split('.')[1]: sameFlg= True

				if not sameFlg:	# Putative fusion pair

					# Donor check module!
					ts1, ts2= [], []	# Transcript strand read1/read2
					for read1exon in re.split(sep, read1Exons):
						for strand in strDict[read1exon.split('.')[0]].split('/'):
							if not strand in ts1: ts1.append(strand)
					for read2exon in re.split(sep, read2Exons):
						for strand in strDict[read2exon.split('.')[0]].split('/'):
							if not strand in ts2: ts2.append(strand)
					donor= check_pDonor(ts1, ts2, strand1, strand2)

					if donor != "NA":
						ctext= "pINTER"
						for pos1 in poss1.split(';'):
							for pos2 in poss2.split(';'):
								if pos1.split(':')[0] == pos2.split(':')[0]: ctext= "pINTRA"
						# No swap!
						fo.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(ctext, poss1, poss2, read1Exons, read2Exons, strand1, strand2, donor))

	fo.close()
	fi.close()
	strDict.clear()

	#---------------------------------------------------
	# Extracting fusion evidences is terminated.
	# Generating gene fusion candidates will start.
	#---------------------------------------------------

	print str(datetime.now())+"\tEstimating homology & distance between genes..."
	txPos, txFASTA= {}, {}
	build_txPos(indir, txPos)
	build_txFASTA(indir, txFASTA)

	# Start reading "outprefix_raw.txt"
	sep= re.compile("[;|]+")	# Set separators
	homologous_pairs= {}
	fi= open(outprefix+"_raw.txt", 'r')
	fusions= []	# Fusion class list
	while 1:
		line= fi.readline()
		if not line: break
		fields= line.rstrip().split("\t")
		hNames, tNames= [], []
		for exon in re.split(sep, fields[3]):
			hNames.append(exon.split('.')[0])
		for exon in re.split(sep, fields[4]):
			tNames.append(exon.split('.')[0])

		# Sequence Homology Detection by bl2seq
		is_homologous= check_homology(hNames, tNames, bl2seqPATH, homologous_pairs, txFASTA, outprefix)

		# Estimate genes' distance
		minDist= estimate_dist(hNames, tNames, txPos)
		if minDist != "NA":
			if minDist < 0: continue	# Overlapping genes

		# Update Fusion list
		if not is_homologous:	# Putative gene fusion
			donor_acceptor= fields[3].split('.')[1]+"\t"+fields[4].split('.')[1]
			if fields[-1] == "read2" or fields[-1] == "tail":
				donor_acceptor= fields[4].split('.')[1]+"\t"+fields[3].split('.')[1]
			type, ctext= fields[0][:1], fields[0][1:]
			does_exist= False
			for fusion in fusions:
				if donor_acceptor == fusion.name:	
					does_exist= True
					if type == 'p': fusion.nPairs+= 1
					else: fusion.nSpans+= 1
					if minDist != "NA":
						if fusion.minDist == "NA": fusion.minDist= str(minDist)
						elif minDist < int(fusion.minDist): fusion.minDist= str(minDist)
					fusion.evidences.append(line.rstrip())
			if not does_exist:
				fusion= Fusion(donor_acceptor, ctext)
				if type == 'p': fusion.nPairs+=1
				else: fusion.nSpans+=1
				fusion.minDist= str(minDist)
				fusion.evidences.append(line.rstrip())
				fusions.append(fusion)

	os.system("rm %s_temp1.fasta %s_temp2.fasta %s_temp.bl2seqout" %(outprefix, outprefix, outprefix))
	fi.close()
	homologous_pairs.clear()
	txFASTA.clear()
	txPos.clear()

	#---------------------------------------------------
	# Generating fusion gene candidates is terminated.
	# Further filtering cascade will be applied.
	#---------------------------------------------------

	print str(datetime.now())+"\tApplying filtering steps & generating output files..."
	fusionNum= 0
	fo_list= open(outprefix+"_fusionList.txt", 'w')
	fo_list.write("ID\tdonor\tacceptor\tcontext\tdist\tnum_pair\tnum_span\n")
	fo_evidence= open(outprefix+"_fusionEvidence.txt", 'w')
	fo_evidence.write("ID\tevidence_type\tdonor_pos\tacceptor_pos\tdonor_exon\tacceptor_exon\n")
	for fusion in fusions:
		if fusion.nPairs < min_pair or fusion.nSpans < min_span: continue
		spans= []
		for evidence in fusion.evidences:
			if evidence.split("\t")[0][:1] == 's': spans.append(evidence)
		spanClusters= cov_filter(spans, min_cov)
		for cluster in spanClusters:
			if len(cluster) < min_span: continue
			shift_pass, fusionNum= shift_filter(cluster, min_shift, fo_evidence, fusionNum)
			if not shift_pass: continue
			for evidence in fusion.evidences:
				evidence_fields= evidence.split("\t")
				if evidence_fields[0][:1] == 'p':
					donor= evidence_fields[-1]
					if donor == "read1":
						fo_evidence.write("GF%d\tread-pair\t%s\t%s\t%s\t%s\n" %(fusionNum, evidence_fields[1], evidence_fields[2], evidence_fields[3], evidence_fields[4]))
					else:
						fo_evidence.write("GF%d\tread-pair\t%s\t%s\t%s\t%s\n" %(fusionNum, evidence_fields[2], evidence_fields[1], evidence_fields[4], evidence_fields[3]))
			fo_list.write("GF%d\t%s\t%s\t%s\t%d\t%d\n" %(fusionNum, fusion.name, fusion.ctext, fusion.minDist, fusion.nPairs, len(cluster)))
	fo_evidence.close()
	fo_list.close()
	print str(datetime.now())+"\tGFP is successfully terminated."


def check_pDonor(ts1, ts2, strand1, strand2):
	donors= []
	dDict= {"+++-":"read1", "++-+":"read2", "---+":"read1", "--+-":"read2",\
		"+-++":"read1", "+---":"read2", "-+--":"read1", "-+++":"read2"}
	for str1 in ts1:
		for str2 in ts2:
			combi= str1+str2+strand1+strand2
			if combi in dDict and not dDict[combi] in donors: donors.append(dDict[combi])
	if len(donors) == 1: return donors[0]
	else: return "NA"


if __name__ == "__main__":
	main()
