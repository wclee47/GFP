#!/usr/local/bin/python
#
# Written by: Won-Chul Lee


#---------------------------------------
# Python default modules
#---------------------------------------
import sys, getopt, string, os


#---------------------------------------
# GFP python modules
#---------------------------------------
from GFP.Bed2pygr import *


def main():

	#----------------------------------------------------------
	# Required parameters
	#----------------------------------------------------------
	indir= ''
	referencefile= ''

	if len(sys.argv) == 1:
		print "\nRequired parameters"
		print "\t-i <string>     Bed file directory."
		print "\t-r <string>     Reference FASTA file."
		print
		sys.exit(1)
	opts, args= getopt.getopt(sys.argv[1:], "i:r:")
	for opt, arg in opts:
		if opt == "-i": indir= arg
		elif opt == "-r": referencefile= arg

	fo_refname= open(os.path.join(indir, "refname.txt"), 'w')
	fo_refname.write(referencefile[referencefile.rfind('/')+1:].split('.')[0]+"\n")
	fo_refname.close()

	bedfiles= []
	for file in os.listdir(indir):
		if file.endswith(".bed") and file != "transcript.bed": bedfiles.append(os.path.join(indir, file))

	for bedfile in bedfiles:
		bed2pygr(bedfile[:bedfile.rfind(".bed")], referencefile, bedfile, indir)


if __name__ == "__main__":
	main()
