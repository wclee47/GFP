#!/usr/local/bin/python

import sys, getopt, string, os
from datetime import datetime


def main():

	#----------------------------------------------------------
	# Required parameters
	#----------------------------------------------------------
	infile= ''
	reffile= ''
	outdir= ''

	if len(sys.argv) == 1:
		print "\nRequired parameters"
		print "\t-i <string>     RefGene file."
		print "\t-r <string>     Reference FASTA file."
		print "\t-o <string>     Output directory."
		print
		sys.exit(1)
	opts, args= getopt.getopt(sys.argv[1:], "i:r:o:")
	for opt, arg in opts:
		if opt == "-i": infile= arg
		elif opt == "-r": reffile= arg
		elif opt == "-o": outdir= arg

	ref_fi= open(reffile, 'r')

	# Spliting the reference FASTA file by chromosome.
	print str(datetime.now())+"\tSpliting reference FASTA file..."
	chroms= []
	os.system("mkdir .Temp")
	fo_temp= ''
	fa_lineLen= 0
	while 1:
		line= ref_fi.readline()
		if not line:
			fo_temp.close()
			break
		if line.startswith('>'):
			chrom= line.rstrip().split()[0][1:]
			print chrom+" is being processed..."
			chroms.append(chrom)
			if fo_temp != '': fo_temp.close()
			fo_temp= open(".Temp/"+chrom+".fasta", 'w')
			fo_temp.write(line)
		else:
			if fa_lineLen == 0: fa_lineLen= len(line.rstrip())
			fo_temp.write(line)

	print str(datetime.now())+"\tSorting RefGene file..."
	os.system("sort -k3,3 %s > refGene_sorted.txt" %(infile))
	fi= open("refGene_sorted.txt", 'r')

	print str(datetime.now())+"\tMaking BED files..."
	os.system("mkdir %s" %(outdir))
	fo_bed= open(os.path.join(outdir, "transcript.bed"), 'w')
	fo_fasta= open(os.path.join(outdir, "transcript.fasta"), 'w')
	# Read the sorted refGene file
	cur_chr, fo= '', ''	# Current chromosome
	while 1:
		line= fi.readline()
		if not line:
			fo.close()
			os.system("sort -k2,2n -k3,3n %s > %s" %(os.path.join(outdir, cur_chr+"_temp.bed"), os.path.join(outdir, cur_chr+".bed")))
			os.system("rm %s" %(os.path.join(outdir, cur_chr+"_temp.bed")))
			if not cur_chr in chroms: os.system("rm %s" %(os.path.join(outdir, cur_chr+".bed")))
			break
		fields= line.rstrip().split("\t")
		chr= fields[2]
		if chr in chroms:
			fo_bed.write("%s\t%d\t%s\t%s\t0\t%s\n" %(chr, int(fields[4])+1, fields[5], fields[1], fields[3]))

		# Retrieve nucleotide sequences
		if chr in chroms:
			chr_fi= open(".Temp/"+chr+".fasta", 'r')
			exonStarts= fields[9].split(',')
			exonEnds= fields[10].split(',')
			del exonStarts[-1]
			del exonEnds[-1]
			seq= ''
			for i in range(len(exonStarts)):
				start= int(exonStarts[i])+1
				end= int(exonEnds[i])
				seq= seq + string.upper(get_seq(chr_fi, start, end, fa_lineLen))
			chr_fi.close()
			fo_fasta.write(">%s\n%s\n" %(fields[1], seq))

		# Switch chromosome if applicable
		if chr != cur_chr:
			if cur_chr != '':
				fo.close()
				os.system("sort -k2,2n -k3,3n %s > %s" %(os.path.join(outdir, cur_chr+"_temp.bed"), os.path.join(outdir, cur_chr+".bed")))
				os.system("rm %s" %(os.path.join(outdir, cur_chr+"_temp.bed")))
				if not cur_chr in chroms: os.system("rm %s" %(os.path.join(outdir, cur_chr+".bed")))
			fo= open(os.path.join(outdir, chr+"_temp.bed"), 'w')
			cur_chr= chr

		strand= fields[3]
		exonNum= int(fields[8])
		exonStarts= fields[9].split(',')
		exonEnds= fields[10].split(',')
		del exonStarts[-1]
		del exonEnds[-1]

		for i in range(len(exonStarts)):
			if strand == '+': fo.write("%s\t%s\t%s\t%s.%s.exon%d\t0\t%s\n"\
				%(fields[2], exonStarts[i], exonEnds[i], fields[1], fields[12], i+1, strand))
			else: fo.write("%s\t%s\t%s\t%s.%s.exon%d\t0\t%s\n"\
				%(fields[2], exonStarts[i], exonEnds[i], fields[1], fields[12], exonNum-i, strand))

	fo_fasta.close()
	fo_bed.close()
	os.system("rm refGene_sorted.txt")
	os.system("rm -rf .Temp")
	print str(datetime.now())+"\tProgram is successfully done."


def get_seq(fi, start, end, fa_lineLen):
	start_line= start/fa_lineLen
	start_char= start%fa_lineLen
	end_line= end/fa_lineLen
	end_char= end%fa_lineLen

	if start_char == 0:
		start_line= start_line-1
		start_char= fa_lineLen
	if end_char == 0:
		end_line= end_line-1
		end_char= fa_lineLen
	line_num= end_line-start_line+1

	fi.seek(0)
	fi.readline()	# for header?
	fi.seek(start_line*(fa_lineLen+1)+start_char-1, 1)

	seq= ''
	if line_num == 1: seq= seq + fi.read(end_char-start_char+1)
	else:
		for i in range(line_num-1):
			seq= seq + fi.readline().strip()
		seq= seq + fi.read(end_char)
	return seq


if __name__ == "__main__":
	main()
