import string, os
from string import join

def extractPos(line):
	fields= line.rstrip().split()
	strand= fields[2][:1]
	chr= fields[2].split(':')[0][1:]
	pos= fields[2].split(':')[1]
	scope= pos.split("..")
	scope.sort(key= lambda x: int(x))
	start, end= int(scope[0]), int(scope[1])
	return chr, start, end, strand

def compStr(str):
	if str == '+': return '-'
	else: return '+'

def check_sDonor(tsH, tsT, strandH, strandT):
	donors= []
	for str1 in tsH:
		for str2 in tsT:
			if str1+str2 == strandH+strandT and not "head" in donors: donors.append("head")
			elif compStr(str1)+compStr(str2) == strandH+strandT and not "tail" in donors: donors.append("tail")
	if len(donors) == 1: return donors[0]
	else: return "NA"

# It does not support more than 2 segments.
# It does not guarantee reads are included exactly within exons.
def exonMapper(aligns, fo, indir, als, als_chrDic, strDict):
	exonsH= []
	chrH, startH, endH, strandH= extractPos(aligns[0]) # sorted start & end
	posH= "%s:%d-%d" %(chrH, startH, endH)
	if chrH in als_chrDic:
		seqH= als[als_chrDic[chrH]].seqDict[open(os.path.join(indir, "refname.txt")).readline().rstrip()+'.'+chrH]
		for anno, edge in als[als_chrDic[chrH]][seqH[startH-1:endH]].items():
			genomeidx= edge.items()[0][0]
			exonsH.append(anno.name)
	if len(aligns) == 1: return posH, join(exonsH, '|'), strandH

	# Another segment exists
	exonsT= []
	chrT, startT, endT, strandT= extractPos(aligns[1])
	posT= "%s:%d-%d" %(chrT, startT, endT)
	if chrT in als_chrDic:
		seqT= als[als_chrDic[chrT]].seqDict[open(os.path.join(indir, "refname.txt")).readline().rstrip()+'.'+chrT]
		for anno, edge in als[als_chrDic[chrT]][seqT[startT-1:endT]].items():
			genomeidx= edge.items()[0][0]
			exonsT.append(anno.name)

	# Check if same genes exist for indel/local-splicing
	sameFlg= False
	for exonH in exonsH:
		for exonT in exonsT:
			if exonH.split('.')[1] == exonT.split('.')[1]: sameFlg= True

	# No swap!
	if len(exonsH) != 0 and len(exonsT) != 0:
		if not sameFlg:   # Putative fusion

			# Donor check module!
			tsH, tsT= [], []  # transcript strand H/T
			for exonH in exonsH:
				for strand in strDict[exonH.split('.')[0]].split('/'):
					if not strand in tsH: tsH.append(strand)
			for exonT in exonsT:
				for strand in strDict[exonT.split('.')[0]].split('/'):
					if not strand in tsT: tsT.append(strand)
			donor= check_sDonor(tsH, tsT, strandH, strandT)
			if donor == "NA": return '', '', '' # return nothing

			ctext= "sINTRA"
			if chrH != chrT: ctext= "sINTER"
			fo.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"\
				%(ctext, posH, posT, join(exonsH, '|'), join(exonsT, '|'), strandH, strandT, donor))
			return '', '', '' # return nothing

		else: # indel/local-splicing
			if strandH == strandT: return posH+';'+posT, join(exonsH, '|')+';'+join(exonsT, '|'), strandH
			else: return '', '', '' # return nothing

	else: return '', '', '' # return nothing
