import string

def sortPos(pos):
	chr= pos.split(':')[0]
	scope= [int(pos.split(':')[1].split('-')[0]), int(pos.split(':')[1].split('-')[1])]
	scope.sort()
	return (chr+':'+str(scope[0])+'-'+str(scope[1]), scope[1]-scope[0])

def intersect(pos1, pos2):
	chr1, start1, end1= pos1.split(':')[0], int(pos1.split(':')[1].split('-')[0]), int(pos1.split(':')[1].split    ('-')[1])
	chr2, start2, end2= pos2.split(':')[0], int(pos2.split(':')[1].split('-')[0]), int(pos2.split(':')[1].split    ('-')[1])
	if chr1 != chr2: return False
	if start2 < start1:
		if start1 <= end2: return True
		else: return False
	elif start1 <= start2 <= end1: return True
	elif end1 < start2: return False

def isOverlap(line1, line2):
	fields1= line1.split("\t")
	fields2= line2.split("\t")
	if (intersect(fields1[1], fields2[1]) or intersect(fields1[2], fields2[1])) and\
		(intersect(fields1[1], fields2[2]) or intersect(fields1[2], fields2[2])):
		return True
	else: return False
	
def cov_filter(spans, min_cov):
	spanClusters= []
	for line in spans:
		fields= line.rstrip().split("\t")
		fields[1], len1= sortPos(fields[1])
		fields[2], len2= sortPos(fields[2])
		new_line= string.join(fields, "\t")
		if len1 < min_cov or len2 < min_cov: continue
		cluster_idx= -1
		for i in range(len(spanClusters)):
			for span in spanClusters[i]:
				if isOverlap(new_line, span): cluster_idx= i
		if cluster_idx == -1: spanClusters.append([new_line])
		else: spanClusters[cluster_idx].append(new_line)
	return spanClusters

def shift_filter(cluster, min_shift, fo_evidence, fusionNum):
	fields0= cluster[0].split("\t")
	start, end= int(fields0[1].split(':')[1].split('-')[0]), int(fields0[1].split(':')[1].split('-')[1])
	minDist= end-start+1
	poss= [start, end]
	for i in range(1, len(cluster)):
		fields= cluster[i].split("\t")
		if intersect(fields0[1], fields[1]):
			cStart, cEnd= int(fields[1].split(':')[1].split('-')[0]), int(fields[1].split(':')[1].split('-')[1])
			if cEnd-cStart+1 < minDist: minDist= cEnd-cStart+1
			poss.append(cStart)
			poss.append(cEnd)
		elif intersect(fields0[1], fields[2]):
			cStart, cEnd= int(fields[2].split(':')[1].split('-')[0]), int(fields[2].split(':')[1].split('-')[1])
			if cEnd-cStart+1 < minDist: minDist= cEnd-cStart+1
			poss.append(cStart)
			poss.append(cEnd)
	poss.sort()
	maxDist= poss[-1]-poss[0]+1
	if maxDist-minDist >= min_shift:
		fusionNum+= 1
		for member in cluster:
			evidence_fields= member.split("\t")
			donor= evidence_fields[-1]
			if donor == "head":
				fo_evidence.write("GF%d\tspanning_read\t%s\t%s\t%s\t%s\n" %(fusionNum, evidence_fields[1], evidence_fields[2], evidence_fields[3], evidence_fields[4]))
			else:
				fo_evidence.write("GF%d\tspanning_read\t%s\t%s\t%s\t%s\n" %(fusionNum, evidence_fields[2], evidence_fields[1], evidence_fields[4], evidence_fields[3]))
		return (True, fusionNum)
	else: return (False, fusionNum)
