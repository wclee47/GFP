import string

def estimate_dist(tx1s, tx2s, txPos):
	poss1, poss2= [], []
	for tx1 in tx1s:
		poss1.append(txPos[tx1])
	for tx2 in tx2s:
		poss2.append(txPos[tx2])
	poss1= string.join(poss1, '|').split('|')
	poss2= string.join(poss2, '|').split('|')
	minDist= "NA"	# Default: INTER
	for pos1 in poss1:
		for pos2 in poss2:
			if pos1.split(':')[0] == pos2.split(':')[0]: # INTRA
				tmp_poss= [pos1.split(':')[1], pos2.split(':')[1]]
				tmp_poss.sort(key= lambda x: int(x.split('-')[0]))
				dist= int(tmp_poss[1].split('-')[0])-int(tmp_poss[0].split('-')[1])
				if minDist == "NA": minDist= dist
				elif dist < minDist: minDist= dist
	return minDist
