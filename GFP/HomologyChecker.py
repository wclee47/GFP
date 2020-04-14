import string, os

def check_homology(tx1s, tx2s, bl2seqPATH, homologous_pairs, txFASTA, outprefix):
	is_homologous= False
	for tx1 in tx1s:
		if is_homologous: break
		for tx2 in tx2s:
			txs= [tx1, tx2]
			txs.sort()
			if txs[0]+'-'+txs[1] in homologous_pairs:
				is_homologous= True
				break
			else:
				open(outprefix+"_temp1.fasta", 'w').write(">%s\n%s" %(tx1, txFASTA[tx1]))
				open(outprefix+"_temp2.fasta", 'w').write(">%s\n%s" %(tx2, txFASTA[tx2]))
				os.system(bl2seqPATH+" -p blastn -e 0.01 -D 1 -i %s_temp1.fasta -j %s_temp2.fasta -o %s_temp.bl2seqout" %(outprefix, outprefix, outprefix))
				if len(open(outprefix+"_temp.bl2seqout").readlines()) > 3:
					is_homologous= True
					homologous_pairs[txs[0]+'-'+txs[1]]= ''
					break
	return is_homologous
