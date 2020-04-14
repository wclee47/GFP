import string, os

def build_txPos(indir, txPos):
	bedfiles= []
	for file in os.listdir(indir):
		if file.endswith(".bed") and file != "transcript.bed" : bedfiles.append(file)
	fi= open(os.path.join(indir, "transcript.bed"), 'r')
	while 1:
		line= fi.readline()
		if not line: break
		fields= line.rstrip().split("\t")
		if not fields[0]+".bed" in bedfiles: continue
		value= "%s:%s-%s" %(fields[0], fields[1], fields[2])
		if fields[3] in txPos: txPos[fields[3]] + '|' + value
		else: txPos[fields[3]]= value

def build_txFASTA(indir, txFASTA):
	name, seq= '', ''
	fi= open(os.path.join(indir, "transcript.fasta"), 'r')
	while 1:
		line= fi.readline()
		if not line:
			txFASTA[name]= seq
			break
		if line.startswith('>'):
			if name != '': txFASTA[name]= seq
			name= line.rstrip()[1:]
			seq= ''
		else: seq= seq + line.rstrip()
