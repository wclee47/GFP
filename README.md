# GFP
## Gene Fusion Program

GFP extracts fusion evidence from RNA-seq alignment results generated from GSNAP (Wu and Watanabe, 2005), which intrinsically reports chimeric spliced reads. A set of filtering steps are followed to remove artificially generated fusion evidence.

### Alignment with GSNAP
GFP reads “GSNAP output format” explained here: http://research-pub.gene.com/gmap/src/README.<br/>When you align read sequence files using GSNAP, please do not forget to specify following 3 options.

1. -s, --use-splices=\<STRING\><br/> 
You might need to prepare required files for <STRING>. Please refer to “Detecting known and novel splice sites in GSNAP” section from here: http://research-pub.gene.com/gmap/src/README.
  
2. -N, --novelsplicing=\<INT\><br/> 
Please set <INT> to 1 so that GSNAP looks for novel splicing.
  
3. -E, --distant-splice-penalty=\<INT\><br/>
Please set <INT> to 0, otherwise GSNAP might not report distance splices naturally generated by fusion events.

### Prerequisites
1. “pygr”, a Python module (can be downloaded from http://code.google.com/p/pygr/).
2. “bl2seq” program in the BLAST package.

### Installation
Install with "setup.py". Users can change the installation directory by specifying “setup.py” options (please refer to the “setup.py” manual).
```
python setup.py install
```

### Preparation (building exon index directory)

1. Preparing preliminary files<br/>
Before building exon index directory, three types of files are necessary. If a user plans to use “refGene,txt” (can be downloaded from UCSC database (http://genome.ucsc.edu/) as a gene annotation file, all the necessary files are automatically created by running “refGene2bed.py” included in GFP package. Otherwise, users should make a directory
containing all the necessary files formatted as below.

	* BED files (e.g. chr1.bed, chr2.bed …) for exon information. Each bed file should be first sorted by
exon’s start position and then by exon’s end position. Tab-delimited column names of the files are:<br/>① chromosome<br/>② exon’s start position<br/>③ exon’s end position<br/>④ transcript’s accession . gene name . exon number<br/>⑤ 0 (always)<br/>⑥ transcribed strand
	* “transcript.bed” for transcript information. Tab-delimited column names of the file are:<br/>① chromosome<br/>② transcript’s start position<br/>③ transcript’s end position<br/>④ transcript’s accession<br/>⑤ 0 (always)<br/>⑥ transcribed strand
	* “transcript.fasta” for transcript nucleotide sequence information (FASTA format).

```
...
chr1 850983 851043 NM_152486.SAMD11.exon1 0 +
chr1 851164 851256 NM_152486.SAMD11.exon2 0 +
chr1 855397 855579 NM_152486.SAMD11.exon3 0 +
...
```
```
chr1 67051161 67163158 NM_024763 0 -
chr1 67075872 67163158 NM_207014 0 -
chr1 8335051 8800286 NM_001042681 0 -
chr1 8335051 8406334 NM_001042682 0 -
```
```
>NM_001005484
ATGGTGACTGAATTCATTTTTCTGGGTCTCTCTGATTCTCAGGAACTCCAGACCTTCCTATTTA
TGTTGTTTTTTGTATTCTATGGAGGAATCGTGTTTGGAAACC…
>NM_001005224
ATGGATGGAGAGAATCACTCAGTGGTATCTGAGTTTTTGTTTCTGGGACTCACTCATT…
…
```

2. Building exon index directory<br/>
Once a directory with all the preliminary files is prepared, building the exon index directory is
straightforward. Please run “build_idxDir.py” also included in the package.

### Running GFP

Once building the exon index directory is completed, you are ready to run GFP.

* Optional parameters
	* --mpair : only gene fusions with ≥ (value) will be reported.
	* --mspan: only gene fusions with ≥ (value) will be reported.
	* --mcov: fusion-spanning reads which cover any of the exons from the two genes by < (value) will
be discarded.
	* --mshift: a fusion point (exon-exon boundary) is considered a genuine fusion point when fusionspanning
reads around the fusion point show at least (value) shifting pattern.
