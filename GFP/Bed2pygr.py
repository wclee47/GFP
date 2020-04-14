#!/usr/bin/env python

from pygr import seqdb, annotation, cnestedlist
import sys, os
import shelve
from collections import defaultdict
from GFP.Bed2pygr1 import *


def parse_bed(filepath):
    for line in open(filepath):
        if line.startswith('browser') or line.startswith('track'):
            continue

        fields = line.split('\t')
        fields[-1] = fields[-1].rstrip()

        yield {
            'sequence': fields[0],
            'start': int(fields[1]),
            'stop': int(fields[2]),
            'name': fields[3],
            'score': int(fields[4]),
            'strand': fields[5],
        }

def load_bed(al, annodb, bedfile, idcounter):
    for e in parse_bed(bedfile):
        idcounter[e['name']] += 1
        nid = '%s|%d' % (e['name'], idcounter[e['name']])

        entry = BEDAnnotation(e['sequence'], e['start'], e['stop'],
                              e['strand'], e['name'], e['score'])
        anno = annodb.new_annotation(nid, entry)
        al.addAnnotation(anno)

def bed2pygr(dbprefix, referencefile, bedfile, indir):

	collision_counter = defaultdict(int)
	chrdb = seqdb.SequenceFileDB(referencefile)
	annodb = annotation.AnnotationDB({}, chrdb)

	al = cnestedlist.NLMSA(dbprefix, 'w', pairwiseMode=True)
	
	load_bed(al, annodb, bedfile, collision_counter)

	al.build(saveSeqDict=True)

	genomeprefix = os.path.basename(referencefile).rsplit('.', 1)[0]
	print >> open(os.path.join(dbprefix) + '.genome', 'w'), genomeprefix
