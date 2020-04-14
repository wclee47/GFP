#!/usr/local/bin/python

class GenomeInterval(object):

    def __init__(self, id, start, stop, strand):
        self.id = id
        self.start = start
        self.stop = stop
        self.strand = strand

    def __repr__(self):
        return '<GenomeInterval id=%s interval=%d-%d(%s)>' % (
                self.id, self.start, self.stop, self.strand)


class BEDAnnotation(GenomeInterval):

    def __init__(self, id, start, stop, strand, name, score):
        GenomeInterval.__init__(self, id, start, stop, strand)
        self.name = name
        self.score = score

    def __repr__(self):
        return '<BEDAnnotation id=%s interval=%d-%d(%s) name=%s score=%d>' % (
                self.id, self.start, self.stop, self.strand,
                self.name, self.score)
