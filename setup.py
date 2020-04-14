#!/usr/bin/python

from distutils.core import setup

setup(
	name= "GFP",
	version= "v0.7",
	description= "A tool to detect fusion genes from RNA-Seq",
	author= "Won-Chul Lee",
	author_email= "wclee47@gmail.com",
	url= "http://www.gmi.ac.kr",
	py_modules= ["gfp", "build_idxDir", "refGene2bed"],
	packages= ["GFP"],
)
