#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Aim: extract sequence(s) from a gzipped fastq file given their full name
# Copyright (C) 2014-2015 Institut National de la Recherche Agronomique
# License: GPL-3+
# Author: Timothée Flutre [cre,aut]
# Version: 1.0.0 # http://semver.org/
# Download: https://github.com/timflutre/quantgen

import sys, gzip
from Bio.SeqIO.QualityIO import FastqGeneralIterator
if len(sys.argv) != 2:
    sys.stderr.write("example: echo -e \"read1\\nread3\" | extract_fastq.py reads.fq.gz\n")
    sys.exit(1)
wanted = [line.strip() for line in sys.stdin]
for name, seq, qual in FastqGeneralIterator(gzip.open(sys.argv[1])):
    if name in wanted:
        sys.stdout.write(("@%s\n%s\n+\n%s\n" % (name,seq,qual)))
