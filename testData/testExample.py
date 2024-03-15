#!/usr/bin/python

import sifi21INTA

pipe = sifi21INTA.SifiPipeline(bowtieDB="Dcitri_genome_index", queryFile="Cs_CHS_mRNA.fas", outputDir=".", mode=0, targetsInRegions=True)
pipe.runPipeline()
