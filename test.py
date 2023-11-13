#!/usr/bin/python

import sifi_pipeline

sifi = sifi_pipeline.SifiPipeline(bowtie_db="/home/sgonzalez/investigador/sadosky/sifi21/testData/Arabidopsis_thaliana.TAIR10.cdna.all.fa",query_sequences="/home/sgonzalez/investigador/sadosky/sifi21/testData/test.fasta", mode=1)
sifi.run_pipeline
