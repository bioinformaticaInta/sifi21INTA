#!/usr/bin/python

#import FileDialog

import sys
import webbrowser
import time
import os

import sifi_pipeline

sifi = sifi_pipeline.SifiPipeline(bowtie_db="/home/sgonzalez/sifi21/testData/Arabidopsis_thaliana.TAIR10.cdna.all.fa",query_sequences="/home/sgonzalez/sifi21/testData/test.fasta", mode=0)
sifi.run_pipeline
