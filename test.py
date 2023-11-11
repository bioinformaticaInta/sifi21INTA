#!/usr/bin/python

#import FileDialog

import sys
import webbrowser
import time
import os

import sifi_pipeline

sifi = sifi_pipeline.SifiPipeline("Arabidopsis_thaliana.TAIR10.cdna.all.fa","/home/sgonzalez/sifi21/testData","testData/test.fasta")
sifi.run_pipeline
