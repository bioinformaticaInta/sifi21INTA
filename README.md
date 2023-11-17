# sifi21INTA

Use example:

import sifi21INTA

Design:

sifiDesign = sifi21INTA.SifiPipeline(bowtie_db=<complete_path>,query_sequences=<complete_path>, mode=0)
sifiDesign.run_pipeline

Off target search:

sifiDesign = sifi21INTA.SifiPipeline(bowtie_db=<complete_path>,query_sequences=<complete_path>, mode=1)
sifiDesign.run_pipeline
