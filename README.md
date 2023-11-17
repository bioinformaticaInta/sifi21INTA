# sifi21INTA Package

A `Python` code library for the prediction of siRNAs on genomes and transcriptomes.

## Dependences
The following programs and packages need to be installed:

- Programs:
  - `ViennaRNA`
  - `bowtie` (version 1)
- `Pytyon` packages:
  - `biopython`
  - `numpy`
  - `json`

Use example:

For fesign:

```
import sifi21INTA
sifiDesign = sifi21INTA.SifiPipeline(bowtie_db=<complete_path>,query_sequences=<complete_path>, mode=0)
sifiDesign.run_pipeline
```

For off target search:

```
import sifi21INTA
sifiOffTargets = sifi21INTA.SifiPipeline(bowtie_db=<complete_path>,query_sequences=<complete_path>, mode=1)
sifiOffTargets.run_pipeline
```
