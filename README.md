# sifi21INTA Package

A `Python` code library for siRNAs prediction in genomes and transcriptomes.

### Dependences
The following programs and packages need to be installed:

- Programs:
  - `ViennaRNA`
  - `bowtie` (version 1)
- `Python` packages:
  - `biopython`
  - `json`
  - `plotly`

### Use example:

For design:

```
import sifi21INTA
sifiDesign = sifi21INTA.SifiPipeline(bowtieDB="<DB index name>", queryFile="<Query file name>", outputDir="<directory name>", mode=0)
sifiDesign.runPipeline()
```

For off target search:

```
import sifi21INTA
sifiOffTargets = sifi21INTA.SifiPipeline(bowtieDB="<DB index name>",queryFile="<Query file name>", outputDir="<directory name>", mode=1)
sifiOffTargets.runPipeline()
```
