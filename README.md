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

For fesign:

```
import sifi21INTA
sifiDesign = sifi21INTA.SifiPipeline(bowtieDB="<complete_path>",queryFile="<complete_path>", outputDir="<complete_path>", mode=0)
sifiDesign.runPipeline()
```

For off target search:

```
import sifi21INTA
sifiOffTargets = sifi21INTA.SifiPipeline(bowtieDB="<complete_path>",queryFile="<complete_path>", outputDir="<complete_path>", mode=1)
sifiOffTargets.runPipeline()
```
