# sifi21INTA Package

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
