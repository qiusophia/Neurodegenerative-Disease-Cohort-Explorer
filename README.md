# Neurodegenerative Disease Cohort Explorer

An interactive dashboard for exploring participant metadata from a multi-cohort genetic study of dementia and related conditions.

The idea is simple: pick a diagnosis, a cohort, or an age range, and instantly see who's in the sample. No digging through spreadsheets.

🔗 **Live site:** [[link]](https://sophiaqiu.shinyapps.io/visualizer/)

## A note on the data

The data in this repo is **not** the real data used on the website. The real participant data is private and can't be shared publicly.

The datasets here are made up. They have the same structure as the real ones (same columns and formats), but the values are fake and don't represent any real participants. They're only here so you can run the code and see how everything works.

## What you can filter by

- **Clinical diagnosis:** Alzheimer's disease, frontotemporal dementia, ALS, early-onset dementia, or healthy participant
- **Cohort:** ReDLat, TANGL, Caribe, Brain-LRS-GNA, BOG
- **Genome quality control:** Keep or Exclude
- **Age at evaluation:** a slider that applies to everyone
- **Age at onset:** a slider that only applies to cases

## The data file

`meta.rds` has one row per participant:

| Column | What it is |
|---|---|
| `vcf_ID` | Sample ID that links to genome data |
| `Quality control of genome` | Keep or Exclude |
| `Diagnosis` | Short diagnosis code (AD, FTD, CN, AFM, ALS, EOD, etc.) |
| `Clinical Diagnosis` | Full diagnosis label |
| `Age at onset` | Age when symptoms started (cases only) |
| `Age at Evaluation` | Age when the participant was assessed |
| `Cohort` | Which study cohort they came from |

## Running it locally

[Adjust this to match your setup.]

```r
install.packages(c("shiny", "tidyverse"))
shiny::runApp()
```

## Built with

[R, Shiny, and whatever else you used]

## Repo contents

- `app.R`: the dashboard code
- `meta.rds`: the made-up participant metadata
