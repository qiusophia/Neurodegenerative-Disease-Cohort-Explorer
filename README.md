# Variant Visualizer

An interactive Shiny app for exploring genetic variants in **PSEN1**, one of the main genes linked to early-onset Alzheimer's disease. https://sophiaqiu.shinyapps.io/visualizer/

The app draws a lollipop plot along the PSEN1 protein. Each variant shows up at its amino acid position, the exons are color-coded underneath, and above each variant you can see how many carriers there are in each diagnosis group (Alzheimer's, FTD, ALS, early-onset dementia, and healthy controls). You can filter both the people and the variants, and the plot updates live.

## A note on the data

The data in this repo is **not** the real data used on the website. The real participant data is private and can't be shared publicly.

The datasets here are made up. They have the same structure as the real ones (same columns and formats), but the sample IDs, genotypes, ages, and cohorts are all fake and don't represent any real participants. They're only here so you can run the app and see how everything works.

## What you can do

**Filter individuals**
- Diagnosis: Alzheimer's disease, frontotemporal dementia, ALS, early-onset dementia, healthy at evaluation
- Cohort
- Age at onset (only applies to cases)
- Age at evaluation
- Include or exclude people with unknown ages

**Filter variants**
- ClinVar classification: pathogenic, likely pathogenic, uncertain, conflicting, benign, or not reported
- Alzforum curation, based on a built-in lookup of known PSEN1 mutations
- Maximum allele frequency in All of Us, on a log10 slider
- Include or exclude variants never seen in All of Us
- CADD and REVEL score ranges

Carrier counts only include people who pass the individual filters. Variant labels are automatically spaced out so they don't overlap when mutations sit close together.

## How it works

1. **Reads the variants** from `lab.txt`. This is a VCF-style table with one row per variant, annotation columns, and one genotype column per sample.
2. **Parses the protein change** (for example, `p.E280A`) from the RefSeq annotation and keeps the canonical PSEN1 transcript, `NM_000021`.
3. **Finds carriers.** Any genotype with a non-reference allele, like `0/1`, counts as a carrier. Missing genotypes (`./.`) don't count.
4. **Joins carriers to participant info** from `filter_individuals.txt` so they can be filtered by diagnosis, cohort, and age.
5. **Draws the plot** with ggplot2: the exon map on the bottom, variant lollipops in the middle, and carrier counts per diagnosis stacked on top.

## Data files

**`lab.txt`**: variants (tab-separated)

| Column | What it is |
|---|---|
| `#CHROM`, `POS`, `REF`, `ALT` | Variant location and alleles |
| `RefSeq_AA_change` | Protein change per transcript, e.g. `PSEN1:NM_000021:exon8:c.E838A:p.E280A` |
| `ClinVar_clinical_significance` | ClinVar classification |
| `AllOfUs_AF_all` | Allele frequency in All of Us |
| `CADD_phred`, `REVEL_score` | Pathogenicity prediction scores |
| `ID1`, `ID2`, ... | One genotype column per sample (`0/0`, `0/1`, `./.`) |

**`filter_individuals.txt`**: participants (tab-separated)

| Column | What it is |
|---|---|
| Sample ID | Matches the genotype column names in `lab.txt` |
| Phenotype | 1 = control, 2 = case |
| Diagnosis | AD, FTD, ALS, EOD, or CONTROL |
| Age at onset | Cases only (NA for controls) |
| Age at evaluation | Everyone |
| Cohort | Which cohort the person belongs to |

## Running it locally

```r
install.packages(c("shiny", "ggplot2", "dplyr", "stringr", "ggrepel"))
shiny::runApp()
```

Keep `app.R`, `lab.txt`, and `filter_individuals.txt` in the same folder.

## Built with

R, Shiny, ggplot2, dplyr, and stringr

## Repo contents

- `app.R`: the full app (UI and server)
- `lab.txt`: made-up PSEN1 variant and genotype data
- `filter_individuals.txt`: made-up participant data for 100 individuals
