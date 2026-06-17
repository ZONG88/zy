# TNF/TNFRSF1A drug-target cis-MR → acute anterior uveitis

Analysis code for: *Genetic target validation of TNF signalling in acute anterior uveitis: a drug-target cis-Mendelian randomization study* (Rheumatology, under review).

## What this does

Two-sample drug-target cis-Mendelian randomization using whole-blood cis-eQTL instruments for TNF (rs1121800, F=599) and TNFRSF1A (rs1800692, F=974). Single cis-SNP Wald ratio estimates for seven outcomes across two channels, following Karhunen et al. (2026) cis-MR framework.

## Data sources

All summary-level data publicly available from IEU OpenGWAS (https://gwas.mrcieu.ac.uk/):

| Data | ID | N |
|------|----|---|
| TNF eQTL (eQTLGen) | eqtl-a-ENSG00000232810 | 14,263 |
| TNFRSF1A eQTL (eQTLGen) | eqtl-a-ENSG00000067182 | 31,684 |
| CRP GWAS | ebi-a-GCST90029070 | 575,531 |
| AAU (iridocyclitis) | finn-b-H7_IRIDOCYCLITIS | 3,622 / 209,287 |
| Ankylosing spondylitis | finn-b-M13_ANKYLOSPON | 1,462 / 164,682 |
| Crohn's disease | ieu-a-12 | 17,897 / 33,977 |
| Rheumatoid arthritis | finn-b-M13_RHEUMA | 6,236 / 147,221 |
| Multiple sclerosis | finn-b-G6_MS | 1,048 / 217,141 |
| Type 2 diabetes | ebi-a-GCST006867 | 655,666 |
| Hypertension | ukb-b-12493 | 54,358 / 408,652 |

## Requirements

- R >= 4.4.1
- TwoSampleMR (>= 0.7.6)
- ieugwasr
- coloc (>= 5.2.3)
- data.table

## Run

```bash
Rscript mr_analysis.R
```

Output: console log with Wald ratio OR, 95% CI, P, and F-statistics for all channel-outcome pairs. Data fetched live from IEU OpenGWAS API.

## Instruments

- **rs1121800** (TNF): chr6:31,535,074, T→A, beta=−0.2915, F=599.4. Within ±10kb of TNF. LD r²=0.025 with HLA-B27 tag (rs4349859).
- **rs1800692** (TNFRSF1A): chr12:6,442,346, C→T, beta=−0.3700, F=974.4. Intronic cis-eQTL inside TNFRSF1A gene body.

## Analysis design

- Primary: TNF cis-eQTL → AAU
- Positive controls: AS, CD, RA (protective), MS (risk-direction)
- Specificity outcomes: T2D, HTN
- Colocalization performed separately (run_coloc.R)

## Scripts

| Script | Contents |
|--------|----------|
| `mr_analysis.R` | Main cis-MR: Wald ratio for 7 outcomes × 2 channels |
| `mr_crp_proxy.R` | CRP proxy channel (rs1800693) — qualitative only |
| `mr_coloc.R` | Bayesian colocalization (TNFRSF1A + TNF/MHC regions) |
| `mr_sensitivity.R` | Steiger directionality + bidirectional MR + HLA-B27 LD |
| `mr_analysis.R` can be run as-is. Other scripts need local eQTL VCF files. |

## Reference

Karhunen V, Woolf B, Bhatnagar P, Gill D, Burgess S. Integrating genetic data with biological insight: a practical guide to cis-Mendelian randomization. Am J Hum Genet. 2026.
