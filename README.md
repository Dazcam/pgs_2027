# Polygenic Scores | Developing Human Brain

> **Polygenic scores for neuropsychiatric disorders and cortical size in prenatal brain donors**  
> Cardiff University · Division of Psychological Medicine and Clinical Neurosciences  
> 2027

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0-brightgreen)](https://snakemake.readthedocs.io)
[![Platform](https://img.shields.io/badge/platform-SLURM%20HPC-blue)](https://slurm.schedmd.com/)

## What this is

This repository computes polygenic scores (PGS) with PRS-CS in 118 donors of European ancestry with prenatal brain single-nucleus RNA sequencing data. Scores are generated for ADHD, anxiety, autism, bipolar disorder, major depressive disorder, schizophrenia and a cortical size factor, plus anxiety and MDD adjusted for ADHD.

## Pipeline overview

| # | Stage | Tools |
| :---: | :--- | :--- |
| 00 | Setup: containers, PRS-CS, LD reference, GWAS downloads | wget, git |
| 01 | Genotypes: European subset, PLINK conversion, genotype PCs | bcftools, PLINK |
| 02 | GWAS preparation: format and QC summary statistics for PRS-CS | Python |
| 03 | PGS: posterior SNP weights and per-donor scores | PRS-CS, PLINK |

## Licence & copyright

© 2027 Cardiff University. See [LICENCE.md](LICENCE.md) for details.
