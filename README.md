---
title: "A Generic project for the Temperton lab"
author: "Ben Temperton"
output: pdf_document
bibliography: biblio.bib
csl: isme.csl
---

## Welcome to A Generic project for the Temperton lab

Many of our projects involve the metagenomic analyses of viromes and cellular genomes. Typically, such projects involve the following steps

1. Read QC
2. Removal of human + vertebrate reads
3. Metagenomic assembly using SPAdes
4. Identification of viruses using VirSorter, VirFinder and CAT
5. Clustering of identified viral genomes into viral populations at 95% ANI
6. Gene calling of population representatives using MetaGeneAnnotator
7. Clustering of population representatives into ICTV-recognized genera using VContact2

This project serves as a template for encoding these analyses in `Snakemake` for running on the ISCA infrastructure. To use it, simply create a new project from the `Use this template` button on this project.
