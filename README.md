# lenga

This repository contains scripts and files aimed at improving and annotating de novo RNA-seq assemblies

- lenga.zip contains an initialization file that the researcher must complete with paths to input files and binaries.
  The structure of the zipped file must be conserved (a common directory containing all files and the scripts/ directory).
  Once the lenga.ini file has been correctly completed, lenga.py can be called in its directory with python2.

- ceg.zip contains a python2 script, a cutoff table and a directory with 248 CEGs (Core Eukaryotic Genes) HMMs (Hidden Markov Models) (Parra, G., Bradnam, K. & Korf, I. CEGMA: a pipeline to accurately annotate core genes in eukaryotic genomes. Bioinformatics 23, 1061–1067 (2007)). 

- SSR_predictions contains 4 tab-delimited tables with SSR predictions for Nothofagus pumilio RNA-seq assemblies.
