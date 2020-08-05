# lenga

This repository contains scripts and files aimed at improving and annotating _de novo_ RNA-seq assemblies

- lenga.py (with its companion folder _scripts_) improves _de novo_ RNA-seq assemblies by eliminating redundant contigs and extending those with more than a certain overlap but considered separate contigs by the assembler.

- annotation.py (with its companion script merge_psl_by_query.php) allows the researcher to annotate an RNA-seq assembly against any sub-set of UniProt peptidic sequences. The script translates each contig into its 6 frames and retains the longest ORF per frame for a 

**Third party software requirements:**

- BLAT: http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/
- Exonerate: http://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate
- Transeq, from the EMBOSS suite: http://emboss.sourceforge.net/download/

**Command line examples:**

python lenga.py name_of_experiment /path/to/FASTA/file /path/to/BLAT/binary /path/to/Exonerate/suite

python annotation.py /path/to/FASTA/file /path/to/UniProt/reference/file /path/to/BLAT/binary /path/to/TranSeq/binary/from/EMBOSS /path/to/merge_psl_by_query.php
