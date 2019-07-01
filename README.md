# Pathogenseq scripts

This repository contains scripts which could be useful to the group.
To be as useable, scripts should:
* Be as agnostic to the organism/input data as possible
* Not contain hard-coded variables which the user must change 

## Whole genome assembly
### Finding/Extracting genes from whole/draft genomes

|  Script name               |  Description   |          Data requirements          | Software requirements |
|----------------------------------|---------|-------------------------|--|
| extract_seq_with_gene_anchors.py    | Extracts the sequence beteewn two genes. This is useful when trying to extract a variable gene where blast doesn't pick up the right coordinates. The user can specify two flanking conserved genes to use as anchors. An alignment is optionally generated using mafft. | Reference genome, Reference GFF, query genome | biopython, blast, samtools, mafft                    |
