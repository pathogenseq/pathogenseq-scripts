# Pathogenseq scripts

This repository contains scripts which could be useful to the group.
To be as useable, scripts should:
* Be as agnostic to the organism/input data as possible
* Not contain hard-coded variables which the user must change
* Work with standard genomic formats such as FASTA, FASTQ, BAM, etc..
* Not be too fussy with versions of software

## Misc scripts

|  Script name |  Description | Data requirements | Software requirements |
|--------------|--------------|-------------------|-----------------------|
| cgrep    | Given a list of keywords as stdin, finds lines in a file with the keyword in a specific column number. |  | |
| revcom    | Given sequence as stdin, prints out the reverse complement |  | |


## Whole genomes
#### Finding/Extracting genes from whole/draft genomes

|  Script name |  Description | Data requirements | Software requirements |
|--------------|--------------|-------------------|-----------------------|
| extract_seq_with_gene_anchors.py    | Extracts the sequence beteewn two genes. This is useful when trying to extract a variable gene where blast doesn't pick up the right coordinates. The user can specify two flanking conserved genes to use as anchors. An alignment is optionally generated using mafft. | Reference genome, Reference GFF, query genome | biopython, blast, samtools, mafft                    |

#### Genome to genome comparison

|  Script name |  Description | Data requirements | Software requirements |
|--------------|--------------|-------------------|-----------------------|
| assembly2vcf.py | Compare a genome to a reference and produce a VCF file with variants | The two genomes in FASTA format | minimap2, paftools.js |

#### Find number of query sequences occurrences in genome
|  Script name |  Description | Data requirements | Software requirements |
|--------------|--------------|-------------------|-----------------------|
| find_primer_matches.py | Finds the number of times a sequence/motif/primer occurs in a reference genome | Reference genome and query sequences in FASTA or CSV format (name first, primer second column) | fuzznuc from the emboss package |

#### Annotating FASTA files
|  Script name |  Description | Data requirements | Software requirements |
|--------------|--------------|-------------------|-----------------------|
| fasta_add_annotations.py | Add text to the names of the sequences (e.g. date or location) | The fasta file. A csv file with an 'id' column containing the exact same IDs as in the tree and other columns with annotation information. |  |

## VCF manipulation
#### Converting a VCF to different formats
|  Script name |  Description | Data requirements | Software requirements |
|--------------|--------------|-------------------|-----------------------|
| vcf2matrix.py | Creates a variant matrix from a VCF | VCF file | bcftools |

## Tree methods
#### Annotating trees
|  Script name |  Description | Data requirements | Software requirements |
|--------------|--------------|-------------------|-----------------------|
| tree_add_annotations.py | Add text to the tip labels of trees (e.g. date or location) | The newick tree. A csv file with an 'id' column containing the exact same IDs as in the tree and other columns with annotation information. | ete3 |
