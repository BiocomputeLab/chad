# Codon Harmonization And Diversification (CHAD)

**Author:** Shivang Joshi (Biocompute Lab, University of Bristol)

## Overview

The CHAD algorithm can be used to generate DNA sequences of highly similar, or identical proteins, such that they can be encoded multiple times within a single DNA contruct or host cell. Diversification helps to reduce the chance of homologous recombination and improves the robustness of the desired function.

CHAD works by taking a set of amino acid sequences in a multi-FASTA format, codon to amino acid mapping and their frequencies to use during the diversification process (typically related to tRNA concentration in the host cell), as well as, parameters such as the threshold of similarity (i.e., longest stretch of homology across all generated sequences). CHAD uses a random sampling approaches to find candidate solutions for each amino acid sequence in turn. The resultant DNA sequences are then written to a multi-FASTA file using the same IDs as their linked amino acid sequence.

## Usage

To see how CHAD might be used, consider a protein called GFP that you need to encode in 4 different places within a genetic construct that will then be expressed in _E. coli_. The input multi-FASTA file will look like the following (all amino acid sequences are identical):
```
>GFP_1
MSKGEELFTGVVPILVELDGDVNGHKFSVS...
>GFP_2
MSKGEELFTGVVPILVELDGDVNGHKFSVS...
>GFP_3
MSKGEELFTGVVPILVELDGDVNGHKFSVS...
>GFP_4
MSKGEELFTGVVPILVELDGDVNGHKFSVS...
```

CHAD can then be run by using the provided _E. coli_ codon to amino acid mapping/frequency table (`ecoli-codon-data.csv`). In this case, we include rare codons during the diversification process, however, we do provide a mapping/frequency table that has these excluded (`ecoli-codon-data-no-rare.csv`), if rare codons should not be included in the output DNA sequences.
```sh
python chad.py 18 ecoli-codon-data.csv gfp_aa.fasta gfp_dna.fasta
```
This will create DNA sequences where there is a maximum homology of 18 bp between all the sequences.

**Note:** You must have the `chad.py` script in your working directory or added to you `PATH`. Also, CHAD will not add stop codons to the DNA sequences. If required, you will need to add a `*` amino acid to the mapping/frequency table manually.
