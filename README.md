# Short Description

This is the code of a project I worked on in 2021 summer, a gene alignment program that calls BLAST in Linux to perform 2 rounds of DNA sequence alignment.

## Python Files

There are 2 python files. IndexCreator.py is for creating the alignment indexes, or potential probes. GENOME_FILTER.py is for the 2 rounds of alignment.

## Round 1

To identify high-repetition DNA fragments in telomere and centromere regions of target chromosone.

## Round 2

The second round of alignment with non-target chromosomes to eliminate non-unique sequences. Thresholds for both rounds of alignment have been optimized for better performance (for identifying most effective chromosome-specific DNA sequences) by quantitatively analyzing the impact of threshold changes on final results.


## Post-processing

Afterwards Data Processing have been are commands in Linux (awk commands) and are written inside the Programs_Instructions.pdf file.
