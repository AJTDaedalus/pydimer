# pydimer - heterodimer analysis tool
Calculates heterodimers for sets of oligonucleotides.

## Introduction:
  The purpose of this tool is to take a set of oligonucleotides in .fasta format, calculate the potential for heterodimers between each possible pair, and report the results (scores and best alignments) in .csv format.  This tool was originally designed with PCR design in mind, to evaluate potential heterodimers between primers and probes.  The only limitation, however, is that for every pair at least one sequence must be less than 60bp.
  
## Technologies/Dependencies: 
Project was created with:<br>
biopython 1.77<br>
pandas 1.1.1<br>
primer3-py 0.6.1<br>
python 3.8.5

## Use:
To run this script, a .fasta file containing all sequences to be analyzed is required.  
