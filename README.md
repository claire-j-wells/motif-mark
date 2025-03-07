# motif-mark

## Author: Claire Wells 

This is a script that runs through a FASTA file and identifies all variations of a given set of motifs. The script then visualizes where these motifs exist along the gene through a visual representation of the gene. The image is rendered to scale where 1 pixel is equivalent to 1 base pair. 

## Inputs: 
The inputs for this script are one FASTA file and one motifs text file. The specifications of the FASTA file are that it can hold a **maximum** of 10 records (genes). The specifications for the motifs text file are that this script can accept a **maximum** of 5 motifs. 

## Outputs:
A .png image that shows all of the genes held in the FASTA file and indicates where the motifs exist on the gene with different colors. This image is to scale. The output image contains a legend that indicates what color corresponds to the respective motif as well as what portions of the gene are exon regions. Motif overlap is accounted for and visible. 

## Command:

This script is adapted to accept the following command where `-f` is for the input FASTA file and `-m` is for the input MOTIF txt file: 

`./motif-mark-oop.py -f [FASTA file] -m [MOTIF txt file]`


