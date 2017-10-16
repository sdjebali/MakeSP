# MakeSP
This is a program to compute segmented projections from annotated transcripts of a genome, as illustrated here https://github.com/sdjebali/MakeSP/blob/master/SP.pdf.
The existence of several transcript isoforms for a given gene indeed allow to define segmented projections (see SP.pdf here), which are maximal exonic segments for wich all bases belong to a commun number of transcripts. This means that if we go from 5' to 3' within a gene, a new segmented projection will be created each time there is a new exon boundary (begining or end of an exon).
This program has 1 mandatory argument which is a gene annotation with at least exon rows in gff version 2 format but where the two first (key,value) pairs are transcript_id and gene_id (in this order).

# Installation
An executable called makeSP is provided for a linux system (64 bit architecture), however if this does not work it is also possible to compile a binary from the source.
For this to work you will need to have ocaml installed, and to follow these steps:

$ git clone https://github.com/sdjebali/MakeSP.git

$ cd MakeSP

$ make

This should create an executable called makeSP which will give you the help when launched with no argument.
