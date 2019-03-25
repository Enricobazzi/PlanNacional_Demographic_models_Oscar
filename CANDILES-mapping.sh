#!/bin/bash

#####################################################################################
## Mapping GENOMA PROJECT FASTQ sequences (Illumina 1.5) to FELIS CATUS reference genome ##
#####################################################################################

# This script will be used to go through the same steps (mapping, adding read groups and merging,
# marking duplicates and realigning) with the remaining sample CANDILES (Lypa23).

# The fastq files were transfered to the genomics-b.ebd.csic.es server from
# from the cesga server in a folder called Enrico_moveLypa with the following path:
# /home/ebazzicalupo/fastqs/Enrico_moveLypa

# tar -zcvf Enrico_moveLypa.tar.gz Enrico_moveLypa
# scp Enrico_moveLypa.tar.gz ebazzicalupo@genomics-b.ebd.csic.es:~
# tar -zxvf Enrico_moveLypa.tar.gz

# As explained in the 0.Mapping_pipeline.Rmd document, this mapping to the cat
# reference genome is necessary in order to generate demographic models through
# machine learning, a step which will be conducted by a collaborator (name ...),
# that will need high coverage sequencing data for at least 2 individuals per species.

# The Read Group Addition step has a starting code which will "extract" the exact run ID
# from the initial fastq list array.
