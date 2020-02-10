#!/bin/bash

## PATH DEFINITION ##
# This will be done on the genomics EBD server

BAM_path="/home/ebazzicalupo/CatRef_bams"
REF="/home/GRUPOS/grupolince/reference_genomes/felis_catus_genome/Felis_catus.Felis_catus_9.0.dna.toplevel.fa"
CONSENSUS_path="/home/ebazzicalupo/PSMC/consensus_fq"
PSMCinput_path="/home/ebazzicalupo/PSMC/inputfiles_psmcfa"
OUTPUT_path="/home/ebazzicalupo/PSMC/output_psmc"
BAM_list=($(ls /home/ebazzicalupo/CatRef_bams/*.bam))

## LOOP ##

for bam in ${BAM_list[@]}
  do
    ID=($( echo $bam | rev | cut -d'/' -f 1 | rev | cut -d '_' -f1,2,3,4 | uniq ))
 	# First step is to generate consensus sequences from BAM files.
	# You generate a consensus with a pipeline of SAMtools mpileup, BCFtools call and a 
	# vcf utility perl script which converts a VCF to FASTQ.
	echo "Generating consensus sequence for ${ID} ... "
	samtools mpileup -Q 30 -q 30 -u -v -f $REF $bam |
	bcftools call -c |
	vcfutils.pl vcf2fq -d 5 -D 34 -Q 30 > $CONSENSUS_path/${ID}.fq
	
	# Then convert the FQ file to the PSMC input format (psmcfa) with a PSCMC utility
	echo "Converting ${ID}.fq to psmc input format ... "
	/opt/psmc/utils/fq2psmcfa $CONSENSUS_path/${ID}.fq > $PSMCinput_path/${ID}.psmcfa
	
	# Then run PSMC
	psmc -p "4+25*2+4+6" -o $OUTPUT_path/${ID}.psmc $PSMCinput_path/${ID}.psmcfa
done
