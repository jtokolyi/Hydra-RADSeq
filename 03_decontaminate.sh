#!/bin/bash
#SBATCH -A lethrus
#SBATCH --job-name=stacks
#SBATCH -n 32
#SBATCH --time=1-00:00:00

module load python/3.6.1

cd $HOME/oligactis_RADSeq/RAD1
bowtie2-build /big/home/h8ftu7/oligactis_RADSeq/RAD1/denovo_default/populations_contaminants.fasta contam

for SampleID in $(sed -r 's/\t.*//' $HOME/oligactis_RADSeq/RAD1/rad1_popmap_final.txt)
do 
bowtie2 --very-sensitive -p 32 -q -x contam \
 -1 $HOME/oligactis_RADSeq/all_samples/$SampleID.1.fq.gz \
 -2 $HOME/oligactis_RADSeq/all_samples/$SampleID.2.fq.gz \
 --al-conc $HOME/oligactis_RADSeq/RAD1/reads_decont/$SampleID.al-conc.fq \
 --un-conc $HOME/oligactis_RADSeq/RAD1/reads_decont/$SampleID.un-conc.fq \
 -S $HOME/oligactis_RADSeq/RAD1/reads_decont/$SampleID.sam \
 2> $HOME/oligactis_RADSeq/RAD1/reads_decont/$SampleID.log
samtools view -bS $HOME/oligactis_RADSeq/RAD1/reads_decont/$SampleID.sam > $HOME/oligactis_RADSeq/RAD1/reads_decont/$SampleID.bam
samtools view -b -f 12 $HOME/oligactis_RADSeq/RAD1/reads_decont/$SampleID.bam > $HOME/oligactis_RADSeq/RAD1/reads_decont/$SampleID.unmapped.bam
samtools bam2fq $HOME/oligactis_RADSeq/RAD1/reads_decont/$SampleID.unmapped.bam > $HOME/oligactis_RADSeq/RAD1/reads_decont/$SampleID.unmapped.fq
cat $HOME/oligactis_RADSeq/RAD1/reads_decont/$SampleID.unmapped.fq | grep '^@.*/1$' -A 3 --no-group-separator > $HOME/oligactis_RADSeq/RAD1/reads_decont/$SampleID.unmapped.1.fq
cat $HOME/oligactis_RADSeq/RAD1/reads_decont/$SampleID.unmapped.fq | grep '^@.*/2$' -A 3 --no-group-separator > $HOME/oligactis_RADSeq/RAD1/reads_decont/$SampleID.unmapped.2.fq

done
