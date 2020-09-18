#!/bin/bash
#SBATCH -A lethrus
#SBATCH --job-name=stacks
#SBATCH -n 160
#SBATCH --time=7-00:00:00

OMP_NUM_THREADS=$SLURM_NTASKS $HOME/local/bin/denovo_map.pl -T 72 -o $HOME/oligactis_RADSeq/RAD1/denovo_default --samples $HOME/oligactis_RADSeq/all_samples --popmap $HOME/oligactis_RADSeq/RAD1/rad1_popmap_final.txt --paired --rm-pcr-duplicates -X "populations: --fasta-loci -R 1.0"

OMP_NUM_THREADS=$SLURM_NTASKS $HOME/local/bin/populations -P $HOME/oligactis_RADSeq/RAD1/denovo_default -O $HOME/oligactis_RADSeq/RAD1/denovo_default --popmap $HOME/oligactis_RADSeq/RAD1/rad1_popmap_final.txt --fasta-loci -R 0.0

blastn -task blastn -db $HOME/oligactis_RADSeq/decontamination/nt -query $HOME/oligactis_RADSeq/RAD1/denovo_default/populations.loci.fa -out $HOME/oligactis_RADSeq/RAD1/oli_ntblast.out -num_threads 160 -outfmt 6 -max_target_seqs 5 -max_hsps 1 -evalue 1e-25
