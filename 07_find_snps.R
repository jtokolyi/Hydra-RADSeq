### Run Stacks populations

if(TRUE){
    setwd("~/mfs/current/oligactis_RADSeq/analysis/RAD1/04_populations")
    system("~/mfs/current/oligactis_RADSeq/soft/stacks-2.5/populations -t 4 -P ~/mfs/current/oligactis_RADSeq/analysis/RAD1/results/denovo_final -O ./ --vcf --fasta-loci --write-single-snp --popmap ~/mfs/current/oligactis_RADSeq/analysis/RAD1/rad1_popmap_final_noreps.txt")

    library(vcfR)
    v <- read.vcfR("populations.snps.vcf")
    quant.depth <- quantile(apply(extract.gt(v,"DP",as.numeric=T),1,mean,na.rm=T),probs=0.975)
    system(paste("vcftools --vcf populations.snps.vcf --max-missing 0.80 --mac 3 --maf 0.05 --minGQ 20 --minDP 3 --out populations_filtered_snps --recode --recode-INFO-all --max-meanDP", quant.depth))
    
    system("vcftools --vcf populations_filtered_snps.recode.vcf --missing-indv")
    miss <- read.csv("out.imiss",sep="\t")
    mean(miss$F_MISS)
    range(miss$F_MISS)

    ### Remove two samples with individual missingness > 0.5 and repeat the whole procedure of SNP selection
    system(paste("vcftools --vcf populations.snps.vcf --out populations_lowmiss --recode --recode-INFO-all", paste(paste("--remove-indv", miss$INDV[miss$F_MISS>0.5], sep=" "),collapse=" ")))
    v <- read.vcfR("populations_lowmiss.recode.vcf")
    quant.depth <- quantile(apply(extract.gt(v,"DP",as.numeric=T),1,mean,na.rm=T),probs=0.975,na.rm=T)
    system(paste("vcftools --vcf populations_lowmiss.recode.vcf --max-missing 0.80 --mac 3 --maf 0.05 --minGQ 20 --minDP 3 --out populations_filtered_snps --recode --recode-INFO-all --max-meanDP", quant.depth))
    
    system("vcftools --vcf populations_filtered_snps.recode.vcf --missing-indv")
    miss <- read.csv("out.imiss",sep="\t")
    mean(miss$F_MISS)
    range(miss$F_MISS)

    v <- read.vcfR("populations_filtered_snps.recode.vcf")
    cat(gsub("_rep.unmapped|.unmapped","",colnames(v@gt))[-1], file="sample_names.txt", sep="\n")
    system("bcftools reheader --samples sample_names.txt -o final.vcf populations_filtered_snps.recode.vcf ")
    file.remove("sample_names.txt")

    r80.loci.names <- paste("CLocus_", v2@fix[,1],sep="")

    system("~/mfs/current/oligactis_RADSeq/soft/bbmap/stats.sh in=populations.loci.fa gc=GC_content.txt")
    gc <- read.csv("GC_content.txt",sep="\t")

    r80.gc <- gc[gc$X.Name%in%r80.loci.names,]
    quantile(r80.gc$GC)
}

