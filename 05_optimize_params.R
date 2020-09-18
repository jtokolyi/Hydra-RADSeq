library(vcfR); library(adegenet); library(ape)

plot.only <- TRUE

#tiff("~/mfs/current/oligactis_RADSeq/analysis/RAD1/03_optimize_params.tif", height=2000, width=3000, res=300, compression="lzw")
#layout(matrix(1:12,ncol=4,byrow=T))
#par(mar=c(4,4,2,1),xpd=NA, oma=c(0,0,4,0))

out <- NULL

setwd("~/mfs/current/oligactis_RADSeq/analysis/RAD1/03_optimize_params")
gr <- as.character(read.csv("~/mfs/current/oligactis_RADSeq/analysis/RAD1/03_optimize_params/replicates_popmap.txt",sep="\t",header=F)[,1])
gr <- paste(paste("[", gr[seq(1,length(gr),2)], ",", gr[seq(2,length(gr),2)], "]",sep=""), collapse=",")
    
m.dirs <- list.dirs(recursive=F)[grep("./m", list.dirs(recursive=F), fixed=T)]
for(i in 1:length(m.dirs)){
    setwd(m.dirs[i])
    M.dirs <- list.dirs(recursive=F)[grep("./M", list.dirs(recursive=F), fixed=T)]

    for(j in 1:length(M.dirs)){
        setwd(M.dirs[j])
        temp <- list()
        temp$m <- gsub("./","",m.dirs[i],fixed=T)
        temp$M <- gsub("./","",M.dirs[j],fixed=T)

        if(!plot.only){
            system(
                paste("denovo_map.pl -T 4 --samples ~/mfs/current/oligactis_RADSeq/analysis/RAD1/03_optimize_params/replicates/ --popmap ~/mfs/current/oligactis_RADSeq/analysis/RAD1/03_optimize_params/replicates_popmap.txt --paired -o ./ --rm-pcr-duplicates",
                      "-m", as.numeric(gsub("m","",temp$m)), "-M", as.numeric(gsub("M","",temp$M)), "-n", as.numeric(gsub("M","",temp$M)))
            )
            system("populations -t 4 -P ./ --popmap ~/mfs/current/oligactis_RADSeq/analysis/RAD1/03_optimize_params/replicates_popmap.txt --vcf -R 0.8")
            v <- read.vcfR("populations.snps.vcf")
            quant.depth <- quantile(apply(extract.gt(v,"DP",as.numeric=T),1,mean,na.rm=T),probs=0.975)
            system(paste("vcftools --vcf populations.snps.vcf --mac 3 --maf 0.05 --minGQ 20 --minDP 3 --out populations_filtered_snps --recode --recode-INFO-all --max-meanDP",quant.depth))
            system(paste("~/mfs/current/oligactis_RADSeq/soft/Tiger/tiger vcf=populations_filtered_snps.recode.vcf task=estimateIndReps outname=temp groups=",gr,sep=""))
        }
        system("~/mfs/current/oligactis_RADSeq/soft/stacks-2.5/scripts/stacks-dist-extract populations.log.distribs snps_per_loc_postfilters > temp.txt")
        x <- read.csv("temp.txt",sep="\t",skip=1)
        temp$polymorphic_loci_r80 <- sum(x[,2][-1])
        y <- read.csv("temp_errorRates.txt",sep="\t")
        temp$ErrorRate <- sum(y$Size*y$ErrorRate)/sum(y$Size)
        temp$ErrorRateHom <- sum(y$Size*y$ErrorRateHom)/sum(y$Size)
        temp$ErrorRateHet <- sum(y$Size*y$ErrorRateHet)/sum(y$Size)
        out <- rbind(out, unlist(temp))

        setwd("../")
    }
    setwd("../")
}

out <- as.data.frame(out,stringsAsFactors=F)
out$ErrorRate <- as.numeric(out$ErrorRate)
out$ErrorRateHom <- as.numeric(out$ErrorRateHom)
out$ErrorRateHet <- as.numeric(out$ErrorRateHet)
out$polymorphic_loci_r80 <- as.numeric(out$polymorphic_loci_r80)

tiff("~/mfs/current/oligactis_RADSeq/analysis/RAD1/04_populations/results/optimize_params.tif",width=7.5, height=7.5, res=300, units="in")
layout(matrix(c(1,2,3,4,4,4),ncol=3,byrow=T))
par(cex.lab=1.5,cex.axis=1.5,font.lab=2,mar=c(4,4,6,2))
matplot(matrix(out$polymorphic_loci_r80, ncol=3, byrow=F),type="b",pch=15:17,ylab="No. polymorphic loci")
text("a",x=-1,y=49000,xpd=NA,font=2,cex=2)
matplot(matrix(out$ErrorRateHom, ncol=3, byrow=F),type="b",pch=15:17,ylab="Homozygote allele error rate",ylim=c(0,0.04))
text("b",x=-1.5,y=0.0425,xpd=NA,font=2,cex=2)
mtext(side=1, line=3, text="No. mismatches allowed between stacks", font=2)
mtext(side=1, line=5, text="within and between individuals (parameters M and n)",font=2)
legend(y=0.049+0.003, x=4, horiz=T, bty="n",legend=c("m=3","m=4","m=5"),col=1:3, lty=1:3, pch=15:17, xpd=NA,cex=1.5)
text(y=0.0455+0.003,x=-2.5,labels="Minimum depth of coverage:",cex=1.5,xpd=NA)
matplot(matrix(out$ErrorRateHet, ncol=3, byrow=F),type="b",pch=15:17,ylab="Heterozygote allele error rate",ylim=c(0,0.04))
text("c",x=-1.5,y=0.0425,xpd=NA,font=2,cex=2)

v=read.vcfR("~/mfs/current/oligactis_RADSeq/analysis/RAD1/03_optimize_params/m3/M3/populations_filtered_snps.recode.vcf")
dist <- dist.gene(genind2df(vcfR2genind(v)), method="perc", pairwise.deletion=T)
nj.rep <- nj(dist); nj.rep$tip.label <- gsub(".unmapped","",nj.rep$tip.label)
plot(nj.rep,direction="upwards",cex=1.25)
box()
mtext(side=1,line=2,text="Neighbor-joining tree of replicates at m=3, M=3 and n=3",cex=1,font=2)
text("d",x=-1.8,y=0.39,xpd=NA,font=2,cex=2)
dev.off()

#system("eog ~/mfs/current/oligactis_RADSeq/analysis/RAD1/04_populations/results/optimize_params.tif")
