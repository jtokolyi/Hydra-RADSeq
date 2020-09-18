### commented lines need to be evaluated only on the first occasion the code is used
setwd("~/mfs/current/oligactis_RADSeq/analysis/RAD1/01_decontamination/")

##system("scp h8ftu7@login.miskolc.hpc.niif.hu:~/oligactis_RADSeq/RAD1/denovo_default/populations.loci.fa ./")
##system("~/mfs/current/oligactis_RADSeq/soft/bbmap/stats.sh in=populations.loci.fa gc=oliGC_content.txt")
system("scp h8ftu7@login.miskolc.hpc.niif.hu:~/oligactis_RADSeq/RAD1/ntblast/oli_ntblast.out1 ./")
system("scp h8ftu7@login.miskolc.hpc.niif.hu:~/oligactis_RADSeq/RAD1/ntblast/oli_ntblast.out2 ./")
system("scp h8ftu7@login.miskolc.hpc.niif.hu:~/oligactis_RADSeq/RAD1/ntblast/oli_ntblast.out3 ./")
system("scp h8ftu7@login.miskolc.hpc.niif.hu:~/oligactis_RADSeq/RAD1/ntblast/oli_ntblast.out4 ./")
system("cat oli_ntblast.out1 oli_ntblast.out2 oli_ntblast.out3 oli_ntblast.out4 > oli_ntblast.full")
y <- read.csv("oli_ntblast.full",sep="\t",header=F)
y2 <- y[which(!duplicated(y$V1)),]

#install.packages("devtools")
#devtools::install_github("sherrillmix/taxonomizr")

library(taxonomizr)

#prepareDatabase('accessionTaxa.sql')
blastAccessions<-y2$V2
ids <- accessionToTaxa(blastAccessions,'accessionTaxa.sql')
res <- getTaxonomy(ids,'accessionTaxa.sql')
y2$blastPhylum <- res[,2]
y2$blastOrder <- res[,4]

gc=read.csv("oliGC_content.txt",sep="\t",header=T)
gc <- gc[-1,]

tiff("decontamination.tif",height=1200,width=1200,res=150,compression="lzw")
layout(matrix(c(1,2,3,3),ncol=2,byrow=T))

##histogram of GC content of all RAD loci with default Stacks parameters
hist(gc$GC, main="", xlab="GC content", ylab="No. loci",ylim=c(0,300000))
mtext(text="GC content of all RAD loci",side=3,line=2,font=2)
text("a",y=325000, x=-0.1, xpd=NA,cex=2)

gc$blastOrder <- NA
gc$blastOrder[match(y2[,1],gc[,1])] <- y2$blastOrder
gc$blastOrder[which(is.na(gc$blastOrder))] <- "NoHit"

gc$blastPhylum <- NA
gc$blastPhylum[match(y2[,1],gc[,1])] <- y2$blastPhylum
gc$blastPhylum[which(is.na(gc$blastPhylum))] <- "NoHit"

gc.decont <- gc[gc$blastPhylum%in%c("Cnidaria","NoHit"),]
hist(gc.decont$GC, main="", xlab="GC content", ylab="No. loci",ylim=c(0,300000))
mtext(text="GC content of RAD loci with",side=3,line=2.5,font=2)
mtext(text="identified contaminants removed",side=3,line=1.5,font=2)
text("b",y=325000, x=-0.1, xpd=NA,cex=2)

pie.data <- rev(sort(table(gc$blastOrder)))
pie.data.short <- c(pie.data[1:6],sum(pie.data[7:length(pie.data)]))
names(pie.data.short)[7] <- "Other"

library(RColorBrewer)
par(mar=c(0,0,2,6))
pie(rev(pie.data.short),col=c(brewer.pal("Dark2",n=6),"grey"),
    cex=1.2,font=4,xpd=NA)
mtext(text="Taxonomic distribution of Blast hits (NT-Blast E-value < 1e-05)",side=3,cex=1,font=2,line=0)
text("c",y=0.9, x=-1.3, xpd=NA,cex=2.5)

dev.off()

system("eog decontamination.tif")

cont <- gc[which(!gc$blastPhylum%in%c("Cnidaria","NoHit")),]
write.table(cont[,1],"contamination_read_IDs.txt",col.names=F,row.names=F,quote=F)
system("~/mfs/current/oligactis_RADSeq/soft/bbmap/filterbyname.sh in=populations.loci.fa out=populations_contaminants.fasta names=contamination_read_IDs.txt include=T overwrite=T")
system("scp ./populations_contaminants.fasta h8ftu7@login.miskolc.hpc.niif.hu:/big/home/h8ftu7/oligactis_RADSeq/RAD1/denovo_default/")
