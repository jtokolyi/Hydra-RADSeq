library(doBy)

all.samples <- c("M25_1_1","M25_2_1","M25_3_1","M25_4_1","M26_1_5","M26_2_15" ,
                 "M26_2_3","M26_2_7","M26_2_9","M26_3_6","M26_4_2","M26_4_3",  
                 "M26_5_1","M28_1_1","M28_10_2","M28_10_3","M28_11_1","M28_11_4", 
                 "M28_11_5","M28_12_10","M28_12_12","M28_12_4","M28_12_7","M28_12_8", 
                 "M28_12_9","M28_2_1","M28_2_10","M28_2_15","M28_2_2","M28_2_20", 
                 "M28_2_3","M28_2_9","M28_3_1","M28_3_2","M28_3_3","M28_3_8",  
                 "M28_3_9","M28_4_1","M28_4_2","M28_5_12","M28_5_15","M28_5_4",  
                 "M28_5_5","M28_5_8","M28_6_10","M28_6_13","M28_6_15","M28_6_2",  
                 "M28_6_5","M28_7_1","M28_7_2","M28_8_3","M28_8_4","M28_8_8",  
                 "M28_9_1","M28_9_11","M28_9_12","M28_9_4","M28_9_8","M38_1_3",  
                 "M38_2_1","M38_3_1","M38_3_2","M38_4_1","M38_5_1","M44_1_1",  
                 "M44_1_13","M44_1_2","M44_1_3","M44_1_6","M44_2_2","M44_2_5",  
                 "M44_3_2","M44_4_3","M44_5_1","M44_6_3","M44_7_3","M47_1_3",  
                 "M47_2_4","M47_3_1","M52_1_1","M52_2_4","M52_3_1","M52_4_4",  
                 "M25_1_1_rep","M26_2_9_rep","M38_3_1_rep","M47_3_1_rep","M52_3_1_rep",
                 "M44_1_13_rep","M28_1_1_rep","M28_6_2_rep","M28_10_2_rep","M28_11_5_rep",
                 "M108_1_5","M108_1_7","M108_1_11","M108_1_12","M108_1_12_rep","M108_2_3","M108_2_5","M67_1_1",         
                 "M67_1_2","M67_2_2_rep","M67_2_2","M67_3_1","M67_3_6","M72_1_3","M72_1_3_rep","M72_1_4",         
                 "M72_2_3","M72_2_6","M72_3_2","M72_4_2","M72_4_5","M72_5_1","M85_1_2","M85_2_2",         
                 "M85_3_2","M85_4_2","M85_4_2_rep","M85_5_2","M90_1_4","M90_1_6","M90_2_1","M90_2_1_rep",     
                 "M90_2_5","M90_3_1","M90_3_5","M90_4_1","M90_4_2","R12_1_1","R12_1_2","R12_1_2_rep",     
                 "R12_2_2","R12_3_4")

### Read in read stats from process_radtags.pl output ##########################################################
l1 <- list.files("~/mfs/current/oligactis_RADSeq/data/lib1/samples_truncated/",full.names=T)
l1 <- l1[grep("process_radtags",l1)]
l2 <- list.files("~/mfs/current/oligactis_RADSeq/data/lib2/samples_truncated/",full.names=T)
l2 <- l2[grep("process_radtags",l2)]
l3 <- list.files("~/mfs/current/oligactis_RADSeq/data/lib3/samples_truncated/",full.names=T)
l3 <- l3[grep("process_radtags",l3)]
l4 <- list.files("~/mfs/current/oligactis_RADSeq/data/lib4/samples_truncated/",full.names=T)
l4 <- l4[grep("process_radtags",l4)]

l <- c(l1, l2, l3, l4)

res <- NULL

for(i in 1:length(l)) res <- append(res,readLines(l[i])[13:17])
res <- res[!duplicated(res)]
cat(res, file="reads.csv", sep="\n")
r <- read.csv("reads.csv",sep="\t")
r <- r[as.character(r$Filename)%in%all.samples,]

sum(r$Total)/2 ## total no. reads
sum(r$Retained)/sum(r$Total) ## prop. reads retained
mean(r$Retained/2) ## mean no. reads / sample
range(r$Retained/2)

### Read in mapping to Hm105 genome ##################################################################
r$Hm105.mappingrate <- NA
for(i in 1:nrow(r)){
    temp <- readLines(paste("~/mfs/current/oligactis_RADSeq/analysis/RAD1/02_map2Hm_genome/logs/",
                            as.character(r$Filename[i]),".bowtielog",sep=""))
    r$Hm105.mappingrate[i] <- as.numeric(strsplit(temp[15], split="%")[[1]][1])
}
tiff("Hm105_mappingrate.tif",compression="lzw",height=750,width=1000,res=150)
par(cex.lab=1.2,cex.axis=1.2,oma=c(0,1,0,0))
hist(r$Hm105.mappingrate,breaks=30, main="",xlab="Mapping percentage")
mtext("Mapping rate to the Hydra 2.0 Genome Assembly",line=2,font=2,cex=1.2)
mtext("(with Bowtie2 --very-sensitive)",line=1)
dev.off()

### Write complete popmap for comparison excluding vulgaris ######################################
write.table(cbind(r$Filename[r$Hm105.mappingrate<40],sapply(strsplit(r$Filename[r$Hm105.mappingrate<40],split="_"),"[",1)),
            file="~/mfs/current/oligactis_RADSeq/analysis/RAD1/rad1_popmap_final.txt",
            col.names=F, row.names=F, quote=F, sep="\t")
system("scp ~/mfs/current/oligactis_RADSeq/analysis/RAD1/rad1_popmap_final.txt h8ftu7@login.miskolc.hpc.niif.hu:~/oligactis_RADSeq/RAD1/")

### Write complete popmap for comparison excluding replicates ######################################
r2 <- r[r$Hm105.mappingrate<40,]
r2$Sample <- gsub("_rep","",r2$Filename)
r3 <- orderBy(~Sample-Retained, r2)
r3 <- r3[!duplicated(r3$Sample),]
r3$Filename <- paste(r3$Filename,".unmapped",sep="")
write.table(cbind(r3$Filename,sapply(strsplit(r3$Filename,split="_"),"[",1)),
            file="~/mfs/current/oligactis_RADSeq/analysis/RAD1/rad1_popmap_final_noreps.txt",
            col.names=F, row.names=F, quote=F, sep="\t")
system("scp ~/mfs/current/oligactis_RADSeq/analysis/RAD1/rad1_popmap_final_noreps.txt h8ftu7@login.miskolc.hpc.niif.hu:~/oligactis_RADSeq/RAD1/")

### Read in coverage stats from denovo_map.pl output #################################################
cov1 <- readLines("~/mfs/current/oligactis_RADSeq/analysis/RAD1/results/denovo_w.rep/denovo_map.log")
cov1 <- cov1[(grep("Depths of Coverage", cov1)+1):(grep("Depths of Coverage", cov1)+100)]
## Mean/range coverage before removal of PCR duplicates
mean(as.numeric(gsub("x$","",sapply(strsplit(cov1, split=": "), "[", 2))))
range(as.numeric(gsub("x$","",sapply(strsplit(cov1, split=": "), "[", 2))))

cov2 <- readLines("~/mfs/current/oligactis_RADSeq/analysis/RAD1/results/denovo_w.rep/gstacks.log")
## PCR duplication rate
cov2[grep("putative PCR duplicates",cov2)]
## Coverage after romoval of PCR duplicates
cov2[grep("effective per-sample coverage",cov2)]

