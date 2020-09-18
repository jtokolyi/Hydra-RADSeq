setwd("~/mfs/current/oligactis_RADSeq/analysis/RAD1/02_map2Hm_genome")

## Enter libraries here: ####################################################################################
ls.lib1 <- list.files("/home/jtokolyi/mfs/current/oligactis_RADSeq/data/lib1/samples_truncated",full.names=T)
ls.lib2 <- list.files("/home/jtokolyi/mfs/current/oligactis_RADSeq/data/lib2/samples_truncated",full.names=T)
ls.lib3 <- list.files("/home/jtokolyi/mfs/current/oligactis_RADSeq/data/lib3/samples_truncated",full.names=T)
ls.lib4 <- list.files("/home/jtokolyi/mfs/current/oligactis_RADSeq/data/lib4/samples_truncated",full.names=T)
ls <- c(ls.lib1, ls.lib2, ls.lib3, ls.lib4)
ls <- ls[-grep("process_radtags",ls)]
ls <- ls[-grep("rem",ls)]
ls <- ls[-grep("tav|osz",ls)]
ls <- ls[-grep("A.1|A.2|N.1|N.2",ls)]
#############################################################################################################

ls1 <- ls[grep(".1.fq.gz$",ls)]
ls2 <- ls[grep(".2.fq.gz$",ls)]

for(i in 1:length(ls1)){
    system(
        paste("bowtie2 --very-sensitive -p 4 -q -x Hm105",
              " -1 ", ls1[i],
              " -2 ", ls2[i],
              " -S ./output/stdout 2> ./logs/",
    gsub(".[1-2].fq.gz$","",tail(strsplit(ls1[i],"/")[[1]],1)), ".bowtielog", sep="")
    )
}


