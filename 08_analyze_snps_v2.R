library(vcfR); library(adegenet); library(poppr); library(ape)
library(readxl); library(ComplexHeatmap); 
library(circlize); library(RColorBrewer); library(RSplitsTree); library(PopGenome)

setwd("~/mfs/current/oligactis_RADSeq/analysis/RAD1/04_populations")

########################### Read in data ####################################################
ids1 <- as.data.frame(read_xls("~/mfs/current/oligactis_RADSeq/data/lib1/Hydra_DNS_Qubit_final.xls"))
ids2 <- as.data.frame(read_xlsx("~/mfs/current/oligactis_RADSeq/data/lib4/Hydra_rad_lib_4_indexed.xlsx"))
ids2 <- ids2[which(!is.na(ids2$"Reproductive Mode")),]
ids2 <- ids2[-grep("_rep",ids2$ID),]

ids1 <- ids1[c("Sex","ID","North","East")]
ids2 <- ids2[c("Reproductive Mode","ID","North","East")]; names(ids2) <- c("Sex","ID","North","East")
ids <- rbind(ids1, ids2)
ids$Sex <- tolower(ids$Sex)
ids$Sex[ids$Sex%in%c("immature","immat","immat+asex")] <- "immat"
ids$Sex <- factor(ids$Sex, levels=c("immat","male","female","nr","asex"))

ids$Type <- ifelse(as.character(ids$Sex)%in%c("female","immature","male","immat","immat+asex"),"sex","asex")
ids$ID2 <- gsub("/", "_", as.character(ids$ID))
ids$Pop_code <- substring(ids$ID2, 1, 3)
v <- read.vcfR("final.vcf")
ids <- ids[which(ids$ID2%in%colnames(v@gt)),]

########################################################################################
######## Clone detection and relatedness ###############################################
### Poppr
vgc <- as.genclone(vcfR2genind(v))
pop(vgc) <- sapply(strsplit(row.names(vgc@tab),split="_"),"[",1)
dist <- dist.gene(genind2df(vcfR2genind(v)), method="perc", pairwise.deletion=T)
poppr.cutoff <- cutoff_predictor(mlg.filter(vgc, distance=dist, stats="THRESHOLDS", threshold=1))
tiff("~/mfs/current/oligactis_RADSeq/analysis/RAD1/04_populations/results/spectrum.tif",compression="lzw")
hist(dist, breaks=100, col="slategrey", xlab="Genetic distance", main="Spectrum of genetic diversity");
abline(v=poppr.cutoff, col=2, lty=2, lwd=2)
dev.off()
mlg.filter(vgc, distance=dist) <- poppr.cutoff
mlgs <- data.frame(polyp.ID=row.names(vgc@tab), mlg.ID=vgc@mlg@mlg[,3])

mlgs <- mlgs[match(ids$ID2, mlgs$polyp.ID),]
identical(ids$ID2, mlgs$polyp.ID)
ids$poppr.mlg.ID <- mlgs$mlg.ID

tiff("~/mfs/current/oligactis_RADSeq/analysis/RAD1/04_populations/results/msn.tif",height=6, width=6, units="in", res=300, compression="lzw")
pmsn <- poppr.msn(gid=vgc, distmat=dist, palette=c("#990F02", "#ED7509", "#A25D34", "#D0A380", "#1CCBEA", "#06126B", "#4B1840", "#000000", "#289085", "#0A5837", "#7BB708", "#F7FC74", "#FFD220"))
plot_poppr_msn(x=vgc, poppr_msn=pmsn, inds="x",nodescale=30,cex=0.5)
dev.off()
system("eog ~/mfs/current/oligactis_RADSeq/analysis/RAD1/04_populations/results/msn.tif")

#########################################################################################
############################# Colony ####################################################
if(FALSE){
    dir.create("relatedness"); setwd("relatedness")
    vdf <- genind2df(vcfR2genind(v), oneColPerAll=T)
    vdf <- as.data.frame(lapply(vdf, as.numeric), row.names=row.names(vdf)) + 1

    er1 <- 0.036; er2 <- 0.010
    colony.base <- "all_snps"
    source("~/mfs/current/oligactis_RADSeq/analysis/RAD1/scripts/write_colony.R")
    system("nohup ~/mfs/current/oligactis_RADSeq/soft/COLONY/colony_v2.0.6.6/colony2s.ifort.out &")
}

setwd("~/mfs/current/oligactis_RADSeq/analysis/RAD1/04_populations/relatedness")
fsd <- read.csv("all_snps.FullSibDyad")
fsd$Pop1 <- substr(fsd$OffspringID1, 1, 3)
fsd$Pop2 <- substr(fsd$OffspringID2, 1, 3)
fsd[fsd$Pop1!=fsd$Pop2 & fsd$Probability>0.9,]

hsd <- read.csv("all_snps.HalfSibDyad")
hsd$Pop1 <- substr(hsd$OffspringID1, 1, 3)
hsd$Pop2 <- substr(hsd$OffspringID2, 1, 3)
hsd[hsd$Pop1!=hsd$Pop2 & hsd$Probability>0.9,]

system("cat all_snps.BestCluster | sed -E 's/([[:blank:]])+/;/g' > all_snps.BestClusterMod")
system("sed -i 's/^;//' all_snps.BestClusterMod")
clust <- read.csv("all_snps.BestClusterMod",sep=";")

system("cat all_snps.BestConfig | sed -E 's/([[:blank:]])+/;/g' > all_snps.BestConfigMod")
system("sed -i 's/^;//' all_snps.BestConfigMod")
col <- read.csv("all_snps.BestConfigMod",sep=";")

system("cat all_snps.BestClone | sed -E 's/([[:blank:]])+/;/g' > all_snps.BestCloneMod")
system("sed -i 's/^;//' all_snps.BestCloneMod")
system("sed -i 's/;$//' all_snps.BestCloneMod")
cl <- read.csv("all_snps.BestCloneMod",sep=";")
names(cl) <- c("cloneID","Clone.prob","Group")

cl$N_polyp <- sapply(strsplit(as.character(cl$Group), split=","), length)
cl$Cluster <- col$ClusterIndex[match(cl$cloneID, col$CloneIndex)]
cl$Cluster.prob <- clust$Probability[match(cl$Cluster, clust$ClusterIndex)]

clone.summary <- cl
clone.summary$Pop <- sapply(strsplit(as.character(cl$Group),split="_"),"[",1)
clone.summary$types <- NULL
for(i in 1:nrow(clone.summary)){
    cc <- unlist(strsplit(as.character(clone.summary$Group[i]), split=","))
    clone.summary$types[i] <- paste(table(ids$Sex[ids$ID2%in%cc]), collapse="/")
}

clone.summary$Group <- NULL
clone.summary <- clone.summary[c("Pop","cloneID","Clone.prob","Cluster","Cluster.prob","N_polyp","types")]
clone.summary$types <- gsub("0","-",clone.summary$types)

clone.summary$Clone.prob <- sprintf("%.2f",round(as.numeric(as.character(clone.summary$Clone.prob)),2))
clone.summary$Cluster.prob <- sprintf("%.2f",round(as.numeric(as.character(clone.summary$Cluster.prob)),2))
cs <- data.frame(Pop=clone.summary$Pop,
                 Clone=paste(clone.summary$cloneID, " (",clone.summary$Clone.prob, ")",sep=""),
                 Cluster=paste(clone.summary$Cluster, " (",clone.summary$Cluster.prob, ")",sep=""),
                 N=paste(clone.summary$N_polyp, " (", clone.summary$types, ")", sep=""))

write.table(cs, file="~/mfs/current/oligactis_RADSeq/analysis/RAD1/04_populations/results/table1.csv",col.names=T, row.names=F, quote=F, sep=";")

col <- col[match(ids$ID2, col$OffspringID),]
ids$Colony.ClusterID <- col$ClusterIndex
ids$Colony.CloneID <- col$CloneIndex

###

clone.summary$sex <- sapply(lapply(lapply(strsplit(gsub("-","0",clone.summary$types), split="/"),"[",1:3),as.numeric), sum)
clone.summary$asex <- sapply(lapply(lapply(strsplit(gsub("-","0",clone.summary$types), split="/"),"[",4:5),as.numeric), sum)

cl.mat <- tapply(clone.summary$N_polyp, list(clone.summary$Pop, clone.summary$Cluster), sum)
cl.mat[is.na(cl.mat)] <- 0

pop.sex.ratio <- prop.table(cbind(tapply(clone.summary$sex, clone.summary$Pop, sum), tapply(clone.summary$asex, clone.summary$Pop, sum)),1)
clust.sex.ratio <- prop.table(cbind(tapply(clone.summary$sex, clone.summary$Cluster, sum), tapply(clone.summary$asex, clone.summary$Cluster, sum)),1)

type.cols <- c("#B7E183", "#E75555"); names(type.cols) <- unique(ids$Type)

hm <- Heatmap(cl.mat, cluster_rows=F, cluster_columns=F, col=colorRamp2(c(0, 2, 5, 40), grey(c(1,2/3,1/3,0))),
              heatmap_legend_param = list(title = "No. polyps"),rect_gp=gpar(col="grey"),
              row_title="Population",row_title_side="left",row_names_side="left",
              column_title="Family cluster",column_title_side="top",column_names_side="top",column_names_rot=0,
              right_annotation=rowAnnotation(psr=anno_barplot(pop.sex.ratio,axis=F,gp=gpar(fill=rev(type.cols))),
                                             show_annotation_name=F,width=unit(2,"cm")),
              bottom_annotation=HeatmapAnnotation(csr=anno_barplot(clust.sex.ratio,axis=F,gp=gpar(fill=rev(type.cols))),
                                                  show_annotation_name=F,height=unit(2,"cm")),
              show_heatmap_legend=F,
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                  if(cl.mat[i,j]>0) grid.text(paste("N=",cl.mat[i,j],sep=""), x, y,
                                              gp=gpar(col="white",font=2))
              }
              )

lgd_list=list(
    Legend(labels=c("Sexual","Asexual"), title="Reproductive mode", legend_gp=gpar(fill=c("#E75555", "#B7E183"))),
    Legend(col_fun=colorRamp2(c(0, 2, 5, 40), grey(c(1,2/3,1/3,0))), title="No. polyps", at=c(0,2,5,40))
)

tiff("~/mfs/current/oligactis_RADSeq/analysis/RAD1/04_populations/results/clmat.tif",height=6,width=9,units="in",res=300,compression="lzw")
draw(hm,heatmap_legend_list=lgd_list)
dev.off()

system("eog ~/mfs/current/oligactis_RADSeq/analysis/RAD1/04_populations/results/clmat.tif")

### Does the Spectrum of Genetic Diversity method agree with Colony v. 2.0.6.6 in inferring clones?
identical(
    lapply(lapply(unique(ids$Colony.CloneID), function(ii) ids$ID[ids$Colony.CloneID==ii]),sort),
    lapply(lapply(unique(ids$poppr.mlg.ID), function(ii) ids$ID[ids$poppr.mlg.ID==ii]),sort)
)
### In how many populations does a clone occur?
apply(table(ids$Colony.CloneID, ids$Pop_code)>0,1,sum)

### Distribution of family clusters by population:
table(ids$Colony.ClusterID, ids$Pop_code)

setwd("../")

#################################################################################################
################################ Basic population genetic stats #################################
system("mkdir popgen")
file.copy("final.vcf","./popgen/final.vcf")

setwd("popgen")
v <- read.vcfR("final.vcf")
popmap <- read.csv("~/mfs/current/oligactis_RADSeq/analysis/RAD1/rad1_popmap_final_noreps.txt",header=F,sep="\t")
popmap[,1] <- gsub("_rep.unmapped|.unmapped","",popmap[,1])
popmap <- popmap[popmap[,1]%in%colnames(v@gt)[-1],]
write.table(popmap, file="popmap_full.txt", sep="\t", col.names=F,row.names=F,quote=F)

### Thin dataset to randomly remove extra individuals from M28 as well as clones from all pops ##

set.seed(1234) ## randomly remove 25 individuals from M28 such that sample size / pop is not very different
sample.set <- colnames(v@gt)[-1]
thinned <- sample(sample.set[grep("M28",sample.set)], 25, replace=F)
red.sample.set <- sample.set[!sample.set%in%thinned]
cat(red.sample.set, file="red.sample.set.txt", sep="\n")
write.table(popmap[popmap[,1]%in%red.sample.set,], file="popmap_red.txt", sep="\t", col.names=F,row.names=F,quote=F)
system(paste("vcftools --vcf final.vcf --keep red.sample.set.txt --out thinned --recode --recode-INFO-all"))

v.red <- read.vcfR("thinned.recode.vcf")
ids.red <- ids[ids$ID2%in%colnames(v.red@gt)[-1],]
red.sample.set.noclones <- ids.red$ID2[!duplicated(ids.red$poppr.mlg.ID)]
cat(red.sample.set.noclones, file="red.sample.set.noclones.txt", sep="\n")
write.table(popmap[popmap[,1]%in%red.sample.set.noclones,], file="popmap_red_noclones.txt", sep="\t", col.names=F,row.names=F,quote=F)
system(paste("vcftools --vcf thinned.recode.vcf --keep red.sample.set.noclones.txt --out thinned_noclones --recode --recode-INFO-all"))
system("rm ./red.sample.set.*")
system("rm ./*.log")

system("populations -V thinned_noclones.recode.vcf -O ./ -M ./popmap_red_noclones.txt --fstats")
temp <- readLines("thinned_noclones.recode.p.sumstats_summary.tsv")
stats <- gsub("# ","",temp[(grep("All positions",temp)+1):length(temp)])
stats.df <- NULL
for(i in 1:30) stats.df <- cbind(stats.df, sapply(sapply(stats, strsplit, split="\t",USE.NAMES=F),"[",i))
colnames(stats.df) <- stats.df[1,]
stats.df <- as.data.frame(stats.df[(-1),][,c(1,2,13,19,28)])
stats.df$N <- table(ids$Pop_code)
stats.df$MLGs <- sapply(tapply(ids$Colony.CloneID, ids$Pop_code, table),length)
write.table(stats.df, file="~/mfs/current/oligactis_RADSeq/analysis/RAD1/04_populations/results/table2.csv",
            col.names=T, row.names=F, quote=F, sep=";")

###########################################################################################
################### Genetic differentiation ###############################################
system("mkdir popgenome")
ids.red2 <- ids[ids$ID2%in%red.sample.set.noclones,]
ids.red2 <- ids.red2[!ids.red2$Pop_code%in%names(which(table(ids.red2$Pop_code)<3)),]
cat(ids.red2$ID2, file="red.sample.set.noclones_pop3.txt", sep="\n")
write.table(popmap[popmap[,1]%in%ids.red2$ID2,], file="popmap_red_noclones_pop3.txt", sep="\t", col.names=F,row.names=F,quote=F)
system(paste("vcftools --vcf thinned_noclones.recode.vcf --keep red.sample.set.noclones_pop3.txt --out ./popgenome/pop3 --recode --recode-INFO-all"))
system("rm ./red.sample.set.*")
system("rm ./popgenome/*.log")

x=readData("./popgenome",format="VCF",include.unknown=T)
inds <- get.individuals(x)[[1]]; inds <- inds[!grepl(".2",inds,fixed=T)]

x <- set.populations(x, tapply(inds, substr(inds,1,3),list), diploid=T)
x <- F_ST.stats(x, mode="nucleotide")

x@nucleotide.F_ST

fst.mat <- as.data.frame(cbind(sapply(strsplit(row.names(x@nuc.F_ST.pairwise),split="/"),"[",1),
                               sapply(strsplit(row.names(x@nuc.F_ST.pairwise),split="/"),"[",2),
                               x@nuc.F_ST.pairwise))
row.names(fst.mat) <- NULL
fst.mat$V1 <- factor(fst.mat$V1,levels=c("pop1","pop2","pop3","pop4","pop5","pop6","pop7","pop8","pop9","pop10"))
fst.mat$V2 <- factor(fst.mat$V2,levels=c("pop2","pop3","pop4","pop5","pop6","pop7","pop8","pop9","pop10","pop11"))
fst.mat$V3 <- as.numeric(fst.mat$V3)

pops <- names(x@populations); pops[1] <- "M108"
pairwise.fst <- xtabs(V3~V1+V2, data=fst.mat)
colnames(pairwise.fst) <- pops[-1]
row.names(pairwise.fst) <- pops[-11]
write.table(pairwise.fst, file="~/mfs/current/oligactis_RADSeq/analysis/RAD1/04_populations/results/table3.csv",
            col.names=T, row.names=F, quote=F, sep=";")

###########################################################################################
#### Relatedness ##########################################################################

library(EcoGenetics)
vgi <- vcfR2genind(v)
veg <- genind2ecogen(vgi)
coef <- eco.kin.loiselle(veg)

ids.spatial <- ids$ID2[order(as.numeric(ids$East))]
coef2 <- coef[match(ids.spatial,colnames(coef)),]
coef3 <- coef2[,match(ids.spatial,colnames(coef2))]

### Create matrix of Colony-inferred relationships
relatedness <- matrix(data=0, nrow=nrow(coef3), ncol=ncol(coef3), dimnames=dimnames(coef3))
diag(relatedness) <- NA
clones <- read.csv("~/mfs/current/oligactis_RADSeq/analysis/RAD1/04_populations/relatedness/all_snps.BestCloneMod",sep=";")[,3]
clones <- strsplit(clones[grep(",",clones)],",")
clone.dyads <- unlist(lapply(clones, combn, 2, simplify=F),recursive=F)
clone.dyads <- data.frame(left=sapply(clone.dyads,"[",1), right=sapply(clone.dyads,"[",2))
fullsib.dyads <- read.csv("~/mfs/current/oligactis_RADSeq/analysis/RAD1/04_populations/relatedness/all_snps.FullSibDyad")[,1:2]
halfsib.dyads <- read.csv("~/mfs/current/oligactis_RADSeq/analysis/RAD1/04_populations/relatedness/all_snps.PairwiseHalfSibDyad")[,1:2]

for(i in 1:nrow(clone.dyads)){
    relatedness[match(clone.dyads[,1][i], colnames(relatedness)), match(clone.dyads[,2][i], row.names(relatedness))] <- 0.5
    relatedness[match(clone.dyads[,2][i], colnames(relatedness)), match(clone.dyads[,1][i], row.names(relatedness))] <- 0.5
}
for(i in 1:nrow(fullsib.dyads)){
    print(fullsib.dyads[i,])
    relatedness[match(fullsib.dyads[,1][i], colnames(relatedness)), match(fullsib.dyads[,2][i], row.names(relatedness))] <- 0.25
    relatedness[match(fullsib.dyads[,2][i], colnames(relatedness)), match(fullsib.dyads[,1][i], row.names(relatedness))] <- 0.25
    print(relatedness[aa,bb])
}
for(i in 1:nrow(halfsib.dyads)){
    relatedness[match(halfsib.dyads[,1][i], colnames(relatedness)), match(halfsib.dyads[,2][i], row.names(relatedness))] <- 0.125
    relatedness[match(halfsib.dyads[,2][i], colnames(relatedness)), match(halfsib.dyads[,1][i], row.names(relatedness))] <- 0.125
}

coef3[upper.tri(coef)] <- relatedness[upper.tri(relatedness)]

pops <- sapply(strsplit(ids$ID2, split="_"),"[",1)
pop.order <- factor(pops, levels=c("M72","M67","M90","M85","M28","M47","M44","M26","M108","M25","M52","M38","R12"))
pop.groups <- c(0, cumsum(table(pop.order))) / length(pop.order)
anno <- rep("",length(ids.spatial)); anno[c(4,11,18,24,45,63,69,79,86,91,95,100,104)] <- levels(pop.order)
topAnno <- HeatmapAnnotation(foo = anno_text(anno, gp=gpar(fontsize=10,font=4,col="black")),which="column")

tiff("~/mfs/current/oligactis_RADSeq/analysis/RAD1/04_populations/results/kinship_matrix.tif", height=5.032, width=6.29, units="in", compression="lzw",res=300)

type.cols <- c("#B7E183", "#E75555"); names(type.cols) <- unique(ids$Type)

col.loiselle <- colorRamp2(c(min(coef3,na.rm=T), 0, 0.125, 0.25, 0.5), grey(c(1,1,2/3,1/3,0)))
col.colony <- colorRamp2(c(0, 0.125, 0.25, 0.5), c("#C4E6FF", "#71ACD6", "#035C96", "#03204C"))

ht <- Heatmap(coef3,cluster_rows=F,cluster_columns=F,
              show_heatmap_legend=F,
        column_names_max_height=unit(100,"mm"),
        bottom_annotation=HeatmapAnnotation(foo = anno_text(anno, gp=gpar(fontsize=10,font=4,col="black")),which="column"),
        left_annotation=HeatmapAnnotation(foo = anno_text(anno, gp=gpar(fontsize=10,font=4,col="black"),
                                                          location=unit(-0.1,"npc")),which="row"),
        show_row_names=F, show_column_names=F,name="Kinship coefficient",rect_gp = gpar(type = "none"),
        heatmap_legend_param = (list(Type=list(title="Reproductive mode",at=c("asex","sex"),labels=c("asex","sex")))),
        cell_fun = function(j, i, x, y, w, h, fill) {
            if(i > j) {
                grid.rect(x, y, w, h, gp = gpar(fill = col.loiselle(coef3[i,j]),col=0))
            } else if(j == i) {
                grid.rect(x, y, w, h, gp = gpar(fill = 0, col=0))
            } else {
                grid.rect(x, y, w, h, gp = gpar(fill = col.colony(coef3[i,j]),col=0))
            }
        }) +
    rowAnnotation(Type=anno_simple(ids$Type[match(ids.spatial, ids$ID2)],col=type.cols,border=T,simple_anno_size=unit(0.5,"cm")))

lgd_list=list(
    Legend(labels=c("Sexual","Asexual"), title="Reproductive mode", legend_gp=gpar(fill=c("#E75555", "#B7E183"))),
    Legend(labels=c("Clones","Full-sibs","Half-sibs","Unrelated"), legend_gp=gpar(fill=col.colony(c(0.5,0.25,0.125,0))),
           title="Kinship class"),
    Legend(col_fun=col.loiselle, title="Kinship coefficient")
)

draw(ht, heatmap_legend_list=lgd_list)

decorate_heatmap_body("Kinship coefficient", {
    for(i in 1:length(pop.groups)) grid.lines(c(pop.groups[i],pop.groups[i]),c(0,1), gp = gpar(lwd = 1, col="black"))
    for(i in 1:length(pop.groups)) grid.lines(c(0,1), 1-c(pop.groups[i],pop.groups[i]), gp = gpar(lwd = 1, col="black"))
})

dev.off()
system("eog ~/mfs/current/oligactis_RADSeq/analysis/RAD1/04_populations/results/kinship_matrix.tif")

#######################################################################################
########################### DAPC ######################################################
v <- read.vcfR("~/mfs/current/oligactis_RADSeq/analysis/RAD1/04_populations/popgen/thinned_noclones.recode.vcf")

vgi <- vcfR2genind(v)
pop(vgi) <- sapply(strsplit(row.names(vgi@tab),split="_"),"[",1)

dapc.vgi <- dapc(vgi, n.pca=150, n.da=100)
temp <- optim.a.score(dapc.vgi)
dapc.vgi <- dapc(vgi, n.pca=temp$best,n.da=100)
tiff("~/mfs/current/oligactis_RADSeq/analysis/RAD1/04_populations/results/dapc.tif",width=1750, height=1750, res=300,compression="lzw")
scatter(dapc.vgi,clab=0,leg=TRUE,posi.da="topleft",txt.leg=sort(unique(pop(vgi))),
        col=c("#990F02", "#ED7509", "#A25D34", "#D0A380", "#1CCBEA", "#06126B", "#4B1840",
              "#000000","#289085", "#0A5837", "#7BB708", "#F7FC74", "#FFD220"),
        pch=19, cex=1.2)
dev.off()
system("eog ~/mfs/current/oligactis_RADSeq/analysis/RAD1/04_populations/results/dapc.tif")
#######################################################################################
####### Reproductive mode vs genetic structure ########################################

table(ids$Colony.ClusterID, ids$Type)
fisher.test(table(ids$Colony.ClusterID, ids$Type))

v <- read.vcfR("~/mfs/current/oligactis_RADSeq/analysis/RAD1/04_populations/final.vcf")
dist <- dist.gene(genind2df(vcfR2genind(v)), method="perc", pairwise.deletion=T)

ids.ordered <- ids[match(attr(dist,"Labels"), ids$ID2),]
identical(attr(dist, "Labels"), ids.ordered$ID2)
type <- as.factor(ids.ordered$Type)

pegas::amova(dist~type,nperm=1000)

