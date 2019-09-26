library(pegas)
library(hierfstat)
library(adegenet)

setwd("~/Documents/WORK/Arabis_resequencing/V5.1/Genotyping/New/zygosity/recalculate4_two_tailed/")
ins=read.table("Genotyped_ins_with_zyg_filtered_MAF0.025_Missing0.1.txt",header=T,stringsAsFactors = F)
abs=read.table("Genotyped_del_with_zyg_filtered_MAF0.025_Missing0.1.txt",header=T,stringsAsFactors = F)
tot=rbind(ins,abs)

cl=read.table("result_reTribe_file.txt",header=T)
tot=unique(tot)
tot1=merge(tot,cl,by="V1",all.x=T)
tot1$V1=paste(tot1$V1,tot1$FTRIBE,sep='_')
tot2=tot1[,c(1:305)]
write.table(tot2,"Genotyped_TOTAL_with_zyg_filtered_MAF0.025_Missing0.1_Tribe.txt",row.names = F,sep='\t',quote = F)
write.table(tot1,"Genotyped_TOTAL_with_zyg_filtered_MAF0.025_Missing0.1_Tribe_FULL.txt",row.names = F,sep='\t',quote = F)
tot3=tot1[tot1$FAMILY!="Clust3416",]
tot2=tot3[,c(1:305)]
write.table(tot2,"Genotyped_TOTAL_with_zyg_filtered_MAF0.025_Missing0.1_Tribe_wo_Clust3416.txt",row.names = F,sep='\t',quote = F)

tot=read.table("Genotyped_TOTAL_with_zyg_filtered_MAF0.025_Missing0.1_Tribe_wo_Clust3416.txt",header=T,stringsAsFactors = F)
###
tot1=tot[ , c("V1","merged.Ma76", "merged.Ma73", "merged.Ma72", "merged.Ma69", "merged.Ma61", "merged.Ma56", "merged.Ma35", "comb.Pi85", "comb.Pi76", "comb.Pi60", "comb.Pi5", "comb.Pi29", "comb.Pi13", "comb.Pi12", "comb.Pi11", "comb.Es9", "comb.Es63", "comb.Es60", "comb.Es36", "comb.Es3", "comb.Es14", "comb.Es13", "comb.Es12", "comb.Pa27", "comb.Pa1", "comb.Pa16", "comb.Pa29", "comb.Pa98", "comb.Pa60", "comb.Pa6")]
tot2=tot1[ tot1$V1 %in% c("A15670_Unknown","A24383_Unknown","A14720_Unknown","I4157_ATGP1_Gypsy"), ]
###
hier=read.table("Plots.txt",header=T,stringsAsFactors = F)
names(hier)=c("ind","pop","plot")
mdat.T <- as.data.frame(t(tot[,2:305]))

mdat.T[mdat.T == 1] <- 12
mdat.T[mdat.T == 0] <- 11
mdat.T[mdat.T == 2] <- 22
mdat.T$pop=1
mdat.T$pop[grepl("2015",row.names(mdat.T))]="Ma"
mdat.T$pop[grepl("merged.Ma",row.names(mdat.T))]="Ma"
mdat.T$pop[grepl("comb.Es",row.names(mdat.T))]="Es"
mdat.T$pop[grepl("comb.Pa",row.names(mdat.T))]="Pa"
mdat.T$pop[grepl("comb.Pi",row.names(mdat.T))]="Pi"
mdat.T$ind=row.names(mdat.T)
mdat.T2=mdat.T[,c(ncol(mdat.T),(ncol(mdat.T)-1),1:(ncol(mdat.T)-2))]
names(mdat.T2)=c("ind","pop",tot$V1)
mdat.T2=mdat.T2[order(mdat.T2$pop),]
df <- mdat.T2[!is.na(names(mdat.T2))]
df1=merge(df,hier,by="ind")
df1$pop=df1$pop.x
df1$pop.x=NULL
df1$pop.y=NULL
write.table(df,"Hierfstat.tot.txt",row.names = F,sep='\t',quote = F)
mdat.T2=df1
mdat1=mdat.T2[mdat.T2$pop=='Ma',]
mdat1.1=mdat1[,apply(mdat1, 2, function(x) !all(is.na(x) | x == 11))]
mdat2=mdat.T2[mdat.T2$pop=='Es',]
mdat2.1=mdat2[,apply(mdat2, 2, function(x) !all(is.na(x) | x == 11))]
mdat3=mdat.T2[mdat.T2$pop=='Pa',]
mdat3.1=mdat3[,apply(mdat3, 2, function(x) !all(is.na(x) | x == 11))]
mdat4=mdat.T2[mdat.T2$pop=='Pi',]
mdat4.1=mdat4[,apply(mdat4, 2, function(x) !all(is.na(x) | x == 11))]
lst=list()


lst[1]=list(names(mdat1.1[,c(-1,-2)]))
lst[2]=list(names(mdat2.1[,c(-1,-2)]))
lst[3]=list(names(mdat3.1[,c(-1,-2)]))
lst[4]=list(names(mdat4.1[,c(-1,-2)]))

library(VennDiagram)
set.seed(1) # For reproducibility of results
xx.1 <- list(A = lst[1][[1]], B = lst[2][[1]], 
             C = lst[3][[1]], D = lst[4][[1]])

venn.diagram(xx.1, filename ="Temp.tiff", height = 2000, width = 2000,category.names=c("Ma","Es","Pa","Pi"),print.mode=(c("raw","percent")),sigdigs = 1)
dat.genind=df2genind(mdat.T2[,c(3:ncol(mdat.T2))],ncode=1,ind.names = mdat.T2$ind,pop = mdat.T2$pop,loc.names = colnames(mdat.T2[,c(3:ncol(mdat.T2))]))
dat.genind=df2genind(mdat.T2[,c(2:ncol(mdat.T2)-2)],ncode=1,ind.names = mdat.T2$ind,pop = mdat.T2$pop,loc.names = colnames(mdat.T2[,c(2:ncol(mdat.T2)-2)]),strata = data.frame(mdat.T2$plot,mdat.T2$pop))
hier(dat.genind) <- ~mdat.T2.plot/mdat.T2.pop


dat.genind1=df2genind(mdat1.1[,c(3:ncol(mdat1.1))],ncode=1,ind.names = mdat1.1$ind,pop = mdat1.1$pop,loc.names = colnames(mdat1.1[,c(3:ncol(mdat1.1))]))
dat.genind2=df2genind(mdat2.1[,c(3:ncol(mdat2.1))],ncode=1,ind.names = mdat2.1$ind,pop = mdat2.1$pop,loc.names = colnames(mdat2.1[,c(3:ncol(mdat2.1))]))
dat.genind3=df2genind(mdat3.1[,c(3:ncol(mdat3.1))],ncode=1,ind.names = mdat3.1$ind,pop = mdat3.1$pop,loc.names = colnames(mdat3.1[,c(3:ncol(mdat3.1))]))
dat.genind4=df2genind(mdat4.1[,c(3:ncol(mdat4.1))],ncode=1,ind.names = mdat4.1$ind,pop = mdat4.1$pop,loc.names = colnames(mdat4.1[,c(3:ncol(mdat4.1))]))

dat.hierfstat=genind2hierfstat(dat.genind)
dat.hierfstat1=genind2hierfstat(dat.genind1)
dat.hierfstat2=genind2hierfstat(dat.genind2)
dat.hierfstat3=genind2hierfstat(dat.genind3)
dat.hierfstat4=genind2hierfstat(dat.genind4)

table(dat.hierfstat[,1])
WC=wc(dat.hierfstat)
bs=basic.stats(dat.hierfstat)
Fis=apply(bs$Fis,2,mean,na.rm=T)
Ho=apply(bs$Ho,2,mean,na.rm=T)
Hs=apply(bs$Hs,2,mean,na.rm=T)


mean(bs$Hs)

##Pairwise Fst
#dat.pFst <- pairwise.fst(dat.genind,res.type="matrix")
dat.pWCFst <- pairwise.WCfst(dat.hierfstat)
geneD=genet.dist(dat.hierfstat,method="WC84")
dat.ppFst <- boot.ppfst(dat.genind,nboot=1000)
dat.ppfis=boot.ppfis(dat.genind,nboot=1000)
main_Ho_by_pop <- colMeans(bs$Ho, na.rm=TRUE)
main_Hs_by_pop <- colMeans(bs$Hs, na.rm=TRUE)
main_Fis_by_pop <- colMeans(bs$Fis, na.rm=TRUE)
main_Fis_by_pop <- colMeans(bs$Fst, na.rm=TRUE)
overall <- colMeans(bs$perloc)
icFis <- colMeans(dat.ppfis$fis.ci)
#datvcomp=varcomp.glob(dat.hierfstat$pop,dat.hierfstat[,-1])$F
datvcomp=varcomp.glob(dat.hierfstat$pop,dat.hierfstat[,-1])
vc=boot.vc(dat.hierfstat$pop,dat.hierfstat[,-1])
#vc1=boot.vc(dat.hierfstat1$pop,dat.hierfstat[,-1])
#vc2=boot.vc(dat.hierfstat2$pop,dat.hierfstat[,-1])
#vc3=boot.vc(dat.hierfstat3$pop,dat.hierfstat[,-1])
#vc4=boot.vc(dat.hierfstat4$pop,dat.hierfstat[,-1])

B=betas(dat.hierfstat,nboot=1000)
fst=rowMeans(B$ci)

#Fis icFis and dat.ppfis
#Fst B
genind_seploc <- seploc(dat.genind)
myfst <- lapply(genind_seploc, function(x) pairwise.fst(x))
genind_perLocusPWFst <- lapply(genind_seploc, pairwise.fst)
##PCA
df1 <- scaleGen(dat.genind, scale = FALSE, NA.method = "mean")
dat.genind@tab=df1
pcaX <- dudi.pca(df1,center=TRUE,scale=FALSE,nf=2,scannf=TRUE)
#select 2 or 3
barplot(pcaX$eig, main = "Eigenvalues")
#pcoX <- dudi.pco(dist(df1), scannf = TRUE, nf = 3)

## scatterplot

myCol <- rainbow(length(levels(pop(dat.genind))))
par(bg = "white")
s.class(pcaX$li, pop(dat.genind), col = myCol)
add.scatter.eig(pcaX$eig, 3, 1, 2)

load=loadingplot(pcaX$c1^2)
