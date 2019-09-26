setwd("~/Documents/WORK/Arabis_resequencing/V5.1/Genotyping/New/zygosity/recalculate4_two_tailed/")

files <- list.files(path="~/Documents/WORK/Arabis_resequencing/V5.1/Genotyping/New/zygosity/output_latest_2to1_2to1_twotailed/", pattern="insertion.txt")
dati1=read.table("~/Documents/WORK/Arabis_resequencing/V5.1/Genotyping/New/zygosity/reformat_genotyped_ins_with_second_pass1.bed",header=T,sep='\t',stringsAsFactors = F)

dati1$V1=gsub( "I", "", as.character(dati1$V1))
dati1$V1=as.numeric(dati1$V1)
DF=as.data.frame(dati1[,c(1)],stringsAsFactors=F)
names(DF)="V1"
setwd("~/Documents/WORK/Arabis_resequencing/V5.1/Genotyping/New/zygosity/output_latest_2to1_2to1_twotailed/")


DF2=NULL
#DF <- NULL
for (f in files) {
 
  
  dat <- read.table(f, header=F, sep="\t")
  fn=strsplit(f,".insertion.txt")
  DF2=merge(DF,dat,by='V1',all.x=T)
  nam=c(names(DF),fn)
  DF =DF2
  names(DF)=nam 
}
setwd("~/Documents/WORK/Arabis_resequencing/V5.1/Genotyping/New/zygosity/recalculate4_two_tailed/")
#dti2=read.table("~/Documents/WORK/Arabis_resequencing/V5.1/Genotyping/New/zygosity/reformat_genotyped_ins_with_second_pass_w_header.bed",header=T,sep='\t')
#heads=names(dti2)
#dput(heads)


#mdf=mapply("*", DF[,c(2:182)], dti[,c(16:319)])

dti=dati1[,c(1,16:319)]
names(dti)=c("V1","20150213.B-Ma28", "20150213.B-Ma29", "20150213.B-Ma30", "20150213.B-Ma31", 
                        "20150213.B-Ma32", "20150213.B-Ma33", "20150213.B-Ma34", "20150213.B-Ma36", 
                        "20150213.B-Ma37", "20150213.B-Ma38", "20150213.B-Ma39", "20150213.B-Ma40", 
                        "20150213.B-Ma41", "20150213.B-Ma42", "20150213.B-Ma43", "20150213.B-Ma44", 
                        "20150213.B-Ma45", "20150213.B-Ma46", "20150213.B-Ma47", "20150213.B-Ma48", 
                        "20150213.B-Ma49", "20150213.B-Ma50", "20150213.B-Ma51", "20150213.B-Ma52", 
                        "20150213.B-Ma53", "20150213.B-Ma54", "20150213.B-Ma55", "20150213.B-Ma57", 
                        "20150213.B-Ma58", "20150213.B-Ma59", "20150213.B-Ma60", "20150213.B-Ma62", 
                        "20150213.B-Ma63", "20150213.B-Ma65", "20150213.B-Ma66", "20150213.B-Ma67", 
                        "20150213.B-Ma68", "20150213.B-Ma70", "20150213.B-Ma71", "20150213.B-Ma74", 
                        "20150213.B-Ma75", "comb.Es100", "comb.Es12", "comb.Es13", "comb.Es14", 
                        "comb.Es16", "comb.Es17", "comb.Es19", "comb.Es20", "comb.Es21", 
                        "comb.Es22", "comb.Es23", "comb.Es24", "comb.Es25", "comb.Es32", 
                        "comb.Es33", "comb.Es34", "comb.Es35", "comb.Es36", "comb.Es38", 
                        "comb.Es39", "comb.Es3", "comb.Es40", "comb.Es41", "comb.Es42", 
                        "comb.Es43", "comb.Es44", "comb.Es46", "comb.Es48", "comb.Es49", 
                        "comb.Es50", "comb.Es52", "comb.Es53", "comb.Es55", "comb.Es56", 
                        "comb.Es57", "comb.Es59", "comb.Es60", "comb.Es63", "comb.Es64", 
                        "comb.Es65-1", "comb.Es66", "comb.Es67", "comb.Es68", "comb.Es6", 
                        "comb.Es71", "comb.Es72", "comb.Es73", "comb.Es75", "comb.Es76", 
                        "comb.Es77", "comb.Es78", "comb.Es79", "comb.Es7", "comb.Es81", 
                        "comb.Es82", "comb.Es83", "comb.Es84", "comb.Es85", "comb.Es86", 
                        "comb.Es88", "comb.Es8", "comb.Es90", "comb.Es91", "comb.Es92", 
                        "comb.Es93", "comb.Es94", "comb.Es95", "comb.Es98", "comb.Es99", 
                        "comb.Es9", "comb.Pa100", "comb.Pa11", "comb.Pa12", "comb.Pa13", 
                        "comb.Pa14", "comb.Pa16", "comb.Pa19", "comb.Pa1", "comb.Pa20", 
                        "comb.Pa23", "comb.Pa24", "comb.Pa26", "comb.Pa27", "comb.Pa28", 
                        "comb.Pa29", "comb.Pa2", "comb.Pa30", "comb.Pa31", "comb.Pa32", 
                        "comb.Pa33", "comb.Pa35", "comb.Pa36", "comb.Pa37", "comb.Pa3", 
                        "comb.Pa40", "comb.Pa42", "comb.Pa43", "comb.Pa45", "comb.Pa46", 
                        "comb.Pa47", "comb.Pa48", "comb.Pa49", "comb.Pa51", "comb.Pa52", 
                        "comb.Pa53", "comb.Pa55", "comb.Pa56", "comb.Pa57", "comb.Pa58", 
                        "comb.Pa5", "comb.Pa60", "comb.Pa61", "comb.Pa63", "comb.Pa64", 
                        "comb.Pa65", "comb.Pa66", "comb.Pa67", "comb.Pa68", "comb.Pa69", 
                        "comb.Pa6", "comb.Pa71", "comb.Pa72", "comb.Pa74", "comb.Pa75", 
                        "comb.Pa76", "comb.Pa77", "comb.Pa78", "comb.Pa84", "comb.Pa85", 
                        "comb.Pa86", "comb.Pa87", "comb.Pa88", "comb.Pa89", "comb.Pa90", 
                        "comb.Pa91", "comb.Pa92", "comb.Pa94", "comb.Pa96", "comb.Pa98", 
                        "comb.Pi11", "comb.Pi12", "comb.Pi13", "comb.Pi14", "comb.Pi15", 
                        "comb.Pi18", "comb.Pi19", "comb.Pi1", "comb.Pi20", "comb.Pi21", 
                        "comb.Pi26", "comb.Pi27", "comb.Pi28", "comb.Pi29", "comb.Pi2", 
                        "comb.Pi31", "comb.Pi33", "comb.Pi34", "comb.Pi36", "comb.Pi39", 
                        "comb.Pi3", "comb.Pi40", "comb.Pi41", "comb.Pi42", "comb.Pi43", 
                        "comb.Pi44", "comb.Pi49", "comb.Pi51", "comb.Pi52", "comb.Pi54", 
                        "comb.Pi55", "comb.Pi56", "comb.Pi57", "comb.Pi58", "comb.Pi59", 
                        "comb.Pi5", "comb.Pi60", "comb.Pi62", "comb.Pi65", "comb.Pi67", 
                        "comb.Pi68", "comb.Pi71", "comb.Pi72", "comb.Pi73", "comb.Pi74", 
                        "comb.Pi75", "comb.Pi76", "comb.Pi77", "comb.Pi78", "comb.Pi79", 
                        "comb.Pi7", "comb.Pi80", "comb.Pi81", "comb.Pi82", "comb.Pi83", 
                        "comb.Pi84", "comb.Pi85", "comb.Pi86", "comb.Pi87", "comb.Pi88", 
                        "comb.Pi89", "comb.Pi90", "comb.Pi92", "comb.Pi93", "comb.Pi94", 
                        "comb.Pi96", "comb.Pi97", "comb.Pi98", "comb.Pi9", "merged.Ma100", 
                        "merged.Ma10", "merged.Ma11", "merged.Ma12", "merged.Ma13", "merged.Ma14", 
                        "merged.Ma15", "merged.Ma16", "merged.Ma17", "merged.Ma18", "merged.Ma19", 
                        "merged.Ma20", "merged.Ma21", "merged.Ma22", "merged.Ma23", "merged.Ma24", 
                        "merged.Ma26", "merged.Ma27", "merged.Ma2", "merged.Ma35", "merged.Ma3", 
                        "merged.Ma4", "merged.Ma56", "merged.Ma5", "merged.Ma61", "merged.Ma69", 
                        "merged.Ma6", "merged.Ma72", "merged.Ma73", "merged.Ma76", "merged.Ma77", 
                        "merged.Ma78", "merged.Ma79", "merged.Ma80", "merged.Ma81", "merged.Ma82", 
                        "merged.Ma83", "merged.Ma84", "merged.Ma85", "merged.Ma86", "merged.Ma87", 
                        "merged.Ma88", "merged.Ma89", "merged.Ma8", "merged.Ma90", "merged.Ma91", 
                        "merged.Ma92", "merged.Ma93", "merged.Ma94", "merged.Ma95", "merged.Ma96", 
                        "merged.Ma97", "merged.Ma98", "merged.Ma99", "merged.Ma9")

mdf=mapply("*", DF[intersect(names(DF), names(dti))],
       dti[intersect(names(dti), names(DF))])
mdat=as.data.frame(mdf)
mdat$V1=sqrt(mdat$V1)
mdat2=mdat*2
mdat2$V1=mdat2$V1/2

mdat2$AF1=rowSums(mdat2[2:305],na.rm =T)/(2*apply(mdat2[2:305], 1, function(x) length(which(!is.na(x)))))
mdat2$AF0=(((2*apply(mdat2[2:305], 1, function(x) length(which(!is.na(x))))))-rowSums(mdat2[2:305],na.rm =T))/((2*apply(mdat2[2:305], 1, function(x) length(which(!is.na(x))))))
mdat2$MAF=ifelse(mdat2$AF0<mdat2$AF1,mdat2$AF0,mdat2$AF1)
mdat2$minimumallele=ifelse(mdat2$AF0<mdat2$AF1,"absence","presence")
mdat2$Missing = (apply(mdat2[,2:305], 1, function(x) sum(is.na(x))))/304
mdat2$V1=paste("I",mdat2$V1,sep='')
write.table(mdat2,"Genotyped_ins_with_zyg.txt",row.names = F,sep='\t',quote = F)

mdat3=mdat2[mdat2$MAF>=0.025,]
mdat4=mdat3[mdat3$Missing<=0.1,]
write.table(mdat4,"Genotyped_ins_with_zyg_filtered_MAF0.025_Missing0.1.txt",row.names = F,sep='\t',quote = F)

