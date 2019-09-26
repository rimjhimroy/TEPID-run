#run in folder where there are tables of split/discordant and properly mapped reads
path = "~/Documents/WORK/Arabis_resequencing/V5.1/Genotyping/New/zygosity/output2/"
outpath="~/Documents/WORK/Arabis_resequencing/V5.1/Genotyping/New/zygosity/output_latest_2to1_2to1_twotailed/"
file.names <- dir(path, pattern ="insertion.txt")
setwd(path)
bt <- function(a, b,p = 0.5) {binom.test(a, b, 0.5, alternative=
                                                "two.sided", conf.level = 0.95)$p.value}

for(i in 1:length(file.names)){
  
outfile=paste(outpath,file.names[i],sep='')
dat=read.table(file.names[i],header=F)
dat$sup=as.integer((dat$V2/2)+1) # 1 is added just to round off to the greater number
dat$nsup=dat$V3
dat$sum=dat$sup+dat$nsup

dat$pval[dat$sum!=0]= mapply(bt, dat$sup[dat$sum!=0],dat$sum[dat$sum!=0])
dat$zyg[dat$sum!=0&dat$V4<0.5&dat$pval<0.05]=0
dat$zyg[dat$sum!=0&dat$V4>0.5&dat$pval<0.05]=1
dat$zyg[dat$sum!=0&dat$pval>0.05]=0.5
dat$zyg[dat$sum==0]=NA
dat$zyg[dat$sum<=5]=NA

dat2=dat[,c(1,9)]
write.table(dat2,outfile,row.names = F,col.names = F,quote = F,sep='\t')

}

file.names <- dir(path, pattern ="absence.txt")
setwd(path)
bt <- function(a, b, alt,p = 0.5) {binom.test(a, b, 0.5, alternative=
                                                "two.sided", conf.level = 0.95)$p.value}

for(i in 1:length(file.names)){
  
  outfile=paste(outpath,file.names[i],sep='')
  dat=read.table(file.names[i],header=F)
  dat$sup=dat$V2
  #dat$nsup=as.integer((dat$V3/2)+1) # 1 is added just to round off to the greater number
  dat$nsup=as.integer((dat$V3/2))
  dat$sum=dat$sup+dat$nsup
  
  dat$pval[dat$sum!=0]= mapply(bt, dat$sup[dat$sum!=0],dat$sum[dat$sum!=0])
  dat$zyg[dat$sum!=0&dat$V4<0.5&dat$pval<0.05]=0
  dat$zyg[dat$sum!=0&dat$V4>0.5&dat$pval<0.05]=1
  dat$zyg[dat$sum!=0&dat$pval>0.05]=0.5
  dat$zyg[dat$sum==0]=NA
  dat$zyg[dat$sum<=5]=NA
  
  dat2=dat[,c(1,9)]
  write.table(dat2,outfile,row.names = F,col.names = F,quote = F,sep='\t')
  
}

