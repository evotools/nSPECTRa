library(plyr)
library(reshape2)
library(ComplexHeatmap)

args <- commandArgs(trailingOnly = TRUE)
popu<-args[1]
samples<-args[2]


temp = list.files(pattern="changes.ALL*")
if(substring(popu, 1, nchar("LBC")) == "LBC")
{
  temp = list.files(pattern="changes.LBC.*")
}
myfiles = lapply(temp, read.table, header=F, sep="\t", stringsAsFactors=FALSE)
dat<-do.call(rbind, myfiles)
#dat<-read.table("/home/jprende3/Scratch/Allan/ParsedCounts/changes.ALL.22.txt", header=F, sep="\t")

colnames(dat)<-c("Chr", "Pos1", "Id1", "Pos2", "Id2", "Dist", "Null", "Cons1", "Cons2", "Anc1", "Anc2", "SameCodon", "StartF", "Mid1F", "Mid2F", "EndF", "StartC", "Mid1C", "Mid2C", "EndC")
dat$Cons1[which((dat$Cons1 == "stop_gained") | (dat$Cons1 == "stop_lost") | (dat$Cons1 == "start_lost") | (dat$Cons1 == "missense_variant") | (dat$Cons1 == "synonymous_variant"))]<-"coding_variant"
dat$Cons2[which((dat$Cons2 == "stop_gained") | (dat$Cons2 == "stop_lost") | (dat$Cons2 == "start_lost") | (dat$Cons2 == "missense_variant") | (dat$Cons2 == "synonymous_variant"))]<-"coding_variant"

#only those in same triplet
dat<-dat[which(dat$SameCodon == 1), ]
#set it so that Mid1 is the more frequent mid change
dat[which(dat$Mid2F > dat$Mid1F), c("Mid1F", "Mid2F", "Mid1C", "Mid2C")] <- dat[which(dat$Mid2F > dat$Mid1F), c("Mid2F", "Mid1F", "Mid2C", "Mid1C")]
#just get MNPs were exactly three allele combinations are observed, one of which is ancestral or just two where where ancestral has a freq of 0
dat$Haps<-rowSums(dat[,c("StartF", "Mid1F", "Mid2F", "EndF")] > 0)
dat<-dat[which(((dat$Haps == 3) & (dat$StartF > 0)) | ((dat$Haps == 2) & (dat$StartF == 0))),]
#only keep sites with confident ancestral alleles
dat<-dat[which((toupper(dat$Anc1) == dat$Anc1) & (toupper(dat$Anc2) == dat$Anc2)),]


#pdf(paste("heatmap_",popu,"_",i,".pdf", sep=""))

#dev.off()


getCounts<-function(datC, pairs, set)
{
  if(set == "first")
  {
    long<-melt(acast(datC, toupper(StartC)~toupper(Mid1C)))
    
  }
  else if(set == "second")
  {
    long<-melt(acast(datC, toupper(Mid1C)~toupper(EndC)))
  }
  long<-merge(pairs, long, all.x=TRUE)
  long[is.na(long)]<-0
  long<-long[with(long, order(Var1, Var2)),]
  long$diffs<-mapply(function(x,y) length(which(x!=y)),strsplit(as.character(long$Var1),""),strsplit(as.character(long$Var2),""))
  long<-long[which(long$diffs == 1),]
  long$diffIndex<-mapply(function(x,y) which(x!=y),strsplit(as.character(long$Var1),""),strsplit(as.character(long$Var2),""))
  long$diffBase<-mapply(function(x,y) y[which(x!=y)],strsplit(as.character(long$Var1),""),strsplit(as.character(long$Var2),""))
  table<-acast(long, Var1~diffIndex+diffBase)
  return(table)
}

getHeatmap<-function(datH)
{
  #types<-names(which(sort(table(datH$Cons1), decreasing=TRUE)> 400))
  #types<-c("intergenic_variant", "intron_variant", "3_prime_UTR_variant", "5_prime_UTR_variant", "coding_variant")
  types<-c("intergenic_variant", "intron_variant", "coding_variant")
  allPairs<-melt(acast(dat, toupper(StartC)~toupper(Mid1C)))
  allPairs$value<-NULL
  
  for(i in 1:length(types))
  {
    print(types[i])
    datSub<-datH[which((datH$Cons1 == types[i]) & (datH$Cons2 == types[i])),]
    
    firstTable<-getCounts(datSub, allPairs, "first")
    secondTable<-getCounts(datSub, allPairs, "second")
    
    ratioTable<-log((secondTable+1)/(firstTable+1), 2)
    
    #pdf(paste("heatmap_",popu,"_",i,".pdf", sep=""))
    if(i == 1)
    {
      ht<-Heatmap(ratioTable, cluster_columns = FALSE, cluster_rows = TRUE, column_title =paste(types[i], dim(datSub)[1], sep=" "))
    }
    else
    {
      ht<-ht+Heatmap(ratioTable, cluster_columns = FALSE, cluster_rows = FALSE, column_title =paste(types[i], dim(datSub)[1], sep=" "))
    }
      #dev.off()
  }
  return(ht)
}

heat<-getHeatmap(dat)
pdf("heatmap.pdf", width=20)
heat=draw(heat)
dev.off()
