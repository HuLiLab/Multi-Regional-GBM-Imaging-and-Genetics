library(lme4)
library(lmerTest)
library(emmeans)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(plyr)
library(matrixStats)
library(RColorBrewer)
library(reshape)
library(ggplot2)

rawdata=read.csv("MultiReg_imaging_genetics.csv")
rawdata=rawdata[rowSums(rawdata[,5:9], na.rm = TRUE)>0,]
rawdata=rawdata[!is.na(rawdata$MRI.contrast.enhancing.annotation),]
garfdata=read.csv("Garafanoclass.csv")
patients=intersect(rawdata[,1],garfdata[,1])
garfdata=garfdata[garfdata[,1]%in%patients,]
rawdata=rawdata[rawdata[,1]%in%patients,]
garfdata=garfdata[order(garfdata[,1]),]
rawdata=rawdata[order(rawdata[,1]),]
rawdata=cbind.data.frame(rawdata,garfdata[,2])
colnames(rawdata)[20]="GarofanoClass"
combos=read.csv("combinations.csv",header=FALSE)

aggres=cbind.data.frame(combos[,c(1,2)],NA)
aggpval=cbind.data.frame(combos[,c(1,2)],NA)
aggmeans=data.frame()
agggarf=data.frame()

for(j in 1:length(combos[,1])){
  data=rawdata[,c("MRI.contrast.enhancing.annotation",as.character(combos[j,2]),"GarofanoClass","Patient_ID")]
  if(combos[j,1] != "Total"){
    data=data[data[,1]==as.character(combos[j,1]),] 
  }
  
  data=data[,-1]
  colnames(data)=c("Measurement","GARF","Patient")
  
  if(combos[j,3]=="Log10"){
    data[,1]=log10(data[,1])
  }
  else if(combos[j,3]=="Sqrt"){
    data[,1]=sqrt(data[,1])
  }
  else if(combos[j,3]=="Recip"){ 
    data[,1]=1/data[,1]
  }
  model<-lmer(Measurement~GARF+(1|Patient),data=data)
  
  qc=plot(model,ylab = "Residual",xlab = "Variable",main = paste("Transform:",combos[j,3]))
  png(file = paste(combos[j,1],combos[j,2],"QC.png",sep="_"),unit="in",width=3,height=3.5,res=600)
  print(qc) 
  dev.off()
  
  aov.res<-anova(model,type=2)
  aov.res$prop.var.explained<-100*aov.res[,1]/((var(data$Measurement,na.rm = TRUE)*((nrow(data)-1))))
  write.csv(aov.res,paste(combos[j,1],combos[j,2],"anova.csv",sep="_"))
  
  em<-summary(emmeans(model,pairwise~GARF,adjust="tukey"))
  if(dim(aggmeans)[1]==0){
    aggmeans=em[["emmeans"]][,1:3]
    colnames(aggmeans)[3]=paste(combos[j,1],combos[j,2],sep="_")
  }
  else {
    aggmeans=cbind.data.frame(aggmeans,em[["emmeans"]][,3])
    colnames(aggmeans)[2+j]=paste(combos[j,1],combos[j,2],sep="_")
  }
  
  split<-strsplit(em$contrasts$contrast," - ")
  em$contrasts$int1<-sapply(split,function(x){x[1]})
  em$contrasts$int2<-sapply(split,function(x){x[2]})
  dims<-unique(unlist(em$contrasts[,7:8]))
  pmat<-matrix(rep(1,length(dims)^2),nrow=length(dims),
               dimnames = list(dims,dims))
  for(i in 1:nrow(em$contrasts)){
    a<-em$contrasts$int1[i]
    b<-em$contrasts$int2[i]
    p<-em$contrasts$p.value[i]
    pmat[a,b]<-pmat[b,a]<-p
  }
  write.csv(pmat,paste(combos[j,1],combos[j,2],"pvalues.csv",sep="_"))
  
  aggres[aggres[,1]==combos[j,1] & aggres[,2]==combos[j,2],3]=aov.res[,7]
  aggpval[aggpval[,1]==aggpval[j,1] & aggpval[,2]==aggpval[j,2],3]=aov.res[,6]
  
  pmat=pmat[rowMeans(pmat,na.rm=TRUE) !=1,colMeans(pmat, na.rm=TRUE) !=1]
  
  col_fun=colorRamp2(seq(from=0.05,to=1,length.out=10),rev(viridis(10)))
  colnames(pmat)=gsub("\\(","",colnames(pmat))
  colnames(pmat)=gsub("\\)","",colnames(pmat))
  rownames(pmat)=gsub("\\(","",rownames(pmat))
  rownames(pmat)=gsub("\\)","",rownames(pmat))
  
  ht=Heatmap(pmat,col=col_fun,heatmap_legend_param = list(title="P-value"))
  png(file = paste(combos[j,1],combos[j,2],"pval_ht.png",sep="_"),unit="in",width=3.5,height=3,res=600)
  print(ht) 
  dev.off()

  print(min(pmat,na.rm=TRUE))
  
  sig_genotypes=melt(pmat) 
  sig_genotypes=sig_genotypes[sig_genotypes[,3]<.1,]
  sig_genotypes=unique(c(sig_genotypes[,1],sig_genotypes[,2]))
  
  em.garf=summary(contrast(emmeans(model,"GARF"),"tukey"))

  if(dim(agggarf)[1]==0){
    em.garf[,2]=paste(combos[j,1],combos[j,2],sep="_")
    agggarf=em.garf[,c(1,2,6)]
  }
  else {
    em.garf[,2]=paste(combos[j,1],combos[j,2],sep="_")
    agggarf=rbind.data.frame(agggarf,em.garf[,c(1,2,6)])
  }
  write.csv(agggarf,paste(combos[j,1],combos[j,2],"garf_contrasts.csv",sep="_"))

  if (length(sig_genotypes)>0){
    sig_genotypes=as.character(sig_genotypes)
    print("Sig")
    meanse=as.data.frame(em[["emmeans"]])

    title=paste(combos[j,1],gsub(".Raw_Mean","",combos[j,2]))
    transform=paste("Transform:",combos[j,3])
    plot=ggplot(meanse,aes(x=GARF,fill=GARF))+geom_boxplot(aes(lower=emmean-SE,upper=emmean+SE,middle=emmean,ymin=emmean-3*SE,ymax=emmean+3*SE),stat="identity")+theme_bw()+xlab("")+ylab(title)+
      coord_flip()+theme(legend.position = "none")
  
    
    png(paste(combos[j,1],gsub(".Raw_Mean","",combos[j,2]),"specifics.png",sep="_"),unit="in",width=3.5,height=1.5,res=600)
    print(plot)
    dev.off()
        
       
  }
}

colnames(aggpval)=c("MRI","Measurment","GARF")
aggres[,2]=gsub(".Raw_Mean","",aggres[,2])
aggres=aggres[order(aggres[,1]),]

write.csv(aggres,"GARF_percentvariance_attributed.csv")
write.csv(aggpval,"GARF_pvalue.csv")

