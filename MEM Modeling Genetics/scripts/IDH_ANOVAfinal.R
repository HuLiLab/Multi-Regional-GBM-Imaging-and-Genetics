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

#Annotate variants with desired labels 
rawdata$IDH.mutation.status[rawdata$IDH.mutation.status=="IDH WT"]="WT"
rawdata$IDH.mutation.status[rawdata$IDH.mutation.status=="IDH MUT"]="MUT"



#keep samples with genomics and imaging data
rawdata=rawdata[rowSums(rawdata[,5:9], na.rm = TRUE)>0,]
rawdata=rawdata[!is.na(rawdata$IDH.mutation.status),]
rawdata=rawdata[!is.na(rawdata$MRI.contrast.enhancing.annotation),]


aggres=cbind.data.frame(combos[,c(1,2)],NA)
aggpval=cbind.data.frame(combos[,c(1,2)],NA)
aggmeans=data.frame()
aggidh=data.frame()

for(j in 1:length(combos[,1])){
  data=rawdata[,c("MRI.contrast.enhancing.annotation",as.character(combos[j,2]),"IDH.mutation.status","Patient_ID")]
  if(combos[j,1] != "Total"){
    data=data[data[,1]==as.character(combos[j,1]),] 
  }
  
  data=data[,-1]
  colnames(data)=c("Measurement","IDH","Patient")
  
  if(combos[j,3]=="Log10"){
    data[,1]=log10(data[,1])
  }
  else if(combos[j,3]=="Sqrt"){
    data[,1]=sqrt(data[,1])
  }
  else if(combos[j,3]=="Recip"){ 
    data[,1]=1/data[,1]
  }
  data=data[!is.na(data[,1]),]
  model<-lmer(Measurement~IDH+(1|Patient),data=data)
  
  qc=plot(model,ylab = "Residual",xlab = "Variable",main = paste("Transform:",combos[j,3]))
  png(file = paste(combos[j,1],combos[j,2],"QC.png",sep="_"),unit="in",width=3,height=3.5,res=600)
  print(qc) 
  dev.off()
  
  aov.res<-anova(model,type=2)
  aov.res$prop.var.explained<-100*aov.res[,1]/((var(data$Measurement,na.rm = TRUE)*((nrow(data)-1))))
  write.csv(aov.res,paste(combos[j,1],combos[j,2],"anova.csv",sep="_"))
  
  em<-summary(emmeans(model,pairwise~IDH,adjust="tukey"))
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
  
  em.idh=summary(contrast(emmeans(model,"IDH"),"tukey"))

  if(dim(aggidh)[1]==0){
    em.idh[,2]=paste(combos[j,1],combos[j,2],sep="_")
    aggidh=em.idh[,c(1,2,6)]
  }
  else {
    em.idh[,2]=paste(combos[j,1],combos[j,2],sep="_")
    aggidh=rbind.data.frame(aggidh,em.idh[,c(1,2,6)])
  }
  write.csv(aggidh,paste(combos[j,1],combos[j,2],"idh_contrasts.csv",sep="_"))

  if (length(sig_genotypes)>0){
    sig_genotypes=as.character(sig_genotypes)
    print("Sig")
    meanse=as.data.frame(em[["emmeans"]])
    if(grepl("EPI",combos[j,2])){
     title=paste(combos[j,1],"EPI+C")
    }
    else{ title=paste(combos[j,1],gsub(".Raw_Mean","",combos[j,2]))}
    transform=paste("Transform:",combos[j,3])
    plot=ggplot(meanse,aes(x=IDH,fill=IDH))+geom_boxplot(aes(lower=emmean-SE,upper=emmean+SE,middle=emmean,ymin=emmean-3*SE,ymax=emmean+3*SE),stat="identity")+theme_classic()+xlab("")+ylab(title)+
      scale_fill_manual(values = c("#CC5500","darkgrey"))+coord_flip()+theme(legend.position = "none")
  
    
    png(paste(combos[j,1],gsub(".Raw_Mean","",combos[j,2]),"specifics.png",sep="_"),unit="in",width=3.5,height=1.5,res=600)
    print(plot)
    dev.off()
    plot=ggplot(meanse,aes(x=IDH,fill=IDH))+geom_boxplot(aes(lower=emmean-SE,upper=emmean+SE,middle=emmean,ymin=emmean-3*SE,ymax=emmean+3*SE),stat="identity")+theme_classic()+xlab("")+ylab(title)+
      scale_fill_manual(values = c("#CC5500","darkgrey"))+theme(legend.position = "none")
    
    
    png(paste(combos[j,1],gsub(".Raw_Mean","",combos[j,2]),"specificsv2.png",sep="_"),unit="in",width=2.3,height=1.5,res=600)
    print(plot)
    dev.off()
        
  }
}

colnames(aggpval)=c("MRI","Measurment","IDH")
aggres[,2]=gsub(".Raw_Mean","",aggres[,2])
aggres=aggres[order(aggres[,1]),]

write.csv(aggres,"IDH_percentvariance_attributed.csv")
write.csv(aggpval,"IDH_pvalue.csv")

