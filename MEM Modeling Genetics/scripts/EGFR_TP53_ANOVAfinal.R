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
rawdata$EGFR_cnv[rawdata$EGFR_cnv==0 ]="WT"
rawdata$EGFR_cnv[rawdata$EGFR_cnv==2 ]="Amp"
rawdata$EGFR_cnv[rawdata$EGFR_cnv==1]="Gain"
rawdata=rawdata[!is.na(rawdata$EGFR_cnv),]
rawdata$EGFR_cnv=factor(rawdata$EGFR_cnv,levels=c("WT","Gain","Amp"))

rawdata$TP53_cnv[rawdata$TP53_cnv==0 ]="WT"
rawdata$TP53_cnv[rawdata$TP53_cnv==-1 ]="Loss-1"
rawdata$TP53_cnv[rawdata$TP53_cnv==-2 ]="Loss-2"
rawdata$TP53_cnv=factor(rawdata$TP53_cnv,levels=c("WT","Loss-1","Loss-2"))

#keep samples with genomics and imaging data
rawdata=rawdata[rowSums(rawdata[,5:9], na.rm = TRUE)>0,]
rawdata=rawdata[!is.na(rawdata$EGFR_cnv),]
rawdata=rawdata[!is.na(rawdata$TP53_cnv),]
rawdata=rawdata[!is.na(rawdata$MRI.contrast.enhancing.annotation),]


combos=read.csv("combinations.csv",header=FALSE)

aggres=cbind.data.frame(combos[,c(1,2)],NA,NA,NA)
aggpval=cbind.data.frame(combos[,c(1,2)],NA,NA,NA)
aggmeans=data.frame()
aggtp53=data.frame()
aggegfr=data.frame()
for(j in 1:length(combos[,1])){
  data=rawdata[,c("MRI.contrast.enhancing.annotation",as.character(combos[j,2]),"EGFR_cnv","TP53_cnv","Patient_ID")]
  if(combos[j,1] != "Total"){
    data=data[data[,1]==as.character(combos[j,1]),] 
  }
  data=data[,-1]
  colnames(data)=c("Measurement","EGFR","TP53","Patient")
  
  if(combos[j,3]=="Log10"){
    data[,1]=log10(data[,1])
  }
  else if(combos[j,3]=="Sqrt"){
    data[,1]=sqrt(data[,1])
  }
  else if(combos[j,3]=="Recip"){ 
    data[,1]=1/data[,1]
  }

  model<-lmer(Measurement~EGFR*TP53+(1|Patient),data=data)
  
  qc=plot(model,ylab = "Residual",xlab = "Variable",main = paste("Transform:",combos[j,3]))
  png(file = paste(combos[j,1],combos[j,2],"EGFRvTP53_QC.png",sep="_"))
  print(qc) 
  dev.off()
  
  aov.res<-anova(model,type=2)
  aov.res$prop.var.explained<-100*aov.res[,1]/((var(data$Measurement,na.rm = TRUE)*((nrow(data)-1))))
  write.csv(aov.res,paste(combos[j,1],combos[j,2],"anova_EGFRvTP53.csv",sep="_"))
  
  em<-summary(emmeans(model,pairwise~EGFR:TP53,adjust="tukey"))
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
  
  aggres[aggres[,1]==combos[j,1] & aggres[,2]==combos[j,2],3:5]=aov.res[,7]
  aggpval[aggpval[,1]==aggpval[j,1] & aggpval[,2]==aggpval[j,2],3:5]=aov.res[,6]
  
  pmat=pmat[rowMeans(pmat,na.rm=TRUE) !=1,colMeans(pmat, na.rm=TRUE) !=1]
  
  col_fun=colorRamp2(seq(from=0.05,to=1,length.out=10),rev(viridis(10)))
  colnames(pmat)=gsub("\\(","",colnames(pmat))
  colnames(pmat)=gsub("\\)","",colnames(pmat))
  rownames(pmat)=gsub("\\(","",rownames(pmat))
  rownames(pmat)=gsub("\\)","",rownames(pmat))
  
  ht=Heatmap(pmat,col=col_fun)
  png(file = paste(combos[j,1],combos[j,2],"pval_ht.png",sep="_"))
  print(ht) 
  dev.off()
  if(min(pmat)<.8){
    pmat=pmat[rowMins(pmat,na.rm=TRUE) <.8,rowMins(pmat, na.rm=TRUE) <.8]
    ht=Heatmap(pmat,col=col_fun)
    png(file = paste(combos[j,1],combos[j,2],"pval_ht_subset.png",sep="_"))
    print(ht) 
    dev.off()
  }
  print(min(pmat,na.rm=TRUE))
  
  sig_genotypes=melt(pmat) 
  sig_genotypes=sig_genotypes[sig_genotypes[,3]<.05,]
  sig_genotypes=unique(c(sig_genotypes[,1],sig_genotypes[,2]))
  
  em.egfr=summary(contrast(emmeans(model,"EGFR"),"tukey"))
  em.tp53=summary(contrast(emmeans(model,"TP53"),"tukey"))
  
  if(dim(aggegfr)[1]==0){
    em.egfr[,2]=paste(combos[j,1],combos[j,2],sep="_")
    aggegfr=em.egfr[,c(1,2,6)]
  }
  else {
    em.egfr[,2]=paste(combos[j,1],combos[j,2],sep="_")
    aggegfr=rbind.data.frame(aggegfr,em.egfr[,c(1,2,6)])
  }
  write.csv(aggegfr,paste(combos[j,1],combos[j,2],"egfr_contrasts.csv",sep="_"))
  
  if(dim(aggtp53)[1]==0){
    em.tp53[,2]=paste(combos[j,1],combos[j,2],sep="_")
    aggtp53=em.tp53[,c(1,2,6)]
  }
  else {
    em.tp53[,2]=paste(combos[j,1],combos[j,2],sep="_")
    aggtp53=rbind.data.frame(aggtp53,em.tp53[,c(1,2,6)])
  }
  write.csv(aggtp53,paste(combos[j,1],combos[j,2],"tp53_contrasts.csv",sep="_"))
  
  
  if (length(sig_genotypes)>0){
    sig_genotypes=t(data.frame(strsplit(as.character(sig_genotypes)," ")))
    print("Sig")
    for(i in 1:length(sig_genotypes[,1])){
      curr=cbind.data.frame(data,0)
      curr[curr[,2]==sig_genotypes[i,1] & curr[,3]==sig_genotypes[i,2],5]="In"
      curr[!(curr[,2]==sig_genotypes[i,1] & curr[,3]==sig_genotypes[i,2]),5]="Out"
      colnames(curr)[5]="Group"
      
      curr[,5]=factor(curr[,5],levels=c("In","Out"))
      contrasts(curr$Group)=c(1,-1)
      model<-lmer(Measurement~Group+EGFR+TP53+(1|Patient),data=curr)
      #model<-lmer(Measurement~Group+(1|Patient),data=curr)
      
      em.group.pval=summary(contrast(emmeans(model,"Group")))
      em.group.sum=summary(emmeans(model,"Group"))
      em.group.pval=cbind.data.frame(em.group.sum[,c(1:3)],em.group.pval[,6])
      em.group.pval[,1]=as.character(em.group.pval[,1])
      em.group.pval[em.group.pval$Group=="In",1]=paste(sig_genotypes[i,1],sig_genotypes[i,2])
      em.group.pval[em.group.pval$Group=="Out",1]="Other"
      title=paste(combos[j,1],gsub(".Raw_Mean","",combos[j,2]),"\n","Transform:",combos[j,3],"\n","P-value:",round(em.group.pval[1,4],4))
      plot=ggplot(em.group.pval,aes(x=Group,y=emmean,fill=Group))+geom_bar(stat="identity")+geom_errorbar(aes(ymin=emmean-SE,ymax=emmean+SE),size=1)+scale_fill_manual(values=c("steelblue1","grey"))+theme_bw()+
        xlab("")+ylab("Value")+ggtitle(title)+theme(legend.position = "none")
      png(paste(combos[j,1],gsub(".Raw_Mean","",combos[j,2]),paste(sig_genotypes[i,1],sig_genotypes[i,2],sep="."),"specifics.png",sep="_"))
      print(plot)
      dev.off()
    }
  }
}

colnames(aggpval)=c("MRI","Measurment","EGFR","TP53","EGFR-TP53")
aggres[,2]=gsub(".Raw_Mean","",aggres[,2])
aggres=aggres[order(aggres[,1]),]

write.csv(aggres,"EGFR_TP53_percentvariance_attributed.csv")
write.csv(aggpval,"EGFR_TP53_pvalue.csv")

