library(lme4)
library(lmerTest)
library(emmeans)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(plyr)
library(matrixStats)
library(RColorBrewer)
rawdata=read.csv("MultiReg_imaging_genetics.csv")
 
rawdata=rawdata[rawdata$IDH.mutation.status=="IDH WT",]

#keep samples with genomics and imaging data
rawdata=rawdata[rowSums(rawdata[,5:9], na.rm = TRUE)>0,]
rawdata=rawdata[!is.na(rawdata$MRI.contrast.enhancing.annotation),]

rawdata[,11]=factor(paste("CNV",rawdata[,11]))
rawdata[,12]=factor(paste("Mut",rawdata[,12]))
rawdata[,13]=factor(paste("CNV",rawdata[,13]))
rawdata[,14]=factor(paste("CNV",rawdata[,14]))
rawdata[,15]=factor(paste("Mut",rawdata[,15]))
rawdata[,16]=factor(paste("CNV",rawdata[,16]))
rawdata[,17]=factor(paste("Mut",rawdata[,17]))
rawdata[,18]=factor(paste("CNV",rawdata[,18]))
rawdata[,19]=factor(paste("Mut",rawdata[,19]))

combos=read.csv("combinations.csv",header=FALSE)



aggres=cbind.data.frame(combos[,c(1,2)],NA,NA,NA,NA,NA,NA,NA,NA,NA)
aggpval=cbind.data.frame(combos[,c(1,2)],NA,NA,NA,NA,NA,NA,NA,NA,NA)
for(j in 1:length(combos[,1])){
  
  data=rawdata[,c(4,grep(combos[j,2],colnames(rawdata)),11:19,2)]
  if(combos[j,1] != "Total"){
    data=data[data[,1]==as.character(combos[j,1]),] 
  }
  data=data[,-1]
  colnames(data)=c("Measurement","EGFR_cnv","EGFR_mut","CDKN2A_cnv","NF1_cnv","NF1_mut","TP53_cnv","TP53_mut","PTEN_cnv","PTEN_mut","Patient")
  if(combos[j,3]=="Log10"){
    data[,1]=log10(data[,1])
  }
  else if(combos[j,3]=="Sqrt"){
    data[,1]=sqrt(data[,1])
  }
  else if(combos[j,3]=="Recip"){
    data[,1]=1/data[,1]
  }
 
  model<-lmer(Measurement~EGFR_cnv+EGFR_mut+CDKN2A_cnv+NF1_cnv+NF1_mut+TP53_cnv+TP53_mut+PTEN_cnv+PTEN_mut+(1|Patient),data=data)
  
  qc=plot(model,ylab = "Residual",xlab = "Variable",main = paste("Transform:",combos[j,3]))
  png(file = paste(combos[j,1],combos[j,2],"QC.png",sep="_"))
  print(qc) 
  dev.off()

  aov.res<-anova(model,type=2)
  aov.res$prop.var.explained<-100*aov.res[,1]/((var(data$Measurement,na.rm = TRUE)*((nrow(data)-1))))

  aggres[aggres[,1]==combos[j,1] & aggres[,2]==combos[j,2],3:11]=aov.res[,7]
  aggpval[aggpval[,1]==combos[j,1] & aggpval[,2]==combos[j,2],3:11]=aov.res[,6]
  
}
colnames(aggpval)=c("MRI","Measurment","EGFR_cnv","EGFR_mut","CDKN2A_cnv","NF1_cnv","NF1_mut","TP53_cnv","TP53_mut","PTEN_cnv","PTEN_mut")
colnames(aggres)=c("MRI","Measurment","EGFR_cnv","EGFR_mut","CDKN2A_cnv","NF1_cnv","NF1_mut","TP53_cnv","TP53_mut","PTEN_cnv","PTEN_mut")

write.csv(aggres,"single_factor_percentvariance_attributed.csv")
write.csv(aggpval,"single_factor_pvalue.csv")


