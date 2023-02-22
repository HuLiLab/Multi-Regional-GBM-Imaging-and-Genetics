library(MuMIn)
library(lme4)
library(lmerTest)
library(reshape2)

rawdata=read.csv("MultiReg_imaging_pathclass.csv")

ceraw=rawdata[rawdata$MRI=="CE",]
neraw=rawdata[rawdata$MRI=="NE",]


rawdata[,8:14]=scale(rawdata[,8:14])
rawdata[,3:6]=scale(rawdata[,3:6])
ceraw[,8:14]=scale(ceraw[,8:14])
ceraw[,3:6]=scale(ceraw[,3:6])
neraw[,8:14]=scale(neraw[,8:14])
neraw[,3:6]=scale(neraw[,3:6])

for(i in c("GPM","MTC","NEU","PPR")){
  result=as.data.frame(matrix(ncol = 15,nrow=7))
  rownames(result)=colnames(rawdata)[8:14]
  colnames(result)=c("All Corr.","All Pval","All MEM R2","All MEM Est.","All MEM SE","CE Corr.","CE Pval","CE MEM R2","CE MEM Est.","CE MEM SE","NE Corr.","NE Pval","NE MEM R2","NE MEM Est.","NE MEM SE")
  
  for (j in colnames(rawdata)[8:14]){
    curr=cor.test(rawdata[,i],rawdata[,j])
    
    result[j,"All Corr."]=curr[["estimate"]]
    result[j,"All Pval"]=curr[["p.value"]]
    
    cecurr=cor.test(ceraw[,i],ceraw[,j])
    result[j,"CE Corr."]=cecurr[["estimate"]]
    result[j,"CE Pval"]=cecurr[["p.value"]]
    
    necurr=cor.test(neraw[,i],neraw[,j])
    result[j,"NE Corr."]=necurr[["estimate"]]
    result[j,"NE Pval"]=necurr[["p.value"]]
    
    model = lmer(paste0(j,' ~ ',i,' + (1 | PatientID.x)'), data=rawdata)
    result[j,"All MEM R2"]=r.squaredGLMM(model)[,'R2c']
    summ = summary(model)
    result[j,"All MEM Est."]=summ[["coefficients"]][2,1]
    result[j,"All MEM SE"]=summ[["coefficients"]][2,2]
    
    
    model = lmer(paste0(j,' ~ ',i,' + (1 | PatientID.x)'), data=ceraw)
    result[j,"CE MEM R2"]=r.squaredGLMM(model)[,'R2c']
    summ = summary(model)
    result[j,"CE MEM Est."]=summ[["coefficients"]][2,1]
    result[j,"CE MEM SE"]=summ[["coefficients"]][2,2]
    
    model = lmer(paste0(j,' ~ ',i,' + (1 | PatientID.x)'), data=neraw)
    result[j,"NE MEM R2"]=r.squaredGLMM(model)[,'R2c']
    summ = summary(model)
    result[j,"NE MEM Est."]=summ[["coefficients"]][2,1]
    result[j,"NE MEM SE"]=summ[["coefficients"]][2,2]
  } 

  write.csv(result,file=paste(i,"statistics.csv",sep="_"))
}

rawdata=read.csv("MultiReg_advimaging_pathclass.csv")

ceraw=rawdata[rawdata$MRI=="CE",]
neraw=rawdata[rawdata$MRI=="NE",]

neraw[,8:15]=scale(neraw[,8:15],)
neraw[,3:6]=scale(neraw[,3:6])
ceraw[,8:15]=scale(ceraw[,8:15],)
ceraw[,3:6]=scale(ceraw[,3:6])
rawdata[,8:15]=scale(rawdata[,8:15],)
rawdata[,3:6]=scale(rawdata[,3:6])


for(i in c("GPM","MTC","NEU","PPR")){
  result=as.data.frame(matrix(ncol = 15,nrow=8))
  rownames(result)=colnames(rawdata)[8:15]
  colnames(result)=c("All Corr.","All Pval","All MEM R2","All MEM Est.","All MEM SE","CE Corr.","CE Pval","CE MEM R2","CE MEM Est.","CE MEM SE","NE Corr.","NE Pval","NE MEM R2","NE MEM Est.","NE MEM SE")
  for (j in colnames(rawdata)[8:15]){
    currdata=rawdata[!is.na(rawdata[,j]),]
    currce=ceraw[!is.na(ceraw[,j]),]
    currne=neraw[!is.na(neraw[,j]),]
    
    curr=cor.test(currdata[,i],currdata[,j])
    
    result[j,"All Corr."]=curr[["estimate"]]
    result[j,"All Pval"]=curr[["p.value"]]
    
    cecurr=cor.test(currce[,i],currce[,j])
    result[j,"CE Corr."]=cecurr[["estimate"]]
    result[j,"CE Pval"]=cecurr[["p.value"]]
    
    necurr=cor.test(currne[,i],currne[,j])
    result[j,"NE Corr."]=necurr[["estimate"]]
    result[j,"NE Pval"]=necurr[["p.value"]]
    print(paste("All",i,j))
    model = lmer(paste0(j,' ~ ',i,' + (1 | PatientID.x)'), data=currdata)
    result[j,"All MEM R2"]=r.squaredGLMM(model)[,'R2c']
    summ = summary(model)
    result[j,"All MEM Est."]=summ[["coefficients"]][2,1]
    result[j,"All MEM SE"]=summ[["coefficients"]][2,2]
    
    print(paste("CE",i,j))
    
    model = lmer(paste0(j,' ~ ',i,' + (1 | PatientID.x)'), data=currce)
    result[j,"CE MEM R2"]=r.squaredGLMM(model)[,'R2c']
    summ = summary(model)
    result[j,"CE MEM Est."]=summ[["coefficients"]][2,1]
    result[j,"CE MEM SE"]=summ[["coefficients"]][2,2]
    print(paste("NE",i,j))
    
    model = lmer(paste0(j,' ~ ',i,' + (1 | PatientID.x)'), data=currne)
    result[j,"NE MEM R2"]=r.squaredGLMM(model)[,'R2c']
    summ = summary(model)
    result[j,"NE MEM Est."]=summ[["coefficients"]][2,1]
    result[j,"NE MEM SE"]=summ[["coefficients"]][2,2]
  } 
  write.csv(result,file=paste(i,"adv_imaging_statistics.csv",sep="_"))
}
