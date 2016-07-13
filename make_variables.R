##source data for input into qPCR_dominance script

make_variables<-function(qPCR)

qPCR=read.table(file.choose(),header=T,sep="\t") ## choose qPCR data for R.txt
qPCR$Rep=as.factor(qPCR$Rep)

AB_only=subset(qPCR,qPCR$Sym_ID=="AB"&qPCR$D_Count==0)
AD_only=subset(qPCR,qPCR$Sym_ID=="AD"&qPCR$B_Count==0)
BD_only=subset(qPCR,qPCR$Sym_ID=="BD"&qPCR$A_Count==0)

AB_onlyT1=subset(qPCR,qPCR$Sym_ID=="AB"&qPCR$D_Count==0&qPCR$Time==1)
AB_onlyT2=subset(qPCR,qPCR$Sym_ID=="AB"&qPCR$D_Count==0&qPCR$Time==2)

AB_onlyT1_Light=subset(qPCR,qPCR$Sym_ID=="AB"&qPCR$D_Count==0&qPCR$Time==1&Env=="Light")
AB_onlyT1_Dark=subset(qPCR,qPCR$Sym_ID=="AB"&qPCR$D_Count==0&qPCR$Time==1&Env=="Dark")
AB_onlyT2_Light=subset(qPCR,qPCR$Sym_ID=="AB"&qPCR$D_Count==0&qPCR$Time==2&Env=="Light")
AB_onlyT2_Dark=subset(qPCR,qPCR$Sym_ID=="AB"&qPCR$D_Count==0&qPCR$Time==2&Env=="Dark")

AD_onlyT1=subset(qPCR,qPCR$Sym_ID=="AD"&qPCR$B_Count==0&qPCR$Time==1)
AD_onlyT2=subset(qPCR,qPCR$Sym_ID=="AD"&qPCR$B_Count==0&qPCR$Time==2)

AD_onlyT1_Light=subset(qPCR,qPCR$Sym_ID=="AD"&qPCR$B_Count==0&qPCR$Time==1&Env=="Light")
AD_onlyT1_Dark=subset(qPCR,qPCR$Sym_ID=="AD"&qPCR$B_Count==0&qPCR$Time==1&Env=="Dark")
AD_onlyT2_Light=subset(qPCR,qPCR$Sym_ID=="AD"&qPCR$B_Count==0&qPCR$Time==2&Env=="Light")
AD_onlyT2_Dark=subset(qPCR,qPCR$Sym_ID=="AD"&qPCR$B_Count==0&qPCR$Time==2&Env=="Dark")

BD_onlyT1=subset(qPCR,qPCR$Sym_ID=="BD"&qPCR$A_Count==0&qPCR$Time==1)
BD_onlyT2=subset(qPCR,qPCR$Sym_ID=="BD"&qPCR$A_Count==0&qPCR$Time==2)

BD_onlyT1_Light=subset(qPCR,qPCR$Sym_ID=="BD"&qPCR$A_Count==0&qPCR$Time==1&Env=="Light")
BD_onlyT1_Dark=subset(qPCR,qPCR$Sym_ID=="BD"&qPCR$A_Count==0&qPCR$Time==1&Env=="Dark")
BD_onlyT2_Light=subset(qPCR,qPCR$Sym_ID=="BD"&qPCR$A_Count==0&qPCR$Time==2&Env=="Light")
BD_onlyT2_Dark=subset(qPCR,qPCR$Sym_ID=="BD"&qPCR$A_Count==0&qPCR$Time==2&Env=="Dark")

ctA_AB_T1=AB_onlyT1$A_Count
ctA_AB_T1[ctA_AB_T1==0]<-1
ctB_AB_T1=AB_onlyT1$B_Count
ctB_AB_T1[ctB_AB_T1==0]<-1
ctA_AB_T2=AB_onlyT2$A_Count
ctA_AB_T2[ctA_AB_T2==0]<-1
ctB_AB_T2=AB_onlyT2$B_Count
ctB_AB_T2[ctB_AB_T2==0]<-1

ctA_AD_T1=AD_onlyT1$A_Count
ctA_AD_T1[ctA_AD_T1==0]<-1
ctD_AD_T1=AD_onlyT1$D_Count
ctD_AD_T1[ctD_AD_T1==0]<-1
ctA_AD_T2=AD_onlyT2$A_Count
ctA_AD_T2[ctA_AD_T2==0]<-1
ctD_AD_T2=AD_onlyT2$D_Count
ctD_AD_T2[ctD_AD_T2==0]<-1

ctB_BD_T1=BD_onlyT1$B_Count
ctB_BD_T1[ctB_BD_T1==0]<-1
ctD_BD_T1=BD_onlyT1$D_Count
ctD_BD_T1[ctD_BD_T1==0]<-1
ctB_BD_T2=BD_onlyT2$B_Count
ctB_BD_T2[ctB_BD_T2==0]<-1
ctD_BD_T2=BD_onlyT2$D_Count
ctD_BD_T2[ctD_BD_T2==0]<-1

ct1ALL_T1_Light=c(AB_onlyT1_Light$A_Count,AD_onlyT1_Light$A_Count,BD_onlyT1_Light$B_Count)

cat("AB_T1","\n")
##qPCR_dominance(ctA_AB_T1,ctB_AB_T1,.95,1000)
qPCR_percentdom(ctA_AB_T1,ctB_AB_T1,.95,1000)
cat("\n")

cat("AD_T1","\n")
##qPCR_dominance(ctA_AD_T1,ctD_AD_T1,.95,1000)
qPCR_percentdom(ctA_AD_T1,ctD_AD_T1,.95,1000)
cat("\n")

cat("BD_T1","\n")
##qPCR_dominance(ctB_BD_T1,ctD_BD_T1,.95,1000)
qPCR_percentdom(ctB_BD_T1,ctD_BD_T1,.95,1000)
cat("\n")

cat("AB_T2","\n")
##qPCR_dominance(ctA_AB_T2,ctB_AB_T2,.95,1000)
qPCR_percentdom(ctA_AB_T2,ctB_AB_T2,.95,1000)
cat("\n")

cat("AD_T2","\n")
##qPCR_dominance(ctA_AD_T2,ctD_AD_T2,.95,1000)
qPCR_percentdom(ctA_AD_T2,ctD_AD_T2,.95,1000)
cat("\n")

cat("BD_T2","\n")
##qPCR_dominance(ctB_BD_T2,ctD_BD_T2,.95,1000)
qPCR_percentdom(ctB_BD_T2,ctD_BD_T2,.95,1000)
cat("\n")
