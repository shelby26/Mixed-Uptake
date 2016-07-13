
##Make variables of cell counts for plact cct1 and cct2 in qPCR dominance scripts. Data includes all Symbiont treatments parced out by Time and Light Environment, then run in both dominance scripts. Results indicate a tendancy for domination but does not specify dominant strain.

##Separate Data for T1 Light
AB_onlyT1_Light=subset(qPCR,qPCR$Sym_ID=="AB"&qPCR$D_Count==0&qPCR$Time==1&Env=="Light")
AD_onlyT1_Light=subset(qPCR,qPCR$Sym_ID=="AD"&qPCR$B_Count==0&qPCR$Time==1&Env=="Light")
BD_onlyT1_Light=subset(qPCR,qPCR$Sym_ID=="BD"&qPCR$A_Count==0&qPCR$Time==1&Env=="Light")
ct1ALL_T1_Light=c(AB_onlyT1_Light$A_Count,AD_onlyT1_Light$A_Count,BD_onlyT1_Light$B_Count)
ct1ALL_T1_Light[ct1ALL_T1_Light==0]<-1
ct2ALL_T1_Light=c(AB_onlyT1_Light$B_Count,AD_onlyT1_Light$D_Count,BD_onlyT1_Light$D_Count)
ct2ALL_T1_Light[ct2ALL_T1_Light==0]<-1

##Separate Data for T1 Dark
AB_onlyT1_Dark=subset(qPCR,qPCR$Sym_ID=="AB"&qPCR$D_Count==0&qPCR$Time==1&Env=="Dark")
AD_onlyT1_Dark=subset(qPCR,qPCR$Sym_ID=="AD"&qPCR$B_Count==0&qPCR$Time==1&Env=="Dark")
BD_onlyT1_Dark=subset(qPCR,qPCR$Sym_ID=="BD"&qPCR$A_Count==0&qPCR$Time==1&Env=="Dark")
ct1ALL_T1_Dark=c(AB_onlyT1_Dark$A_Count,AD_onlyT1_Dark$A_Count,BD_onlyT1_Dark$B_Count)
ct1ALL_T1_Dark[ct1ALL_T1_Dark==0]<-1
ct2ALL_T1_Dark=c(AB_onlyT1_Dark$B_Count,AD_onlyT1_Dark$D_Count,BD_onlyT1_Dark$D_Count)
ct2ALL_T1_Dark[ct2ALL_T1_Dark==0]<-1

##Separate Data for T2 Light
AB_onlyT2_Light=subset(qPCR,qPCR$Sym_ID=="AB"&qPCR$D_Count==0&qPCR$Time==2&Env=="Light")
AD_onlyT2_Light=subset(qPCR,qPCR$Sym_ID=="AD"&qPCR$B_Count==0&qPCR$Time==2&Env=="Light")
BD_onlyT2_Light=subset(qPCR,qPCR$Sym_ID=="BD"&qPCR$A_Count==0&qPCR$Time==2&Env=="Light")
ct1ALL_T2_Light=c(AB_onlyT2_Light$A_Count,AD_onlyT2_Light$A_Count,BD_onlyT2_Light$B_Count)
ct1ALL_T2_Light[ct1ALL_T2_Light==0]<-1
ct2ALL_T2_Light=c(AB_onlyT2_Light$B_Count,AD_onlyT2_Light$D_Count,BD_onlyT2_Light$D_Count)
ct2ALL_T2_Light[ct2ALL_T2_Light==0]<-1

##Separate Data for T2 Dark
AB_onlyT2_Dark=subset(qPCR,qPCR$Sym_ID=="AB"&qPCR$D_Count==0&qPCR$Time==2&Env=="Dark")
AD_onlyT2_Dark=subset(qPCR,qPCR$Sym_ID=="AD"&qPCR$B_Count==0&qPCR$Time==2&Env=="Dark")
BD_onlyT2_Dark=subset(qPCR,qPCR$Sym_ID=="BD"&qPCR$A_Count==0&qPCR$Time==2&Env=="Dark")
ct1ALL_T2_Dark=c(AB_onlyT2_Dark$A_Count,AD_onlyT2_Dark$A_Count,BD_onlyT2_Dark$B_Count)
ct1ALL_T2_Dark[ct1ALL_T2_Dark==0]<-1
ct2ALL_T2_Dark=c(AB_onlyT2_Dark$B_Count,AD_onlyT2_Dark$D_Count,BD_onlyT2_Dark$D_Count)
ct2ALL_T2_Dark[ct2ALL_T2_Dark==0]<-1

##Run each subset through both qPCR dominance scripts (log based and percent based)
cat("***ALL T1 Light***","\n")
qPCR_percentdom(ct1ALL_T1_Light,ct2ALL_T1_Light,.95,100)
##qPCR_dominance(ct1ALL_T1_Light,ct2ALL_T1_Light,.95,100)
cat("AB")
qPCR_percentdom(AB_onlyT1_Light$A_Count,AB_onlyT1_Light$B_Count,.95,100)
cat("AD")
qPCR_percentdom(AD_onlyT1_Light$A_Count,AD_onlyT1_Light$D_Count,.95,100)
cat("BD")
qPCR_percentdom(BD_onlyT1_Light$B_Count,BD_onlyT1_Light$D_Count,.95,100)

cat("\n","***ALL T1 Dark***","\n")
qPCR_percentdom(ct1ALL_T1_Dark,ct2ALL_T1_Dark,.95,100)
##qPCR_dominance(ct1ALL_T1_Dark,ct2ALL_T1_Dark,.95,100)
cat("AB")
qPCR_percentdom(AB_onlyT1_Dark$A_Count,AB_onlyT1_Dark$B_Count,.95,100)
cat("AD")
qPCR_percentdom(AD_onlyT1_Dark$A_Count,AD_onlyT1_Dark$D_Count,.95,100)
cat("BD")
qPCR_percentdom(BD_onlyT1_Dark$B_Count,BD_onlyT1_Dark$D_Count,.95,100)

cat("\n","***ALL T2 Light***","\n")
qPCR_percentdom(ct1ALL_T2_Light,ct2ALL_T2_Light,.95,100)
##qPCR_dominance(ct1ALL_T2_Light,ct2ALL_T2_Light,.95,100)
cat("AB")
qPCR_percentdom(AB_onlyT2_Light$A_Count,AB_onlyT2_Light$B_Count,.95,100)
cat("AD")
qPCR_percentdom(AD_onlyT2_Light$A_Count,AD_onlyT2_Light$D_Count,.95,100)
cat("BD")
qPCR_percentdom(BD_onlyT2_Light$B_Count,BD_onlyT2_Light$D_Count,.95,100)

cat("\n","***ALL T2 Dark***","\n")
qPCR_percentdom(ct1ALL_T2_Dark,ct2ALL_T2_Dark,.95,100)
##qPCR_dominance(ct1ALL_T2_Dark,ct2ALL_T2_Dark,.95,100)
cat("AB")
qPCR_percentdom(AB_onlyT2_Dark$A_Count,AB_onlyT2_Dark$B_Count,.95,100)
cat("AD")
qPCR_percentdom(AD_onlyT2_Dark$A_Count,AD_onlyT2_Dark$D_Count,.95,100)
cat("BD")
qPCR_percentdom(BD_onlyT2_Dark$B_Count,BD_onlyT2_Dark$D_Count,.95,100)



##Make sweet, sweet figure to visualize percentage dominance
quartz()
par(mfrow=c(2,3))
hist(AB_onlyT1_Dark$Per_A,col="red",main="",xlab="Percent A in AB-Dark-T1")
hist(AD_onlyT1_Dark$Per_A,col="red",main="",xlab="Percent A in AD-Dark")
hist(BD_onlyT1_Dark$Per_B,col="blue",main="",xlab="Percent B in BD-Dark")
hist(AB_onlyT1_Dark$Per_B,col="blue",main="",xlab="Percent B in AB-Dark")
hist(AD_onlyT1_Dark$Per_D,col="green",main="",xlab="Percent D in AD-Dark")
hist(BD_onlyT1_Dark$Per_D,col="green",main="",xlab="Percent D in BD-Dark")

##Make sweet, sweet figure to visualize percentage dominance
quartz()
par(mfrow=c(2,3))
hist(AB_onlyT1_Light$Per_A,col="red",main="",xlab="Percent A in AB-Light-T1")
hist(AD_onlyT1_Light$Per_A,col="red",main="",xlab="Percent A in AD-Light")
hist(BD_onlyT1_Light$Per_B,col="blue",main="",xlab="Percent B in BD-Light")
hist(AB_onlyT1_Light$Per_B,col="blue",main="",xlab="Percent B in AB-Light")
hist(AD_onlyT1_Light$Per_D,col="green",main="",xlab="Percent D in AD-Light")
hist(BD_onlyT1_Light$Per_D,col="green",main="",xlab="Percent D in BD-Light")

##Make sweet, sweet figure to visualize percentage dominance
quartz()
par(mfrow=c(2,3))
hist(AB_onlyT2_Dark$Per_A,col="red",main="",xlab="Percent A in AB-Dark-T2")
hist(AD_onlyT2_Dark$Per_A,col="red",main="",xlab="Percent A in AD-Dark")
hist(BD_onlyT2_Dark$Per_B,col="blue",main="",xlab="Percent B in BD-Dark")
hist(AB_onlyT2_Dark$Per_B,col="blue",main="",xlab="Percent B in AB-Dark")
hist(AD_onlyT2_Dark$Per_D,col="green",main="",xlab="Percent D in AD-Dark")
hist(BD_onlyT2_Dark$Per_D,col="green",main="",xlab="Percent D in BD-Dark")

##Make sweet, sweet figure to visualize percentage dominance
quartz()
par(mfrow=c(2,3))
hist(AB_onlyT2_Light$Per_A,col="red",main="",xlab="Percent A in AB-Light-T2")
hist(AD_onlyT2_Light$Per_A,col="red",main="",xlab="Percent A in AD-Light")
hist(BD_onlyT2_Light$Per_B,col="blue",main="",xlab="Percent B in BD-Light")
hist(AB_onlyT2_Light$Per_B,col="blue",main="",xlab="Percent B in AB-Light")
hist(AD_onlyT2_Light$Per_D,col="green",main="",xlab="Percent D in AD-Light")
hist(BD_onlyT2_Light$Per_D,col="green",main="",xlab="Percent D in BD-Light")

