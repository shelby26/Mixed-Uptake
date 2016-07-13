##Make variables of cells counts for cct1 and cct2 in qPCR dominance scripts pooling all times and environments into Sym_ID categories

##Separate Data for AB_only
AB_only=subset(qPCR,qPCR$Sym_ID=="AB"&qPCR$D_Count==0)
AD_only=subset(qPCR,qPCR$Sym_ID=="AD"&qPCR$B_Count==0)
BD_only=subset(qPCR,qPCR$Sym_ID=="BD"&qPCR$A_Count==0)

ct1AB_A_Count=AB_only$A_Count
ct1AB_A_Count[ct1AB_A_Count==0]<-1
ct2AB_B_Count=AB_only$B_Count
ct2AB_B_Count[ct2AB_B_Count==0]<-1

ct1AD_A_Count=AD_only$A_Count
ct1AD_A_Count[ct1AD_A_Count==0]<-1
ct2AD_D_Count=AD_only$D_Count
ct2AD_D_Count[ct2AD_D_Count==0]<-1

ct1BD_B_Count=BD_only$B_Count
ct1BD_B_Count[ct1BD_B_Count==0]<-1
ct2BD_D_Count=BD_only$D_Count
ct2BD_D_Count[ct2BD_D_Count==0]<-1

##Run each subset through both qPCR dominance scripts (log based and percent based)
cat("***ALL AB ONLY***")
qPCR_percentdom(ct1AB_A_Count,ct2AB_B_Count,.95,100)
qPCR_dominance(ct1AB_A_Count,ct2AB_B_Count,.95,100)

cat("***ALL AD ONLY***")
qPCR_percentdom(ct1AD_A_Count,ct2AD_D_Count,.95,100)
qPCR_dominance(ct1AD_A_Count,ct2AD_D_Count,.95,100)

cat("***ALL BD ONLY***")
qPCR_percentdom(ct1BD_B_Count,ct2BD_D_Count,.95,100)
qPCR_dominance(ct1BD_B_Count,ct2BD_D_Count,.95,100)

quartz()
par(mfrow=c(2,3))
hist(AB_only$Per_A,col="red",main="",xlab="Percent A in ALL AB")
hist(AD_only$Per_A,col="red",main="",xlab="Percent A in ALL AD")
hist(BD_only$Per_B,col="blue",main="",xlab="Percent B in ALL BD")
hist(AB_only$Per_B,col="blue",main="",xlab="Percent B in ALL AB")
hist(AD_only$Per_D,col="green",main="",xlab="Percent D in ALL AD")
hist(BD_only$Per_D,col="green",main="",xlab="Percent D in ALL BD")
