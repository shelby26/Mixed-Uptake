library(plyr)

# Load data
data <- read.table("qPCR Data_cont rem_R.txt", header=T)
# Subset into time 1 and time 2
data1 <- subset(data, Time==1)
data2 <- subset(data, Time==2)


# This analysis takes a two part ("delta") approach -- first analyzing the probability of an observation
# being zero, and the second part analyzing non-zero observations coming from a log-normal distribution.

# Calculate time 1 probability of zero -----
Azero <- aggregate(data.frame(A.Zero=data1$A_Count), 
                   by=list(data1$Sym_ID, data1$Env), FUN=function(x) table(x==0)[2]/length(x))
Bzero <- aggregate(data.frame(B.Zero=data1$B_Count), 
                   by=list(data1$Sym_ID, data1$Env), FUN=function(x) table(x==0)[2]/length(x))
Dzero <- aggregate(data.frame(D.Zero=data1$D_Count), 
                   by=list(data1$Sym_ID, data1$Env), FUN=function(x) table(x==0)[2]/length(x))
Zeros <- join_all(list(Azero, Bzero, Dzero))
colnames(Zeros) <- c("Sym_ID", "Env", "A.Zero", "B.Zero", "D.Zero")
Zeros[is.na(Zeros)] <- 0

# Calculate time 1 means (not including zeros) -----
Ares <- aggregate(data.frame(A.Count=data1$A_Count), 
                  by=list(data1$Sym_ID, data1$Env), FUN=function(x) exp(mean(log(x[x!=0]))))
Bres <- aggregate(data.frame(B.Count=data1$B_Count), 
                  by=list(data1$Sym_ID, data1$Env), FUN=function(x) exp(mean(log(x[x!=0]))))
Dres <- aggregate(data.frame(D.Count=data1$D_Count), 
                  by=list(data1$Sym_ID, data1$Env), FUN=function(x) exp(mean(log(x[x!=0]))))
res <- join_all(list(Ares, Bres, Dres), type="full")
colnames(res) <- c("Sym_ID", "Env", "A.Count", "B.Count", "D.Count")
res[is.na(res)] <- 0


# Plot time 1 -----
par(mar=c(4,3,1,1))
layout(mat=matrix(c(1,2,3,4,5,6), nrow=3, byrow=F))
# Plot Dark
ABdf <- subset(res, Sym_ID %in% c("A","AB","B") & Env=="Dark")
plot(A.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf), type="b", col="red", ylim=c(0,5000),
     main="AB - Low Light", xaxt="n", xlab="")
axis(side=1, at=c(1,2,3), labels=c("A only", "AB mix", "B only"))
lines(A.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf)[-2,], lty=2, col="red")
points(B.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf), type="b", col="blue")
lines(B.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf)[-2,], lty=2, col="blue")
points(A.Count + B.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf), type="b", col="darkgray")
lines(A.Count + B.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf)[-2,], lty=2, col="darkgray")
mtext(side=1, line=2, text=paste(with(Zeros, round(1-Zeros[Sym_ID=="AB" & Env=="Dark", "A.Zero"],2)), "/",
                                   with(Zeros, round(1-Zeros[Sym_ID=="AB" & Env=="Dark", "B.Zero"],2))),
      cex=0.75)


ADdf <- subset(res, Sym_ID %in% c("A","AD","D") & Env=="Dark")
plot(A.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf), type="b", col="red", ylim=c(0,5000),
     main="AD - Low Light", xaxt="n", xlab="")
axis(side=1, at=c(1,2,3), labels=c("A only", "AD mix", "D only"))
lines(A.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf)[-2,], lty=2, col="red")
points(D.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf), type="b", col="purple")
lines(D.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf)[-2,], lty=2, col="purple")
points(A.Count + D.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf), type="b", col="darkgray")
lines(A.Count + D.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf)[-2,], lty=2, col="darkgray")
mtext(side=1, line=2, text=paste(with(Zeros, round(1-Zeros[Sym_ID=="AD" & Env=="Dark", "A.Zero"],2)), "/",
                                 with(Zeros, round(1-Zeros[Sym_ID=="AD" & Env=="Dark", "D.Zero"],2))),
      cex=0.75)


BDdf <- subset(res, Sym_ID %in% c("B","BD","D") & Env=="Dark")
plot(B.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf), type="b", col="blue", ylim=c(0,5000),
     main="BD - Low Light", xaxt="n", xlab="")
axis(side=1, at=c(1,2,3), labels=c("B only", "BD mix", "D only"))
lines(B.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf)[-2,], lty=2, col="blue")
points(D.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf), type="b", col="purple")
lines(D.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf)[-2,], lty=2, col="purple")
points(B.Count + D.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf), type="b", col="darkgray")
lines(B.Count + D.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf)[-2,], lty=2, col="darkgray")
mtext(side=1, line=2, font=2, text=paste(with(Zeros, round(1-Zeros[Sym_ID=="BD" & Env=="Dark", "B.Zero"],2)), "/",
                                 with(Zeros, round(1-Zeros[Sym_ID=="BD" & Env=="Dark", "D.Zero"],2))),
      cex=0.75)

# Plot Light
ABdf <- subset(res, Sym_ID %in% c("A","AB","B") & Env=="Light")
plot(A.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf), type="b", col="red", ylim=c(0,15000),
     main="AB - High Light", xaxt="n", xlab="")
axis(side=1, at=c(1,2,3), labels=c("A only", "AB mix", "B only"))
lines(A.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf)[-2,], lty=2, col="red")
points(B.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf), type="b", col="blue")
lines(B.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf)[-2,], lty=2, col="blue")
points(A.Count + B.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf), type="b", col="darkgray")
lines(A.Count + B.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf)[-2,], lty=2, col="darkgray")
mtext(side=1, line=2, font=2, text=paste(with(Zeros, round(1-Zeros[Sym_ID=="AB" & Env=="Light", "A.Zero"],2)), "/",
                                 with(Zeros, round(1-Zeros[Sym_ID=="AB" & Env=="Light", "B.Zero"],2))),
      cex=0.75)

ADdf <- subset(res, Sym_ID %in% c("A","AD","D") & Env=="Light")
plot(A.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf), type="b", col="red", ylim=c(0,15000),
     main="AD - High Light", xaxt="n", xlab="")
axis(side=1, at=c(1,2,3), labels=c("A only", "AD mix", "D only"))
lines(A.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf)[-2,], lty=2, col="red")
points(D.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf), type="b", col="purple")
lines(D.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf)[-2,], lty=2, col="purple")
points(A.Count + D.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf), type="b", col="darkgray")
lines(A.Count + D.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf)[-2,], lty=2, col="darkgray")
mtext(side=1, line=2, text=paste(with(Zeros, round(1-Zeros[Sym_ID=="AD" & Env=="Light", "A.Zero"],2)), "/",
                                 with(Zeros, round(1-Zeros[Sym_ID=="AD" & Env=="Light", "D.Zero"],2))),
      cex=0.75)

BDdf <- subset(res, Sym_ID %in% c("B","BD","D") & Env=="Light")
plot(B.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf), type="b", col="blue", ylim=c(0,15000),
     main="BD - High Light", xaxt="n", xlab="")
axis(side=1, at=c(1,2,3), labels=c("B only", "BD mix", "D only"))
lines(B.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf)[-2,], lty=2, col="blue")
points(D.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf), type="b", col="purple")
lines(D.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf)[-2,], lty=2, col="purple")
points(B.Count + D.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf), type="b", col="darkgray")
lines(B.Count + D.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf)[-2,], lty=2, col="darkgray")
mtext(side=1, line=2, text=paste(with(Zeros, round(1-Zeros[Sym_ID=="BD" & Env=="Light", "B.Zero"]), 2), "/",
                                 with(Zeros, round(1-Zeros[Sym_ID=="BD" & Env=="Light", "D.Zero"]), 2)), cex=0.75)



# Calculate time 2 probability of zero -----
Azero <- aggregate(data.frame(A.Zero=data2$A_Count), 
                   by=list(data2$Sym_ID, data2$Env), FUN=function(x) table(x==0)[2]/length(x))
Bzero <- aggregate(data.frame(B.Zero=data2$B_Count), 
                   by=list(data2$Sym_ID, data2$Env), FUN=function(x) table(x==0)[2]/length(x))
Dzero <- aggregate(data.frame(D.Zero=data2$D_Count), 
                   by=list(data2$Sym_ID, data2$Env), FUN=function(x) table(x==0)[2]/length(x))
Zeros <- join_all(list(Azero, Bzero, Dzero))
colnames(Zeros) <- c("Sym_ID", "Env", "A.Zero", "B.Zero", "D.Zero")
Zeros[is.na(Zeros)] <- 0

# Calculate time 2 means (not including zeros) -----
Ares <- aggregate(data.frame(A.Count=data2$A_Count), 
                  by=list(data2$Sym_ID, data2$Env), FUN=function(x) exp(mean(log(x[x!=0]))))
Bres <- aggregate(data.frame(B.Count=data2$B_Count), 
                  by=list(data2$Sym_ID, data2$Env), FUN=function(x) exp(mean(log(x[x!=0]))))
Dres <- aggregate(data.frame(D.Count=data2$D_Count), 
                  by=list(data2$Sym_ID, data2$Env), FUN=function(x) exp(mean(log(x[x!=0]))))
res <- join_all(list(Ares, Bres, Dres), type="full")
colnames(res) <- c("Sym_ID", "Env", "A.Count", "B.Count", "D.Count")
res[is.na(res)] <- 0


# Plot time 2 -----
par(mar=c(4,3,1,1))
layout(mat=matrix(c(1,2,3,4,5,6), nrow=3, byrow=F))
# Plot Dark
ABdf <- subset(res, Sym_ID %in% c("A","AB","B") & Env=="Dark")
plot(A.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf), type="b", col="red", ylim=c(0,40000),
     main="AB - Low Light", xaxt="n", xlab="")
axis(side=1, at=c(1,2,3), labels=c("A only", "AB mix", "B only"))
lines(A.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf)[-2,], lty=2, col="red")
points(B.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf), type="b", col="blue")
lines(B.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf)[-2,], lty=2, col="blue")
points(A.Count + B.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf), type="b", col="darkgray")
lines(A.Count + B.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf)[-2,], lty=2, col="darkgray")
mtext(side=1, line=2, font=2, text=paste(with(Zeros, round(1-Zeros[Sym_ID=="AB" & Env=="Dark", "A.Zero"],2)), "/",
                                 with(Zeros, round(1-Zeros[Sym_ID=="AB" & Env=="Dark", "B.Zero"],2))),
      cex=0.75)


ADdf <- subset(res, Sym_ID %in% c("A","AD","D") & Env=="Dark")
plot(A.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf), type="b", col="red", ylim=c(0,40000),
     main="AD - Low Light", xaxt="n", xlab="")
axis(side=1, at=c(1,2,3), labels=c("A only", "AD mix", "D only"))
lines(A.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf)[-2,], lty=2, col="red")
points(D.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf), type="b", col="purple")
lines(D.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf)[-2,], lty=2, col="purple")
points(A.Count + D.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf), type="b", col="darkgray")
lines(A.Count + D.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf)[-2,], lty=2, col="darkgray")
mtext(side=1, line=2, text=paste(with(Zeros, round(1-Zeros[Sym_ID=="AD" & Env=="Dark", "A.Zero"],2)), "/",
                                 with(Zeros, round(1-Zeros[Sym_ID=="AD" & Env=="Dark", "D.Zero"],2))),
      cex=0.75)


BDdf <- subset(res, Sym_ID %in% c("B","BD","D") & Env=="Dark")
plot(B.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf), type="b", col="blue", ylim=c(0,40000),
     main="BD - Low Light", xaxt="n", xlab="")
axis(side=1, at=c(1,2,3), labels=c("B only", "BD mix", "D only"))
lines(B.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf)[-2,], lty=2, col="blue")
points(D.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf), type="b", col="purple")
lines(D.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf)[-2,], lty=2, col="purple")
points(B.Count + D.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf), type="b", col="darkgray")
lines(B.Count + D.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf)[-2,], lty=2, col="darkgray")
text(c(-10000,10000,-10000),"*",cex=2,col="purple")
text(c(-10000,0,-10000),"*",cex=2,col="blue")
mtext(side=1, line=2, text=paste(with(Zeros, round(1-Zeros[Sym_ID=="BD" & Env=="Dark", "B.Zero"],2)), "/",
                                 with(Zeros, round(1-Zeros[Sym_ID=="BD" & Env=="Dark", "D.Zero"],))),
      cex=0.75)

# Plot Light
ABdf <- subset(res, Sym_ID %in% c("A","AB","B") & Env=="Light")
plot(A.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf), type="b", col="red", ylim=c(0,55000),
     main="AB - High Light", xaxt="n", xlab="")
axis(side=1, at=c(1,2,3), labels=c("A only", "AB mix", "B only"))
lines(A.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf)[-2,], lty=2, col="red")
points(B.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf), type="b", col="blue")
lines(B.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf)[-2,], lty=2, col="blue")
points(A.Count + B.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf), type="b", col="darkgray")
lines(A.Count + B.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf)[-2,], lty=2, col="darkgray")
text(c(-10000,0,-10000),"*",cex=2,col="red")
mtext(side=1, line=2, text=paste(with(Zeros, round(1-Zeros[Sym_ID=="AB" & Env=="Light", "A.Zero"],2)), "/",
                                 with(Zeros, round(1-Zeros[Sym_ID=="AB" & Env=="Light", "B.Zero"],2))),
      cex=0.75)

ADdf <- subset(res, Sym_ID %in% c("A","AD","D") & Env=="Light")
plot(A.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf), type="b", col="red", ylim=c(0,55000),
     main="AD - High Light", xaxt="n", xlab="")
axis(side=1, at=c(1,2,3), labels=c("A only", "AD mix", "D only"))
lines(A.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf)[-2,], lty=2, col="red")
points(D.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf), type="b", col="purple")
lines(D.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf)[-2,], lty=2, col="purple")
points(A.Count + D.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf), type="b", col="darkgray")
lines(A.Count + D.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf)[-2,], lty=2, col="darkgray")
text(c(-10000,0,-10000),"*",cex=2,col="purple")
mtext(side=1, line=2, text=paste(with(Zeros, round(1-Zeros[Sym_ID=="AD" & Env=="Light", "A.Zero"],2)), "/",
                                 with(Zeros, round(1-Zeros[Sym_ID=="AD" & Env=="Light", "D.Zero"],2))),
      cex=0.75)

BDdf <- subset(res, Sym_ID %in% c("B","BD","D") & Env=="Light")
plot(B.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf), type="b", col="blue", ylim=c(0,55000),
     main="BD - High Light", xaxt="n", xlab="")
axis(side=1, at=c(1,2,3), labels=c("B only", "BD mix", "D only"))
lines(B.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf)[-2,], lty=2, col="blue")
points(D.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf), type="b", col="purple")
lines(D.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf)[-2,], lty=2, col="purple")
points(B.Count + D.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf), type="b", col="darkgray")
lines(B.Count + D.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf)[-2,], lty=2, col="darkgray")
mtext(side=1, line=2, text=paste(with(Zeros, round(1-Zeros[Sym_ID=="BD" & Env=="Light", "B.Zero"]), 2), "/",
                                 with(Zeros, round(1-Zeros[Sym_ID=="BD" & Env=="Light", "D.Zero"]), 2)), cex=0.75)

