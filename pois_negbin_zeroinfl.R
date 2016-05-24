# Website with info on poisson/negbin models
#http://datavoreconsulting.com/programming-tips/count-data-glms-choosing-poisson-negative-binomial-zero-inflated-poisson/

library(bbmle)
library(fishMod)
library(MASS)
library(pscl)
library(plyr)

# Load data
data <- read.table("qPCR Data_cont rem_R.txt", header=T)
data1 <- subset(data, Time==1)
data2 <- subset(data, Time==2)
# Create newdata frame for predicting
nd <- expand.grid(Sym_ID=c("A","AB","AD","B","BD","D"), Env=c("Dark","Light"))
# Calculate raw means of each group...
Ameans <- aggregate(data1$A_Count, by=list(data1$Sym_ID, data1$Env), FUN=mean)
Ameans
aggregate(data1$A_Count, by=list(data1$Sym_ID, data1$Env), FUN=c)


# TIME POINT 1
# CLADE A COUNTS -----
# Test for Poisson data
Adf <- subset(data1, Sym_ID %in% c("A","AB","AD"))
Andf <- subset(nd, Sym_ID %in% c("A","AB","AD"))
model.pois <- glm(A_Count ~ Sym_ID * Env, data=Adf, family="poisson")
summary(model.pois)
# Goodness of fit test
1 - pchisq(summary(model.pois)$deviance, 
           summary(model.pois)$df.residual)  # p=0, reject poisson distribution

# Test for negative binomial data
model.nb = glm.nb(A_Count ~ Sym_ID * Env, data = Adf)
summary(model.nb)
# Goodness of fit test
1 - pchisq(summary(model.nb)$deviance,
           summary(model.nb)$df.residual)  # p=0, negative binomial model does not fit the data

# Test zero-inflated poisson distribution
model.zip = zeroinfl(A_Count ~ Sym_ID*Env|Sym_ID*Env, data = Adf)
summary(model.zip)
# Fitted values
Amod <- cbind(Andf, 
      Count = predict(model.zip, newdata = Andf, type = "count"),
      Zero = predict(model.zip, newdata = Andf, type = "zero"))
plot(residuals(model.zip, "response"))

# Test zero-inflated negative binomial (in the case of overdispersion)
model.zip.3 = zeroinfl(A_Count ~ Sym_ID*Env|Sym_ID*Env, data = Adf, dist = "negbin")
summary(model.zip.3)  # Theta parameter is not significant --> poisson is appropriate

# CLADE B COUNTS -----
# Test for Poisson data
Bdf <- subset(data1, Sym_ID %in% c("AB","B","BD"))
Bndf <- subset(nd, Sym_ID %in% c("AB","B","BD"))
model.pois <- glm(B_Count ~ Sym_ID * Env, data=Bdf, family="poisson")
summary(model.pois)
# Goodness of fit test
1 - pchisq(summary(model.pois)$deviance, 
           summary(model.pois)$df.residual)  # p=0, reject poisson distribution

# Test for negative binomial data
model.nb = glm.nb(B_Count ~ Sym_ID * Env, data = Bdf)
summary(model.nb)
# Goodness of fit test
1 - pchisq(summary(model.nb)$deviance,
           summary(model.nb)$df.residual)  # p=0, negative binomial model does not fit the data

# Test zero-inflated poisson distribution
model.zip = zeroinfl(B_Count ~ Sym_ID*Env|Sym_ID*Env, data = Bdf)
summary(model.zip)
# Fitted values
Bmod <- cbind(Bndf, 
              Count = predict(model.zip, newdata = Bndf, type = "count"),
              Zero = predict(model.zip, newdata = Bndf, type = "zero"))
plot(residuals(model.zip, "response"))

# Test zero-inflated negative binomial (in the case of overdispersion)
model.zip.3 = zeroinfl(B_Count ~ Sym_ID*Env|Sym_ID*Env, data = Bdf, dist = "negbin")
summary(model.zip.3)  # Theta parameter is not significant --> poisson is appropriate

# CLADE D COUNTS -----
# Test for Poisson data
Ddf <- subset(data1, Sym_ID %in% c("AD","BD","D"))
Dndf <- subset(nd, Sym_ID %in% c("AD","BD","D"))
model.pois <- glm(D_Count ~ Sym_ID * Env, data=Ddf, family="poisson")
summary(model.pois)
# Goodness of fit test
1 - pchisq(summary(model.pois)$deviance, 
           summary(model.pois)$df.residual)  # p=0, reject poisson distribution

# Test for negative binomial data
model.nb = glm.nb(D_Count ~ Sym_ID * Env, data = Ddf)
summary(model.nb)
# Goodness of fit test
1 - pchisq(summary(model.nb)$deviance,
           summary(model.nb)$df.residual)  # p=0.023, negative binomial model does not fit the data

# Test zero-inflated poisson distribution
model.zip = zeroinfl(D_Count ~ Sym_ID*Env|Sym_ID*Env, data = Ddf)
summary(model.zip)
# Fitted values
Dmod <- cbind(Dndf, 
              Count = predict(model.zip, newdata = Dndf, type = "count"),
              Zero = predict(model.zip, newdata = Dndf, type = "zero"))
plot(residuals(model.zip, "response"))

# Test zero-inflated negative binomial (in the case of overdispersion)
model.zip.3 = zeroinfl(D_Count ~ Sym_ID*Env|Sym_ID*Env, data = Ddf, dist = "negbin")
summary(model.zip.3)  # Theta parameter is not significant --> poisson is appropriate


# RESULTS -----
colnames(Amod) <- c("Sym_ID", "Env", "A.Count", "A.Zero")
colnames(Bmod) <- c("Sym_ID", "Env", "B.Count", "B.Zero")
colnames(Dmod) <- c("Sym_ID", "Env", "D.Count", "D.Zero")
res <- join_all(dfs = list(Amod, Bmod, Dmod), type="full")


ABdf <- subset(res, Sym_ID %in% c("A","AB","B") & Env=="Light")
ABdf[is.na(ABdf)] <- 0

plot(A.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf), type="b", col="red")
lines(A.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf)[-2,], lty=2, col="red")
points(B.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf), type="b", col="blue")
lines(B.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf)[-2,], lty=2, col="blue")
points(A.Count + B.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf), type="b", col="darkgray")
lines(A.Count + B.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf)[-2,], lty=2, col="darkgray")


ADdf <- subset(res, Sym_ID %in% c("A","AD","D") & Env=="Light")
ADdf[is.na(ADdf)] <- 0

plot(A.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf), type="b", col="red", ylim=c(0,30000))
lines(A.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf)[-2,], lty=2, col="red")
points(D.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf), type="b", col="purple")
lines(D.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf)[-2,], lty=2, col="purple")
points(A.Count + D.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf), type="b", col="darkgray")
lines(A.Count + D.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf)[-2,], lty=2, col="darkgray")


BDdf <- subset(res, Sym_ID %in% c("B","BD","D") & Env=="Light")
BDdf[is.na(BDdf)] <- 0

plot(B.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf), type="b", col="blue", ylim=c(0,30000))
lines(B.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf)[-2,], lty=2, col="blue")
points(D.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf), type="b", col="purple")
lines(D.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf)[-2,], lty=2, col="purple")
points(B.Count + D.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf), type="b", col="darkgray")
lines(B.Count + D.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf)[-2,], lty=2, col="darkgray")



# LOW LIGHT
ABdf <- subset(res, Sym_ID %in% c("A","AB","B") & Env=="Dark")
ABdf[is.na(ABdf)] <- 0

plot(A.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf), type="b", col="red", ylim=c(0,4000))
lines(A.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf)[-2,], lty=2, col="red")
points(B.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf), type="b", col="blue")
lines(B.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf)[-2,], lty=2, col="blue")
points(A.Count + B.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf), type="b", col="darkgray")
lines(A.Count + B.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf)[-2,], lty=2, col="darkgray")


ADdf <- subset(res, Sym_ID %in% c("A","AD","D") & Env=="Dark")
ADdf[is.na(ADdf)] <- 0

plot(A.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf), type="b", col="red", ylim=c(0,10000))
lines(A.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf)[-2,], lty=2, col="red")
points(D.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf), type="b", col="purple")
lines(D.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf)[-2,], lty=2, col="purple")
points(A.Count + D.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf), type="b", col="darkgray")
lines(A.Count + D.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf)[-2,], lty=2, col="darkgray")


BDdf <- subset(res, Sym_ID %in% c("B","BD","D") & Env=="Dark")
BDdf[is.na(BDdf)] <- 0

plot(B.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf), type="b", col="blue", ylim=c(0,10000))
lines(B.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf)[-2,], lty=2, col="blue")
points(D.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf), type="b", col="purple")
lines(D.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf)[-2,], lty=2, col="purple")
points(B.Count + D.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf), type="b", col="darkgray")
lines(B.Count + D.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf)[-2,], lty=2, col="darkgray")



# TIME POINT 2
# CLADE A COUNTS -----
# Test for Poisson data
Adf <- subset(data2, Sym_ID %in% c("A","AB","AD"))
Andf <- subset(nd, Sym_ID %in% c("A","AB","AD"))
model.pois <- glm(A_Count ~ Sym_ID * Env, data=Adf, family="poisson")
summary(model.pois)
# Goodness of fit test
1 - pchisq(summary(model.pois)$deviance, 
           summary(model.pois)$df.residual)  # p=0, reject poisson distribution

# Test for negative binomial data
model.nb = glm.nb(A_Count ~ Sym_ID * Env, data = Adf)
summary(model.nb)
# Goodness of fit test
1 - pchisq(summary(model.nb)$deviance,
           summary(model.nb)$df.residual)  # p=0, negative binomial model does not fit the data

# Test zero-inflated poisson distribution   # DOES NOT FIT.....maybe because of overdispersion?
model.zip = zeroinfl(A_Count ~ Sym_ID*Env|Sym_ID*Env, data = Adf)
summary(model.zip)
# Fitted values

# Test zero-inflated negative binomial (in the case of overdispersion)
model.zip.3 = zeroinfl(A_Count ~ Sym_ID*Env|Sym_ID*Env, data = Adf, dist = "negbin")
summary(model.zip.3)  # Theta parameter is significant --> NEGBIN is appropriate
Amod <- cbind(Andf, 
              Count = predict(model.zip.3, newdata = Andf, type = "count"),
              Zero = predict(model.zip.3, newdata = Andf, type = "zero"))
plot(residuals(model.zip.3, "response"))

# CLADE B COUNTS -----
# Test for Poisson data
Bdf <- subset(data2, Sym_ID %in% c("AB","B","BD"))
Bndf <- subset(nd, Sym_ID %in% c("AB","B","BD"))
model.pois <- glm(B_Count ~ Sym_ID * Env, data=Bdf, family="poisson")
summary(model.pois)
# Goodness of fit test
1 - pchisq(summary(model.pois)$deviance, 
           summary(model.pois)$df.residual)  # p=0, reject poisson distribution

# Test for negative binomial data
model.nb = glm.nb(B_Count ~ Sym_ID * Env, data = Bdf)
summary(model.nb)
# Goodness of fit test
1 - pchisq(summary(model.nb)$deviance,
           summary(model.nb)$df.residual)  # p=0, negative binomial model does not fit the data

# Test zero-inflated poisson distribution
model.zip = zeroinfl(B_Count ~ Sym_ID*Env|Sym_ID*Env, data = Bdf)
summary(model.zip)

# Test zero-inflated negative binomial (in the case of overdispersion)
model.zip.3 = zeroinfl(B_Count ~ Sym_ID*Env|Sym_ID*Env, data = Bdf, dist = "negbin")
summary(model.zip.3)  # Theta parameter is significant --> NEGBIN is appropriate
# Fitted values
Bmod <- cbind(Bndf, 
              Count = predict(model.zip.3, newdata = Bndf, type = "count"),
              Zero = predict(model.zip.3, newdata = Bndf, type = "zero"))
plot(residuals(model.zip.3, "response"))

# CLADE D COUNTS -----
# Test for Poisson data
Ddf <- subset(data2, Sym_ID %in% c("AD","BD","D"))
Dndf <- subset(nd, Sym_ID %in% c("AD","BD","D"))
model.pois <- glm(D_Count ~ Sym_ID * Env, data=Ddf, family="poisson")
summary(model.pois)
# Goodness of fit test
1 - pchisq(summary(model.pois)$deviance, 
           summary(model.pois)$df.residual)  # p=0, reject poisson distribution

# Test for negative binomial data
model.nb = glm.nb(D_Count ~ Sym_ID * Env, data = Ddf)
summary(model.nb)
# Goodness of fit test
1 - pchisq(summary(model.nb)$deviance,
           summary(model.nb)$df.residual)  # p=0.023, negative binomial model does not fit the data

# Test zero-inflated poisson distribution  # DOES NOT FIT -- overdispersion?
model.zip = zeroinfl(D_Count ~ Sym_ID*Env|Sym_ID*Env, data = Ddf)
summary(model.zip)

# Test zero-inflated negative binomial (in the case of overdispersion)
model.zip.3 = zeroinfl(D_Count ~ Sym_ID*Env|Sym_ID*Env, data = Ddf, dist = "negbin")
summary(model.zip.3)  # Theta parameter is significant --> NEGBIN is appropriate
# Fitted values
Dmod <- cbind(Dndf, 
              Count = predict(model.zip.3, newdata = Dndf, type = "count"),
              Zero = predict(model.zip.3, newdata = Dndf, type = "zero"))
plot(residuals(model.zip.3, "response"))



# RESULTS -----
colnames(Amod) <- c("Sym_ID", "Env", "A.Count", "A.Zero")
colnames(Bmod) <- c("Sym_ID", "Env", "B.Count", "B.Zero")
colnames(Dmod) <- c("Sym_ID", "Env", "D.Count", "D.Zero")
res <- join_all(dfs = list(Amod, Bmod, Dmod), type="full")
res

ABdf <- subset(res, Sym_ID %in% c("A","AB","B") & Env=="Light")
ABdf[is.na(ABdf)] <- 0

plot(A.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf), type="b", col="red")
lines(A.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf)[-2,], lty=2, col="red")
points(B.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf), type="b", col="blue")
lines(B.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf)[-2,], lty=2, col="blue")
points(A.Count + B.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf), type="b", col="darkgray")
lines(A.Count + B.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf)[-2,], lty=2, col="darkgray")


ADdf <- subset(res, Sym_ID %in% c("A","AD","D") & Env=="Light")
ADdf[is.na(ADdf)] <- 0

plot(A.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf), type="b", col="red", ylim=c(0,160000))
lines(A.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf)[-2,], lty=2, col="red")
points(D.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf), type="b", col="purple")
lines(D.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf)[-2,], lty=2, col="purple")
points(A.Count + D.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf), type="b", col="darkgray")
lines(A.Count + D.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf)[-2,], lty=2, col="darkgray")


BDdf <- subset(res, Sym_ID %in% c("B","BD","D") & Env=="Light")
BDdf[is.na(BDdf)] <- 0

plot(B.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf), type="b", col="blue", ylim=c(0,160000))
lines(B.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf)[-2,], lty=2, col="blue")
points(D.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf), type="b", col="purple")
lines(D.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf)[-2,], lty=2, col="purple")
points(B.Count + D.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf), type="b", col="darkgray")
lines(B.Count + D.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf)[-2,], lty=2, col="darkgray")



# LOW LIGHT
ABdf <- subset(res, Sym_ID %in% c("A","AB","B") & Env=="Dark")
ABdf[is.na(ABdf)] <- 0

plot(A.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf), type="b", col="red", ylim=c(0,90000))
lines(A.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf)[-2,], lty=2, col="red")
points(B.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf), type="b", col="blue")
lines(B.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf)[-2,], lty=2, col="blue")
points(A.Count + B.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf), type="b", col="darkgray")
lines(A.Count + B.Count ~ as.numeric(Sym_ID), data=droplevels(ABdf)[-2,], lty=2, col="darkgray")


ADdf <- subset(res, Sym_ID %in% c("A","AD","D") & Env=="Dark")
ADdf[is.na(ADdf)] <- 0

plot(A.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf), type="b", col="red", ylim=c(0,100000))
lines(A.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf)[-2,], lty=2, col="red")
points(D.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf), type="b", col="purple")
lines(D.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf)[-2,], lty=2, col="purple")
points(A.Count + D.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf), type="b", col="darkgray")
lines(A.Count + D.Count ~ as.numeric(Sym_ID), data=droplevels(ADdf)[-2,], lty=2, col="darkgray")


BDdf <- subset(res, Sym_ID %in% c("B","BD","D") & Env=="Dark")
BDdf[is.na(BDdf)] <- 0

plot(B.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf), type="b", col="blue", ylim=c(0,100000))
lines(B.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf)[-2,], lty=2, col="blue")
points(D.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf), type="b", col="purple")
lines(D.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf)[-2,], lty=2, col="purple")
points(B.Count + D.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf), type="b", col="darkgray")
lines(B.Count + D.Count ~ as.numeric(Sym_ID), data=droplevels(BDdf)[-2,], lty=2, col="darkgray")








# TRY DELTA LN MODELS -----
Adf <- subset(data1, Sym_ID %in% c("A","AB","AD"))
Andf <- subset(nd, Sym_ID %in% c("A","AB","AD"))
fm.dln <- deltaLN(ln.form=A_Count~Sym_ID*Env, binary.form=~Sym_ID*Env, data=droplevels(Adf), residuals=T)
Adf$fitted <- fm.dln$fitted
aggregate(Adf$fitted, by=list(Adf$Sym_ID, Adf$Env), FUN=mean)
head(Adf)
plot(fm.dln$fitted, fm.dln$residuals[,"quantile"], pch=20, main="Delta Log-Normal quantile residuals")
model.frame(fm.dln$lnMod)
summary(fm.dln$lnMod)
predict(fm.dln$lnMod, Andf)
aggregate(Adf$A_Count, by=list(Adf$Sym_ID, Adf$Env), FUN=mean)
aggregate(Adf$A_Count, by=list(Adf$Sym_ID, Adf$Env), FUN=function(x) mean(x[x!=0]))
aggregate(Adf$A_Count, by=list(Adf$Sym_ID, Adf$Env), FUN=function(x) exp(mean(log(x[x!=0]))))

my.coef <- c(0.6, 1.2, 0, -0.3, 0, -0.5, 0.85)
sim.dat <- simReg( n=250, lambda.tau=my.coef[1:3], mu.Z.tau=my.coef[4:6], alpha=my.coef[7])
fm <- pgm( p.form=y~1+x1+x2, g.form=~1+x1+x2, data=sim.dat)
tmp <- matrix( c( my.coef, fm$coef, sqrt( diag( fm$vcov))), ncol=3)
tmp[nrow( tmp),1] <- log( tmp[nrow( tmp),1])  #putting values on same scale
colnames( tmp) <- c("actual","estiated","SE")
rownames( tmp) <- names( fm$coef)
print( tmp)








fm.Tweedie <- tglm(A_Count~Sym_ID*Env, data=droplevels(Adf))  #estimate power param too!
fm.Tweedie$fitted
Adf$Tfitted <- fm.Tweedie$fitted
aggregate(Adf$Tfitted, by=list(Adf$Sym_ID, Adf$Env), FUN=mean)
plot( fm.Tweedie$fitted, fm.Tweedie$residuals[,"random"], pch=20, main="Tweedie GLM Randomised quantile residuals")
abline( h=0, col="red")




# Do my own delta model on log data
Adf <- subset(data1, Sym_ID %in% c("A","AB","AD"))
Andf <- subset(nd, Sym_ID %in% c("A","AB","AD"))
A0 <- subset(Adf, A_Count == 0)
A1 <- subset(Adf, A_Count > 0)

binom <- glm(A_Count==0 ~ Sym_ID * Env, data=Adf, family="binomial")
summary(binom)

lognorm <- lm(Log_A ~ Sym_ID * Env, data=A1)
summary(lognorm)
anova(lognorm)
cbind(Andf, predict(lognorm, Andf))
aggregate(Adf$A_Count, by=list(Adf$Sym_ID, Adf$Env), FUN=function(x) mean(log10(x[x!=0])))




# -----
Adf <- subset(data1, Sym_ID %in% c("A","AB","AD"))
Bdf <- subset(data1, Sym_ID %in% c("AB","B","BD"))
Ddf <- subset(data1, Sym_ID %in% c("AD","BD","D"))
Ares <- aggregate(data.frame(A.Count=Adf$A_Count), by=list(Adf$Sym_ID, Adf$Env), FUN=function(x) exp(mean(log(x[x!=0]))))
Bres <- aggregate(data.frame(B.Count=Bdf$B_Count), by=list(Bdf$Sym_ID, Bdf$Env), FUN=function(x) exp(mean(log(x[x!=0]))))
Dres <- aggregate(data.frame(D.Count=Ddf$D_Count), by=list(Ddf$Sym_ID, Ddf$Env), FUN=function(x) exp(mean(log(x[x!=0]))))
res <- join_all(list(Ares, Bres, Dres), type="full")
colnames(res) <- c("Sym_ID", "Env", "A.Count", "B.Count", "D.Count")
res[is.na(res)] <- 0
res



