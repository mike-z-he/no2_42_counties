## Standard Errors for Cumulative Effect EStimates - DLM

######DATA CLEANING (DEFAULT)##########
setwd("D:/Users/profu/Documents/Schoolwork/PhD/3rd Year Summer/ISEE Analysis")
options(mc.cores=parallel::detectCores())
dta <- read.csv("databydate-PM2.5update0404.csv", header=TRUE, stringsAsFactors = FALSE)
str(dta)
summary(dta)
names(dta)

library(mgcv)
library(splines)
library(dplyr)
#library(plyr)
library(dlnm)
library(mvmeta)

## Code for aggregating into monthly data (exposure, outcome, and covariates)
head(dta$mcode)

dta$newdate <- as.Date(dta$date, format="%m/%d/%Y")
dta$month   <- as.POSIXlt(dta$newdate)$mon + 1    
dta$year    <- as.POSIXlt(dta$newdate)$year + 1900
dta$dow     <- as.POSIXlt(dta$newdate)$wday

dta$pm25b <- ifelse(dta$pm2.5new %in% c(0,3), NA, dta$pm2.5new)

mort.dta <- aggregate(cbind(all, allg1, allg2, cir, cirg1, cirg2, res, resg1, resg2) ~ month + year + mcode + ecounty, data = dta, FUN=sum, na.rm=TRUE)
head(mort.dta, 20)

exp.dta <- aggregate(cbind(pm25b, no2, tmean, rh) ~ month + year + mcode + ecounty, data = dta, FUN=mean, na.rm=TRUE)
head(exp.dta, 20)

monthly.dta <- merge(mort.dta, exp.dta, by=c("month", "year", "mcode", "ecounty"), all=TRUE)
dim(monthly.dta)[1] == length(unique(dta$mcode))*12*3
summary(monthly.dta)

months <- unique(monthly.dta[c("month", "year")])
months <- na.omit(months[order(c(months$year, months$month)),])
months$monthvar <- 1:dim(months)[1] #months is 1-12, monthvar is 1-36

monthly.dta <- merge(monthly.dta, months, by=c("month", "year"))

table(monthly.dta$monthvar)
head(monthly.dta)
tail(monthly.dta)

## Adding other covariates
names(dta)
dta2 <- unique(dta[,c(5,33:48)])
dim(dta2)
head(dta2)

monthly.dta <- merge(monthly.dta, dta2, by="mcode", all=TRUE)
dim(monthly.dta)
monthly.dta$unemp.rate <- monthly.dta$unemployed/monthly.dta$pop

## Creating lags
monthly.dta <- monthly.dta %>% group_by(mcode) %>% mutate(lag.tmean = lag(tmean, 1))
monthly.dta <- monthly.dta %>% group_by(mcode) %>% mutate(lag.pm25b = lag(pm25b, 1))
monthly.dta <- monthly.dta %>% group_by(mcode) %>% mutate(lag.no2 = lag(no2, 1))

dta <- dta %>% group_by(mcode) %>% mutate(lag1.tmean = lag(tmean, 1))
dta <- dta %>% group_by(mcode) %>% mutate(lag2.tmean = lag(tmean, 2))
dta <- dta %>% group_by(mcode) %>% mutate(lag3.tmean = lag(tmean, 3))

dta <- dta %>% group_by(mcode) %>% mutate(lag.pm25b = lag(pm25b, 1))
dta <- dta %>% group_by(mcode) %>% mutate(lag2.pm25b = lag(pm25b, 2))
dta <- dta %>% group_by(mcode) %>% mutate(lag3.pm25b = lag(pm25b, 3))
dta <- dta %>% group_by(mcode) %>% mutate(lag4.pm25b = lag(pm25b, 4))
dta <- dta %>% group_by(mcode) %>% mutate(lag5.pm25b = lag(pm25b, 5))
dta <- dta %>% group_by(mcode) %>% mutate(lag6.pm25b = lag(pm25b, 6))
dta <- dta %>% group_by(mcode) %>% mutate(lag7.pm25b = lag(pm25b, 7))

dta <- dta %>% group_by(mcode) %>% mutate(lag.no2 = lag(no2, 1))
dta <- dta %>% group_by(mcode) %>% mutate(lag2.no2 = lag(no2, 2))
dta <- dta %>% group_by(mcode) %>% mutate(lag3.no2 = lag(no2, 3))
dta <- dta %>% group_by(mcode) %>% mutate(lag4.no2 = lag(no2, 4))
dta <- dta %>% group_by(mcode) %>% mutate(lag5.no2 = lag(no2, 5))
dta <- dta %>% group_by(mcode) %>% mutate(lag6.no2 = lag(no2, 6))
dta <- dta %>% group_by(mcode) %>% mutate(lag7.no2 = lag(no2, 7))


##Lags are created above (including monthly lag)
head(dta$lag1.tmean, 40)
head(dta$lag2.tmean, 40)
head(dta$lag3.tmean, 40)

dta$lag1to3.tmean <- (dta$lag1.tmean + dta$lag2.tmean + dta$lag3.tmean) / 3
head(dta$lag1to3.tmean)
##This is the average lag from day 1 to day 3

dta$lag0to1.no2 <- (dta$no2 + dta$lag.no2) / 2
dta$lag0to2.no2 <- (dta$no2 + dta$lag.no2 + dta$lag2.no2) / 3
dta$lag0to7.no2 <- (dta$no2 + dta$lag.no2 + dta$lag2.no2 + dta$lag3.no2 + dta$lag4.no2 +
                      dta$lag5.no2 + dta$lag6.no2 + dta$lag7.no2) / 8

dta$lag0to1.pm25b <- (dta$pm25b + dta$lag.pm25b) / 2
dta$lag0to2.pm25b <- (dta$pm25b + dta$lag.pm25b + dta$lag2.pm25b) / 3
dta$lag0to7.pm25b <- (dta$pm25b + dta$lag.pm25b + dta$lag2.pm25b + dta$lag3.pm25b + dta$lag4.pm25b +
                        dta$lag5.pm25b + dta$lag6.pm25b + dta$lag7.pm25b) / 8
##Lag 0-2 and lag 0-7 exposure

setwd("D:/Users/profu/Documents/Schoolwork/PhD/3rd Year Spring/Research/NO2 Manuscript")
numbers <- read.csv("County Numbers.csv", header=TRUE)

dta <- merge(dta, numbers)
monthly.dta <- merge(monthly.dta, numbers)

dta <- dta[order(dta$ecounty, dta$newdate),]

setwd("D:/Users/profu/Documents/Schoolwork/PhD/Research Projects/no2_42_counties/environmental_pollution_submission/round_2_revision/data")
holidays <- read.csv("holiday_13_15.csv", header=TRUE)

dta <- merge(dta, holidays, by = "date", all = TRUE)
dta$holiday[is.na(dta$holiday)] <- 0
dta <- dta[order(dta$ecounty, dta$newdate),]


######PICK ONE COUNTY (CAIDIAN)######
dta.c <- dta[which(dta$ecounty=="Caidian"),]

cb1.pm <- crossbasis(dta.c$pm25b, lag=30, argvar=list(fun="lin"),
                     arglag=list(df=3))
cb2.no2 <- crossbasis(dta.c$no2, lag=30, argvar=list(fun="lin"),
                      arglag=list(df=3))
cb1.temp <- crossbasis(dta.c$tmean, lag=30, argvar=list(df=3),
                       arglag=list(df=3))
cb1.rh <- crossbasis(dta.c$rh, lag=30, argvar=list(df=3),
                     arglag=list(df=3))
all <- glm(all ~ offset(log(pop)) + cb2.no2 + cb1.pm + cb1.temp + cb1.rh + 
             ns(as.numeric(newdate), df = 12) + as.factor(dow) + as.factor(holiday), data = dta.c, family = quasipoisson())
pred.all <- crosspred(cb2.no2, all, at=0:330, bylag=0.2, cumul=TRUE)


##Lag-Specific Plot
dtaX <- data.frame(matrix(NA, length(pred.all$matfit[10,]), 5))
names(dtaX) <- c("fit", "se", "pe", "lci", "uci")
dtaX$fit <- pred.all$matfit[2,]*10
dtaX$se  <- pred.all$matse[2,]*10
dtaX$pe  <- 100*(exp(dtaX$fit)-1)
dtaX$lci <- 100*(exp(dtaX$fit-1.96*dtaX$se)-1)
dtaX$uci <- 100*(exp(dtaX$fit+1.96*dtaX$se)-1)
rownames(dtaX) <- colnames(pred.all$matfit)
par(cex=1.3, las=1)
plot(c(1,dim(dtaX)[1]), c(min(dtaX$lci), max(dtaX$uci)), type="n",xaxt="n",
     xlab="Lag", ylab="% Change", main="Lag Specific Plot from Matfit, Caidian (Non-Accidental Mortality)")
axis(1, at = seq(0,200,25), labels = seq(0,40,5))
polygon(c(rev(1:dim(dtaX)[1]), 1:dim(dtaX)[1]), c(rev(dtaX$lci), dtaX$uci),
        col = 'grey80', border = NA, lty=0)
lines(1:dim(dtaX)[1], dtaX$pe, lwd=2, col="black")
abline(h=0, lty=2)


##Cumulative Plot
dtaX <- data.frame(matrix(NA, length(pred.all$cumfit[10,]), 5))
names(dtaX) <- c("fit", "se", "pe", "lci", "uci")
dtaX$fit <- pred.all$cumfit[2,]*10
dtaX$se  <- pred.all$cumse[2,]*10
dtaX$pe  <- 100*(exp(dtaX$fit)-1)
dtaX$lci <- 100*(exp(dtaX$fit-1.96*dtaX$se)-1)
dtaX$uci <- 100*(exp(dtaX$fit+1.96*dtaX$se)-1)
rownames(dtaX) <- colnames(pred.all$cumfit)
plot(c(1,dim(dtaX)[1]), c(min(dtaX$lci), max(dtaX$uci)), type="n",xaxt="n",
     xlab="Lag", ylab="% Change")
par(cex=1.3, las=1)
axis(1, at = seq(0,30,5))
polygon(c(rev(1:dim(dtaX)[1]), 1:dim(dtaX)[1]), c(rev(dtaX$lci), dtaX$uci),
        col = 'grey80', border = NA, lty=0)
lines(1:dim(dtaX)[1], dtaX$pe, lwd=2, col="black")
abline(h=0, lty=2)


## Read in coeff and vcov Matrix
setwd("D:/Users/profu/Documents/Schoolwork/PhD/Research Projects/no2_42_counties/meta_analysis_holidays/final_model_3df_meta/all/coeff")
coeff <- read.csv("Caidian_coeff.csv", header=TRUE, stringsAsFactors = FALSE)
coeff <- as.matrix(coeff$x)

setwd("D:/Users/profu/Documents/Schoolwork/PhD/Research Projects/no2_42_counties/meta_analysis_holidays/final_model_3df_meta/all/vcov")
vcov <- read.csv("Caidian_vcov.csv", header=TRUE, stringsAsFactors = FALSE)
vcov <- cbind(vcov$cb2.no2v1.l1, vcov$cb2.no2v1.l2, vcov$cb2.no2v1.l3)

B <- ns(0:30, df=3, intercept=TRUE)
DLp     <- B%*%coeff
DLvcovp <- B%*%vcov%*%t(B)
DLsep   <- sqrt(diag(DLvcovp))
lcp <- DLp - 1.96*DLsep
ucp <- DLp + 1.96*DLsep

DLp1 <- 100*(exp(DLp*10)-1)
lcp1 <- 100*(exp((DLp - 1.96* DLsep)*10)-1)
ucp1 <- 100*(exp((DLp + 1.96* DLsep)*10)-1)


## Lag-Specific Plot (Reconstructed)
par(cex=1.3, las=1)
plot(c(0,30), c(min(lcp1), max(ucp1)), type="n",xaxt="n",
     xlab="Lag", ylab="% Change", main="Reconstructed Lag-Specific Plot, Caidian (Non-Accidental Mortality)")
axis(1, at = seq(0,30,5))
polygon(c(rev(0:30), 0:30), c((lcp1[31:1]), ucp1[1:31]), 
        col="gray90", border = NA, lty=0)
lines(0:30, DLp1, lty=1, lwd=2, col="black")
abline(h=0, lty=2)


## Loop for Cumulative Values
cumulative <- matrix(nrow = 31, ncol = 2)
cumulative[1,1] <- DLp[1,1]
cumulative[1,2] <- DLvcovp[1,1]

for (c in 2:dim(DLvcovp)[1]){
  
  cumulative[c,1] <- sum(DLp[1:c])
  cumulative[c,2] <- sum(DLvcovp[1:c, 1:c])
} 

DLcum <- as.matrix(cumulative[,1])
DLsecum <- sqrt(as.matrix(cumulative[,2]))

DLcum1 <- 100*(exp(DLcum*10)-1)
lccum1 <- 100*(exp((DLcum - 1.96* DLsecum)*10)-1)
uccum1 <- 100*(exp((DLcum + 1.96* DLsecum)*10)-1)


## Cumulative Plot (Reconstructed)
par(cex=1.3, las=1)
plot(c(0,30), c(min(lccum1), max(uccum1)), type="n",xaxt="n",
     xlab="Lag", ylab="% Change", main="Reconstructed Cumulative Plot, Caidian (Non-Accidental Mortality)")
axis(1, at = seq(0,30,5))
polygon(c(rev(0:30), 0:30), c((lccum1[31:1]), uccum1[1:31]), 
        col="gray90", border = NA, lty=0)
lines(0:30, DLcum1, lty=1, lwd=2, col="black")
abline(h=0, lty=2)



######META-ANALYSIS ALL CAUSES######
setwd("D:/Users/profu/Documents/Schoolwork/PhD/Research Projects/no2_42_counties/meta_analysis_holidays/final_model_3df_meta/all")

#coeff
coef.files <- (Sys.glob("./coeff/*.csv"))
coeffs <- lapply(coef.files, function(x) read.table(x, header = TRUE)) 
length(coeffs)
ourY <- matrix(NA, length(coeffs), nrow(coeffs[[1]]))
for (i in 1:dim(ourY)[1]){ ## i = 1
  ourY[i,] <- coeffs[[i]]$x
}

#vcov
vcov.files <- (Sys.glob("./vcov/*.csv"))
vcovs <- lapply(vcov.files, function(x) read.csv(x, header = TRUE)) 
length(vcovs)
ourV <- array(NA, dim = c(nrow(vcovs[[1]]), ncol(vcovs[[1]]), length(vcovs)))
for (i in 1:dim(ourV)[3]){ ## i = 1
  ourV[,,i] <- as.matrix(vcovs[[i]])
}
ourV <- lapply(seq(dim(ourV)[3]), function(x) ourV[ , , x])

metaAll <- mvmeta(ourY~1,ourV,method="reml")
summary(metaAll)

coeff <- t(as.matrix(metaAll$coefficients))
vcov <- as.matrix(metaAll$vcov)

#Lag-Specific
B <- ns(0:30, df=3, intercept=TRUE)
DLp     <- B%*%coeff
DLvcovp <- B%*%vcov%*%t(B)
DLsep   <- sqrt(diag(DLvcovp))
lcp <- DLp - 1.96*DLsep
ucp <- DLp + 1.96*DLsep

DLp1 <- 100*(exp(DLp*10)-1)
lcp1 <- 100*(exp((DLp - 1.96* DLsep)*10)-1)
ucp1 <- 100*(exp((DLp + 1.96* DLsep)*10)-1)

# (For sex-stratified only)
DLp2 <- 100*(exp(DLp*10)-1)
lcp2 <- 100*(exp((DLp - 1.96* DLsep)*10)-1)
ucp2 <- 100*(exp((DLp + 1.96* DLsep)*10)-1)

#Cumulative
cumulative <- matrix(nrow = 31, ncol = 2)
cumulative[1,1] <- DLp[1,1]
cumulative[1,2] <- DLvcovp[1,1]

for (c in 2:dim(DLvcovp)[1]){
  
  cumulative[c,1] <- sum(DLp[1:c])
  cumulative[c,2] <- sum(DLvcovp[1:c, 1:c])
} 

DLcum <- as.matrix(cumulative[,1])
DLsecum <- sqrt(as.matrix(cumulative[,2]))

DLcum1 <- 100*(exp(DLcum*10)-1)
lccum1 <- 100*(exp((DLcum - 1.96* DLsecum)*10)-1)
uccum1 <- 100*(exp((DLcum + 1.96* DLsecum)*10)-1)

# (For sex-stratified only)
DLcum2 <- 100*(exp(DLcum*10)-1)
lccum2 <- 100*(exp((DLcum - 1.96* DLsecum)*10)-1)
uccum2 <- 100*(exp((DLcum + 1.96* DLsecum)*10)-1)

View(cbind(DLp1, lcp1, ucp1))
View(cbind(DLcum1, lccum1, uccum1))

View(cbind(DLp2, lcp2, ucp2))
View(cbind(DLcum2, lccum2, uccum2))

##Plots
par(cex=1.3, las=1)
plot(c(0,30), c(min(lcp1), max(ucp1)), type="n",xaxt="n",
     xlab="Lag", ylab="% Change", ylim=c(-0.3, 0.33))
axis(1, at = seq(0,30,5))
polygon(c(rev(0:30), 0:30), c((lcp1[31:1]), ucp1[1:31]), 
        col="gray90", border = NA, lty=0)
lines(0:30, DLp1, lty=1, lwd=2, col="black")
abline(h=0, lty=2)
legend("topleft", "(c)", bty = "n")

par(cex=1.3, las=1)
plot(c(0,30), c(min(lccum1), max(uccum1)), type="n",xaxt="n",
     xlab="Lag", ylab="% Change", ylim=c(-3.2, 5.5))
axis(1, at = seq(0,30,5))
polygon(c(rev(0:30), 0:30), c((lccum1[31:1]), uccum1[1:31]), 
        col="gray90", border = NA, lty=0)
lines(0:30, DLcum1, lty=1, lwd=2, col="black")
abline(h=0, lty=2)
legend("topleft", "(c)", bty = "n")


## Plots Sex Stratified
par(cex=1.3, las=1)
plot(c(0,30), c(min(lcp1), max(lcp1)), type="n",xaxt="n",
     xlab="Lag", ylab="% Change",
     ylim=c(-0.5, 0.6))
axis(1, at = seq(0,30,5))
polygon(c(rev(0:30), 0:30), c(rev(lcp1), ucp1),
        col = 'grey80', border = NA, lty=0)
polygon(c(rev(0:30), 0:30), c(rev(lcp2), ucp2),
        col = 'grey40', border = NA, density = 25)
lines(0:30, DLp1, lty=1, lwd=2, col="black")
lines(0:30, DLp2, lty=2, lwd=2, col="black")
abline(h=0, lty=2)
legend("topleft", "(b)", bty = "n")
#legend("topright", inset = 0.02, legend=c("Men", "Women"), lty=1:2)

par(cex=1.3, las=1)
plot(c(0,30), c(min(lccum1), max(lccum1)), type="n",xaxt="n",
     xlab="Lag", ylab="% Change",
     ylim=c(-4.1, 8.2))
axis(1, at = seq(0,30,5))
polygon(c(rev(0:30), 0:30), c(rev(lccum1), uccum1),
        col = 'grey80', border = NA, lty=0)
polygon(c(rev(0:30), 0:30), c(rev(lccum2), uccum2),
        col = 'grey40', border = NA, density = 25)
lines(0:30, DLcum1, lty=1, lwd=2, col="black")
lines(0:30, DLcum2, lty=2, lwd=2, col="black")
abline(h=0, lty=2)
legend("topleft", "(a)", bty = "n")
legend("topright", inset = 0.02, legend=c("Men", "Women"), lty=1:2)


######LAG0-1, LAG0-7, LAG0-30 FOREST PLOTS######
setwd("D:/Users/profu/Documents/Schoolwork/PhD/Research Projects/no2_42_counties/environmental_pollution_submission/round_2_revision/forest_plot")
no2all <- matrix(nrow = 42, ncol = 7)

my.counties <- unique(dta$ecounty)

for (c in 1:length(my.counties)){
  
  dta.c <- dta[which(dta$ecounty == my.counties[c]),]
  cb1.pm <- crossbasis(dta.c$pm25b, lag=30, argvar=list(fun="lin"),
                       arglag=list(df=3))
  cb2.no2 <- crossbasis(dta.c$no2, lag=30, argvar=list(fun="lin"),
                        arglag=list(df=3))
  cb1.temp <- crossbasis(dta.c$tmean, lag=30, argvar=list(df=3),
                         arglag=list(df=3))
  cb1.rh <- crossbasis(dta.c$rh, lag=30, argvar=list(df=3),
                       arglag=list(df=3))
  all <- glm(all ~ offset(log(pop)) + cb2.no2 + cb1.pm + cb1.temp + cb1.rh + 
               ns(as.numeric(newdate), df = 12) + as.factor(dow) + as.factor(holiday), 
             data = dta.c, family = quasipoisson())
  pred.all <- crosspred(cb2.no2, all, at=0:330, bylag=0.2, cumul=TRUE)
  dtaX <- data.frame(matrix(NA, length(pred.all$matfit[10,]), 5))
  
  names(dtaX) <- c("fit", "se", "pe", "lci", "uci")
  dtaX <- data.frame(matrix(NA, length(pred.all$cumfit[10,]), 5))
  names(dtaX) <- c("fit", "se", "pe", "lci", "uci")
  dtaX$fit <- pred.all$cumfit[2,]*10
  dtaX$se  <- pred.all$cumse[2,]*10
  dtaX$pe  <- 100*(exp(dtaX$fit)-1)
  dtaX$lci <- 100*(exp(dtaX$fit-1.96*dtaX$se)-1)
  dtaX$uci <- 100*(exp(dtaX$fit+1.96*dtaX$se)-1)
  rownames(dtaX) <- colnames(pred.all$cumfit)
  
  no2all[c,1] <- unique(dta.c$ecounty)
  no2all[c,2] <- unique(dta.c$mcode)
  no2all[c,3] <- dtaX$pe[31]
  no2all[c,4] <- dtaX$lci[31]
  no2all[c,5] <- dtaX$uci[31]
  no2all[c,6] <- dtaX$fit[31]
  no2all[c,7] <- dtaX$se[31]
  
} 

no2all <- as.data.frame(no2all)
no2all2 <- rename(no2all, c("V1"="county", "V2"="mcode", "V3"="pe", "V4"="lci", "V5"="uci", "V6"="coeff", "V7"="se"))
write.csv(no2all2, "D:/Users/profu/Documents/Schoolwork/PhD/Research Projects/no2_42_counties/environmental_pollution_submission/round_2_revision/forest_plot/no2all_lag030.csv")

## Meta-Analysis and Forest Plot
library(forestplot)
library(mvmeta)
setwd("D:/Users/profu/Documents/Schoolwork/PhD/Research Projects/no2_42_counties/environmental_pollution_submission/round_2_revision/forest_plot")

#lag01
all <- read.csv("no2all_lag01.csv", header=TRUE, stringsAsFactors = FALSE)
model.all <- mvmeta(I(all$coeff)~1, S=(all$se)^2, method="reml", na.action="na.omit")
summary(model.all)
100*(exp(model.all$coefficients[1])-1)
100*(exp((summary(model.all)$coef[1] - 1.96* summary(model.all)$coef[2]))-1)
100*(exp((summary(model.all)$coef[1] + 1.96* summary(model.all)$coef[2]))-1)

forestplot(labeltext = c("County", all$county, NA, "Summary"), mean = c(NA, all$pe, NA, 0.24), 
           hrzl_lines = gpar(col="#444444"), lower = c(NA, all$lci, NA, 0.05),
           upper = c(NA, all$uci, NA, 0.43), is.summary=c(TRUE, rep(FALSE, 42), TRUE), 
           col=fpColors(box="black", line="black", summary = "black"), vertices=T)

#lag07
all <- read.csv("no2all_lag07.csv", header=TRUE, stringsAsFactors = FALSE)
model.all <- mvmeta(I(all$coeff)~1, S=(all$se)^2, method="reml", na.action="na.omit")
summary(model.all)
100*(exp(model.all$coefficients[1])-1)
100*(exp((summary(model.all)$coef[1] - 1.96* summary(model.all)$coef[2]))-1)
100*(exp((summary(model.all)$coef[1] + 1.96* summary(model.all)$coef[2]))-1)

forestplot(labeltext = c("County", all$county, NA, "Summary"), mean = c(NA, all$pe, NA, 0.57), 
           hrzl_lines = gpar(col="#444444"), lower = c(NA, all$lci, NA, -0.04),
           upper = c(NA, all$uci, NA, 1.18), is.summary=c(TRUE, rep(FALSE, 42), TRUE), 
           col=fpColors(box="black", line="black", summary = "black"), vertices=T, clip = c(-10, 15))

#lag030
all <- read.csv("no2all_lag030.csv", header=TRUE, stringsAsFactors = FALSE)
model.all <- mvmeta(I(all$coeff)~1, S=(all$se)^2, method="reml", na.action="na.omit")
summary(model.all)
100*(exp(model.all$coefficients[1])-1)
100*(exp((summary(model.all)$coef[1] - 1.96* summary(model.all)$coef[2]))-1)
100*(exp((summary(model.all)$coef[1] + 1.96* summary(model.all)$coef[2]))-1)

forestplot(labeltext = c("County", all$county, NA, "Summary"), mean = c(NA, all$pe, NA, -0.14), 
           hrzl_lines = gpar(col="#444444"), lower = c(NA, all$lci, NA, -1.63),
           upper = c(NA, all$uci, NA, 1.37), is.summary=c(TRUE, rep(FALSE, 42), TRUE), 
           col=fpColors(box="black", line="black", summary = "black"), vertices=T, clip = c(-20,30))

#Plotting (700 x 750)
pushViewport(viewport(layout = grid.layout(1, 2)))
pushViewport(viewport(layout.pos.col = 1))
forestplot(labeltext = c("County", all$county, NA, "Summary"), mean = c(NA, all$pe, NA, 0.24), 
           hrzl_lines = gpar(col="#444444"), lower = c(NA, all$lci, NA, 0.05),
           upper = c(NA, all$uci, NA, 0.43), is.summary=c(TRUE, rep(FALSE, 42), TRUE), 
           col=fpColors(box="black", line="black", summary = "black"), vertices=T, new_page = FALSE)

pushViewport(viewport(layout = grid.layout(1, 2)))
pushViewport(viewport(layout.pos.col = 1))
forestplot(labeltext = c("County", all$county, NA, "Summary"), mean = c(NA, all$pe, NA, 0.57), 
           hrzl_lines = gpar(col="#444444"), lower = c(NA, all$lci, NA, -0.04),
           upper = c(NA, all$uci, NA, 1.18), is.summary=c(TRUE, rep(FALSE, 42), TRUE), 
           col=fpColors(box="black", line="black", summary = "black"), vertices=T, new_page = FALSE, clip = c(-10, 15))

pushViewport(viewport(layout = grid.layout(1, 2)))
pushViewport(viewport(layout.pos.col = 1))
forestplot(labeltext = c("County", all$county, NA, "Summary"), mean = c(NA, all$pe, NA, -0.14), 
           hrzl_lines = gpar(col="#444444"), lower = c(NA, all$lci, NA, -1.63),
           upper = c(NA, all$uci, NA, 1.37), is.summary=c(TRUE, rep(FALSE, 42), TRUE), 
           col=fpColors(box="black", line="black", summary = "black"), vertices=T, new_page = FALSE, clip = c(-20,25))
