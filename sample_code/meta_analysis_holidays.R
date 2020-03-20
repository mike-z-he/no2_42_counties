## Meta-Analysis for NO2-Mortality

setwd("D:/Users/profu/Documents/Schoolwork/PhD/3rd Year Summer/ISEE Analysis")
options(mc.cores=parallel::detectCores())

##########DATA CLEANING##########
dta <- read.csv("databydate-PM2.5update0404.csv", header=TRUE, stringsAsFactors = FALSE)
str(dta)
summary(dta)
names(dta)

library(mgcv)
library(splines)
library(dplyr)
library(MuMIn)
library(gplots)
library(reshape)

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


########ALL-CAUSE MORTALITY DLM MODEL#########
library(dlnm)
library(splines)

## Crossbasis Matrix
cb1.pm <- crossbasis(dta$pm25b, lag=30, argvar=list(fun="lin"),
                     arglag=list(df=3), group = dta$ecounty)
cb2.no2 <- crossbasis(dta$no2, lag=30, argvar=list(fun="lin"),
                      arglag=list(df=3), group = dta$ecounty)
cb1.temp <- crossbasis(dta$tmean, lag=30, argvar=list(df=3),
                       arglag=list(df=3), group = dta$ecounty)
cb1.rh <- crossbasis(dta$rh, lag=30, argvar=list(df=3),
                     arglag=list(df=3), group = dta$ecounty)

all <- glm(all ~ offset(log(pop)) + cb2.no2 + cb1.pm + cb1.temp + cb1.rh + as.factor(mcode) + 
             ns(as.numeric(newdate), df = 12) + as.factor(dow) + as.factor(holiday), data = dta, family = quasipoisson())
pred.all <- crosspred(cb2.no2, all, at=0:330, bylag=0.2, cumul=TRUE)


######SUBSET MODEL (ONLY COUNTIES WITH MONITORS)######

dta2 <- dta[which(dta$ecounty=="Pukou" | dta$ecounty=="Zhangjiagang" | dta$ecounty=="Caidian"
                  | dta$ecounty=="Jiang'an" | dta$ecounty=="Daowai" | dta$ecounty=="Panyu"
                  | dta$ecounty=="Liwan" | dta$ecounty=="Yuexiu" | dta$ecounty=="Lianhu"
                  | dta$ecounty=="Lintong" | dta$ecounty=="Changping" | dta$ecounty=="Dongcheng"
                  | dta$ecounty=="Fengtai" | dta$ecounty=="Mentougou" | dta$ecounty=="Miyun"
                  | dta$ecounty=="Tongzhou" | dta$ecounty=="Yanqing" | dta$ecounty=="Xinghualing"),]

## Further removes the rural counties (Changping, Mentougou, Miyun, Yanqing)
dta2 <- dta[which(dta$ecounty=="Pukou" | dta$ecounty=="Zhangjiagang" | dta$ecounty=="Caidian"
                  | dta$ecounty=="Jiang'an" | dta$ecounty=="Daowai" | dta$ecounty=="Panyu"
                  | dta$ecounty=="Liwan" | dta$ecounty=="Yuexiu" | dta$ecounty=="Lianhu"
                  | dta$ecounty=="Lintong" | dta$ecounty=="Dongcheng" | dta$ecounty=="Fengtai" 
                  | dta$ecounty=="Tongzhou" | dta$ecounty=="Xinghualing"),]

## Subset model that removes those with a lot of missing data
dta2 <- dta[which(dta$ecounty!="Changning" & dta$ecounty != "Fengxian" & dta$ecounty != "Huangpu"
                  & dta$ecounty != "Jinshan" & dta$ecounty != "Lintong" & dta$ecounty != "Minhang"
                  & dta$ecounty != "Putuo" & dta$ecounty != "Qiaokou" & dta$ecounty != "Songjiang"
                  & dta$ecounty != "Xinghualing"),]

dta2 <- dta2[order(dta2$ecounty, dta2$newdate),]

cb1.pm <- crossbasis(dta2$pm25b, lag=30, argvar=list(fun="lin"),
                     arglag=list(df=3), group = dta2$ecounty)
cb2.no2 <- crossbasis(dta2$no2, lag=30, argvar=list(fun="lin"),
                      arglag=list(df=3), group = dta2$ecounty)
cb1.temp <- crossbasis(dta2$tmean, lag=30, argvar=list(df=3),
                       arglag=list(df=3), group = dta2$ecounty)
cb1.rh <- crossbasis(dta2$rh, lag=30, argvar=list(df=3),
                     arglag=list(df=3), group = dta2$ecounty)

all_subset <- glm(all ~ offset(log(pop)) + cb2.no2 + cb1.pm + cb1.temp + cb1.rh + as.factor(mcode) + 
             ns(as.numeric(newdate), df = 12) + as.factor(dow) + as.factor(holiday), data = dta2, family = quasipoisson())
pred.all_subset <- crosspred(cb2.no2, all_subset, at=0:330, bylag=0.2, cumul=TRUE)


######COMPARING FULL VS SUBSET MODELS######
dtaX <- data.frame(matrix(NA, length(pred.all$matfit[10,]), 5))
names(dtaX) <- c("fit", "se", "pe", "lci", "uci")
dtaX$fit <- pred.all$matfit[2,]*10
dtaX$se  <- pred.all$matse[2,]*10
dtaX$pe  <- 100*(exp(dtaX$fit)-1)
dtaX$lci <- 100*(exp(dtaX$fit-1.96*dtaX$se)-1)
dtaX$uci <- 100*(exp(dtaX$fit+1.96*dtaX$se)-1)
rownames(dtaX) <- colnames(pred.all$matfit)
dtaY <- data.frame(matrix(NA, length(pred.all_subset$matfit[10,]), 5))
names(dtaY) <- c("fit", "se", "pe", "lci", "uci")
dtaY$fit <- pred.all_subset$matfit[2,]*10
dtaY$se  <- pred.all_subset$matse[2,]*10
dtaY$pe  <- 100*(exp(dtaY$fit)-1)
dtaY$lci <- 100*(exp(dtaY$fit-1.96*dtaY$se)-1)
dtaY$uci <- 100*(exp(dtaY$fit+1.96*dtaY$se)-1)
rownames(dtaY) <- colnames(pred.all_subset$matfit)
par(cex=1.3, las=1)
plot(c(1,dim(dtaX)[1]), c(min(dtaX$lci), max(dtaX$uci)), type="n",xaxt="n",
     xlab="Lag", ylab="% Change",
     ylim=c(-0.15, 0.35))
axis(1, at = seq(0,200,25), labels = seq(0,40,5))
polygon(c(rev(1:dim(dtaX)[1]), 1:dim(dtaX)[1]), c(rev(dtaX$lci), dtaX$uci),
        col = 'grey80', border = NA, lty=0)
polygon(c(rev(1:dim(dtaY)[1]), 1:dim(dtaY)[1]), c(rev(dtaY$lci), dtaY$uci),
        col = 'grey40', border = NA, density = 25)
lines(1:dim(dtaX)[1], dtaX$pe, lty=1, lwd=2, col="black")
lines(1:dim(dtaY)[1], dtaY$pe, lty=2, lwd=2, col="black")
abline(h=0, lty=2)
legend("topright", inset = 0.02, legend=c("Full Analysis", "Subset w/ Least Missing Data"), lty=1:2)

######COUNTY-SPECIFIC MODELS######
lin.coeff <- summary(all)$coefficients[2:4,1]
lin.vcov <- vcov(all)[2:4,2:4]

setwd("D:/Users/profu/Documents/Schoolwork/PhD/Research Projects/no2_42_counties/meta_analysis_holidays/final_model_3df_meta")

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
  all <- glm(resg2 ~ offset(log(pop)) + cb2.no2 + cb1.pm + cb1.temp + cb1.rh + 
               ns(as.numeric(newdate), df = 12) + as.factor(dow) + as.factor(holiday), data = dta.c, family = quasipoisson())
  lin.coeff <- summary(all)$coefficients[2:4,1]
  lin.vcov <- vcov(all)[2:4,2:4]
  #nm1 <- paste(my.counties[c], "_coeff.csv", sep="")
  nm2 <- paste(my.counties[c], "_vcov.csv", sep="")
  #write.csv(lin.coeff, file = nm1, row.names = FALSE)
  write.csv(lin.vcov, file = nm2, row.names = FALSE)
}


######COUNTY-SPECIFIC MATFIT/CUMFIT OUTPUTS######
setwd("D:/Users/profu/Documents/Schoolwork/PhD/Research Projects/no2_42_counties/meta_analysis_holidays")

#matfit
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
               ns(as.numeric(newdate), df = 12) + as.factor(dow) + as.factor(holiday), data = dta.c, family = quasipoisson())
  pred.all <- crosspred(cb2.no2, all, at=0:330, bylag=0.2, cumul=TRUE)
  
  dtaX <- data.frame(matrix(NA, length(pred.all$matfit[10,]), 5))
  names(dtaX) <- c("fit", "se", "pe", "lci", "uci")
  dtaX$fit <- pred.all$matfit[2,]*10
  dtaX$se  <- pred.all$matse[2,]*10
  dtaX$pe  <- 100*(exp(dtaX$fit)-1)
  dtaX$lci <- 100*(exp(dtaX$fit-1.96*dtaX$se)-1)
  dtaX$uci <- 100*(exp(dtaX$fit+1.96*dtaX$se)-1)
  rownames(dtaX) <- colnames(pred.all$matfit)
  nm1 <- paste(my.counties[c], "_matfit.csv", sep="")
  write.csv(dtaX, file = nm1, row.names = FALSE)
}

#cumfit
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
               ns(as.numeric(newdate), df = 12) + as.factor(dow) + as.factor(holiday), data = dta.c, family = quasipoisson())
  pred.all <- crosspred(cb2.no2, all, at=0:330, bylag=0.2, cumul=TRUE)
  
  dtaX <- data.frame(matrix(NA, length(pred.all$cumfit[10,]), 5))
  names(dtaX) <- c("fit", "se", "pe", "lci", "uci")
  dtaX$fit <- pred.all$cumfit[2,]*10
  dtaX$se  <- pred.all$cumse[2,]*10
  dtaX$pe  <- 100*(exp(dtaX$fit)-1)
  dtaX$lci <- 100*(exp(dtaX$fit-1.96*dtaX$se)-1)
  dtaX$uci <- 100*(exp(dtaX$fit+1.96*dtaX$se)-1)
  rownames(dtaX) <- colnames(pred.all$cumfit)
  nm1 <- paste(my.counties[c], "_cumfit.csv", sep="")
  write.csv(dtaX, file = nm1, row.names = FALSE)
}

######FOREST PLOT######
## Extracting lag0-1 % changes for all counties
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
  no2all[c,3] <- dtaX$pe[1]
  no2all[c,4] <- dtaX$lci[1]
  no2all[c,5] <- dtaX$uci[1]
  no2all[c,6] <- dtaX$fit[1]
  no2all[c,7] <- dtaX$se[1]
  
} 

no2all <- as.data.frame(no2all)
no2all2 <- rename(no2all, c("V1"="county", "V2"="mcode", "V3"="pe", "V4"="lci", "V5"="uci", "V6"="coeff", "V7"="se"))
write.csv(no2all2, "D:/Users/profu/Documents/Schoolwork/PhD/Research Projects/no2_42_counties/environmental_pollution_submission/round_2_revision/forest_plot/no2all_lag0.csv")

## Meta-Analysis and Forest Plot
library(forestplot)
library(mvmeta)
setwd("D:/Users/profu/Documents/Schoolwork/PhD/Research Projects/no2_42_counties/environmental_pollution_submission/round_2_revision/forest_plot")

#lag0
all <- read.csv("no2all_lag0.csv", header=TRUE, stringsAsFactors = FALSE)
model.all <- mvmeta(I(all$coeff)~1, S=(all$se)^2, method="reml", na.action="na.omit")
summary(model.all)
100*(exp(model.all$coefficients[1])-1)
100*(exp((summary(model.all)$coef[1] - 1.96* summary(model.all)$coef[2]))-1)
100*(exp((summary(model.all)$coef[1] + 1.96* summary(model.all)$coef[2]))-1)

forestplot(labeltext = c("County", all$county, NA, "Summary"), mean = c(NA, all$pe, NA, 0.129), 
           hrzl_lines = gpar(col="#444444"), lower = c(NA, all$lci, NA, 0.026),
           upper = c(NA, all$uci, NA, 0.231), is.summary=c(TRUE, rep(FALSE, 42), TRUE), 
           col=fpColors(box="royalblue", line="darkblue", summary = "royalblue"), vertices=T)

#lag01
all <- read.csv("no2all_lag01.csv", header=TRUE, stringsAsFactors = FALSE)
model.all <- mvmeta(I(all$coeff)~1, S=(all$se)^2, method="reml", na.action="na.omit")
summary(model.all)
100*(exp(model.all$coefficients[1])-1)
100*(exp((summary(model.all)$coef[1] - 1.96* summary(model.all)$coef[2]))-1)
100*(exp((summary(model.all)$coef[1] + 1.96* summary(model.all)$coef[2]))-1)

forestplot(labeltext = c("County", all$county, NA, "Summary"), mean = c(NA, all$pe, NA, 0.239), 
           hrzl_lines = gpar(col="#444444"), lower = c(NA, all$lci, NA, 0.041),
           upper = c(NA, all$uci, NA, 0.436), is.summary=c(TRUE, rep(FALSE, 42), TRUE), 
           col=fpColors(box="royalblue", line="darkblue", summary = "royalblue"), vertices=T)

######POOLED ANALYSIS######
setwd("D:/Users/profu/Documents/Schoolwork/PhD/Research Projects/no2_42_counties/meta_analysis_holidays")

library(MASS)
library(devtools) 
library(tlnise)
library(ggplot2)
library(splines)


coef.files <- (Sys.glob("./coeff/*.csv"))

coeffs <- lapply(coef.files, function(x) read.table(x, header = TRUE)) 
length(coeffs)

ourY <- matrix(NA, length(coeffs), nrow(coeffs[[1]]))

for (i in 1:dim(ourY)[1]){ ## i = 1
  
  ourY[i,] <- coeffs[[i]]$x
  
}

vcov.files <- (Sys.glob("./vcov/*.csv"))

vcovs <- lapply(vcov.files, function(x) read.csv(x, header = TRUE)) 
length(vcovs)

ourV <- array(NA, dim = c(nrow(vcovs[[1]]), ncol(vcovs[[1]]), length(vcovs)))

for (i in 1:dim(ourV)[3]){ ## i = 1
  
  ourV[,,i] <- as.matrix(vcovs[[i]])
  
}


### vcov_test ###
vcov.files <- (Sys.glob("./vcov_test/*.csv"))

vcovs <- lapply(vcov.files, function(x) read.csv(x, header = TRUE)) 
length(vcovs)

ourV <- array(NA, dim = c(nrow(vcovs[[1]]), ncol(vcovs[[1]]), length(vcovs)))

for (i in 1:dim(ourV)[3]){ ## i = 1
  
  ourV[,,i] <- as.matrix(vcovs[[i]])
  
}

###

pooled <- tlnise(Y=ourY, V=ourV, intercept=FALSE, seed=68460, maxiter=10000)

B <- ns(0:30, df=3, intercept=TRUE)

DLp     <- B%*%pooled$gamma[,1]
DLvcovp <- B%*%pooled$Dgamma%*%t(B)
DLsep   <- sqrt(diag(DLvcovp))

lcp <- DLp - 1.96*DLsep
ucp <- DLp + 1.96*DLsep


oneAnal.coef <- read.csv("lin.coeff.csv", header = TRUE)
oneAnal.vcov <- read.csv("lin.vcov.csv", header = TRUE)

DL1     <- B%*%oneAnal.coef[,2]
DLvcov1 <- B%*%as.matrix(oneAnal.vcov[,2:4])%*%t(B)
DLse1   <- sqrt(diag(DLvcov1))

lc1 <- DL1 - 1.96*DLse1
uc1 <- DL1 + 1.96*DLse1


######PLOT COMPARISON######

par(cex=1.3, las=1)
plot(c(0,30), c(min(lcp), max(ucp)), type="n",xaxt="n",
     xlab="Lag", ylab=expression(beta), main="Pooled vs. Fixed Effects (Non-Accidental Mortality)", ylim=c(-0.00015, 0.0004))

axis(1, at = seq(0,30,5))

polygon(c(rev(0:30), 0:30), c((lcp[31:1]), ucp[1:31]), 
         col="gray90", border = NA, lty=0)

polygon(c(rev(0:30), 0:30), c((lc1[nrow(lc1):1]), uc1[1:nrow(uc1)]), 
         col="gray40", border = NA, density=25)

lines(0:30, DLp, lty=1, lwd=2, col="black")
lines(0:30, DL1, lty=2, lwd=2, col="black")

abline(h=0, lty=2)

legend("topleft", inset = 0.02, legend=c("Meta-Analysis", "Fixed Effects"), lty=1:2)


######CUMULATIVE PLOT######
setwd("D:/Users/profu/Documents/Schoolwork/PhD/Research Projects/no2_42_counties/meta_analysis_holidays")
dta <- read.csv("DLp.csv", header=TRUE, stringsAsFactors = FALSE)

par(cex=1.3, las=1)
plot(c(0,30), c(min(lcp), max(ucp)), type="n",xaxt="n",
     xlab="Lag", ylab=expression(beta), main="Cumulative Pooled vs. Fixed Effects (Non-Accidental Mortality)", ylim=c(-0.0002, 0.0014))

axis(1, at = seq(0,30,5))
lines(0:30, dta$cumDLp, lty=1, lwd=2, col="black")
lines(0:30, dta$cumDL1, lty=2, lwd=2, col="black")
abline(h=0, lty=2)
legend("topleft", inset = 0.02, legend=c("Meta-Analysis", "Fixed Effects"), lty=1:2)



######META-ANALYSIS USING MVMETA######
setwd("D:/Users/profu/Documents/Schoolwork/PhD/Research Projects/no2_42_counties/meta_analysis_holidays/final_model_3df_meta/all")

library(mvmeta)

coef.files <- (Sys.glob("./coeff/*.csv"))

coeffs <- lapply(coef.files, function(x) read.table(x, header = TRUE)) 
length(coeffs)

ourY <- matrix(NA, length(coeffs), nrow(coeffs[[1]]))

for (i in 1:dim(ourY)[1]){ ## i = 1
  
  ourY[i,] <- coeffs[[i]]$x
  
}

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

B <- ns(0:30, df=3, intercept=TRUE)

DLp1     <- B%*%t(as.matrix(metaAll$coefficients))
DLvcovp1 <- B%*%metaAll$vcov%*%t(B)
DLsep1   <- sqrt(diag(DLvcovp1))

#lcp1 <- DLp1 - 1.96*DLsep1
#ucp1 <- DLp1 + 1.96*DLsep1


lcp1 <- 100*(exp((DLp1 - 1.96* DLsep1)*10)-1)
ucp1 <- 100*(exp((DLp1 + 1.96* DLsep1)*10)-1)
DLp1     <- 100*(exp(DLp1*10)-1)

a <- as.data.frame(cbind(DLp1, lcp1, ucp1))
write.csv(a, "D:/Users/profu/Documents/Schoolwork/PhD/Research Projects/no2_42_counties/meta_analysis_holidays/final_model_3df_meta/all.csv")



######TLNISE VS MVMETA PLOTS######
par(cex=1.3, las=1)
plot(c(0,30), c(min(lcp), max(ucp)), type="n",xaxt="n",
     xlab="Lag", ylab=expression(beta), main="Tlnise vs Mvmeta (Non-Accidental Mortality)", ylim=c(-0.00015, 0.0004))

axis(1, at = seq(0,30,5))

polygon(c(rev(0:30), 0:30), c((lcp[31:1]), ucp[1:31]), 
        col="gray90", border = NA, lty=0)

polygon(c(rev(0:30), 0:30), c((lcp1[31:1]), ucp1[1:31]), 
        col="gray40", border = NA, density=25)

lines(0:30, DLp, lty=1, lwd=2, col="black")
lines(0:30, DLp1, lty=2, lwd=2, col="black")

abline(h=0, lty=2)

legend("topleft", inset = 0.02, legend=c("Tlnise", "Mvmeta"), lty=1:2)



######NEW: MVMETA VS FIXED EFFECTS######
setwd("D:/Users/profu/Documents/Schoolwork/PhD/Research Projects/no2_42_counties/meta_analysis_holidays")
oneAnal.coef <- read.csv("lin.coeff.csv", header = TRUE)
oneAnal.vcov <- read.csv("lin.vcov.csv", header = TRUE)

B <- ns(0:30, df=3, intercept=TRUE)
DL1     <- B%*%oneAnal.coef[,2]
DLvcov1 <- B%*%as.matrix(oneAnal.vcov[,2:4])%*%t(B)
DLse1   <- sqrt(diag(DLvcov1))

lc1 <- 100*(exp((DL1 - 1.96* DLse1)*10)-1)
uc1 <- 100*(exp((DL1 + 1.96* DLse1)*10)-1)
DL1     <- 100*(exp(DL1*10)-1)

par(cex=1.3, las=1)
plot(c(0,30), c(min(lcp1), max(ucp1)), type="n",xaxt="n",
     xlab="Lag", ylab="% Change", ylim = c(-0.15, 0.35))

axis(1, at = seq(0,30,5))

polygon(c(rev(0:30), 0:30), c((lcp1[31:1]), ucp1[1:31]), 
        col="gray90", border = NA, lty=0)

polygon(c(rev(0:30), 0:30), c((lc1[31:1]), uc1[1:31]), 
        col="gray40", border = NA, density=25)

lines(0:30, DLp1, lty=1, lwd=2, col="black")
lines(0:30, DL1, lty=2, lwd=2, col="black")

abline(h=0, lty=2)

legend("topleft", inset = 0.02, legend=c("Meta-Analysis", "Fixed Effects"), lty=1:2)



######NEW: FULL VS SUBSET MODELS######
setwd("D:/Users/profu/Documents/Schoolwork/PhD/Research Projects/no2_42_counties/meta_analysis_holidays/final_model_3df_meta/all/subset")
library(mvmeta)
coef.files <- (Sys.glob("./coeff/*.csv"))
coeffs <- lapply(coef.files, function(x) read.table(x, header = TRUE)) 
length(coeffs)
ourY <- matrix(NA, length(coeffs), nrow(coeffs[[1]]))
for (i in 1:dim(ourY)[1]){ ## i = 1
  ourY[i,] <- coeffs[[i]]$x
}

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

B <- ns(0:30, df=3, intercept=TRUE)

DL1     <- B%*%t(as.matrix(metaAll$coefficients))
DLvcov1 <- B%*%metaAll$vcov%*%t(B)
DLse1   <- sqrt(diag(DLvcov1))


lc1 <- 100*(exp((DL1 - 1.96* DLse1)*10)-1)
uc1 <- 100*(exp((DL1 + 1.96* DLse1)*10)-1)
DL1     <- 100*(exp(DL1*10)-1)


par(cex=1.3, las=1)
plot(c(0,30), c(min(lcp1), max(ucp1)), type="n",xaxt="n",
     xlab="Lag", ylab="% Change", ylim = c(-0.15, 0.45))
axis(1, at = seq(0,30,5))
polygon(c(rev(0:30), 0:30), c((lcp1[31:1]), ucp1[1:31]), 
        col="gray90", border = NA, lty=0)
polygon(c(rev(0:30), 0:30), c((lc1[31:1]), uc1[1:31]), 
        col="gray40", border = NA, density=25)
lines(0:30, DLp1, lty=1, lwd=2, col="black")
lines(0:30, DL1, lty=2, lwd=2, col="black")
abline(h=0, lty=2)
legend("topleft", inset = 0.02, legend=c("Full Analysis", "Subset w/ Monitors"), lty=1:2)




######42 DLMs on the same plot######
setwd("D:/Users/profu/Documents/Schoolwork/PhD/Research Projects/no2_42_counties/meta_analysis_holidays")
matfit.files <- (Sys.glob("./matfit/*.csv"))


matfit <- lapply(matfit.files, function(x) read.csv(x, header = TRUE)) 
length(matfit)
par(cex=1.3, las=1)
plot(c(1,151), c(min(dtaX$lci), max(dtaX$uci)), type="n",xaxt="n",
     main="42 Counties, Matfit", xlab="Lag", ylab="% Change",
     ylim=c(-1.5, 2))
axis(1, at = seq(0,200,25), labels = seq(0,40,5))
lines(1:151, matfit[[1]]$pe, lty=1, lwd=2, col="black")
lines(1:151, matfit[[2]]$pe, lty=1, lwd=2, col="black")
lines(1:151, matfit[[3]]$pe, lty=1, lwd=2, col="black")
lines(1:151, matfit[[4]]$pe, lty=1, lwd=2, col="black")
lines(1:151, matfit[[5]]$pe, lty=1, lwd=2, col="black")
lines(1:151, matfit[[6]]$pe, lty=1, lwd=2, col="black")
lines(1:151, matfit[[7]]$pe, lty=1, lwd=2, col="black")
lines(1:151, matfit[[8]]$pe, lty=1, lwd=2, col="black")
lines(1:151, matfit[[9]]$pe, lty=1, lwd=2, col="black")
lines(1:151, matfit[[10]]$pe, lty=1, lwd=2, col="black")
lines(1:151, matfit[[11]]$pe, lty=1, lwd=2, col="black")
lines(1:151, matfit[[12]]$pe, lty=1, lwd=2, col="black")
lines(1:151, matfit[[13]]$pe, lty=1, lwd=2, col="black")
lines(1:151, matfit[[14]]$pe, lty=1, lwd=2, col="black")
lines(1:151, matfit[[15]]$pe, lty=1, lwd=2, col="black")
lines(1:151, matfit[[16]]$pe, lty=1, lwd=2, col="black")
lines(1:151, matfit[[17]]$pe, lty=1, lwd=2, col="black")
lines(1:151, matfit[[18]]$pe, lty=1, lwd=2, col="black")
lines(1:151, matfit[[19]]$pe, lty=1, lwd=2, col="black")
lines(1:151, matfit[[20]]$pe, lty=1, lwd=2, col="black")
lines(1:151, matfit[[21]]$pe, lty=1, lwd=2, col="black")
lines(1:151, matfit[[22]]$pe, lty=1, lwd=2, col="black")
lines(1:151, matfit[[23]]$pe, lty=1, lwd=2, col="black")
lines(1:151, matfit[[24]]$pe, lty=1, lwd=2, col="black")
lines(1:151, matfit[[25]]$pe, lty=1, lwd=2, col="black")
lines(1:151, matfit[[26]]$pe, lty=1, lwd=2, col="black")
lines(1:151, matfit[[27]]$pe, lty=1, lwd=2, col="black")
lines(1:151, matfit[[28]]$pe, lty=1, lwd=2, col="black")
lines(1:151, matfit[[29]]$pe, lty=1, lwd=2, col="black")
lines(1:151, matfit[[30]]$pe, lty=1, lwd=2, col="black")
lines(1:151, matfit[[31]]$pe, lty=1, lwd=2, col="black")
lines(1:151, matfit[[32]]$pe, lty=1, lwd=2, col="black")
lines(1:151, matfit[[33]]$pe, lty=1, lwd=2, col="black")
lines(1:151, matfit[[34]]$pe, lty=1, lwd=2, col="black")
lines(1:151, matfit[[35]]$pe, lty=1, lwd=2, col="black")
lines(1:151, matfit[[36]]$pe, lty=1, lwd=2, col="black")
lines(1:151, matfit[[37]]$pe, lty=1, lwd=2, col="black")
lines(1:151, matfit[[38]]$pe, lty=1, lwd=2, col="black")
lines(1:151, matfit[[39]]$pe, lty=1, lwd=2, col="black")
lines(1:151, matfit[[40]]$pe, lty=1, lwd=2, col="black")
lines(1:151, matfit[[41]]$pe, lty=1, lwd=2, col="black")
lines(1:151, matfit[[42]]$pe, lty=1, lwd=2, col="black")
abline(h=0, lty=2)


setwd("D:/Users/profu/Documents/Schoolwork/PhD/Research Projects/no2_42_counties/meta_analysis_holidays")
cumfit.files <- (Sys.glob("./cumfit/*.csv"))

cumfit <- lapply(cumfit.files, function(x) read.csv(x, header = TRUE)) 
length(cumfit)
par(cex=1.3, las=1)
plot(c(1,31), c(min(dtaX$lci), max(dtaX$uci)), type="n",xaxt="n",
     main="42 Counties, Cumulative Fit", xlab="Lag", ylab="% Change",
     ylim=c(-13, 18))
axis(1, at = seq(0,30,5))
lines(1:31, cumfit[[1]]$pe, lty=1, lwd=2, col="black")
lines(1:31, cumfit[[2]]$pe, lty=1, lwd=2, col="black")
lines(1:31, cumfit[[3]]$pe, lty=1, lwd=2, col="black")
lines(1:31, cumfit[[4]]$pe, lty=1, lwd=2, col="black")
lines(1:31, cumfit[[5]]$pe, lty=1, lwd=2, col="black")
lines(1:31, cumfit[[6]]$pe, lty=1, lwd=2, col="black")
lines(1:31, cumfit[[7]]$pe, lty=1, lwd=2, col="black")
lines(1:31, cumfit[[8]]$pe, lty=1, lwd=2, col="black")
lines(1:31, cumfit[[9]]$pe, lty=1, lwd=2, col="black")
lines(1:31, cumfit[[10]]$pe, lty=1, lwd=2, col="black")
lines(1:31, cumfit[[11]]$pe, lty=1, lwd=2, col="black")
lines(1:31, cumfit[[12]]$pe, lty=1, lwd=2, col="black")
lines(1:31, cumfit[[13]]$pe, lty=1, lwd=2, col="black")
lines(1:31, cumfit[[14]]$pe, lty=1, lwd=2, col="black")
lines(1:31, cumfit[[15]]$pe, lty=1, lwd=2, col="black")
lines(1:31, cumfit[[16]]$pe, lty=1, lwd=2, col="black")
lines(1:31, cumfit[[17]]$pe, lty=1, lwd=2, col="black")
lines(1:31, cumfit[[18]]$pe, lty=1, lwd=2, col="black")
lines(1:31, cumfit[[19]]$pe, lty=1, lwd=2, col="black")
lines(1:31, cumfit[[20]]$pe, lty=1, lwd=2, col="black")
lines(1:31, cumfit[[21]]$pe, lty=1, lwd=2, col="black")
lines(1:31, cumfit[[22]]$pe, lty=1, lwd=2, col="black")
lines(1:31, cumfit[[23]]$pe, lty=1, lwd=2, col="black")
lines(1:31, cumfit[[24]]$pe, lty=1, lwd=2, col="black")
lines(1:31, cumfit[[25]]$pe, lty=1, lwd=2, col="black")
lines(1:31, cumfit[[26]]$pe, lty=1, lwd=2, col="black")
lines(1:31, cumfit[[27]]$pe, lty=1, lwd=2, col="black")
lines(1:31, cumfit[[28]]$pe, lty=1, lwd=2, col="black")
lines(1:31, cumfit[[29]]$pe, lty=1, lwd=2, col="black")
lines(1:31, cumfit[[30]]$pe, lty=1, lwd=2, col="black")
lines(1:31, cumfit[[31]]$pe, lty=1, lwd=2, col="black")
lines(1:31, cumfit[[32]]$pe, lty=1, lwd=2, col="black")
lines(1:31, cumfit[[33]]$pe, lty=1, lwd=2, col="black")
lines(1:31, cumfit[[34]]$pe, lty=1, lwd=2, col="black")
lines(1:31, cumfit[[35]]$pe, lty=1, lwd=2, col="black")
lines(1:31, cumfit[[36]]$pe, lty=1, lwd=2, col="black")
lines(1:31, cumfit[[37]]$pe, lty=1, lwd=2, col="black")
lines(1:31, cumfit[[38]]$pe, lty=1, lwd=2, col="black")
lines(1:31, cumfit[[39]]$pe, lty=1, lwd=2, col="black")
lines(1:31, cumfit[[40]]$pe, lty=1, lwd=2, col="black")
lines(1:31, cumfit[[41]]$pe, lty=1, lwd=2, col="black")
lines(1:31, cumfit[[42]]$pe, lty=1, lwd=2, col="black")
abline(h=0, lty=2)



######42 DLMs GGPLOT######
library(ggplot2)

setwd("D:/Users/profu/Documents/Schoolwork/PhD/Research Projects/no2_42_counties/meta_analysis_holidays/matfit")
matfit <- read.csv("all_matfit.csv", header=TRUE, stringsAsFactors = FALSE)

setwd("D:/Users/profu/Documents/Schoolwork/PhD/Research Projects/no2_42_counties/meta_analysis_holidays/cumfit")
cumfit <- read.csv("all_cumfit.csv", header=TRUE, stringsAsFactors = FALSE)

ggplot(matfit, aes(x=lag, y=estimate)) +
  geom_line(aes(color = county)) +
  labs(
    title = "NO2-Mortality DLM Plots",
    x = "Lag",
    y = "% Change"
  )


ggplot(cumfit, aes(x=lag, y=estimate)) +
  geom_line(aes(color = county)) +
  labs(
    title = "NO2-Mortality DLM Plots, Cumulative",
    x = "Lag",
    y = "% Change"
  )



######FIXED EFFECTS AND INTERACTIONS COMPARISON######
## Meta Coefficients
setwd("D:/Users/profu/Documents/Schoolwork/PhD/Research Projects/no2_42_counties/meta_analysis_holidays")

no2all <- matrix(nrow = 42, ncol = 12)

for (c in 1:42){
  
  subset <- dta[which(dta$countyNumber==c), ]
  daily <- glm(all ~ offset(log(pop)) + lag0to1.no2 + lag0to1.pm25b + as.factor(dow) + ns(as.numeric(newdate), df = 12) + 
                 ns(tmean, df=3) + ns(lag1to3.tmean, df=3) + ns(rh, df=3) + as.factor(holiday), data = subset, family=quasipoisson())
  summary(daily)
  
  ##NO2 Matrix
  # MCode, Coeff, Exp(Coeff), SE, Exp(SE), IQR, Estimate/Lower CI/Upper CI (IQR), Estimate/Lower CI/Upper CI (per 10)
  no2all[c,1] <- unique(subset$ecounty)
  no2all[c,2] <- daily$coefficients[2]
  no2all[c,3] <- exp(daily$coefficients[2])
  no2all[c,4] <- summary(daily)$coef[2,2]
  no2all[c,5] <- exp(summary(daily)$coef[2,2])
  no2all[c,6] <- IQR(subset$no2, na.rm=TRUE)
  no2all[c,7] <- 100*(exp(daily$coefficients[2]*IQR(subset$no2, na.rm=TRUE))-1)
  no2all[c,8] <- 100*(exp((summary(daily)$coef[2,1] - 1.96* summary(daily)$coef[2,2])*IQR(subset$no2, na.rm=TRUE))-1)
  no2all[c,9] <- 100*(exp((summary(daily)$coef[2,1] + 1.96* summary(daily)$coef[2,2])*IQR(subset$no2, na.rm=TRUE))-1)
  no2all[c,10] <- 100*(exp(daily$coefficients[2]*10)-1)
  no2all[c,11] <- 100*(exp((summary(daily)$coef[2,1] - 1.96* summary(daily)$coef[2,2])*10)-1)
  no2all[c,12] <- 100*(exp((summary(daily)$coef[2,1] + 1.96* summary(daily)$coef[2,2])*10)-1)
  
} 

no2all2 <- as.data.frame(no2all)
no2all2 <- rename(no2all2, c("V1"="county", "V2"="coeff", "V3"="expCoeff", "V4"="se", "V5"="expSE", "V6"="IQR", 
                             "V7"="estimate (IQR)", "V8"="lowerCI (IQR)", "V9"="upperCI (IQR)", 
                             "V10"="estimate (per 10)", "V11"="lowerCI (per 10)", "V12"="upperCI (per 10"))

write.csv(no2all2, "D:/Users/profu/Documents/Schoolwork/PhD/Research Projects/no2_42_counties/meta_analysis_holidays/all.csv")


## Fixed Effects Coefficients (+ Interaction Term)
daily <- glm(all ~ offset(log(pop)) + lag0to1.no2 + lag0to1.pm25b + as.factor(dow) + 
               ns(as.numeric(newdate), df = 12) + ns(tmean, df=3) + ns(lag1to3.tmean, df=3) + 
               ns(rh, df=3) + as.factor(holiday) + as.factor(ecounty) + 
               lag0to1.no2*as.factor(ecounty), data = dta, family=quasipoisson())
summary(daily)
table <- as.data.frame(daily)
write.csv(table, "D:/Users/profu/Documents/Schoolwork/PhD/Research Projects/no2_42_counties/meta_analysis_holidays/all_fixed.csv")


## Plotting
library(ggplot2)

setwd("D:/Users/profu/Documents/Schoolwork/PhD/Research Projects/no2_42_counties/meta_analysis_holidays/meta_vs_fixed_lag0")
dta <- read.csv("meta_vs_fixed_lag0.csv", header=TRUE, stringsAsFactors = FALSE)

ggplot(dta, aes(x=coeff_meta, y=coeff_fixed)) +
  geom_point(aes(color = county), size = 4) +
  geom_abline(intercept = 0, slope = 1) +
  geom_text(aes(coeff_meta, coeff_fixed, label=county), size = 2.5) +
  theme(legend.position = "none") +
  labs(
    title = "Fixed Effects vs Meta Analysis Coefficient for NO2",
    x = "Coefficients from Meta Analysis",
    y = "Coefficients from Fixed Effects Model"
  )
  

######Meta-Analysis vs Fixed Effect Estimates
library(mvmeta)
setwd("D:/Users/profu/Documents/Schoolwork/PhD/Research Projects/no2_42_counties/meta_analysis_holidays/meta_vs_fixed_lag0")
all <- read.csv("all_meta.csv", header=TRUE, stringsAsFactors = FALSE)

model.all <- mvmeta(I(all$coeff)~1, S=(all$se)^2, method="reml", na.action="na.omit")
summary(model.all)
100*(exp(model.all$coefficients[1]*10)-1)
100*(exp((summary(model.all)$coef[1] - 1.96* summary(model.all)$coef[2])*10)-1)
100*(exp((summary(model.all)$coef[1] + 1.96* summary(model.all)$coef[2])*10)-1)

daily <- glm(all ~ offset(log(pop)) + lag0to1.no2 + lag0to1.pm25b + as.factor(dow) + 
               ns(as.numeric(newdate), df = 12) + ns(tmean, df=3) + ns(lag1to3.tmean, df=3) + 
               ns(rh, df=3) + as.factor(holiday) + as.factor(ecounty), data = dta, family=quasipoisson())
100*(exp(daily$coefficients[2]*10)-1)
100*(exp((summary(daily)$coef[2] - 1.96* summary(daily)$coef[2,2])*10)-1)
100*(exp((summary(daily)$coef[2] + 1.96* summary(daily)$coef[2,2])*10)-1)



######COUNTY-SPECIFIC ARGLAG TESTING######

##QAIC Function
fqaic <- function(model) {
  loglik <- sum(dpois(model$y,model$fitted.values,log=TRUE))
  phi <- summary(model)$dispersion
  qaic <- -2*loglik + 2*summary(model)$df[3]*phi
  return(qaic)
}

##Loop for arglag 3-7

setwd("D:/Users/profu/Documents/Schoolwork/PhD/Research Projects/no2_42_counties/meta_analysis_holidays/county_specific_df")
no2all <- matrix(nrow = 42, ncol = 6)

my.counties <- unique(dta$ecounty)
for (c in 1:length(my.counties)){ 
  dta.c <- dta[which(dta$ecounty == my.counties[c]),]
  cb1.pm <- crossbasis(dta.c$pm25b, lag=30, argvar=list(fun="lin"),
                       arglag=list(df=3))
  cb2.no2 <- crossbasis(dta.c$no2, lag=30, argvar=list(fun="lin"),
                        arglag=list(df=7))
  cb1.temp <- crossbasis(dta.c$tmean, lag=30, argvar=list(df=3),
                         arglag=list(df=3))
  cb1.rh <- crossbasis(dta.c$rh, lag=30, argvar=list(df=3),
                       arglag=list(df=3))
  all <- glm(all ~ offset(log(pop)) + cb2.no2 + cb1.pm + cb1.temp + cb1.rh + 
               ns(as.numeric(newdate), df = 12) + as.factor(dow) + as.factor(holiday), data = dta.c, family = quasipoisson())
  pred.all <- crosspred(cb2.no2, all, at=0:330, bylag=0.2, cumul=TRUE)
  
  #no2all[c,1] <- unique(dta.c$ecounty)
  #no2all[c,2] <- fqaic(all)
  #no2all[c,3] <- fqaic(all)
  #no2all[c,4] <- fqaic(all)
  #no2all[c,5] <- fqaic(all)
}

write.csv(no2all, "D:/Users/profu/Documents/Schoolwork/PhD/Research Projects/no2_42_counties/meta_analysis_holidays/county_specific_df/3to6.csv")


##County-Specific coeff/vcov with optimal DFs
setwd("D:/Users/profu/Documents/Schoolwork/PhD/Research Projects/no2_42_counties/meta_analysis_holidays/county_specific_df")
bestDF <- read.csv("bestDF.csv", header=TRUE, stringsAsFactors = FALSE)

dta <- merge(dta, bestDF, by="ecounty")

my.counties <- unique(dta$ecounty)
for (c in 1:length(my.counties)){ 
  dta.c <- dta[which(dta$ecounty == my.counties[c]),]
  cb1.pm <- crossbasis(dta.c$pm25b, lag=30, argvar=list(fun="lin"),
                       arglag=list(df=3))
  cb2.no2 <- crossbasis(dta.c$no2, lag=30, argvar=list(fun="lin"),
                        arglag=list(df=dta.c$lowest_df))
  cb1.temp <- crossbasis(dta.c$tmean, lag=30, argvar=list(df=3),
                         arglag=list(df=3))
  cb1.rh <- crossbasis(dta.c$rh, lag=30, argvar=list(df=3),
                       arglag=list(df=3))
  all <- glm(all ~ offset(log(pop)) + cb2.no2 + cb1.pm + cb1.temp + cb1.rh + 
               ns(as.numeric(newdate), df = 12) + as.factor(dow) + as.factor(holiday), data = dta.c, family = quasipoisson())
  pred.all <- crosspred(cb2.no2, all, at=0:330, bylag=0.2, cumul=TRUE)
  
  a <- dta.c$lowest_df + 1
  #lin.coeff <- summary(all)$coefficients[2:a,1]
  lin.vcov <- vcov(all)[2:a,2:a]
  #nm1 <- paste(my.counties[c], "_coeff.csv", sep="")
  nm2 <- paste(my.counties[c], "_vcov.csv", sep="")
  #write.csv(lin.coeff, file = nm1, row.names = FALSE)
  write.csv(lin.vcov, file = nm2, row.names = FALSE)
}



######META-ANALYSIS USING OPTIMAL DF######
setwd("D:/Users/profu/Documents/Schoolwork/PhD/Research Projects/no2_42_counties/meta_analysis_holidays/county_specific_df")

library(mvmeta)
library(splines)

coef.files <- (Sys.glob("./coeff/*.csv"))

coeffs <- lapply(coef.files, function(x) read.table(x, header = TRUE)) 
length(coeffs)


vcov.files <- (Sys.glob("./vcov/*.csv"))

vcovs <- lapply(vcov.files, function(x) read.csv(x, header = TRUE)) 
length(vcovs)

B3 <- ns(0:30, df=3, intercept=TRUE)
B4 <- ns(0:30, df=4, intercept=TRUE)
B5 <- ns(0:30, df=5, intercept=TRUE)
B6 <- ns(0:30, df=6, intercept=TRUE)

#####Creating correctly specified matrix (by hand because I can't think how to automate at the moment)#####
coeffs[[1]] <- B5%*%as.matrix(coeffs[[1]])
coeffs[[2]] <- B6%*%as.matrix(coeffs[[2]])
coeffs[[3]] <- B6%*%as.matrix(coeffs[[3]])
coeffs[[4]] <- B3%*%as.matrix(coeffs[[4]])
coeffs[[5]] <- B3%*%as.matrix(coeffs[[5]])
coeffs[[6]] <- B3%*%as.matrix(coeffs[[6]])
coeffs[[7]] <- B6%*%as.matrix(coeffs[[7]])
coeffs[[8]] <- B3%*%as.matrix(coeffs[[8]])
coeffs[[9]] <- B5%*%as.matrix(coeffs[[9]])
coeffs[[10]] <- B4%*%as.matrix(coeffs[[10]])

coeffs[[11]] <- B5%*%as.matrix(coeffs[[11]])
coeffs[[12]] <- B4%*%as.matrix(coeffs[[12]])
coeffs[[13]] <- B3%*%as.matrix(coeffs[[13]])
coeffs[[14]] <- B4%*%as.matrix(coeffs[[14]])
coeffs[[15]] <- B3%*%as.matrix(coeffs[[15]])
coeffs[[16]] <- B3%*%as.matrix(coeffs[[16]])
coeffs[[17]] <- B4%*%as.matrix(coeffs[[17]])
coeffs[[18]] <- B3%*%as.matrix(coeffs[[18]])
coeffs[[19]] <- B3%*%as.matrix(coeffs[[19]])
coeffs[[20]] <- B3%*%as.matrix(coeffs[[20]])

coeffs[[21]] <- B3%*%as.matrix(coeffs[[21]])
coeffs[[22]] <- B3%*%as.matrix(coeffs[[22]])
coeffs[[23]] <- B3%*%as.matrix(coeffs[[23]])
coeffs[[24]] <- B6%*%as.matrix(coeffs[[24]])
coeffs[[25]] <- B3%*%as.matrix(coeffs[[25]])
coeffs[[26]] <- B3%*%as.matrix(coeffs[[26]])
coeffs[[27]] <- B4%*%as.matrix(coeffs[[27]])
coeffs[[28]] <- B3%*%as.matrix(coeffs[[28]])
coeffs[[29]] <- B5%*%as.matrix(coeffs[[29]])
coeffs[[30]] <- B3%*%as.matrix(coeffs[[30]])

coeffs[[31]] <- B3%*%as.matrix(coeffs[[31]])
coeffs[[32]] <- B3%*%as.matrix(coeffs[[32]])
coeffs[[33]] <- B3%*%as.matrix(coeffs[[33]])
coeffs[[34]] <- B6%*%as.matrix(coeffs[[34]])
coeffs[[35]] <- B3%*%as.matrix(coeffs[[35]])
coeffs[[36]] <- B3%*%as.matrix(coeffs[[36]])
coeffs[[37]] <- B3%*%as.matrix(coeffs[[37]])
coeffs[[38]] <- B6%*%as.matrix(coeffs[[38]])
coeffs[[39]] <- B3%*%as.matrix(coeffs[[39]])
coeffs[[40]] <- B3%*%as.matrix(coeffs[[40]])
coeffs[[41]] <- B4%*%as.matrix(coeffs[[41]])
coeffs[[42]] <- B5%*%as.matrix(coeffs[[42]])


##vcovs
vcovs[[1]] <- B5%*%as.matrix(vcovs[[1]])%*%t(B5)
vcovs[[2]] <- B6%*%as.matrix(vcovs[[2]])%*%t(B6)
vcovs[[3]] <- B6%*%as.matrix(vcovs[[3]])%*%t(B6)
vcovs[[4]] <- B3%*%as.matrix(vcovs[[4]])%*%t(B3)
vcovs[[5]] <- B3%*%as.matrix(vcovs[[5]])%*%t(B3)
vcovs[[6]] <- B3%*%as.matrix(vcovs[[6]])%*%t(B3)
vcovs[[7]] <- B6%*%as.matrix(vcovs[[7]])%*%t(B6)
vcovs[[8]] <- B3%*%as.matrix(vcovs[[8]])%*%t(B3)
vcovs[[9]] <- B5%*%as.matrix(vcovs[[9]])%*%t(B5)
vcovs[[10]] <- B4%*%as.matrix(vcovs[[10]])%*%t(B4)

vcovs[[11]] <- B5%*%as.matrix(vcovs[[11]])%*%t(B5)
vcovs[[12]] <- B4%*%as.matrix(vcovs[[12]])%*%t(B4)
vcovs[[13]] <- B3%*%as.matrix(vcovs[[13]])%*%t(B3)
vcovs[[14]] <- B4%*%as.matrix(vcovs[[14]])%*%t(B4)
vcovs[[15]] <- B3%*%as.matrix(vcovs[[15]])%*%t(B3)
vcovs[[16]] <- B3%*%as.matrix(vcovs[[16]])%*%t(B3)
vcovs[[17]] <- B4%*%as.matrix(vcovs[[17]])%*%t(B4)
vcovs[[18]] <- B3%*%as.matrix(vcovs[[18]])%*%t(B3)
vcovs[[19]] <- B3%*%as.matrix(vcovs[[19]])%*%t(B3)
vcovs[[20]] <- B3%*%as.matrix(vcovs[[20]])%*%t(B3)

vcovs[[21]] <- B3%*%as.matrix(vcovs[[21]])%*%t(B3)
vcovs[[22]] <- B3%*%as.matrix(vcovs[[22]])%*%t(B3)
vcovs[[23]] <- B3%*%as.matrix(vcovs[[23]])%*%t(B3)
vcovs[[24]] <- B6%*%as.matrix(vcovs[[24]])%*%t(B6)
vcovs[[25]] <- B3%*%as.matrix(vcovs[[25]])%*%t(B3)
vcovs[[26]] <- B3%*%as.matrix(vcovs[[26]])%*%t(B3)
vcovs[[27]] <- B4%*%as.matrix(vcovs[[27]])%*%t(B4)
vcovs[[28]] <- B3%*%as.matrix(vcovs[[28]])%*%t(B3)
vcovs[[29]] <- B5%*%as.matrix(vcovs[[29]])%*%t(B5)
vcovs[[30]] <- B3%*%as.matrix(vcovs[[30]])%*%t(B3)

vcovs[[31]] <- B3%*%as.matrix(vcovs[[31]])%*%t(B3)
vcovs[[32]] <- B3%*%as.matrix(vcovs[[32]])%*%t(B3)
vcovs[[33]] <- B3%*%as.matrix(vcovs[[33]])%*%t(B3)
vcovs[[34]] <- B6%*%as.matrix(vcovs[[34]])%*%t(B6)
vcovs[[35]] <- B3%*%as.matrix(vcovs[[35]])%*%t(B3)
vcovs[[36]] <- B3%*%as.matrix(vcovs[[36]])%*%t(B3)
vcovs[[37]] <- B3%*%as.matrix(vcovs[[37]])%*%t(B3)
vcovs[[38]] <- B6%*%as.matrix(vcovs[[38]])%*%t(B6)
vcovs[[39]] <- B3%*%as.matrix(vcovs[[39]])%*%t(B3)
vcovs[[40]] <- B3%*%as.matrix(vcovs[[40]])%*%t(B3)
vcovs[[41]] <- B4%*%as.matrix(vcovs[[41]])%*%t(B4)
vcovs[[42]] <- B5%*%as.matrix(vcovs[[42]])%*%t(B5)

ourY <- matrix(NA, length(coeffs), nrow(coeffs[[1]]))

for (i in 1:dim(ourY)[1]){ ## i = 1
  
  ourY[i,] <- coeffs[[i]]
  
}


ourV <- array(NA, dim = c(nrow(vcovs[[1]]), ncol(vcovs[[1]]), length(vcovs)))

for (i in 1:dim(ourV)[3]){ ## i = 1
  
  ourV[,,i] <- as.matrix(vcovs[[i]])
  
}


for (i in 1:dim(ourV)[3]){
  
  print(na.omit(dim(ourV[,,i])))
  
}


######TLNISE and MVMETA######
library(tlnise)
library(mvmeta)

##tlnise: doesn't work
pooled <- tlnise(Y=ourY, V=ourV, intercept=FALSE, seed=68460, maxiter=1000, Tol = 1e-30)

DLp     <- B%*%pooled$gamma[,1]
DLvcovp <- B%*%pooled$Dgamma%*%t(B)
DLsep   <- sqrt(diag(DLvcovp))
lcp <- DLp - 1.96*DLsep
ucp <- DLp + 1.96*DLsep


##mvmeta: it runs, but the output is all wrong
ourV <- lapply(seq(dim(ourV)[3]), function(x) ourV[ , , x])
metaAll <- mvmeta(ourY~1,ourV,method="reml")
summary(metaAll)

DLp1     <- B%*%t(as.matrix(metaAll$coefficients))
DLvcovp1 <- B%*%metaAll$vcov%*%t(B)
DLsep1   <- sqrt(diag(DLvcovp1))
lcp1 <- DLp1 - 1.96*DLsep1
ucp1 <- DLp1 + 1.96*DLsep1



######SATURATED MODELS######
## Fixed effects + NO2*County only
daily <- glm(all ~ offset(log(pop)) + lag0to1.no2*as.factor(ecounty) + lag0to1.pm25b + as.factor(dow) + 
               ns(as.numeric(newdate), df = 12) + ns(tmean, df=3) + ns(lag1to3.tmean, df=3) + 
               ns(rh, df=3) + as.factor(holiday), data = dta, family=quasipoisson())
table <- as.matrix(daily$coefficients)
table2 <- as.matrix(summary(daily)$coefficients[,2])
table3 <- as.matrix(vcov(daily)[,1])
table4 <- as.data.frame(cbind(table, table2, table3))
write.csv(table4, "D:/Users/profu/Documents/Schoolwork/PhD/Research Projects/no2_42_counties/meta_analysis_holidays/fixed.csv")



## Saturated model
saturated <- glm(all ~ offset(log(pop)) + lag0to1.no2*as.factor(ecounty) + lag0to1.pm25b*as.factor(ecounty) + 
                   as.factor(dow)*as.factor(ecounty) + ns(as.numeric(newdate), df = 12)*as.factor(ecounty) + 
                   ns(tmean, df=3)*as.factor(ecounty) + ns(lag1to3.tmean, df=3)*as.factor(ecounty) + 
                   ns(rh, df=3)*as.factor(ecounty) + as.factor(holiday)*as.factor(ecounty), data = dta, family=quasipoisson())

table <- as.matrix(saturated$coefficients)
table2 <- as.matrix(summary(saturated)$coefficients[,2])
table3 <- as.matrix(vcov(saturated)[,1])
table4 <- as.data.frame(cbind(table, table2, table3))
write.csv(table4, "D:/Users/profu/Documents/Schoolwork/PhD/Research Projects/no2_42_counties/meta_analysis_holidays/saturated.csv")



######POINT PLOTS: INTERACTION AND SATURATED VS COUNTY-SPECIFIC######
library(ggplot2)

setwd("D:/Users/profu/Documents/Schoolwork/PhD/Research Projects/no2_42_counties/meta_analysis_holidays/saturated_model")
plot <- read.csv("meta_fixed_saturated.csv", header=TRUE, stringsAsFactors = FALSE)

ggplot(plot, aes(x=coeff_meta, y=coeff_fixed)) +
  geom_point(aes(color = county), size = 4) +
  #coord_cartesian(xlim = c(-0.002, 0.005), ylim = c(-0.002, 0.004)) +
  coord_cartesian(xlim = c(-0.0035, 0.0065), ylim = c(-0.01, 0.01)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_text(aes(coeff_meta, coeff_fixed, label=county), size = 2.5) +
  geom_errorbarh(aes(xmin = coeff_meta-se_meta, xmax = coeff_meta+se_meta)) +
  geom_errorbar(aes(ymin = coeff_fixed-se_fixed, ymax = coeff_fixed+se_fixed)) +
  theme(legend.position = "none") +
  labs(
    title = "Fixed Effects vs Meta Analysis Coefficient for NO2",
    x = "Coefficients from Meta Analysis",
    y = "Coefficients from Fixed Effects Model"
  )

ggplot(plot, aes(x=coeff_meta, y=coeff_saturated)) +
  geom_point(aes(color = county), size = 4) +
  #coord_cartesian(xlim = c(-0.002, 0.005), ylim = c(-0.002, 0.004)) +
  coord_cartesian(xlim = c(-0.0035, 0.0065), ylim = c(-0.01, 0.01)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_text(aes(coeff_meta, coeff_saturated, label=county), size = 2.5) +
  geom_errorbarh(aes(xmin = coeff_meta-se_meta, xmax = coeff_meta+se_meta)) +
  geom_errorbar(aes(ymin = coeff_saturated-se_saturated, ymax = coeff_saturated+se_saturated)) +
  theme(legend.position = "none") +
  labs(
    title = "(Saturated) Effects vs Meta Analysis Coefficient for NO2",
    x = "Coefficients from Meta Analysis",
    y = "Coefficients from (Saturated) Model"
  )


######SANITY CHECK: TESTING IF ORIGINAL AND RECONSTRUCTED DLM LOOKS THE SAME######
dta.c <- dta[which(dta$ecounty=="Caidian"),]

cb1.pm <- crossbasis(dta.c$pm25b, lag=30, argvar=list(fun="lin"),
                     arglag=list(df=3))
cb2.no2 <- crossbasis(dta.c$no2, lag=30, argvar=list(fun="lin"),
                      arglag=list(df=5))
cb1.temp <- crossbasis(dta.c$tmean, lag=30, argvar=list(df=3),
                       arglag=list(df=3))
cb1.rh <- crossbasis(dta.c$rh, lag=30, argvar=list(df=3),
                     arglag=list(df=3))
all <- glm(all ~ offset(log(pop)) + cb2.no2 + cb1.pm + cb1.temp + cb1.rh + 
             ns(as.numeric(newdate), df = 12) + as.factor(dow) + as.factor(holiday), data = dta.c, family = quasipoisson())
pred.all <- crosspred(cb2.no2, all, at=0:330, bylag=0.2, cumul=TRUE)

plot(pred.all, "slices", var=10, col=3, ylab="RR",  ci.arg=list(density=15,lwd=2),
     main="Association with a 10-unit increase in NO2")


setwd("D:/Users/profu/Documents/Schoolwork/PhD/Research Projects/no2_42_counties/meta_analysis_holidays/county_specific_df/coeff")
oneAnal.coef <- read.csv("Caidian_coeff.csv", header = TRUE)

setwd("D:/Users/profu/Documents/Schoolwork/PhD/Research Projects/no2_42_counties/meta_analysis_holidays/county_specific_df/vcov")
oneAnal.vcov <- read.csv("Caidian_vcov.csv", header = TRUE)

B <- ns(0:30, df=5, intercept=TRUE)

DL1     <- B%*%oneAnal.coef[,1]
DLvcov1 <- B%*%as.matrix(oneAnal.vcov[,1:5])%*%t(B)
DLse1   <- sqrt(diag(DLvcov1))
lc1 <- DL1 - 1.96*DLse1
uc1 <- DL1 + 1.96*DLse1

par(cex=1.3, las=1)
plot(c(0,30), c(min(lc1), max(uc1)), type="n",xaxt="n",
     xlab="Lag", ylab=expression(beta), main="Caidian, 5df")
axis(1, at = seq(0,30,5))
polygon(c(rev(0:30), 0:30), c((lc1[nrow(lc1):1]), uc1[1:nrow(uc1)]), 
        col="gray40", border = NA, density=25)
lines(0:30, DL1, lty=2, lwd=2, col="black")
abline(h=0, lty=2)



######3, 4, 5, AND 6 DF TLNISE (AND MVMETA) PLOTS######
library(MASS)
library(devtools) 
library(tlnise)
library(mvmeta)
library(ggplot2)
library(splines)
library(dlnm)

setwd("D:/Users/profu/Documents/Schoolwork/PhD/Research Projects/no2_42_counties/meta_analysis_holidays/county_specific_df/3df")
setwd("D:/Users/profu/Documents/Schoolwork/PhD/Research Projects/no2_42_counties/meta_analysis_holidays/county_specific_df/4df")
setwd("D:/Users/profu/Documents/Schoolwork/PhD/Research Projects/no2_42_counties/meta_analysis_holidays/county_specific_df/5df")
setwd("D:/Users/profu/Documents/Schoolwork/PhD/Research Projects/no2_42_counties/meta_analysis_holidays/county_specific_df/6df")


## Run this separately for each df
coef.files <- (Sys.glob("./coeff/*.csv"))

coeffs <- lapply(coef.files, function(x) read.table(x, header = TRUE)) 
length(coeffs)

ourY <- matrix(NA, length(coeffs), nrow(coeffs[[1]]))

for (i in 1:dim(ourY)[1]){ ## i = 1
  
  ourY[i,] <- coeffs[[i]]$x
  
}

vcov.files <- (Sys.glob("./vcov/*.csv"))

vcovs <- lapply(vcov.files, function(x) read.csv(x, header = TRUE)) 
length(vcovs)

ourV <- array(NA, dim = c(nrow(vcovs[[1]]), ncol(vcovs[[1]]), length(vcovs)))

for (i in 1:dim(ourV)[3]){ ## i = 1
  
  ourV[,,i] <- as.matrix(vcovs[[i]])
  
}
ourV <- lapply(seq(dim(ourV)[3]), function(x) ourV[ , , x])

#3df
pooled <- tlnise(Y=ourY, V=ourV, intercept=FALSE, seed=68460, maxiter=10000)
B <- ns(0:30, df=3, intercept=TRUE)
DL3     <- B%*%pooled$gamma[,1]
DLvcovp <- B%*%pooled$Dgamma%*%t(B)
DLsep   <- sqrt(diag(DLvcovp))
lc3 <- DL3 - 1.96*DLsep
uc3 <- DL3 + 1.96*DLsep

metaAll <- mvmeta(ourY~1,ourV,method="reml")
xlag <- 0:300/10
blag <- do.call("onebasis",c(list(x=xlag),attr(cb2.no2,"arglag")))
cpall <- crosspred(blag,coef=coef(metaAll),vcov=vcov(metaAll),
                   model.link="log",at=0:30)
pe3  <- 100*(exp(cpall$matfit*10)-1)
lci3 <- 100*(exp((cpall$matfit - 1.96* cpall$matse)*10)-1)
uci3 <- 100*(exp((cpall$matfit + 1.96* cpall$matse)*10)-1)


#4df
pooled <- tlnise(Y=ourY, V=ourV, intercept=FALSE, seed=68460, maxiter=10000)
B <- ns(0:30, df=4, intercept=TRUE)
DL4     <- B%*%pooled$gamma[,1]
DLvcovp <- B%*%pooled$Dgamma%*%t(B)
DLsep   <- sqrt(diag(DLvcovp))
lc4 <- DL4 - 1.96*DLsep
uc4 <- DL4 + 1.96*DLsep

metaAll <- mvmeta(ourY~1,ourV,method="reml")
xlag <- 0:300/10
blag <- do.call("onebasis",c(list(x=xlag),attr(cb2.no2,"arglag")))
cpall <- crosspred(blag,coef=coef(metaAll),vcov=vcov(metaAll),
                   model.link="log",at=0:30)
pe4  <- 100*(exp(cpall$matfit*10)-1)
lci4 <- 100*(exp((cpall$matfit - 1.96* cpall$matse)*10)-1)
uci4 <- 100*(exp((cpall$matfit + 1.96* cpall$matse)*10)-1)

#5df
pooled <- tlnise(Y=ourY, V=ourV, intercept=FALSE, seed=68460, maxiter=10000)
B <- ns(0:30, df=5, intercept=TRUE)
DL5     <- B%*%pooled$gamma[,1]
DLvcovp <- B%*%pooled$Dgamma%*%t(B)
DLsep   <- sqrt(diag(DLvcovp))
lc5 <- DL5 - 1.96*DLsep
uc5 <- DL5 + 1.96*DLsep

metaAll <- mvmeta(ourY~1,ourV,method="reml")
xlag <- 0:300/10
blag <- do.call("onebasis",c(list(x=xlag),attr(cb2.no2,"arglag")))
cpall <- crosspred(blag,coef=coef(metaAll),vcov=vcov(metaAll),
                   model.link="log",at=0:30)
pe5  <- 100*(exp(cpall$matfit*10)-1)
lci5 <- 100*(exp((cpall$matfit - 1.96* cpall$matse)*10)-1)
uci5 <- 100*(exp((cpall$matfit + 1.96* cpall$matse)*10)-1)

#6df
pooled <- tlnise(Y=ourY, V=ourV, intercept=FALSE, seed=68460, maxiter=10000)
B <- ns(0:30, df=6, intercept=TRUE)
DL6     <- B%*%pooled$gamma[,1]
DLvcovp <- B%*%pooled$Dgamma%*%t(B)
DLsep   <- sqrt(diag(DLvcovp))
lc6 <- DL6 - 1.96*DLsep
uc6 <- DL6 + 1.96*DLsep

metaAll <- mvmeta(ourY~1,ourV,method="reml")
xlag <- 0:300/10
blag <- do.call("onebasis",c(list(x=xlag),attr(cb2.no2,"arglag")))
cpall <- crosspred(blag,coef=coef(metaAll),vcov=vcov(metaAll),
                   model.link="log",at=0:30)
pe6  <- 100*(exp(cpall$matfit*10)-1)
lci6 <- 100*(exp((cpall$matfit - 1.96* cpall$matse)*10)-1)
uci6 <- 100*(exp((cpall$matfit + 1.96* cpall$matse)*10)-1)

allDL <- as.data.frame(cbind(DL3, DL4, DL5, DL6))
write.csv(allDL, "D:/Users/profu/Documents/Schoolwork/PhD/Research Projects/no2_42_counties/meta_analysis_holidays/dl3to6.csv")

## Plotting
par(cex=1.3, las=1)
plot(c(0,30), c(min(lc6), max(uc6)), type="n",xaxt="n",
     xlab="Lag", ylab=expression(beta), main="Tlnise, 6df")
#ylim=c(-0.00015, 0.0004)
axis(1, at = seq(0,30,5))
polygon(c(rev(0:30), 0:30), c((lc3[31:1]), uc3[1:31]), 
        col="gray90", border = NA, lty=0)
polygon(c(rev(0:30), 0:30), c((lc4[31:1]), uc4[1:31]), 
        col="gray90", border = NA, lty=0)
polygon(c(rev(0:30), 0:30), c((lc5[31:1]), uc5[1:31]), 
        col="gray90", border = NA, lty=0)
polygon(c(rev(0:30), 0:30), c((lc6[31:1]), uc6[1:31]), 
        col="gray90", border = NA, lty=0)
lines(0:30, DL3, lty=1, lwd=2, col="black")
lines(0:30, DL4, lty=2, lwd=2, col="black")
lines(0:30, DL5, lty=3, lwd=2, col="black")
lines(0:30, DL6, lty=4, lwd=2, col="black")
abline(h=0, lty=2)
legend("topleft", inset = 0.02, legend=c("3df", "4df", "5df", "6df"), lty=1:4, ncol=2)

par(cex=1.3, las=1)
plot(c(0,30), c(min(lci5), max(uci6)), type="n",xaxt="n",
     xlab="Lag", ylab="% Change")
#ylim=c(-0.8, 1.0)
axis(1, at = seq(0,30,5))
polygon(c(rev(0:30), 0:30), c((lci3[31:1]), uci3[1:31]), 
        col="gray90", border = NA, lty=0)
polygon(c(rev(0:30), 0:30), c((lci4[31:1]), uci4[1:31]), 
        col="gray90", border = NA, lty=0)
polygon(c(rev(0:30), 0:30), c((lci5[31:1]), uci5[1:31]), 
        col="gray90", border = NA, lty=0)
polygon(c(rev(0:30), 0:30), c((lci6[31:1]), uci6[1:31]), 
        col="gray90", border = NA, lty=0)
lines(0:30, pe3, lty=1, lwd=2, col="black")
lines(0:30, pe4, lty=1, lwd=2, col="black")
lines(0:30, pe5, lty=1, lwd=2, col="black")
lines(0:30, pe6, lty=1, lwd=2, col="black")
abline(h=0, lty=2)
legend("topleft", "(d)", bty = "n")
#legend("topleft", inset = 0.02, legend=c("3df", "4df", "5df", "6df"), lty=1:4, ncol=2)


## Cumulative Plot
setwd("D:/Users/profu/Documents/Schoolwork/PhD/Research Projects/no2_42_counties/meta_analysis_holidays")
allDL <- read.csv("dl3to6.csv", header=TRUE, stringsAsFactors = FALSE)

par(cex=1.3, las=1)
plot(c(0,30), c(-0.002, 0.003), type="n",xaxt="n",
     xlab="Lag", ylab=expression(beta), main="Tlnise Cumulative, 3-6df")
axis(1, at = seq(0,30,5))

lines(0:30, allDL$X3df_cum, lty=1, lwd=2, col="black")
lines(0:30, allDL$X4df_cum, lty=2, lwd=2, col="black")
lines(0:30, allDL$X5df_cum, lty=3, lwd=2, col="black")
lines(0:30, allDL$X6df_cum, lty=4, lwd=2, col="black")

abline(h=0, lty=2)

legend("topleft", inset = 0.02, legend=c("3df", "4df", "5df", "6df"), lty=1:4, ncol=2)



######NEW MVMETA PLOTS USING GGPLOT (THESE LOOK TERRIBLE)######
library(ggplot2)
setwd("D:/Users/profu/Documents/Schoolwork/PhD/Research Projects/no2_42_counties/meta_analysis_holidays/final_model_3df_meta")
plot <- read.csv("model_estimates.csv", header=TRUE, stringsAsFactors = FALSE)

## All
ggplot(plot, aes(x=lag, y=pe_all)) +
  geom_line() +
  geom_ribbon(data=plot, aes(ymin=lci_all, ymax=uci_all), alpha=0.3) +
  geom_hline(yintercept=0, linetype="dashed") +
  theme_classic(base_size = 15) +
  labs(
    x = "Lag",
    y = "% Change"
  )

## All Men and Women
#line_types <- c("LINE1" = 1, "LINE2" = 5)
ggplot(plot, aes(x=lag)) +
  geom_line(data=plot, aes(x=lag, y=pe_allm), linetype = 1) +
  geom_line(data=plot, aes(x=lag, y=pe_allf), linetype = 5) +
  scale_linetype_manual(values=line_types) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_ribbon(data=plot, aes(ymin=lci_allm, ymax=uci_allm), alpha=0.3) +
  geom_ribbon(data=plot, aes(ymin=lci_allf, ymax=uci_allf), alpha=0.1) +
  theme_classic(base_size = 15) +
  labs(
    x = "Lag",
    y = "% Change"
  )

## CVD
ggplot(plot, aes(x=lag, y=pe_cvd)) +
  geom_line() +
  geom_ribbon(data=plot, aes(ymin=lci_cvd, ymax=uci_cvd), alpha=0.3) +
  geom_hline(yintercept=0, linetype="dashed") +
  theme_classic(base_size = 15) +
  labs(
    x = "Lag",
    y = "% Change"
  )

## CVD Men and Women
ggplot(plot, aes(x=lag)) +
  geom_line(data=plot, aes(x=lag, y=pe_cvdm), linetype = 1) +
  geom_line(data=plot, aes(x=lag, y=pe_cvdf), linetype = 5) +
  scale_linetype_manual(values=line_types) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_ribbon(data=plot, aes(ymin=lci_cvdm, ymax=uci_cvdm), alpha=0.3) +
  geom_ribbon(data=plot, aes(ymin=lci_cvdf, ymax=uci_cvdf), alpha=0.1) +
  theme_classic(base_size = 15) +
  labs(
    x = "Lag",
    y = "% Change"
  )

## Resp
ggplot(plot, aes(x=lag, y=pe_res)) +
  geom_line() +
  geom_ribbon(data=plot, aes(ymin=lci_res, ymax=uci_res), alpha=0.3) +
  geom_hline(yintercept=0, linetype="dashed") +
  theme_classic(base_size = 15) +
  labs(
    x = "Lag",
    y = "% Change"
  )


## Resp Men and Women
ggplot(plot, aes(x=lag)) +
  geom_line(data=plot, aes(x=lag, y=pe_resm), linetype = 1) +
  geom_line(data=plot, aes(x=lag, y=pe_resf), linetype = 5) +
  scale_linetype_manual(values=line_types) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_ribbon(data=plot, aes(ymin=lci_resm, ymax=uci_resm), alpha=0.3) +
  geom_ribbon(data=plot, aes(ymin=lci_resf, ymax=uci_resf), alpha=0.1) +
  theme_classic(base_size = 15) +
  labs(
    x = "Lag",
    y = "% Change"
  )



## New Cumulative Plot (Mvmeta, 3df)
plot <- read.csv("model_estimates_cumulative.csv", header=TRUE, stringsAsFactors = FALSE)

## All Men and Women
ggplot(plot, aes(x=lag)) +
  geom_line(data=plot, aes(x=lag, y=pe_all_cum)) +
  geom_line(data=plot, aes(x=lag, y=pe_allm_cum), linetype = "longdash") +
  geom_line(data=plot, aes(x=lag, y=pe_allf_cum), linetype = "dotdash") +
  geom_hline(yintercept=0, linetype="dashed") +
  theme_classic(base_size = 15) +
  labs(
    x = "Lag",
    y = "% Change"
  )

## Cvd Men and Women
ggplot(plot, aes(x=lag)) +
  geom_line(data=plot, aes(x=lag, y=pe_cvd_cum)) +
  geom_line(data=plot, aes(x=lag, y=pe_cvdm_cum), linetype = "longdash") +
  geom_line(data=plot, aes(x=lag, y=pe_cvdf_cum), linetype = "dotdash") +
  geom_hline(yintercept=0, linetype="dashed") +
  theme_classic(base_size = 15) +
  labs(
    x = "Lag",
    y = "% Change"
  )

## Resp Men and Women
ggplot(plot, aes(x=lag)) +
  geom_line(data=plot, aes(x=lag, y=pe_res_cum)) +
  geom_line(data=plot, aes(x=lag, y=pe_resm_cum), linetype = "longdash") +
  geom_line(data=plot, aes(x=lag, y=pe_resf_cum), linetype = "dotdash") +
  geom_hline(yintercept=0, linetype="dashed") +
  theme_classic(base_size = 15) +
  labs(
    x = "Lag",
    y = "% Change"
  )


####Subsetting to figure out final sample size####
myvars <- c("pm25b", "no2", "tmean", "rh", "mcode", "ecounty", "newdate", "dow", "holiday", "pop", "all", 
            "allg2", "cir", "cirg1", "cirg2", "res", "resg1", "resg2")
newdata <- dta[myvars]
newdata2 <- na.omit(newdata)
table(newdata2$ecounty)
