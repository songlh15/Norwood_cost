# lib2 <- c('stringi','stringr','reshape','plyr','dplyr','tidyr',
#           'gridExtra','data.table','lattice','latticeExtra',
#           'ggplot2','sqldf','haven','readr','shiny','DT','tools','cellranger',
#           'broom','cellranger','clipr','curl','dbplyr','fs','hexbin','httr',
#           'modelr','lubridate','plotly','selectr','tidyverse','whisker','xml2',
#           'readxl','reprex','rvest','sqldf','DT','stringr','tools','data.table',
#           'haven','ggplot2','leaflet','shinyWidgets','gdtools','maps','systemfonts',
#           'uuid','ggiraph','shinyjs','shinythemes','googlesheets4','lqmm','qrLMM','table1','gtsummary')
# 
# 
# if (length(setdiff(lib2, rownames(installed.packages()))) > 0) {
#   suppressMessages(suppressPackageStartupMessages((install.packages(setdiff(lib2, rownames(installed.packages())),repos = "http://cran.us.r-project.org"))))
# }
# install.packages('gtsummary')

library(data.table)
library(table1)
library(lqmm)
library(qrLMM)
library(quantreg)
library(gtsummary)

cost <- fread("M:\\projects\\Lihai\\Jing\\Norwood\\data\\costdata.csv")
attach(cost)
metadata <- list(
  labels=list(
    SHUNTRND = "Cohort",
    sext = "Sex",
    racet = "Race",
    norwage = "Age at Norwood",
    payort = "Insurance"
  ),
  units=list(
    age = "years",
    wgt = "kg",
    norwage = "Days"),
  categoricals=list(
    sex = list(
      `1` = "Female",
      `2` = "Male"),
    racet = list(
      `1` = "White",
      `2` = "Black",
      `3` = 'Other'),
    payort = list(
      '1' = 'Government',
      '2' = 'Private',
      '4' = 'Other'
    )
    
  ))

data <- t1read(cost, metadata)


rndr <- function(x, name, ...) {
  if (length(x) == 0) {
    y <- data[[name]]
    s <- rep("", length(render.default(x=y, name=name, ...)))
    if (is.numeric(y)) {
      p <- t.test(y ~ data$SHUNTRND)$p.value
    } else {
      p <- chisq.test(table(y, droplevels(data$SHUNTRND)))$p.value
    }
    s[2] <- sub("<", "&lt;", format.pval(p, digits=3, eps=0.001))
    s
  } else {
    render.default(x=x, name=name, ...)
  }
}

rndr.strat <- function(label, n, ...) {
  ifelse(n==0, label, render.strat.default(label, n, ...))
}
table1(~ sex + racet + payort + norwage| SHUNTRND, data=data,droplevels=F, render=rndr,
       render.start=rndr.strat)

# ?tbl_summary()
# start prepare matrixmatrix
attach(cost)
treat = c()
treat[SHUNTRND=="1:MBTS"]=1
treat[!SHUNTRND=="1:MBTS"]=2

y = cum_cost #response
x = cbind(1,treat,stage) #design matrix for fixed effects
z = cbind(1,stage) #design matrix for random effects
#A median regression
# median_reg = QRLMM(y,x,z,blind_id,MaxIter = 500)
cost_qt <- QRLMM(y,x,z,blind_id,p=c(0.5,0.9),precision=0.001,MaxIter=100,M=10,cp=0.25,
      beta=NA,sigma=NA,Psi=NA,show.convergence=TRUE,CI=95)
summary(cost_qt)




y=cum_cost
x=cbind(1,SHUNTRND)
z=cbind(1,stage)
costqt <- QRLMM(y,x,z,blind_id,MaxIter = 10)

# Random intercept model
library(lqmm)

cost$y <- scale(cost$cum_cost, center = T, scale = T)
cost$stage <- as.factor(cost$stage)

# lqmm(y ~ aoi, random =  ~ 1, group = vpName, data = stackoverflow, tau = 0.15, control = lqmmControl(method = "df", UP_max_iter = 200))

cost.lqmm1 <- lqmm( y ~ SHUNTRND + stage , random = ~ stage, group = blind_id,
                    tau = c(0.5,0.75,0.9), data = cost,  control = lqmmControl(method = "df", UP_max_iter = 800))

coef(cost.lqmm1)

# Random slope model
cost.lqmm2 <- lqmm( y ~ SHUNTRND + stage + SHUNTRND*stage, random = ~ 1, group = blind_id,
                    tau = c(0.25,0.5,0.75,0.9), data = cost,  control = lqmmControl(method = "df", UP_max_iter = 800))
coef(cost.lqmm2)


# using quantreg::AIC.nlrq
qtrg <- nlrq(cum_cost~, data=cost, start, tau=0.5, 
             control, trace=FALSE,method="L-BFGS-B")
## S3 method for class 'nlrq'
summary(object, ...)
## S3 method for class 'summary.nlrq'
print(x, digits = max(5, .Options$digits - 2), ...)


Dat.nlrq <- nlrq(cum_cost ~ SHUNTRND + stage + sex +racet, data=cost, tau=0.5, trace=TRUE)
lines(1:25, predict(Dat.nlrq, newdata=list(x=1:25)), col=2)


# the 1st and 3rd quartiles regressions
Dat.nlrq <- nlrq(y ~ SSlogis(x, Asym, mid, scal), data=Dat, tau=0.25, trace=TRUE)
lines(1:25, predict(Dat.nlrq, newdata=list(x=1:25)), col=3)
Dat.nlrq <- nlrq(y ~ SSlogis(x, Asym, mid, scal), data=Dat, tau=0.75, trace=TRUE)
lines(1:25, predict(Dat.nlrq, newdata=list(x=1:25)), col=3)
# and finally "external envelopes" holding 95 percent of the data
Dat.nlrq <- nlrq(y ~ SSlogis(x, Asym, mid, scal), data=Dat, tau=0.025, trace=TRUE)
lines(1:25, predict(Dat.nlrq, newdata=list(x=1:25)), col=4)
Dat.nlrq <- nlrq(y ~ SSlogis(x, Asym, mid, scal), data=Dat, tau=0.975, trace=TRUE)
lines(1:25, predict(Dat.nlrq, newdata=list(x=1:25)), col=4)
leg <- c("least squares","median (0.5)","quartiles (0.25/0.75)",".95 band (0.025/0.975)")
legend(1, 12.5, legend=leg, lty=1, col=1:4)


## rq
 # library(quantreg)
 # data(engel)
 # attach(engel)
 # plot(income,foodexp,cex=.25,type="n",xlab="Household Income", ylab="Food Expenditure")
 # points(income,foodexp,cex=.5,col="blue")
 # abline(rq(foodexp~income,tau=.5),col="blue")
 # abline(lm(foodexp~income),lty=2,col="red") #the dreaded ols line
 # taus <- c(.05,.1,.25,.75,.90,.95)
 # 
 # for( i in 1:length(taus)){
 #  abline(rq(foodexp~income,tau=taus[i]),col="gray")
 # }
 
 library(quantreg)
 plot(stage,cum_cost,cex=.25,type="n",xlab="Stage", ylab="Cost")
 points(stage,cum_cost,cex=.5,col="blue")
 abline(rq(cum_cost~stage,tau=.5),col="blue")
 abline(lm(cum_cost~stage),lty=2,col="red") #the dreaded ols line
 taus <- c(.05,.1,.25,.75,.90,.95)
 
 for( i in 1:length(taus)){
   abline(rq(foodexp~income,tau=taus[i]),col="gray")
 }
 

