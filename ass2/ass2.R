library(rpart)
library(foreign)
library(dplyr)
library(psych)
library(survival)
library(ggExtra)
library(ggplot2)
library(splines)

####################################################

# acclft   = age of the patient at transplantation (years)
# dgf      = delayed graft function (yes/no)
# gsurv    = number of months the patient was followed until rejection or censoring
# aantalre = number of acute rejection-treatments (0,1,2,3,......)
# creat    = serum creatinine (mumol/L) (marker of kidney-function)
# predias  = diastolic bloodpressure pre-transplantation
# prac     = panel reactive antibodies (marker of inflammation)
# uprotein = protein in urine (0,1,2) (marker of kidney function)
# cregsh   = number of HLA-genes matched between kidney-donor and acceptor

# gstatus  = graft-rejection (=1) or not (=0)

df <- read.spss("renaltx.sav", use.value.label=TRUE,to.data.frame=TRUE) # Importing the SPSS file
df$id = 1:698

quant = c("acclft", "gsurv", "creat", "predias")

event = filter(df, df$gstatus == 1)
noEvent = filter(df, df$gstatus == 0)
describe(event)
describe(noEvent)

cols <- c('gstatus', 'dgf', 'uprotein', 'cregsh')  #
df[cols] <- lapply(df[cols], as.factor)  # converting 'dgf', 'gstatus', 'uprotein' to factors

str(df) 
summary(df)
class(df)

#### Creatinine Level ##### <<< Skewed and a lot of outliers

summary(creat)
boxplot(creat)
boxplot(log(creat))
hist(log(creat))
plot(density(creat), main = "Creatinine Density Spread")

pie(table(uprotein), main = "protein in urine (0,1,2) (marker of kidney function)")
plot(predias)

########

reject <- filter(df, gstatus == "1")

plot(reject$uprotein)
scatterHist(df$gsurv, log(df$creat))
str(reject$gsurv)

bestShape = function(param){
  mod1 = coxph(Surv(df$gsurv, df$gstatus, type = "right") ~ log(param), data = df, id = df$id)
  print("log")
  print(AIC(mod1))
  
  mod2 = coxph(Surv(df$gsurv, df$gstatus, type = "right") ~ sqrt(param), data = df, id = df$id)
  print("sqrt")
  print(AIC(mod2))
  
  AICs = c()
  for (n in 1:10) {
    mod = coxph(Surv(df$gsurv, df$gstatus, type = "right") ~ poly(param, n), data = df, id = df$id)
    AICs[n] = AIC(mod)
  }
  AICs = data.frame(degree = 1:10, aic = AICs)
  
  print("polynomial")
  print(AICs)
  
  print(min(AICs$aic))
  plot(AICs)
  
  splineAICs = c()
  for (i in 1:12) {
    splineAICs[i] = AIC(coxph(Surv(df$gsurv, df$gstatus, type = "right") ~ ns(param, df=i), data = df, id = df$id))
  }
  splineAICs = data.frame(degree = 1:12, aic = splineAICs)
  print("splines")
  print(splineAICs)
  #plot(splineAICs)
}

############################# Cumulative Hazard
survival = c()
for (y in 1:150) {
  d = filter(df, df$gsurv < y & df$gstatus == 1)
  survival[y] = dim(d)[1]/698
}

plot(1:150, survival)

#############################

event = c()
for (re in 0:3) {
  d = filter(df, df$aantalre == re)
  dEvent = filter(d, gstatus == 1)
  event[re+1] = dim(dEvent)[1]/dim(d)[1]
}
plot(1:4, event, type = "b", xgap.axis = 1)

######

min = min(log(df$creat))
stepsize = (max(log(df$creat)) - min) / 10
logCreatinine = 1:10*stepsize + min

eventProportion = c()
for (c in 1:10) {
  cremin = (c - 1) * stepsize + min
  cremax = c * stepsize + min
  d = filter(df, log(df$creat) >= cremin & log(df$creat) < cremax)
  dEvent = filter(d, gstatus == 1)
  eventProportion[c] = dim(dEvent)[1]/dim(d)[1]
}
plot(logCreatinine, eventProportion, type = "b", main = "Graft Rejection Proportion for different levels of Serum Creatinine", xlab = "log(Creatinine)", ylab = "Graft Rejection Proportion")

####

min = min(df$predias)
stepsize = 20

event = c()
for (c in 1:4) {
  cremin = (c - 1) * stepsize + min
  cremax = c * stepsize + min
  d = filter(df, predias > cremin & predias <= cremax)
  dEvent = filter(d, gstatus == 1)
  event[c] = dim(dEvent)[1]/dim(d)[1]
}
plot((1:4*stepsize + (min - 10)), event, xlim = c(50, 130), main = "Graft Rejection Proportion versus Diastolic Blood Pressure", xlab = "Diastolic Blood Pressure Before Transplantation [mmHg]", ylab = "Graft Rejection Proportion", type = "b")

####

event = c()
for (c in 1:8) {
  d = filter(df, cregsh == c)
  dEvent = filter(d, gstatus == 1)
  event[c] = dim(dEvent)[1]/dim(d)[1]
}
plot(1:8, event)

####

steps = 10
min = min(df$acclft)
stepsize = (max(df$acclft) - min) / steps

event = c()
for (c in 1:steps) {
  cremin = (c - 1) * stepsize + min
  cremax = c * stepsize + min
  d = filter(df, acclft > cremin & acclft <= cremax)
  dEvent = filter(d, gstatus == 1)
  event[c] = dim(dEvent)[1]/dim(d)[1]
}
plot((1:steps*stepsize + (min - (stepsize/2))), event, main = "Graft Rejection Proportion versus age at Transplantation", xlab = "Age at transplantation [years]", ylab = "Graft Rejection Proportion", type = "b")

####

event = c()
for (c in 1:6) {
  d = filter(df, df$pracCat == c-1)
  dEvent = filter(d, gstatus == 1)
  event[c] = dim(dEvent)[1]/dim(d)[1]
}
plot((0:5), event)

####################### Creat

p2 <- df %>% ggplot(aes(x=creat, gsurv, color=gstatus, groupFill = F)) + geom_point() + theme(legend.position="left")
ggMarginal(p2, type="histogram", groupColour = TRUE, groupFill = TRUE)

#################### Acclft

p3 <- df %>%
  ggplot(aes(x=acclft, gsurv, color=gstatus))+
  geom_point() +
  theme(legend.position="left")
ggMarginal(p3, type="histogram", groupColour = TRUE, groupFill = TRUE)

####################   cregsh

p4 <- df %>%
  ggplot(aes(x=cregsh, gsurv, color=gstatus))+
  geom_point() +
  theme(legend.position="left")
ggMarginal(p4, type="histogram", groupColour = TRUE, groupFill = TRUE)

#################   Prac

p5 <- df %>%
  ggplot(aes(x=prac, gsurv, color=gstatus))+
  geom_point() +
  theme(legend.position="left")
ggMarginal(p5, type="histogram", groupColour = TRUE, groupFill = TRUE)

################   predias

p6 <- df %>%
  ggplot(aes(x=predias, gsurv, color=gstatus))+
  geom_point() +
  theme(legend.position="left")
ggMarginal(p6, type="histogram", groupColour = TRUE, groupFill = TRUE)

#################   uprotein

p7 <- df %>%
  ggplot(aes(x=uprotein, gsurv, color=gstatus))+
  geom_point() +
  theme(legend.position="left")
ggMarginal(p7, type="histogram", groupColour = TRUE, groupFill = TRUE)

#################   dgf

p8 <- df %>%
  ggplot(aes(x=dgf, gsurv, color=gstatus))+
  geom_point() +
  theme(legend.position="left")
ggMarginal(p8, type="histogram", groupColour = TRUE, groupFill = TRUE)

##################  

p9 <- df %>%
  ggplot(aes(x=aantalre, gsurv, color=gstatus))+
  geom_point() +
  theme(legend.position="left")
ggMarginal(p9, type="histogram", groupColour = TRUE, groupFill = TRUE)

################ survival analysis

model = survfit( Surv(df$gsurv,df$gstatus) ~ df$uprotein )
summary(model)
model

plot(model, conf.int = F, xlab = "Time (Months)", ylab = "Proportion With Event", main = "Kaplan-Meier Curve for the Urinary Protein", las = 1,
     col = c("green", "blue", "red"), lwd = 2)
legend(145,0.5, legend = c("uprotein 0", "uprotein 1", "uprotein 2"), lty = 1, lwd = 2,
       col = c("green", "blue", "red"), bty = "", cex = 0.7 )

model = survfit( Surv(df$gsurv,df$gstatus) ~ df$dgf )
plot(model, conf.int = F, xlab = "Time (Months)", ylab = "Proportion With Event", main = "Kaplan-Meier Curve for the Delayed Graft Function", las = 1,
     col = c("green", "red"), lwd = 2)
legend(125,0.1, legend = c("Not Delayed", "Delayed"), lty = 1, lwd = 2,
       col = c("green", "red"), bty = "", cex = 0.7 )

model = survfit( Surv(df$gsurv,df$gstatus) ~ df$aantalre )
plot(model, conf.int = F, xlab = "Time (Months)", ylab = "Proportion With Event", main = "Kaplan-Meier Curve for the Number of Rejection Treatments", las = 1,
     col = c("green", "orange", "red", "purple"), lwd = 2)
legend(0,0.5, legend = 0:3, lty = 1, lwd = 2,
       col = c("green", "orange", "red", "purple"), bty = "", cex = 0.7 )


df$pracCat = round((df$prac)^(1/3))
df$pracCat = 
model = survfit( Surv(df$gsurv,df$gstatus) ~ df$pracCat )
plot(model, conf.int = F, xlab = "Time (Months)", ylab = "Proportion With Event", main = "Kaplan-Meier Curve for the Panel Reactive Antibodies", las = 1,
     col = c("green", "orange", "red", "purple", "brown", "black"), lwd = 2)
legend(0,0.4, legend = 0:5, lty = 1, lwd = 2,
       col = c("green", "orange", "red", "purple", "brown", "black"), bty = "", cex = 0.7 )

df$hla = ceiling(df$cregsh / 2)

tr = c(1, 1, 1, 2, 3, 4, 4, 4)

df$hla2 = tr[df$cregsh]

model = survfit( Surv(df$gsurv,df$gstatus) ~ df$hla2)
plot(model, conf.int = F, xlab = "Time (Months)", ylab = "% With Event", main = "Kaplan-Meier Curve for the HLA genes", las = 1,
     col = c("green", "orange", "red", "purple"), lwd = 2)
legend(0,0.3, legend = c("1", "2", "3", "4"), lty = 1, lwd = 2,
       col = c("green", "orange", "red", "purple"), bty = "", cex = 0.7 )

# Do the LOG-RANK Test 
# H0 : survival in three groups is the same
# Ha : survival in three groups is not the same

df$gstatus = as.numeric(df$gstatus)

df$gstatus <- as.factor(df$gstatus)
survdiff( Surv(df$gsurv,df$gstatus) ~ df$uprotein)



params = c(acclft, gsurv, aantalre, creat, predias, prac)
df$gstatus = as.factor(as.numeric(df$gstatus) + 1)

mod = lm(gstatus ~ poly(creat, 3) + poly(acclft, 5) + pracCat + log(predias) + aantalre + uprotein, data = df, na.action = na.exclude)
mod$coefficients
preds = predict(mod, type = "response")
plot(df$creat, mod$residuals)

df$c = 1

cx = coxph( Surv(df$gsurv, df$gstatus, type = "right") ~ poly(df$creat, 5) + df$acclft + df$pracCat + df$predias + df$aantalre + df$uprotein + df$cregsh + df$dgf, df, id = df$id)
summary(cx)
ROCR::performance(predict(cx), "auc")

coxx = coxph( Surv(df$gsurv, df$gstatus, type = "right") ~ df$acclft, df, id = df$id)
plot(cox.zph(coxx))

predict(cx, newdata = df, type = "risk")

summary = summary(cx)
summary$coefficients

write.csv(summary$coefficients, file = "coef.csv")


