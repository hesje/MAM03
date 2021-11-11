library(rpart)
library(foreign)
library(dplyr)
library(psych)
library(survival)
library(ggExtra)
library(ggplot2)

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
describe(df)

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

bestShape = function(param, y = gstatus){
  
  mod1 = lm(y ~ log(param), data = df, na.action = na.exclude)
  print("log")
  print(AIC(mod1))
  
  mod2 = lm(y ~ sqrt(param), data = df, na.action = na.exclude)
  print("sqrt")
  print(AIC(mod2))
  
  AICs = c()
  for (n in 1:10) {
    mod = lm(y ~ poly(param, n), data = df, na.action = na.exclude)
    AICs[n] = AIC(mod)
  }
  AICs = data.frame(degree = 1:10, aic = AICs)
  
  print("polynomial")
  print(AICs)
  
  print(min(AICs$aic))
  plot(AICs)
}

params = c(acclft, gsurv, aantalre, creat, predias, prac)

mod = lm(gsurv ~ poly(creat, 6) + poly(acclft, 5) + poly(prac, 5) + log(predias) + gstatus + aantalre + uprotein, data = df, na.action = na.exclude)
plot(creat, mod$residuals)


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
plot(1:4, event)

######

min = min(log(df$creat))
stepsize = (max(log(df$creat)) - min) / 10

event = c()
for (c in 1:10) {
  cremin = (c - 1) * stepsize + min
  cremax = c * stepsize + min
  d = filter(df, log(df$creat) >= cremin & log(df$creat) < cremax)
  dEvent = filter(d, gstatus == 1)
  event[c] = dim(dEvent)[1]/dim(d)[1]
}
plot((1:10*stepsize + min), event)

####

min = min(df$predias)
stepsize = 10

event = c()
for (c in 1:10) {
  cremin = (c - 1) * stepsize + min
  cremax = c * stepsize + min
  d = filter(df, predias >= cremin & predias < cremax)
  dEvent = filter(d, gstatus == 1)
  event[c] = dim(dEvent)[1]/dim(d)[1]
}
plot((1:10*stepsize + min), event)

####

event = c()
for (c in 1:8) {
  d = filter(df, cregsh == c)
  dEvent = filter(d, gstatus == 1)
  event[c] = dim(dEvent)[1]/dim(d)[1]
}
plot(1:8, event)

####################### Creat

p2 <- df %>% ggplot(aes(x=creat, gsurv, color=gstatus)) + geom_point() + theme(legend.position="left")
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

plot(model, conf.int = F, xlab = "Time (Months)", ylab = "%Alive", main = "KM-Model", las = 1,
     col = c("green", "blue", "red"), lwd = 2)
legend(145,0.2, legend = c("uprotein 0", "uprotein 1", "uprotein 2"), lty = 1, lwd = 2,
       col = c("green", "blue", "red"), bty = "", cex = 0.7 )

# Do the LOG-RANK Test 
# H0 : survival in three groups is the same
# Ha : survival in three groups is not the same

df$gstatus <- as.factor(df$gstatus)
survdiff( Surv(df$gsurv,df$gstatus) ~ df$uprotein)












