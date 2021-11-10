
d = read.sav()


library(rpart)
library(foreign)
library(dplyr)

df <- read.spss("renaltx.sav", use.value.label=TRUE,to.data.frame=TRUE)

head(df)

str(df)

hist(df$creat, col = 'blue', main = "Histpgram of Creatinine Levels",
     xlab = "Creatinine", ylab = "Frequency" )

help(hist)
table(df$gstatus)

attach(df)

df_subset <- subset(df, select = c(creat,prac))
boxplot(df_subset, col = c('blue','red'))

table(gstatus)

library(psych)
describe(df)

pairs.panels(df)
hist(cregsh)

cols <- c('gstatus', 'dgf', 'uprotein')
df[cols] <- lapply(df[cols], as.factor)

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

table(prac)

table(cregsh)
hist(cregsh)

cols <- c('gstatus', 'dgf', 'uprotein', 'cregsh')  #
df[cols] <- lapply(df[cols], as.factor)  # converting 'dgf', 'gstatus', 'uprotein' to factors

str(df) 
summary(df)
class(df)

#### Creatinine Level ##### <<< Skewed and a lot of outliers

summary(creat)
boxplot(creat)
hist(creat)
plot(density(creat), main = "Creatinine Density Spread")

pie(table(uprotein), main = "protein in urine (0,1,2) (marker of kidney function)")
plot(predias)

########

reject <- filter(df, gstatus == "1")

table(reject$cregsh)

plot(reject$uprotein)
scatterHist(reject$gsurv)
