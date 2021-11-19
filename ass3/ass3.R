library(dplyr)
library(ggplot2)
library(nlme)
library(numbers)
library(table1)

###################### Getting the data

d = read.csv("marfan.csv")

dm = dplyr::filter(d, d$sexe == 1)
df = dplyr::filter(d, d$sexe == 0)

d$sexe = recode_factor(d$sexe, `1` = "male", `0` = "female")

######################## 
plot(d$age, d$diameter)

####################### Q1 #################
dunique = unique(select(d, c("patnr", "sexe")))
dunique$age = d$age[d$metingnr == 1]

label(dunique$age) = "Age at first measurement [years]"
dunique$sexe = factor(dunique$sexe, levels=c("male","female"), labels=c("Male", "Female"))
numMeas = c()
ageL = c()
diffDia = c()
for (i in 1:159) {
  dat = dplyr::filter(d, d$patnr == i)
  numMeas[i] = dim(dat)[1]
  ageL[i] = dat$age[order(-dat$age)[1]]
  diffDia[i] = dat$diameter[order(-dat$metingnr)[1]] - dat$diameter[order(dat$metingnr)[1]]
}
dunique$numMeas = numMeas
dunique$ageLast = ageL
dunique$diffDiameter = diffDia
label(dunique$ageLast) = "Age at last measurement [years]"
label(dunique$numMeas) = "Number of aorta diameter measurements"
label(dunique$diffDiameter) = "Δ Ø aorta between first and last measurement [mm]"

rndr <- function(x, name, ...) {
  if (!is.numeric(x)) return(render.categorical.default(x))
  what <- switch(name,
                 age = "Median [Min, Max]",
                 ageLast = "Median [Min, Max]",
                 numMeas  = "Median [Min, Max]",
                 diffDiameter = "Median [Min, Max]",
                 diameter = "Median [Min, Max]")
  parse.abbrev.render.code(c("", what))(x)
}

table1(~ age + ageLast + numMeas + diffDiameter | sexe, data=dunique, render = rndr)

dd = d
dd$age = floor(dd$age/10)*10
dd$ageF = factor(dd$age, levels=c(10, 20, 30, 40, 50, 60, 70), labels=c("10-19", "20-29", "30-39", "40-49", "50-79", "50-79", "50-79"))

label(dd$ageF) = "Age at measurement [years]"
label(dd$diameter) = "Aorta Ø [mm]"
strata = c(list(Total = dd), split(dd, dd$ageF))
labels = list(variables = list(sexe="Sexe", diameter="Aorta Ø [mm]"), groups = list("", "Age at measurement [years]"))
table1(strata, labels, data=dd, render = rndr, groupspan=c(1, 7))

####################### Q2 ################# (qq plots van residuals, average = 0)

mod1 = lme(diameter ~ age + sexe, random=(~1|patnr), data=d, method = "ML")
mod2 = lme(diameter ~ age + sexe, random=(~0 + age|patnr), data=d, method = "ML")
mod3 = lme(diameter ~ age + sexe, random=(~1 + age|patnr), data=d, method = "ML")
a = anova(mod1, mod2, mod3) # kan niet voor REML --> method = "ML" (alleen voor anova)
a
write.csv(a, "anova.csv")

mod = lme(diameter ~ age + sexe, random=(~1 + age|patnr), data=d, method = "REML")
summary(mod)

par(mfrow = c(3,2))

# scatter of residuals
plot(d$age, resid(mod), main = "Residuals of Model Fitting", xlab = "Age [years]", ylab = "Residual")
abline(h=0, lty=2)

# QQ of fixed coeff
qqnorm(resid(mod), main = "Q-Q Plot of the Residuals of the Coefficients")
qqline(resid(mod))

# hist + QQ of random intercepts
hist(mod$coefficients$random$patnr[,1], main = "Distribution of the Random Effects Intercept", xlab = "Random Effect", breaks = 12)
qqnorm(mod$coefficients$random$patnr[,1], main = "Q-Q Plot of Random Effects Intercept Coefficients")
qqline(mod$coefficients$random$patnr[,1])

# hist + QQ of random slopes
hist(mod$coefficients$random$patnr[,2], main = "Distribution of the Random Effects Slope", xlab = "Random Effect", breaks = 12)
qqnorm(mod$coefficients$random$patnr[,2], main = "Q-Q Plot of the Random Effects Slope Coefficients")
qqline(mod$coefficients$random$patnr[,2])

####################### Q3 #################
sexe = c(rep('male', 21), rep('female', 21))

newdata = data.frame(age = rep(seq(20, 40, 1), 2), patnr = 0:41, sexe = sexe)
designmatrix = model.matrix(~ age + sexe, newdata, contrasts.arg = list(sexe = contr.treatment(c("male", "female"), base = 2)))
predvar = diag(designmatrix %*% vcov(mod) %*% t(designmatrix))
newdata$SE = sqrt(predvar)
newdata$diameter = designmatrix %*% mod$coefficients$fixed

newdata$seLower = newdata$diameter - 1.96 * newdata$SE
newdata$seUpper = newdata$diameter + 1.96 * newdata$SE

#p + geom_line(aes(x = newdata$age, y = newdata$diameter, colour = newdata$sexe, group = newdata$sexe)) + geom_line(aes(x = newdata$age, y = newdata$seLower, group = newdata$sexe)) + geom_line(aes(x = newdata$age, y = newdata$seUpper, group = newdata$sexe))
p = ggplot() + geom_line(aes(x = newdata$age, y = newdata$diameter, colour = newdata$sexe, group = newdata$sexe)) + geom_ribbon(data = dplyr::filter(newdata, newdata$sexe == "male"), aes(x = age, ymin = seLower, ymax = seUpper), alpha = 0.2)+ geom_line(aes(x = newdata$age, y = newdata$diameter, colour = newdata$sexe, group = newdata$sexe)) + geom_ribbon(data = dplyr::filter(newdata, newdata$sexe == "female"), aes(x = age, ymin = seLower, ymax = seUpper), alpha = 0.2)
p + labs(colour = "Sexe", x = "Age [years]", y = "Aorta Ø [mm]", title = "Aorta diameter prediction for ages 20-40")

####################### Q4 #################
dex = select(d, c("patnr", "age", "sexe", "diameter"))
exd = data.frame()

for (ptnr in unique(dex$patnr)) {
  measurements = dplyr::filter(dex, dex$patnr == ptnr)
  meas = measurements[order(-measurements$age)[1],]
  exd = rbind(exd, data.frame(patnr = ptnr, age = meas$age + 1, sexe = meas$sexe, diameter = 0))
}

randomEff = mod$coefficients$random$patnr
randomEff = cbind(randomEff, rep(0, 159))

designmatrix = model.matrix(~ age + sexe, exd, contrasts.arg = list(sexe = contr.treatment(c("male", "female"), base = 1)))
predvar = diag(designmatrix %*% vcov(mod) %*% t(designmatrix))
SE = sqrt(predvar)

coeff = randomEff+rep(mod$coefficients$fixed,each=nrow(randomEff))

dia = c()
for (i in 1:159) {
  dia[i] = designmatrix[i,] %*% coeff[i,]
}
exd$diameter = dia

dex = rbind(dex, exd)

plotPatient = function(ptnr) {
  measurements = dplyr::filter(dex, dex$patnr == ptnr)
  m = measurements[order(-measurements$age)[1],]
  #ggplot() + geom_point(aes(x=measurements$age, y=measurements$diameter)) + geom_errorbar(aes(x=m$age, ymin=m$diameter - (1.96 * SE[ptnr]), ymax=m$diameter + (1.96 * SE[ptnr]), width=.1))
  ggplot() + geom_point(aes(x=measurements$age, y=measurements$diameter)) + geom_abline(intercept = coeff[ptnr, 1] + (m$sexe == "female") * coeff[ptnr, 3], slope = coeff[ptnr, 2]) + geom_errorbar(aes(x=m$age, ymin=m$diameter - (1.96 * SE[ptnr]), ymax=m$diameter + (1.96 * SE[ptnr]), width=.1))
}











