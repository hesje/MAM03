library(dplyr)
library(ggplot2)
library(nlme)
library(numbers)

###################### Getting the data

d = read.csv("marfan.csv")

dm = dplyr::filter(d, d$sexe == 1)
df = dplyr::filter(d, d$sexe == 0)

d$sexe = recode_factor(d$sexe, `1` = "male", `0` = "female")

######################## 
plot(d$age, d$diameter)

####################### Q1 #################


####################### Q2 ################# (qq plots van residuals, average = 0)

mod1 = lme(diameter ~ age + sexe, random=(~1|patnr), data=d, method = "ML")
mod2 = lme(diameter ~ age + sexe, random=(~0 + age|patnr), data=d, method = "ML")
mod3 = lme(diameter ~ age + sexe, random=(~1 + age|patnr), data=d, method = "ML")
anova(mod1, mod2, mod3) # kan niet voor REML --> method = "ML" (alleen voor anova)

mod = lme(diameter ~ age + sexe, random=(~1 + age|patnr), data=d, method = "REML")
summary(mod)

par(mfrow = c(2,2))
plot(d$age, resid(mod), main = "Residuals of Model Fitting", xlab = "Age [years]", ylab = "Residual")
abline(h=0, lty=2)

qqnorm(resid(mod), main = "Q-Q Plot of the Residuals of the Coefficients")
qqline(resid(mod))

hist(mod$coefficients$random$patnr, main = "Distribution of the Random Effect per Patient", xlab = "Random Effect")

qqnorm(mod$coefficients$random$patnr, main = "Q-Q Plot of the Residuals of the Random Effects")
qqline(mod$coefficients$random$patnr)

####################### Q3 #################
sexe = c(rep('male', 21), rep('female', 21))

newdata = data.frame(age = rep(seq(20, 40, 1), 2), patnr = 0:41, sexe = sexe)
designmatrix = model.matrix(~ age + sexe, newdata, contrasts.arg = list(sexe = contr.treatment(c("male", "female"), base = 2)))
predvar = diag(designmatrix %*% vcov(mod) %*% t(designmatrix))
newdata$SE = sqrt(predvar)
newdata$diameter = designmatrix %*% mod$coefficients$fixed

newdata$seLower = newdata$diameter - 1.96 * newdata$SE
newdata$seUpper = newdata$diameter + 1.96 * newdata$SE

p = ggplot()
p + geom_line(aes(x = newdata$age, y = newdata$diameter, colour = newdata$sexe, group = newdata$sexe)) + geom_line(aes(x = newdata$age, y = newdata$seLower, group = newdata$sexe)) + geom_line(aes(x = newdata$age, y = newdata$seUpper, group = newdata$sexe))
p + geom_line(aes(x = newdata$age, y = newdata$diameter, colour = newdata$sexe, group = newdata$sexe)) + geom_ribbon(data = dplyr::filter(newdata, newdata$sexe == "male"), aes(x = age, ymin = seLower, ymax = seUpper), alpha = 0.2)+ geom_line(aes(x = newdata$age, y = newdata$diameter, colour = newdata$sexe, group = newdata$sexe)) + geom_ribbon(data = dplyr::filter(newdata, newdata$sexe == "female"), aes(x = age, ymin = seLower, ymax = seUpper), alpha = 0.2)

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
# predvar = diag(designmatrix %*% vcov(mod) %*% t(designmatrix))
# newdata$SE = sqrt(predvar)

coeff = randomEff+rep(mod$coefficients$fixed,each=nrow(randomEff))

dia = c()
for (i in 1:159) {
  dia[i] = designmatrix[i,] %*% coeff[i,]
}
exd$diameter = dia

dex = rbind(dex, exd)

plotPatient = function(ptnr) {
  measurements = dplyr::filter(dex, dex$patnr == ptnr)
  plot(measurements$age, measurements$diameter)
}











