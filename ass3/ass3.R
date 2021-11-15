library(dplyr)
library(ggplot2)
library(nlme)
library(numbers)

###################### Getting the data

d = read.csv("marfan.csv")

dm = dplyr::filter(d, d$sexe == 0)
df = dplyr::filter(d, d$sexe == 1)

d$sexe = as.factor(d$sexe)

######################## Some random statistics shit and copy-paste from slides.
dia = c()
for (i in 1:22) {
  dia[i] = mean(dplyr::filter(d, metingnr == i)$diameter)
}
plot(1:22, dia)

plot(d$metingnr, d$diameter)

p = ggplot()
p + geom_boxplot(aes(x = dm$metingnr, y = dm$diameter, group = dm$metingnr))

mod = lme(diameter ~ age, random=(~1 + age|patnr), data=d, method = "REML")
summary(mod)

par(mfrow = c(2,2))
plot(d$age, resid(mod))
abline(h=0, lty=2)

lines(unique(d$age), tapply(resid(mod), d$age, mean), col=2, lwd=3)

qqnorm(resid(mod))
qqline(resid(mod))

qqnorm(mod$coefficients$random$patnr)
qqline(mod$coefficients$random$patnr)

mod1 = lme(diameter ~ age + sexe, random=(~1|patnr), data=d, method = "REML")
mod2 = lme(diameter ~ age + sexe, random=(~0 + age|patnr), data=d, method = "REML")
mod3 = lme(diameter ~ age + sexe + (age * sexe), random=(~1 + age|patnr), data=d, method = "REML")

anova(mod1, mod2, mod3)

mod4 = lme(diameter ~ age, random=(~1 + age|patnr), data=d, method = "REML")
newdata <- data.frame(age = 20:40, patnr = 1:21)
designmatix <- model.matrix(~age, newdata)
predvar <- diag(designmatix %*% vcov(mod4) %*% t(designmatix))
newdata$SE <- sqrt(predvar)

ggplot() + geom_line(aes(x = x, y = ym, colour = "red")) + geom_line(aes(x = x, y = yf, colour = "blue"))

####################### Q2 of assignment I think? #################
mod5 = lme(diameter ~ age + sexe + (age * sexe), random=(~1 + age|patnr), data=d, method = "REML")

####################### Q3 of assignment I think? #################
newdata <- data.frame(age = seq(20, 40, 0.5), patnr = 0:40, sexe = mod(0:40, 2))
designmatix <- model.matrix(~ age + sexe + (age * sexe), newdata)
predvar <- diag(designmatix %*% vcov(mod5) %*% t(designmatix))
newdata$SE <- sqrt(predvar)
newdata$diameter = designmatix %*% mod5$coefficients$fixed


newdata$sexe = recode_factor(newdata$sexe, `1` = "male", `0` = "female")

newdata$seLower = newdata$diameter - 1.96 * newdata$SE
newdata$seUpper = newdata$diameter + 1.96 * newdata$SE

p = ggplot()
p + geom_line(aes(x = newdata$age, y = newdata$diameter, colour = newdata$sexe, group = newdata$sexe)) + geom_ribbon(data = newdata, aes(x = age, ymin = seLower, ymax = seUpper), alpha = 0.2, group = newdata$sexe)
############################








