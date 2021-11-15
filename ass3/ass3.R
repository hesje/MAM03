library(dplyr)
library(ggplot2)
library(nlme)

d = read.csv("marfan.csv")

dm = dplyr::filter(d, d$sexe == 0)
df = dplyr::filter(d, d$sexe == 1)

d$sexe = as.factor(d$sexe)

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
mod3 = lme(diameter ~ age + sexe, random=(~1 + age|patnr), data=d, method = "REML")

anova(mod1, mod2, mod3)

x = 20:40
ym = 27.18 + x * 0.397
yf = 27.18 + x * 0.397 + 8.46

ggplot() + geom_line(aes(x = x, y = ym, colour = "red")) + geom_line(aes(x = x, y = yf, colour = "blue"))

mod4 = lme(diameter ~ age, random=(~1 + age|patnr), data=d, method = "REML")
newdata <- data.frame(age = 20:40, patnr = 1:21)
designmatix <- model.matrix(~age, newdata)
predvar <- diag(designmatix %*% vcov(mod4) %*% t(designmatix))
newdata$SE <- sqrt(predvar)


ggplot() + geom_line(aes(x = x, y = ym, colour = "red")) + geom_line(aes(x = x, y = yf, colour = "blue"))












