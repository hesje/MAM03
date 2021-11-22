library(nlme)
library(survival)
library(JM)
library(tidyverse)

load("kidney transplant.RData")

d$years = d$time_to_death
dlongg = dlong %>% drop_na(gfr)

l = lme(gfr ~ sex_pat + age_at_tx + bmi + retrans + years + sexdon, random = ~1 + years|ID, data = dlongg, na.action = na.exclude)

preds = predict(l, level = 1, data = dlongg)

plot(dlongg$years[dlongg$ID==2], dlongg$gfr[dlongg$ID==2])
lines(dlongg$years[dlongg$ID==2], preds[dlongg$ID==2])

cox = coxph(Surv(d$time_to_death, d$stat_pat) ~ d$sex_pat + d$bmi + d$gfr1 + d$screat1 + d$map1 + cluster(d$ID), data=d, x=TRUE, id=d$ID, model = TRUE, na.action = na.pass)

years = dlongg$years

j = jointModel(l, cox, timeVar="years")

predss = predict(j, level = 1, newdata = dlongg)

plotPrediction = function(patientId){
  plot(dlongg$years[dlongg$ID==patientId], dlongg$gfr[dlongg$ID==patientId])
  lines(dlongg$years[dlongg$ID==patientId], predss[dlongg$ID==patientId])
}
