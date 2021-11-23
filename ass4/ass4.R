library(nlme)
library(survival)
library(JM)
library(tidyverse)
library(ggplot2)
library(survminer)

############################################## Data
load("kidney transplant.RData")

d$sexDonEqual = d$sexdon == d$sex_pat
dlong$sexDonEqual = dlong$sexdon == dlong$sex_pat

dcox = dlong
dcox$start = dcox$years
dcox$end[dcox$years == 0.25] = 0.5
dcox$end[dcox$years == 0.5] = 1
dcox$end[dcox$years == 1] = 2
dcox$end[dcox$years == 2] = 5
dcox$end[dcox$years == 5] = 10
dcox$end[dcox$years == 10] = 15
dcox$end[dcox$years == 15] = 20



for (i in unique(dcox$ID)) {
  dcox$end[dcox$ID == i & dcox$start < dcox$time_to_death & dcox$end > dcox$time_to_death] = dcox$time_to_death[dcox$ID == i][1]
}

for (j in 1:dim(dcox)[1]) {
  dd = dcox[j,]
  dcox$gEvent[j] = dd$time_to_graft_failure > dd$start & dd$time_to_graft_failure <= dd$end & dd$stat_gra == "Graft loss"
  dcox$dEvent[j] = dd$time_to_death > dd$start & dd$time_to_death <= dd$end & dd$stat_pat == "overleden"
}

dcox = subset(dcox, start < time_to_death)

gcox = subset(dcox, start < time_to_graft_failure)

for (i in unique(gcox$ID)) {
  gcox$end[gcox$ID == i & gcox$start < gcox$time_to_graft_failure & gcox$end > gcox$time_to_graft_failure] = gcox$time_to_graft_failure[gcox$ID == i][1]
}

######################################################## models
## COX Baseline models
cbld = coxph(Surv(time_to_death, stat_pat) ~ sex_pat + sexDonEqual + bmi + age_at_tx + type_dia + retrans + duur_dia, data=d, x=TRUE, model = TRUE, na.action = na.exclude, id = d$ID)
cblg = coxph(Surv(time_to_graft_failure, stat_gra) ~ sex_pat + sexDonEqual + bmi + age_at_tx + type_dia + retrans + duur_dia, data=d, x=TRUE, model = TRUE, na.action = na.exclude, id = d$ID)

ggforest(cbld, data = d)

## LME
l = lme(gfr ~ sex_pat + age_at_tx + bmi + retrans + years + sexdon + type_dia , random = ~1 + years|ID, data = dcox, na.action = na.exclude)

preds = predict(l, level = 1, data = dcox)

plotPredLme = function(pt){
  plot(dcox$years[dcox$ID==pt], dcox$gfr[dcox$ID==pt])
  lines(dcox$years[dcox$ID==pt], preds[dcox$ID==pt])
}

## COX for time to death
cd = coxph(Surv(start, end, dEvent) ~ sex_pat + bmi + age_at_tx + type_dia + sexDonEqual + preds + cluster(dcox$ID), data=dcox, x=TRUE, model = TRUE, na.action = na.exclude)
ggforest(cd, data = dcox)

## COX for time to graft failure
cg = coxph(Surv(start, end, gEvent) ~ sex_pat + bmi + gfr + creat + map + age_at_tx + type_dia + sexdon + cluster(gcox$ID), data=gcox, x=TRUE, id=gcox$ID, model = TRUE, na.action = na.exclude)
ggforest(cg, data = gcox)

mod = survfit(Surv(start, end, dEvent) ~ type_dia, data = dcox)
plot(mod)

plot(cox.zph(cd))
cox.zph(cd)

## Joint model
years = dcox$years

j = jointModel(l, cd, timeVar="years")

predss = predict(j, level = 1, newdata = dlongg)

plotPrediction = function(patientId){
  plot(dlongg$years[dlongg$ID==patientId], dlongg$gfr[dlongg$ID==patientId])
  lines(dlongg$years[dlongg$ID==patientId], predss[dlongg$ID==patientId])
}
