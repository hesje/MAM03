library(nlme)
library(survival)
library(JM)
library(tidyverse)
library(ggplot2)
library(survminer)
library(gtsummary)
library(ggfortify)
library(table1)
library(ROCR)
library(splines)

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


########################### Space assignment 4 #############################

# (!) Baseline Cox model: 

mod1 = coxph(Surv(time_to_death,as.numeric(stat_pat)-1)~1,data=d)
mod2graft = coxph(Surv(time_to_graft_failure ,as.numeric(stat_gra)-1)~+age_at_tx+agedon+type_dia+
                    duur_dia,data=d,x=TRUE)


summary(mod2graft)

modtest = survfit(Surv(time_to_death,as.numeric(stat_pat)-1)~1,data=d)
modtest2 = survfit(Surv(time_to_graft_failure,as.numeric(stat_gra)-1)~1,data=d)

summary(modtest)

fit <- list(PFS = modtest, OS = modtest2)

ggsurvplot_combine(fit, d, risk.table = TRUE,                  # Add risk table
                   conf.int = TRUE,                    # Add confidence interval
                   conf.int.style = "step",            # CI style, use "step" or "ribbon"
                   censor = FALSE,                     # Remove censor points
                   tables.theme = theme_cleantable(),  # Clean risk table
                   palette = "jco")

plot(table(d$stat_pat, d$stat_gra))


model1death = coxph(Surv(time_to_death,as.numeric(stat_pat)-1)~+age_at_tx+agedon+type_dia+ retrans +bmi+
                         duur_dia,data=d,x=TRUE)

model1graft = coxph(Surv(time_to_graft_failure,as.numeric(stat_gra)-1)~+age_at_tx+agedon+type_dia+ retrans +bmi+
                      duur_dia,data=d,x=TRUE)

model1death


install.packages("gtsummary")
install.packages("ggfortify")
install.packages("table1")

############################### tables descriptive ##########
t1 = tbl_regression(model1death)
t2 = tbl_regression(model1graft)

t2

ggforest(model1death, data = d)

tbl_merge_ex1 <-
  tbl_merge(
    tbls = list(t1, t2),
    tab_spanner = c("**Time to Death**", "**Time to Graft loss**")
  )

tbl_merge_ex1

aa_fit <-aareg(Surv(time_to_death,as.numeric(stat_pat)-1)~+age_at_tx+agedon+type_dia+ retrans +bmi+
                 duur_dia,data=d,x=TRUE)
aa_fit

autoplot(aa_fit)

table1::label(d$time_to_death) <- "Time to event (Death/Censored)"
table1::label(d$agedon) <- "Age donnor"
table1::label(d$age_at_tx) <- "Age at Treatment"
table1::label(d$sex_pat) <- "Sex Patient"
table1::label(d$stat_gra) <- "Graft Status"
table1::label(d$duur_dia) <- "Duration of dialysis"
table1::label(d$type_dia) <- "Type of Dialysis"
table1::label(d$retrans) <- "First-re-Transplantation"



table1::table1(~time_to_death + type_dia + duur_dia + stat_gra + retrans + sex_pat + age_at_tx + agedon | stat_pat, data = d)

#################### Assumption Cox model #######

### Proportionality #### 

par(mfrow = c(2,1))
plot(cox.zph(model1death))

testprop = cox.zph(model1death)

testprop

ggcoxzph(testprop)

######## Linearity ##########

ggcoxfunctional(Surv(time_to_death, as.numeric(stat_pat)-1) ~ log(age_at_tx) + (agedon) + sqrt(duur_dia) , data = d)

linearmodel = coxph(formula = Surv(time_to_death, as.numeric(stat_pat) - 1,) ~ 
        ns(age_at_tx, df = 4) + ns(agedon, df = 4) + ns(duur_dia, df = 4), 
      data = d, x = TRUE)

linearmodel2 = coxph(formula = Surv(time_to_death, as.numeric(stat_pat) - 1,) ~ 
                      (age_at_tx) + poly(agedon, 4) + (duur_dia), 
                    data = d, x = TRUE)


linearmodel3 = coxph(formula = Surv(time_to_death, as.numeric(stat_pat) - 1,) ~ 
                       pspline(age_at_tx) + pspline(agedon) + pspline(duur_dia), 
                     data = d, x = TRUE)

(AIC(linearmodel))
(AIC(linearmodel2))
(AIC(linearmodel3))


(AIC(model1death))






linearmodel2 = coxph(formula = Surv(time_to_death, as.numeric(stat_pat) - 1,) ~ 
                       pspline(duur_dia), 
                    data = d, x = TRUE)


summary(coxph(formula = Surv(time_to_death, as.numeric(stat_pat) - 1,) ~ 
                pspline(age_at_tx) + pspline(agedon) + pspline(duur_dia), 
              data = d, x = TRUE))

termplot(linearmodel3, term=2, se=TRUE, col.term=1, col.se=1)

summary(linearmodel)


AIC(model1death)
AIC(linearmodel)

linearmodel2 = coxph(formula = Surv(time_to_death, as.numeric(stat_pat) - 1,) ~ 
                       spline(duur_dia), 
                     data = d, x = TRUE)




############# repeatedly measured GFR ######


dlong1=subset(dlong,is.na(gfr)==FALSE)
dim(dlong1)

dlong1=subset(dlong1,(dlong1$years>dlong1$time_to_death)==FALSE)
dim(dlong1)

#### copy-paste code###
dlong1$gfr[dlong1$gfr > 200]= NA
dlong1=subset(dlong,is.na(gfr)==FALSE)


uniekepats=unique(dlong1$ID)
i=1
plot(dlong1$years[dlong1$ID==uniekepats[i]],dlong1$gfr[dlong1$ID==uniekepats[i]],type="l",
     xlab="years since Tx",ylab="GFR",   xlim=c(0,max(dlong1$years,na.rm=T)), 
     ylim=c(min(dlong1$gfr,na.rm=T),max(dlong1$gfr,na.rm=T)),
     col=mean(as.numeric(dlong1$stat_pat[dlong1$ID==uniekepats[i]],na.rm=T)),lwd=1)
for (i in 2:length(uniekepats)) {
  lines(dlong1$years[dlong1$ID==uniekepats[i]],dlong1$gfr[dlong1$ID==uniekepats[i]],
        col=mean(as.numeric(dlong1$stat_pat[dlong1$ID==uniekepats[i]],na.rm=T)),lwd=1)
}

########### lme models ###########

p1=lme(gfr~ns(years,df=3),random=~1|ID,data=dlong1,method="ML")
p2=lme(gfr~ns(years,df=3),random=~1+years|ID,data=dlong1,method="ML")
p3=lme(gfr~ns(years,df=3),random=~1+ns(years,df=2)|ID,data=dlong1,method="ML")
anova(p1,p2,p3)

p3a=lme(gfr~ns(years,df=2),random=~1+ns(years,df=2)|ID,data=dlong1,method="ML")
p3a1= lme(gfr~years,random=~1+ns(years,df=2)|ID,data=dlong1,method="ML")
anova(p3,p3a,p3a1)

p3a=lme(gfr~ns(years,df=2),random=~1+ns(years,df=2)|ID,data=dlong1,method="REML")

p3b=lme(gfr~ns(years,df=2)+sex_pat+age_at_tx+bmi+sexdon+agedon+type_dia+duur_dia,
        random=~1+ns(years,df=2)|ID,data=dlong1,method="REML")

summary(p3a)
summary(p3b)

predgfr = predict(p3a,newdata= data.frame(years=c(0,0.25,0.5,1,2,5,10,15)), level = 0)

lines(c(0,0.25,0.5,1,2,5,10,15),predgfr,col=7,lwd=5)







##### Test ######


Testlandmark = coxph(Surv(time_to_death,as.numeric(stat_pat)-1)~ gfr + map + creat +sex_pat +age_at_tx + bmi + sexdon +agedon + type_dia + duur_dia + retrans + strata(years) + cluster(ID), data = dlong)

summary(Testlandmark)





