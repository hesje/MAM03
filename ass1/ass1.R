library(dplyr)
library(boot)

CEdata = read.csv(file="costefficacydata.csv", header=TRUE, sep=",")
CEdata$event = as.factor(CEdata$event)


trtOne = filter(CEdata, CEdata$trt == 1)
trtTwo = filter(CEdata, CEdata$trt == 2)

eventProp <- function(d){
  alive = sum(d$event == 1)
  return(alive/dim(d)[1])
}

CE <- function(data, i){
  
  data = data[i,]
  
  trtOne = filter(data, data$trt == 1)
  trtTwo = filter(data, data$trt == 2)
  
  p1 <- eventProp(trtOne)
  p2 <- eventProp(trtTwo)
  E <- p1 - p2
  
  C = mean(trtOne$costs) - mean(trtTwo$costs)

  return(C/E)
}

res <- boot(CEdata, CE, R=10000, stype = "i")
boot.ci(res, type = "all")

mean(res$t)
res$t0

sd = sd(res$t)
se = sd/sqrt(10000)

res$t0 + 1.96*se
res$t0 - 1.96*se

hist(log(res$t))
hist(res$t, breaks = 100)

dev.new()


n1 = dim(trtOne)[1]
n2 = dim(trtTwo)[1]
m1 = mean(trtOne$costs) 
m2 = mean(trtTwo$costs)
p1 = eventProp(trtOne)
p2 = eventProp(trtTwo)
s1 = sd(trtOne$costs)
s2 = sd(trtTwo$costs)

ratio = c()
for (j in 1:10000) {
  dd = data.frame(rbind(
    cbind(rbinom(n1, 1, p1), rnorm(n1, m1, s1), rep(1, n1)),
    cbind(rbinom(n2, 1, p2), rnorm(n2, m2, s2), rep(2, n2))
  ))
  names(dd) = c("event", "costs", "trt")
  ratio[j] = (mean(dd$costs[dd$trt==1]) - mean(dd$costs[dd$trt==2])) / (mean(dd$event[dd$trt==1]) - mean(dd$event[dd$trt==2]))
}

hist(log(ratio))
mean(ratio)












sep1 <- sqrt(p1*(1-p1)/ dim(data)[1])
SE_C <- sqrt(SEmeanOne^2 + SEmeanTwo^2)


SDone <- sd(trtOne$costs)  # standard deviation cost one
SDTwo <- sd(trtTwo$costs)  # standard deviation cost two

SEmeanOne <- SDone / sqrt(dim(trtOne)[1])
SEmeanTwo <- SDTwo / sqrt(dim(trtTwo)[1])

CIsemean1Plus <- mean(trtOne$costs) + 1.96 * SEmeanOne
CIsemean1Min <- mean(trtOne$costs) - 1.96 * SEmeanOne

# CI(E) = (1024 - 966)

CIsemean2Plus <- mean(trtTwo$costs) + 1.96 * SEmeanTwo
CIsemean2Min <- mean(trtTwo$costs) - 1.96 * SEmeanTwo

# CI(E) = (40 - 38)



seE <- sqrt(sep1^2 + sep2^2)

# 95% CI for seE
CIeffPlus <- E + 1.96 * seE 
CieffMin  <- E - 1.96 * seE



sep1 <- sqrt(p1*(1-p1)/ dim(CEdata)[1])

CIp1Plus <- p1 + 1.96 * sep1
CIp1Min <- p1 - 1.96 * sep1

# 95% CI P1 = (0.423 - 0.329)

sep2 <- sqrt(p2*(1-p2)/409)

CIp2Plus <- p2 + 1.96 * sep2
CIp2Min <- p2 - 1.96 * sep2

# 95% CI P2 = (0.36 - 0.27)