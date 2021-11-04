

costefficacydata$event = as.factor(costefficacydata$event)


aliveTrtOne <- sum(costefficacydata$trt == 1 & costefficacydata$event ==1)
aliveTrtTwo <- sum(costefficacydata$trt == 2 & costefficacydata$event == 1)
deadTrtOne <- sum(costefficacydata$trt == 1 & costefficacydata$event == 0)
deadTrtTwo <- sum(costefficacydata$trt == 2 & costefficacydata$event == 0)

p1 <- aliveTrtOne / 409
p2 <- aliveTrtTwo / 409

sep1 <- sqrt(p1*(1-p1)/409)

CIp1Plus <- p1 + 1.96 * sep1
CIp1Min <- p1 - 1.96 * sep1

# 95% CI P1 = (0.423 - 0.329)

sep2 <- sqrt(p2*(1-p2)/409)

CIp2Plus <- p2 + 1.96 * sep2
CIp2Min <- p2 - 1.96 * sep2

# 95% CI P2 = (0.36 - 0.27)

### treatment difference = efficacy E
E <- p1 - p2

seE <- sqrt(sep1^2 + sep2^2)

# 95% CI for seE
CIeffPlus <- E + 1.96 * seE 
CieffMin  <- E - 1.96 * seE

######################
#                    #
#      COST          #
#                    #
######################

costTrtone <- subset(costefficacydata$costs, costefficacydata$trt == 1)
costTrtTwo <- subset(costefficacydata$costs, costefficacydata$trt == 2)

meanCostone <- mean(costTrtone)
meanCosttwo <- mean(costTrtTwo)

SDone <- sd(costTrtone)  # standard deviation cost one
SDTwo <- sd(costTrtTwo)  # standard deviation cost two

SEmeanOne <- SDone / sqrt(length(costTrtone))
SEmeanTwo <- SDTwo / sqrt(length(costTrtTwo))

CIsemean1Plus <- meanCostone + 1.96 * (SDone / 206) 
CIsemean1Min <- meanCostone - 1.96 * (SDone / 206) 

# CI(E) = (1024 - 966)

CIsemean2Plus <- meanCosttwo + 1.96 * (SDTwo / 203)
CIsemean2Min <- meanCosttwo - 1.96 * (SDTwo / 203)

# CI(E) = (40 - 38)

### Cost difference = C
C = meanCostone - meanCosttwo

SE_C <- sqrt(SEmeanOne^2 + SEmeanTwo^2)

## CE-ratio

CE_ratio = C/E # 15638 (euro per)
lnCE_ratio = log(C/E)
##############

p1_minus_p2 = rnorm(409, 0.0611, 0.0331)
meancosts1_meancosts2 = rnorm(409, 955.9, 211)
plot(p1_minus_p2, meancosts1_meancosts2, ylab="difference of the mean costs", xlab="difference of the proportions alive")

CE_Ratios <- meancosts1_meancosts2 / p1_minus_p2

hist(CE_Ratios)

hist(log(CE_Ratios))

help(boot)

fn <- function(costefficacydata,e){
  
  
  
}



boot(costefficacydata, fn, R=1000)
