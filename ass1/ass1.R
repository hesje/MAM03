library(dplyr)
library(boot)

CEdata = read.csv(file="costefficacydata.csv", header=TRUE, sep=",")
CEdata$event = as.factor(CEdata$event)

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
  
  sep1 <- sqrt(p1*(1-p1)/ dim(data)[1])
  
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
  
  meanCostone <- mean(trtOne$costs)
  meanCosttwo <- mean(trtTwo$costs)
  
  SDone <- sd(trtOne$costs)  # standard deviation cost one
  SDTwo <- sd(trtTwo$costs)  # standard deviation cost two
  
  SEmeanOne <- SDone / sqrt(dim(trtOne)[1])
  SEmeanTwo <- SDTwo / sqrt(dim(trtTwo)[1])
  
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

  return(CE_ratio)
}


res <- boot(CEdata, CE, R=1000, stype = "i")
boot.ci(res, type = "all")

hist(log(res$t))

sep1 <- sqrt(p1*(1-p1)/ dim(data)[1])



