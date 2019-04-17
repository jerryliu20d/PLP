############################Analysis of coal data##########################################
rm(list=ls())
set.seed(123)
require(R2jags)
setwd("C:/Users/jerry/Dropbox/paper_PLP/coal_Data")
for(i in 0:3){
  nchange <- i
  if(nchange == 0){
    var_names <- c("a1","b1")
    jags.data <- c("t1","n1","TT","zeros","n")
    jags.par <- c("alpha[1]","beta[1]")
    model.jags <- function()
    {
      for(i in 1:n){
        for(j in 1:n1[i])
        {
          phi[i,j] <- -l[i,j]
          zeros[i,j] ~ dpois(phi[i,j])
          
          k[i,j] <- 1
          llambda[i,j] <- log(alpha[k[i,j]])-log(beta[k[i,j]])+(alpha[k[i,j]]-1)*(log(t1[i,j])-log(beta[k[i,j]]))
          
          l[i,j] <- llambda[i,j]+mT[i]/n1[i]
        }
        
        mT[i] <- -pow(TT[i]/beta[1],alpha[1])
      }
      alpha[1] ~ dunif(0,100)
      beta[1] ~ dunif(0,80)
    }
  }else if(nchange == 1){
    var_names <- c("a1","a2","b1","b2","tau")
    jags.data <- c("t1","n1","TT","zeros","n")
    jags.par <- c("alpha[1]","alpha[2]","beta[1]","beta[2]","tau")
    model.jags <- function()
    {
      for(i in 1:n){
        for(j in 1:n1[i])
        {
          phi[i,j] <- -l[i,j]
          zeros[i,j] ~ dpois(phi[i,j])
          
          k[i,j] <- 2-step(tau-t1[i,j])
          llambda[i,j] <- log(alpha[k[i,j]])-log(beta[k[i,j]])+(alpha[k[i,j]]-1)*(log(t1[i,j])-log(beta[k[i,j]]))
          
          l[i,j] <- llambda[i,j]+mT[i]/n1[i]
        }
        
        mT[i] <- -pow(tau/beta[1],alpha[1])-pow(TT[i]/beta[2],alpha[2])+pow(tau/beta[2],alpha[2])
      }
      tau ~ dunif(0,max(TT))
      alpha[1] ~ dunif(0,100)
      alpha[2] ~ dunif(0,100)
      beta[1] ~ dunif(0,80)
      beta[2] ~ dunif(0,80)
    }
  }else if(nchange == 2){
    var_names <- c("a1","a2","a3","b1","b2","b3","tau1","tau2")
    jags.data <- c("t1","n1","TT","zeros","n")
    jags.par <- c("alpha[1]","alpha[2]","alpha[3]","beta[1]","beta[2]","beta[3]","tau[1]","tau[2]")
    model.jags <- function()
    {
      for(i in 1:n){
        for(j in 1:n1[i])
        {
          phi[i,j] <- -l[i,j]
          zeros[i,j] ~ dpois(phi[i,j])
          
          k[i,j] <- 3-step(tau[1]-t1[i,j])-step(tau[2]-t1[i,j])
          llambda[i,j] <- log(alpha[k[i,j]])-log(beta[k[i,j]])+(alpha[k[i,j]]-1)*(log(t1[i,j])-log(beta[k[i,j]]))
          
          l[i,j] <- llambda[i,j]+mT[i]/n1[i]
        }
        
        mT[i] <- -pow(tau[1]/beta[1],alpha[1])-
          pow(tau[2]/beta[2],alpha[2])+pow(tau[1]/beta[2],alpha[2])-
          pow(TT[i]/beta[3],alpha[3])+pow(tau[2]/beta[3],alpha[3])
      }
      tau[1] ~ dunif(0,max(TT))
      tau[2] ~ dunif(tau[1],max(TT))
      alpha[1] ~ dunif(0,100)
      alpha[2] ~ dunif(0,100)
      alpha[3] ~ dunif(0,100)
      beta[1] ~ dunif(0,80)
      beta[2] ~ dunif(0,80)
      beta[3] ~ dunif(0,80)
    }
  }else if(nchange == 3){
    var_names <- c("a1","a2","a3","a4","b1","b2","b3","b4","tau1","tau2","tau3")
    jags.data <- c("t1","n1","TT","zeros","n")
    jags.par <- c("alpha[1]","alpha[2]","alpha[3]","alpha[4]",
                  "beta[1]","beta[2]","beta[3]","beta[4]",
                  "tau[1]","tau[2]","tau[3]")
    model.jags <- function()
    {
      for(i in 1:n){
        for(j in 1:n1[i])
        {
          phi[i,j] <- -l[i,j]
          zeros[i,j] ~ dpois(phi[i,j])
          
          k[i,j] <- 4-step(tau[1]-t1[i,j])-step(tau[2]-t1[i,j])-step(tau[3]-t1[i,j])
          llambda[i,j] <- log(alpha[k[i,j]])-log(beta[k[i,j]])+
            (alpha[k[i,j]]-1)*(log(t1[i,j])-log(beta[k[i,j]]))
          
          l[i,j] <- llambda[i,j]+mT[i]/n1[i]
        }
        
        mT[i] <- -pow(tau[1]/beta[1],alpha[1])-
          pow(tau[2]/beta[2],alpha[2])+pow(tau[1]/beta[2],alpha[2])-
          pow(tau[3]/beta[3],alpha[3])+pow(tau[2]/beta[3],alpha[3])-
          pow(TT[i]/beta[4],alpha[4])+pow(tau[3]/beta[4],alpha[4])
      }
      tau[1] ~ dunif(0,max(TT))
      tau[2] ~ dunif(tau[1],max(TT))
      tau[3] ~ dunif(tau[2],max(TT))
      alpha[1] ~ dunif(0,100)
      alpha[2] ~ dunif(0,100)
      alpha[3] ~ dunif(0,100)
      alpha[4] ~ dunif(0,100)
      beta[1] ~ dunif(0,80)
      beta[2] ~ dunif(0,80)
      beta[3] ~ dunif(0,80)
      beta[4] ~ dunif(0,80)
    }
  }
  if(nchange != 0){
    bayes.mod.inits <- function(){
      list("alpha" = runif(nchange+1,0,5), "beta" = runif(nchange+1,0,100),
           "tau" = sort(runif(nchange,0,max(TT))))
    }
  }else{
    bayes.mod.inits <- function(){
      list("alpha" = runif(nchange+1,0,5), "beta" = runif(nchange+1,0,100))
    }
  }
  iter <- 1
  C = read.table(paste0("censortm",iter,".csv"),head=F)/4/3
  evtm = read.csv(paste0("eventtm",iter,".csv"),head=F)/4/3
  nt = read.csv(paste0("Nt",iter,".csv"),head=F)
  
  
  t1 = evtm
  TT = C[[1]]
  n1 = sapply(1:nrow(t1), FUN = function(x) sum(t1[x,]!=0))
  zeros = matrix(0,length(TT),max(n1))
  n = length(TT)
  estimativas <- jags(              data = jags.data,
                                    #     inits = jags.inits,
                                    parameters.to.save = jags.par,
                                    model.file = model.jags,
                                    n.iter = 20*1000,
                                    n.burnin = 10000,
                                    n.chain = 1,
                                    DIC = TRUE
                                    )
  print(estimativas$BUGSoutput$DIC)
  output <- estimativas$BUGSoutput$sims.matrix
  seq <- 1:dim(output)[2]
  seq <- seq[-((nchange+1)*2+1)]
  library(ggplot2)
  p=list()
  count = 1
  png(file=paste0(nchange,"change-post.png"), width = 1080, height = 603)
  for(i in seq){
    name <- colnames(output)[i]
    ldata <- data.frame(simulation = output[,i])
    p[[count]] = ggplot(ldata,aes(x=simulation))+geom_density()+ggtitle(name)+
      theme(plot.title = element_text(hjust = 0.5))
    count = count + 1
  }
  pic <- Rmisc::multiplot(plotlist=p, cols=2+nchange)
  dev.off()
  start <- as.Date("1850-01-01")
  if(nchange != 0){
    mean <- sapply(1:nchange, FUN = function(x) mean(output[,(1+nchange)*2+1+x]))*12
    mode <- sapply(1:nchange, FUN = function(x) median(output[,(1+nchange)*2+1+x]))*12
    sd <- sapply(1:nchange, FUN = function(x) sd(output[,(1+nchange)*2+1+x]))/12
    CI_l <- mean - 1.96*sd
    CI_u <- mean + 1.96*sd
    mean <- start + round(mean)
    mode <- start + round(mode)
    CI_l <- CI_l + start
    CI_u <- CI_u + start
    ana <- data.frame("Change Points" = 1:nchange ,Mean = mean,Mode = mode , sd = sd,"CI lower" = CI_l,
                      "CI upper" = CI_u)
    write.table(ana, paste0(nchange, "change-analysis.csv"),sep=",", col.names = 
                  c("Change Point","Mean", "Mode","sd","CI lower","CI upper"),row.names = F)
    names(output) <- var_names
    write.table(output[,seq], paste0("simulation_data2_change",nchange,".csv"), col.names = var_names, row.names = F, sep=",")
    print(sd)
  }else{
    names(output) <- var_names
    write.table(output, paste0("simulation_data2_change",nchange,".csv"), col.names = T, row.names = F,sep=",")
  }
}
if(nchange == 1){
  seq=c(1:4,6)
  mean = round(sapply(seq, FUN = function(x) mean(output[,x])),2)
  mode = round(sapply(seq, FUN= function(x) median(output[,x])),2)
  sd = round(sapply(seq, FUN = function(x) sd(output[,x])),2)
  CI_low = sapply(seq, FUN = function(x) quantile(output[,x],0.025))
  CI_up = sapply(seq, FUN = function(x) quantile(output[,x],0.975))
  CI <- paste0("[",CI_low,",",CI_up,"]")
  write.table( file = "1change_par_ana2.csv",data.frame(mean = mean, mode = mode, sd = sd, CI_low = CI_low, CI_up = CI_up), 
               row.names = F,sep=",")
}

event = read.csv(paste0("eventtm",iter,".csv"),head=F)
sum = sapply(1:24, FUN = function(x) sum((event[x,]!=0)))
mean(sum);sd(sum);min(sum);max(sum);sum(sum)
