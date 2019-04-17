#############################Setting hyperparameters####################################
args=commandArgs(TRUE)
datasets <- read.table("dataset",sep=",",header=F,fill=T,col.names = paste0("V",seq_len(11)))
print(args)
print(datasets)
print(typeof(args))
args = as.numeric(args)
dataset=as.numeric(datasets[args,1])
nchange=as.numeric(datasets[args,2])
m=as.numeric(datasets[args,3])
args=datasets[args,]
args=args[!is.na(args)]
nset=30
dir=paste('data',dataset,sep="")
write(dir,"dir")
setwd(dir)
library(R2jags)
evtm_all=data.frame()
#############################build the model and do the simulation for different number of change points#############################
if(nchange==1){
  var_names <- c("a1","a2","b1","b2","tau")
  
  for (iter in 1:nset){
    C = read.table(paste0("censortm",iter,".csv"),head=F)
    evtm = read.csv(paste0("eventtm",iter,".csv"),head=F)
    for (i in 1:2){
      len = sum(evtm[i,]!=0)
      evtm[i,1:len]=sort(evtm[i,1:len])
    }
    nt = read.csv(paste0("Nt",iter,".csv"),head=F)
    evtm_all=merge(evtm_all,evtm,all=T)
    C<-C[1:m,]
    evtm<-evtm[1:m,]
    nt=nt[1:m]
    
    set.seed(123)
    t1 = evtm
    TT = C
    n1 = sapply(1:nrow(t1), FUN = function(x) sum(t1[x,]!=0))
    zeros = matrix(0,length(TT),max(n1))+1
    n = length(TT)
    jags.data <- c("t1","n1","TT","zeros","n")
    
    jags.par <- c("alpha[1]","alpha[2]","sigma[1]","sigma[2]","tau")
    #############################build model########################################
    model.jags <- function()
    {
      for(i in 1:n){
        for(j in 1:n1[i])
        {
          phi[i,j] <- -l[i,j]
          zeros[i,j] ~ dpois(phi[i,j])
          
          k[i,j] <- 2-step(tau-t1[i,j])
          llambda[i,j] <- log(alpha[k[i,j]])-log(sigma[k[i,j]])+(alpha[k[i,j]]-1)*(log(t1[i,j])-log(sigma[k[i,j]]))
          
          l[i,j] <- llambda[i,j]+mT[i]/n1[i]
        }
        
        mT[i] <- -pow(tau/sigma[1],alpha[1])-pow(TT[i]/sigma[2],alpha[2])+pow(tau/sigma[2],alpha[2])
      }
      tau ~ dunif(0,max(TT))
      alpha[1] ~ dunif(0,100)
      alpha[2] ~ dunif(0,100)
      sigma[1] ~ dunif(0,100)
      sigma[2] ~ dunif(0,100)
    }
    ###############################The output of the simulation###########################################
    estimativas <- jags(              data = jags.data,
                                      #     inits = jags.inits,
                                      parameters.to.save = jags.par,
                                      model.file = model.jags,
                                      n.iter = 15000,
                                      n.burnin = 5000,
                                      n.chain = 1,
				      DIC = TRUE)
    
    alpha1 <- estimativas$BUGSoutput$sims.list$alpha[,1]
    alpha2 <- estimativas$BUGSoutput$sims.list$alpha[,2]
    sigma1 <- estimativas$BUGSoutput$sims.list$sigma[,1]
    sigma2 <- estimativas$BUGSoutput$sims.list$sigma[,2]
    tau <- estimativas$BUGSoutput$sims.list$tau
    
    # print(sapply(1:5, FUN = function(x) quantile(data.frame(alpha1,alpha2,sigma1,sigma2,tau)[,x],c(0.05,0.95))))
    # print(sapply(1:5, FUN = function(x) mean(data.frame(alpha1,alpha2,sigma1,sigma2,tau)[,x],na.rm = T)))
    write.table(alpha1,paste0("a1_",nchange,"_",iter,".csv"),row.names = F, col.names = F)
    write.table(alpha2,paste0("a2_",nchange,"_",iter,".csv"),row.names = F, col.names = F)
    write.table(sigma1,paste0("b1_",nchange,"_",iter,".csv"),row.names = F, col.names = F)
    write.table(sigma2,paste0("b2_",nchange,"_",iter,".csv"),row.names = F, col.names = F)
    write.table(tau,paste0("tau_",nchange,"_",iter,".csv"),row.names = F, col.names = F)
    
    npars=(length(args)-3)%/%3
    true_value <- as.numeric(c(args[(4+npars):(4+npars*3+1)],args[4:(4+npars-1)]))
    file_names <- paste(var_names,"_",nchange,"_",iter,".csv",sep="")
    a1 <- read.csv(file_names[1],header = F)[[1]]
    a2 <- read.csv(file_names[2],header = F)[[1]]
    b1 <- read.csv(file_names[3],header = F)[[1]]
    b2 <- read.csv(file_names[4],header = F)[[1]]
    tau <- read.csv(file_names[5],header = F)[[1]]
    # dt <- data.frame(k = 2, a1 = mean(a1), a2 = mean(a2), a3 = mean(a3), b1 = mean(b1), b2 = mean(b2), b3 = mean(b3) ,tau1 = mean(tau1), tau2 = mean(tau2))
    dt <- c(1,mean(a1),mean(a2),mean(b1),mean(b2),mean(tau))
    sim_data=data.frame(a1=a1,a2=a2,b1=b1,b2=b2,tau=tau)
    ###########################################Save the summary table of the output############################################
    setwd("..")
    if(iter==1){
      write.table(x=sim_data,file=paste0("simulation_data",dataset,"_nchange",nchange,".csv"),row.names = F,sep=",")
    }else{
      write.table(x=sim_data,file=paste0("simulation_data",dataset,"_nchange",nchange,".csv"),col.names = F,row.names = F,append = T,sep=",")
    }
    #names(true_value)=var_names
    write.table(t(true_value),paste0("truevalue_data",dataset,".csv"),sep=",",row.names = F)
    write.table(estimativas$BUGSoutput$DIC,file=paste0("DIC_data",dataset,"_nchange",nchange,".csv"),col.names = F,row.names = F,append = T,sep=",")
    setwd(dir)
    
    
    # write(t(dt),paste("data",dataset,"_change",nchange,".csv",sep = ""),sep=',',ncolumns = length(dt))
    # # write.table(t(dt),paste("data",dataset,"_change",nchange,".csv",sep = ""),col.names = F,row.names = F)
    # C <- read.csv(paste("censortm",iter,".csv",sep=""),header = F)[[1]]
    # N <- read.csv(paste("Nt",iter,".csv",sep=""),header = F)
    # t <- read.csv(paste("eventtm",iter,".csv",sep = ""),header = F)
    # f = file(paste("data",dataset,"_change",nchange,".csv",sep = ""),"a")
    # for(i in 1:length(C)){
    #   t_tmp <- t[i,]
    #   t_tmp <- t_tmp[t_tmp!=0]
    #   tmp <- c(N[[i]], t_tmp, C[i])
    #   cat(t(tmp),file =A f, sep = ",")
    #   cat("\n",file = f)
    # write(t(tmp),file = paste0("data",dataset,"_change",nchange,".csv",sep = ""),sep=",", append = T,ncolumns = length(tmp))
    # write(t(tmp),file = paste0("data",dataset,"_change",nchange,".csv",sep = ""),append = T)
  }
  data <- data.frame()
  for ( i in 1:nset){
    file_names <- paste(var_names,"_",nchange,"_",i,".csv",sep="")
    a1 <- read.csv(file_names[1],header = F)[[1]]
    a2 <- read.csv(file_names[2],header = F)[[1]]
    b1 <- read.csv(file_names[3],header = F)[[1]]
    b2 <- read.csv(file_names[4],header = F)[[1]]
    tau <- read.csv(file_names[5],header = F)[[1]]
    data_fm <- data.frame(a1 = a1,a2 = a2, b1 = b1, b2 = b2, tau = tau)
    
    data <- rbind(data,data_fm)
  }
  coverate_cal <- function(dtf,true,prob){
    coverage <- c()
    for(i in 1:dim(dtf)[2]){
      count = 0
      for(j in seq(0,dim(dtf)[1]-1,1000)){
        tmp <- dtf[j:(j+999)+1,i]
        ran <- quantile(tmp,c(prob/2,1-prob/2))
        if(true[i]<=ran[2] & true[i] >= ran[1]){count = count + 1}
      }
      coverage <- c(coverage,count/(dim(dtf)[1]/1000))
    }
    return(coverage)
  }
  write.csv(data[,1],'data.csv')
  write.csv(true_value,'value.csv')
  if(nchange*3+2 != length(true_value)){stop('Incorrect value for analysis')}
  mean <- colMeans(data)
  RMSE <- sqrt(sapply(1:length(true_value), FUN = function(x) mean((data[,x]-true_value[x])^2)))
  bias <- abs(mean-true_value)/true_value*100
  coverage <- coverate_cal(data,true_value,0.05)*100
  output <- data.frame(par_name = var_names, true_value = true_value,mean = mean, RMSE =RMSE, BIAS = bias, coverage = coverage)
  setwd('..')
  evtm_all[evtm_all==0]=NA
  write.csv(output,paste0("analysis_data",dataset,"_nchange",nchange,".csv"))
  write.csv(evtm_all,paste0("evetm_data",dataset,".csv"))
  # close(f)
}else if(nchange == 2){
  var_names <- c("a1","a2","a3","b1","b2","b3","tau1","tau2")
  for (iter in 1:nset){
    C = read.table(paste0("censortm",iter,".csv"),head=F)
    evtm = read.csv(paste0("eventtm",iter,".csv"),head=F)
    for (i in 1:2){
      len = sum(evtm[i,]!=0)
      evtm[i,1:len]=sort(evtm[i,1:len])
    }
    nt = read.csv(paste0("Nt",iter,".csv"),head=F)
    evtm_all=merge(evtm_all,evtm,all=T)
    C<-C[1:m,]
    evtm<-evtm[1:m,]
    nt=nt[1:m]
    set.seed(123)
    
    t1 = evtm
    TT = C
    n1 = sapply(1:nrow(t1), FUN = function(x) sum(t1[x,]!=0))
    zeros = matrix(0,length(TT),max(n1))+1
    n = length(TT)
    jags.data <- c("t1","n1","TT","zeros","n")
    
    jags.par <- c("alpha[1]","alpha[2]","alpha[3]","sigma[1]","sigma[2]","sigma[3]","tau[1]","tau[2]")
    
    model.jags <- function()
    {
      for(i in 1:n){
        for(j in 1:n1[i])
        {
          phi[i,j] <- -l[i,j]
          zeros[i,j] ~ dpois(phi[i,j])
          
          k[i,j] <- 3-step(tau[1]-t1[i,j])-step(tau[2]-t1[i,j])
          llambda[i,j] <- log(alpha[k[i,j]])-log(sigma[k[i,j]])+(alpha[k[i,j]]-1)*(log(t1[i,j])-log(sigma[k[i,j]]))
          
          l[i,j] <- llambda[i,j]+mT[i]/n1[i]
        }
        
        mT[i] <- -pow(tau[1]/sigma[1],alpha[1])-
          pow(tau[2]/sigma[2],alpha[2])+pow(tau[1]/sigma[2],alpha[2])-
          pow(TT[i]/sigma[3],alpha[3])+pow(tau[2]/sigma[3],alpha[3])
      }
      tau[1] ~ dunif(0,max(TT))
      tau[2] ~ dunif(tau[1],max(TT))
      alpha[1] ~ dunif(0,100)
      alpha[2] ~ dunif(0,100)
      alpha[3] ~ dunif(0,100)
      sigma[1] ~ dunif(0,100)
      sigma[2] ~ dunif(0,100)
      sigma[3] ~ dunif(0,100)
    }
    
    estimativas <- jags(              data = jags.data,
                                      #     inits = jags.inits,
                                      parameters.to.save = jags.par,
                                      model.file = model.jags,
                                      n.iter = 15000,
                                      n.burnin = 5000,
                                      n.chain = 1,
				      DIC = TRUE)
    
    alpha1 <- estimativas$BUGSoutput$sims.list$alpha[,1]
    alpha2 <- estimativas$BUGSoutput$sims.list$alpha[,2]
    alpha3 <- estimativas$BUGSoutput$sims.list$alpha[,3]
    sigma1 <- estimativas$BUGSoutput$sims.list$sigma[,1]
    sigma2 <- estimativas$BUGSoutput$sims.list$sigma[,2]
    sigma3 <- estimativas$BUGSoutput$sims.list$sigma[,3]
    tau1 <- estimativas$BUGSoutput$sims.list$tau[,1]
    tau2 <- estimativas$BUGSoutput$sims.list$tau[,2]
    # sapply(1:8, FUN = function(x) quantile(data.frame(alpha1,alpha2,alpha3,sigma1,sigma2,sigma3,tau1,tau2)[,x],c(0.05,0.95)))
    # sapply(1:8, FUN = function(x) mean(data.frame(alpha1,alpha2,alpha3,sigma1,sigma2,sigma3,tau1,tau2)[,x],na.rm = T))
    write.table(alpha1,paste0("a1_",nchange,"_",iter,".csv"),row.names = F, col.names = F)
    write.table(alpha2,paste0("a2_",nchange,"_",iter,".csv"),row.names = F, col.names = F)
    write.table(alpha3,paste0("a3_",nchange,"_",iter,".csv"),row.names = F, col.names = F)
    write.table(sigma1,paste0("b1_",nchange,"_",iter,".csv"),row.names = F, col.names = F)
    write.table(sigma2,paste0("b2_",nchange,"_",iter,".csv"),row.names = F, col.names = F)
    write.table(sigma3,paste0("b3_",nchange,"_",iter,".csv"),row.names = F, col.names = F)
    write.table(tau1,paste0("tau1_",nchange,"_",iter,".csv"),row.names = F, col.names = F)
    write.table(tau2,paste0("tau2_",nchange,"_",iter,".csv"),row.names = F, col.names = F)
    
    npars=(length(args)-3)%/%3
    true_value <- as.numeric(c(args[(4+npars):(4+npars*3+1)],args[4:(4+npars-1)]))
    file_names <- paste(var_names,"_",nchange,"_",iter,".csv",sep="")
    a1 <- read.csv(file_names[1],header = F)[[1]]
    a2 <- read.csv(file_names[2],header = F)[[1]]
    a3 <- read.csv(file_names[3],header = F)[[1]]
    b1 <- read.csv(file_names[4],header = F)[[1]]
    b2 <- read.csv(file_names[5],header = F)[[1]]
    b3 <- read.csv(file_names[6],header = F)[[1]]
    tau1 <- read.csv(file_names[7],header = F)[[1]]
    tau2 <- read.csv(file_names[8],header = F)[[1]]
    dt <- data.frame(k = 2, a1 = mean(a1), a2 = mean(a2), a3 = mean(a3), b1 = mean(b1), b2 = mean(b2), b3 = mean(b3) ,tau1 = mean(tau1), tau2 = mean(tau2))
    
    sim_data=data.frame(a1=a1,a2=a2,a3=a3,b1=b1,b2=b2,b3=b3,tau1=tau1,tau2=tau2)
    
    setwd("..")
    if(iter==1){
      write.table(x=sim_data,file=paste0("simulation_data",dataset,"_nchange",nchange,".csv"),row.names = F,sep=",")
    }else{
      write.table(x=sim_data,file=paste0("simulation_data",dataset,"_nchange",nchange,".csv"),col.names = F,row.names = F,append = T,sep=",")
    }
    #names(true_value)=var_names
    write.table(t(true_value),paste0("truevalue_data",dataset,".csv"),sep=",",row.names = F)
    write.table(estimativas$BUGSoutput$DIC,file=paste0("DIC_data",dataset,"_nchange",nchange,".csv"),col.names = F,row.names = F,append = T,sep=",")
    setwd(dir)
  }
  data <- data.frame()
  for ( i in 1:nset){
    file_names <- paste(var_names,"_",nchange,"_",i,".csv",sep="")
    a1 <- read.csv(file_names[1],header = F)[[1]]
    a2 <- read.csv(file_names[2],header = F)[[1]]
    a3 <- read.csv(file_names[3],header = F)[[1]]
    b1 <- read.csv(file_names[4],header = F)[[1]]
    b2 <- read.csv(file_names[5],header = F)[[1]]
    b3 <- read.csv(file_names[6],header = F)[[1]]
    tau1 <- read.csv(file_names[7],header = F)[[1]]
    tau2 <- read.csv(file_names[8],header = F)[[1]]
    data_fm <- data.frame(a1 = a1,a2 = a2, a3=a3,b1 = b1, b2 = b2,b3=b3,tau1 = tau1,tau2=tau2)
    
    data <- rbind(data,data_fm)
  }
  coverate_cal <- function(dtf,true,prob){
    coverage <- c()
    for(i in 1:dim(dtf)[2]){
      count = 0
      for(j in seq(0,dim(dtf)[1]-1,1000)){
        tmp <- dtf[j:(j+999)+1,i]
        ran <- quantile(tmp,c(prob/2,1-prob/2))
        if(true[i]<=ran[2] & true[i] >= ran[1]){count = count + 1}
      }
      coverage <- c(coverage,count/(dim(dtf)[1]/1000))
    }
    return(coverage)
  }
  write.csv(data[,1],'data.csv')
  write.csv(true_value,'value.csv')
  if(nchange*3+2 != length(true_value)){stop('Incorrect value for analysis')}
  mean <- colMeans(data)
  RMSE <- sqrt(sapply(1:length(true_value), FUN = function(x) mean((data[,x]-true_value[x])^2)))
  bias <- abs(mean-true_value)/true_value*100
  coverage <- coverate_cal(data,true_value,0.05)*100
  output <- data.frame(par_name = var_names, true_value = true_value,mean = mean, RMSE =RMSE, BIAS = bias, coverage = coverage)
  setwd('..')
  evtm_all[evtm_all==0]=NA
  write.csv(output,paste0("analysis_data",dataset,"_nchange",nchange,".csv"))
  write.csv(evtm_all,paste0("evetm_data",dataset,".csv"))
  # close(f)
  # dt <- c(2,mean(a1),mean(a2),mean(a3),mean(b1),mean(b2),mean(b3),mean(tau1),mean(tau2))
  # write(t(dt),paste("data",dataset,"_change",nchange,".csv",sep = ""),sep=',',ncolumns = length(dt))
  # # write.table(t(dt),paste("data",dataset,"_change",nchange,".csv",sep = ""),col.names = F,row.names = F)
  # C <- read.csv(paste("censortm",iter,".csv",sep=""),header = F)[[1]]
  # N <- read.csv(paste("Nt",iter,".csv",sep=""),header = F)
  # t <- read.csv(paste("eventtm",iter,".csv",sep = ""),header = F)
  # f = file(paste("data",dataset,"_change",nchange,".csv",sep = ""),"a")
  # for(i in 1:length(C)){
  #   t_tmp <- t[i,]
  #   t_tmp <- t_tmp[t_tmp!=0]
  #   tmp <- c(N[[i]], t_tmp, C[i])
  #   cat(t(tmp),file = f, sep = ",")
  #   cat("\n",file = f)
  #   # write(t(tmp),file = paste0("data",dataset,"_change",nchange,".csv",sep = ""),sep=",", append = T,ncolumns = length(tmp))
  #   # write(t(tmp),file = paste0("data",dataset,"_change",nchange,".csv",sep = ""),append = T)
  # }
  # close(f) 
}else if(nchange == 3){
  var_names <- c("a1","a2","a3","a4","b1","b2","b3","b4","tau1","tau2","tau3")
  for (iter in 1:nset){
    C = read.table(paste0("censortm",iter,".csv"),head=F)
    evtm = read.csv(paste0("eventtm",iter,".csv"),head=F)
    for (i in 1:2){
      len = sum(evtm[i,]!=0)
      evtm[i,1:len]=sort(evtm[i,1:len])
    }
    nt = read.csv(paste0("Nt",iter,".csv"),head=F)
    evtm_all=merge(evtm_all,evtm,all=T)
    C<-C[1:m,]
    evtm<-evtm[1:m,]
    nt=nt[1:m]
    set.seed(123)
    
    t1 = evtm
    TT = C
    n1 = sapply(1:nrow(t1), FUN = function(x) sum(t1[x,]!=0))
    zeros = matrix(0,length(TT),max(n1))+1
    n = length(TT)
    jags.data <- c("t1","n1","TT","zeros","n")
    
    jags.par <- c("alpha[1]","alpha[2]","alpha[3]","alpha[4]",
                  "sigma[1]","sigma[2]","sigma[3]","sigma[4]",
                  "tau[1]","tau[2]","tau[3]")
    
    model.jags <- function()
    {
      for(i in 1:n){
        for(j in 1:n1[i])
        {
          phi[i,j] <- -l[i,j]
          zeros[i,j] ~ dpois(phi[i,j])
          
          k[i,j] <- 4-step(tau[1]-t1[i,j])-step(tau[2]-t1[i,j])-step(tau[3]-t1[i,j])
          llambda[i,j] <- log(alpha[k[i,j]])-log(sigma[k[i,j]])+
            (alpha[k[i,j]]-1)*(log(t1[i,j])-log(sigma[k[i,j]]))
          
          l[i,j] <- llambda[i,j]+mT[i]/n1[i]
        }
        
        mT[i] <- -pow(tau[1]/sigma[1],alpha[1])-
          pow(tau[2]/sigma[2],alpha[2])+pow(tau[1]/sigma[2],alpha[2])-
          pow(tau[3]/sigma[3],alpha[3])+pow(tau[2]/sigma[3],alpha[3])-
          pow(TT[i]/sigma[4],alpha[4])+pow(tau[3]/sigma[4],alpha[4])
      }
      tau[1] ~ dunif(0,max(TT))
      tau[2] ~ dunif(tau[1],max(TT))
      tau[3] ~ dunif(tau[2],max(TT))
      alpha[1] ~ dunif(0,100)
      alpha[2] ~ dunif(0,100)
      alpha[3] ~ dunif(0,100)
      alpha[4] ~ dunif(0,100)
      sigma[1] ~ dunif(0,100)
      sigma[2] ~ dunif(0,100)
      sigma[3] ~ dunif(0,100)
      sigma[4] ~ dunif(0,100)
    }
    
    estimativas <- jags(              data = jags.data,
                                      #     inits = jags.inits,
                                      parameters.to.save = jags.par,
                                      model.file = model.jags,
                                      n.iter = 15000,
                                      n.burnin = 5000,
                                      n.chain = 1,
                                      DIC = TRUE)
    
    alpha1 <- estimativas$BUGSoutput$sims.list$alpha[,1]
    alpha2 <- estimativas$BUGSoutput$sims.list$alpha[,2]
    alpha3 <- estimativas$BUGSoutput$sims.list$alpha[,3]
    alpha4 <- estimativas$BUGSoutput$sims.list$alpha[,4]
    sigma1 <- estimativas$BUGSoutput$sims.list$sigma[,1]
    sigma2 <- estimativas$BUGSoutput$sims.list$sigma[,2]
    sigma3 <- estimativas$BUGSoutput$sims.list$sigma[,3]
    sigma4 <- estimativas$BUGSoutput$sims.list$sigma[,4]
    tau1 <- estimativas$BUGSoutput$sims.list$tau[,1]
    tau2 <- estimativas$BUGSoutput$sims.list$tau[,2]
    tau3 <- estimativas$BUGSoutput$sims.list$tau[,3]
    # sapply(1:8, FUN = function(x) quantile(data.frame(alpha1,alpha2,alpha3,sigma1,sigma2,sigma3,tau1,tau2)[,x],c(0.05,0.95)))
    # sapply(1:8, FUN = function(x) mean(data.frame(alpha1,alpha2,alpha3,sigma1,sigma2,sigma3,tau1,tau2)[,x],na.rm = T))
    write.table(alpha1,paste0("a1_",nchange,"_",iter,".csv"),row.names = F, col.names = F)
    write.table(alpha2,paste0("a2_",nchange,"_",iter,".csv"),row.names = F, col.names = F)
    write.table(alpha3,paste0("a3_",nchange,"_",iter,".csv"),row.names = F, col.names = F)
    write.table(alpha4,paste0("a4_",nchange,"_",iter,".csv"),row.names = F, col.names = F)
    write.table(sigma1,paste0("b1_",nchange,"_",iter,".csv"),row.names = F, col.names = F)
    write.table(sigma2,paste0("b2_",nchange,"_",iter,".csv"),row.names = F, col.names = F)
    write.table(sigma3,paste0("b3_",nchange,"_",iter,".csv"),row.names = F, col.names = F)
    write.table(sigma4,paste0("b4_",nchange,"_",iter,".csv"),row.names = F, col.names = F)
    write.table(tau1,paste0("tau1_",nchange,"_",iter,".csv"),row.names = F, col.names = F)
    write.table(tau2,paste0("tau2_",nchange,"_",iter,".csv"),row.names = F, col.names = F)
    write.table(tau3,paste0("tau3_",nchange,"_",iter,".csv"),row.names = F, col.names = F)
    
    npars=(length(args)-3)%/%3
    true_value <- as.numeric(c(args[(4+npars):(4+npars*3+1)],args[4:(4+npars-1)]))
    file_names <- paste(var_names,"_",nchange,"_",iter,".csv",sep="")
    a1 <- read.csv(file_names[1],header = F)[[1]]
    a2 <- read.csv(file_names[2],header = F)[[1]]
    a3 <- read.csv(file_names[3],header = F)[[1]]
    a4 <- read.csv(file_names[4],header = F)[[1]]
    b1 <- read.csv(file_names[5],header = F)[[1]]
    b2 <- read.csv(file_names[6],header = F)[[1]]
    b3 <- read.csv(file_names[7],header = F)[[1]]
    b4 <- read.csv(file_names[8],header = F)[[1]]
    tau1 <- read.csv(file_names[9],header = F)[[1]]
    tau2 <- read.csv(file_names[10],header = F)[[1]]
    tau3 <- read.csv(file_names[11],header = F)[[1]]
    dt <- data.frame(k = 2, a1 = mean(a1), a2 = mean(a2), a3 = mean(a3), a4 = mean(a4),
                     b1 = mean(b1), b2 = mean(b2), b3 = mean(b3), b4 = mean(b4),
                     tau1 = mean(tau1), tau2 = mean(tau2), tau3 = mean(tau3))
    
    sim_data=data.frame(a1=a1,a2=a2,a3=a3,a4 = a4,
                        b1=b1,b2=b2,b3=b3,b4=b4,
                        tau1=tau1,tau2=tau2,tau3=tau3)
    
    setwd("..")
    if(iter==1){
      write.table(x=sim_data,file=paste0("simulation_data",dataset,"_nchange",nchange,".csv"),row.names = F,sep=",")
    }else{
      write.table(x=sim_data,file=paste0("simulation_data",dataset,"_nchange",nchange,".csv"),col.names = F,row.names = F,append = T,sep=",")
    }
    #names(true_value)=var_names
    write.table(t(true_value),paste0("truevalue_data",dataset,".csv"),sep=",",row.names = F)
    write.table(estimativas$BUGSoutput$DIC,file=paste0("DIC_data",dataset,"_nchange",nchange,".csv"),col.names = F,row.names = F,append = T,sep=",")
    setwd(dir)
  }
  data <- data.frame()
  for ( i in 1:nset){
    file_names <- paste(var_names,"_",nchange,"_",i,".csv",sep="")
    a1 <- read.csv(file_names[1],header = F)[[1]]
    a2 <- read.csv(file_names[2],header = F)[[1]]
    a3 <- read.csv(file_names[3],header = F)[[1]]
    a4 <- read.csv(file_names[4],header = F)[[1]]
    b1 <- read.csv(file_names[5],header = F)[[1]]
    b2 <- read.csv(file_names[6],header = F)[[1]]
    b3 <- read.csv(file_names[7],header = F)[[1]]
    b4 <- read.csv(file_names[8],header = F)[[1]]
    tau1 <- read.csv(file_names[9],header = F)[[1]]
    tau2 <- read.csv(file_names[10],header = F)[[1]]
    tau3 <- read.csv(file_names[11],header = F)[[1]]
    data_fm <- data.frame(a1 = a1,a2 = a2, a3=a3,a4=a4,
                          b1 = b1, b2 = b2,b3=b3,b4=b4,
                          tau1 = tau1,tau2=tau2,tau3=tau3)
    
    data <- rbind(data,data_fm)
  }
  coverate_cal <- function(dtf,true,prob){
    coverage <- c()
    for(i in 1:dim(dtf)[2]){
      count = 0
      for(j in seq(0,dim(dtf)[1]-1,1000)){
        tmp <- dtf[j:(j+999)+1,i]
        ran <- quantile(tmp,c(prob/2,1-prob/2))
        if(true[i]<=ran[2] & true[i] >= ran[1]){count = count + 1}
      }
      coverage <- c(coverage,count/(dim(dtf)[1]/1000))
    }
    return(coverage)
  }
  write.csv(data[,1],'data.csv')
  write.csv(true_value,'value.csv')
  if(nchange*3+2 != length(true_value)){stop('Incorrect value for analysis')}
  mean <- colMeans(data)
  RMSE <- sqrt(sapply(1:length(true_value), FUN = function(x) mean((data[,x]-true_value[x])^2)))
  bias <- abs(mean-true_value)/true_value*100
  coverage <- coverate_cal(data,true_value,0.05)*100
  output <- data.frame(par_name = var_names, true_value = true_value,mean = mean, RMSE =RMSE, BIAS = bias, coverage = coverage)
  setwd('..')
  evtm_all[evtm_all==0]=NA
  write.csv(output,paste0("analysis_data",dataset,"_nchange",nchange,".csv"))
  write.csv(evtm_all,paste0("evetm_data",dataset,".csv"))
  # close(f)
  # dt <- c(2,mean(a1),mean(a2),mean(a3),mean(b1),mean(b2),mean(b3),mean(tau1),mean(tau2))
  # write(t(dt),paste("data",dataset,"_change",nchange,".csv",sep = ""),sep=',',ncolumns = length(dt))
  # # write.table(t(dt),paste("data",dataset,"_change",nchange,".csv",sep = ""),col.names = F,row.names = F)
  # C <- read.csv(paste("censortm",iter,".csv",sep=""),header = F)[[1]]
  # N <- read.csv(paste("Nt",iter,".csv",sep=""),header = F)
  # t <- read.csv(paste("eventtm",iter,".csv",sep = ""),header = F)
  # f = file(paste("data",dataset,"_change",nchange,".csv",sep = ""),"a")
  # for(i in 1:length(C)){
  #   t_tmp <- t[i,]
  #   t_tmp <- t_tmp[t_tmp!=0]
  #   tmp <- c(N[[i]], t_tmp, C[i])
  #   cat(t(tmp),file = f, sep = ",")
  #   cat("\n",file = f)
  #   # write(t(tmp),file = paste0("data",dataset,"_change",nchange,".csv",sep = ""),sep=",", append = T,ncolumns = length(tmp))
  #   # write(t(tmp),file = paste0("data",dataset,"_change",nchange,".csv",sep = ""),append = T)
  # }
  # close(f) 
}else if(nchange==0){
  var_names <- c("a1","b1")
  
  for (iter in 1:nset){
    C = read.table(paste0("censortm",iter,".csv"),head=F)
    evtm = read.csv(paste0("eventtm",iter,".csv"),head=F)
    for (i in 1:2){
      len = sum(evtm[i,]!=0)
      evtm[i,1:len]=sort(evtm[i,1:len])
    }
    nt = read.csv(paste0("Nt",iter,".csv"),head=F)
    evtm_all=merge(evtm_all,evtm,all=T)
    C<-C[1:m,]
    evtm<-evtm[1:m,]
    nt=nt[1:m]
    
    set.seed(123)
    t1 = evtm
    TT = C
    n1 = sapply(1:nrow(t1), FUN = function(x) sum(t1[x,]!=0))
    zeros = matrix(0,length(TT),max(n1))+1
    n = length(TT)
    jags.data <- c("t1","n1","TT","zeros","n")
    
    jags.par <- c("alpha[1]","sigma[1]")
    
    model.jags <- function()
    {
      for(i in 1:n){
        for(j in 1:n1[i])
        {
          phi[i,j] <- -l[i,j]
          zeros[i,j] ~ dpois(phi[i,j])
          
          k[i,j] <- 1
          llambda[i,j] <- log(alpha[k[i,j]])-log(sigma[k[i,j]])+(alpha[k[i,j]]-1)*(log(t1[i,j])-log(sigma[k[i,j]]))
          
          l[i,j] <- llambda[i,j]+mT[i]/n1[i]
        }
        
        mT[i] <- -pow(TT[i]/sigma[1],alpha[1])
      }
      alpha[1] ~ dunif(0,100)
      sigma[1] ~ dunif(0,100)
    }
    
    estimativas <- jags(              data = jags.data,
                                      #     inits = jags.inits,
                                      parameters.to.save = jags.par,
                                      model.file = model.jags,
                                      n.iter = 15000,
                                      n.burnin = 5000,
                                      n.chain = 1,
                                      DIC = TRUE)
    
    alpha1 <- estimativas$BUGSoutput$sims.list$alpha[,1]
    sigma1 <- estimativas$BUGSoutput$sims.list$sigma[,1]

    
    # print(sapply(1:5, FUN = function(x) quantile(data.frame(alpha1,alpha2,sigma1,sigma2,tau)[,x],c(0.05,0.95))))
    # print(sapply(1:5, FUN = function(x) mean(data.frame(alpha1,alpha2,sigma1,sigma2,tau)[,x],na.rm = T)))
    write.table(alpha1,paste0("a1_",nchange,"_",iter,".csv"),row.names = F, col.names = F)
    write.table(sigma1,paste0("b1_",nchange,"_",iter,".csv"),row.names = F, col.names = F)

    npars=(length(args)-3)%/%3
    true_value <- as.numeric(c(args[(4+npars):(4+npars*3+1)],args[4:(4+npars-1)]))
    file_names <- paste(var_names,"_",nchange,"_",iter,".csv",sep="")
    a1 <- read.csv(file_names[1],header = F)[[1]]
    b1 <- read.csv(file_names[2],header = F)[[1]]
    # dt <- data.frame(k = 2, a1 = mean(a1), a2 = mean(a2), a3 = mean(a3), b1 = mean(b1), b2 = mean(b2), b3 = mean(b3) ,tau1 = mean(tau1), tau2 = mean(tau2))
    dt <- c(1,mean(a1),mean(b1))
    sim_data=data.frame(a1=a1,b1=b1)
    
    setwd("..")
    if(iter==1){
      write.table(x=sim_data,file=paste0("simulation_data",dataset,"_nchange",nchange,".csv"),row.names = F,sep=",")
    }else{
      write.table(x=sim_data,file=paste0("simulation_data",dataset,"_nchange",nchange,".csv"),col.names = F,row.names = F,append = T,sep=",")
    }
    #names(true_value)=var_names
    write.table(t(true_value),paste0("truevalue_data",dataset,".csv"),sep=",",row.names = F)
    write.table(estimativas$BUGSoutput$DIC,file=paste0("DIC_data",dataset,"_nchange",nchange,".csv"),col.names = F,row.names = F,append = T,sep=",")
    setwd(dir)
    
    
    # write(t(dt),paste("data",dataset,"_change",nchange,".csv",sep = ""),sep=',',ncolumns = length(dt))
    # # write.table(t(dt),paste("data",dataset,"_change",nchange,".csv",sep = ""),col.names = F,row.names = F)
    # C <- read.csv(paste("censortm",iter,".csv",sep=""),header = F)[[1]]
    # N <- read.csv(paste("Nt",iter,".csv",sep=""),header = F)
    # t <- read.csv(paste("eventtm",iter,".csv",sep = ""),header = F)
    # f = file(paste("data",dataset,"_change",nchange,".csv",sep = ""),"a")
    # for(i in 1:length(C)){
    #   t_tmp <- t[i,]
    #   t_tmp <- t_tmp[t_tmp!=0]
    #   tmp <- c(N[[i]], t_tmp, C[i])
    #   cat(t(tmp),file =A f, sep = ",")
    #   cat("\n",file = f)
    # write(t(tmp),file = paste0("data",dataset,"_change",nchange,".csv",sep = ""),sep=",", append = T,ncolumns = length(tmp))
    # write(t(tmp),file = paste0("data",dataset,"_change",nchange,".csv",sep = ""),append = T)
  }
  data <- data.frame()
  for ( i in 1:nset){
    file_names <- paste(var_names,"_",nchange,"_",i,".csv",sep="")
    a1 <- read.csv(file_names[1],header = F)[[1]]
    b1 <- read.csv(file_names[2],header = F)[[1]]
    data_fm <- data.frame(a1 = a1, b1 = b1)
    
    data <- rbind(data,data_fm)
  }
  coverate_cal <- function(dtf,true,prob){
    coverage <- c()
    for(i in 1:dim(dtf)[2]){
      count = 0
      for(j in seq(0,dim(dtf)[1]-1,1000)){
        tmp <- dtf[j:(j+999)+1,i]
        ran <- quantile(tmp,c(prob/2,1-prob/2))
        if(true[i]<=ran[2] & true[i] >= ran[1]){count = count + 1}
      }
      coverage <- c(coverage,count/(dim(dtf)[1]/1000))
    }
    return(coverage)
  }
  write.csv(data[,1],'data.csv')
  write.csv(true_value,'value.csv')
  if(nchange*3+2 != length(true_value)){stop('Incorrect value for analysis')}
  mean <- colMeans(data)
  RMSE <- sqrt(sapply(1:length(true_value), FUN = function(x) mean((data[,x]-true_value[x])^2)))
  bias <- abs(mean-true_value)/true_value*100
  coverage <- coverate_cal(data,true_value,0.05)*100
  output <- data.frame(par_name = var_names, true_value = true_value,mean = mean, RMSE =RMSE, BIAS = bias, coverage = coverage)
  setwd('..')
  evtm_all[evtm_all==0]=NA
  write.csv(output,paste0("analysis_data",dataset,"_nchange",nchange,".csv"))
  write.csv(evtm_all,paste0("evetm_data",dataset,".csv"))
  # close(f)
}
