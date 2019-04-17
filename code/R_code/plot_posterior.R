############################### Reading data ####################################
setwd("C:/Users/jerry/Dropbox/paper_PLP/coal_Data")
library(ggplot2)
post_data <- read.csv("simulation_data2_change1.csv",header = T)
exp_name=c(expression(paste(alpha[1])),expression(paste(alpha[2])),
           expression(paste(beta[1])),expression(paste(beta[2])),
           expression(paste(tau)))

############# Plot the posterior distribution of intensity parameter###############
a_min <- round(min(post_data$a1,post_data$a2)-0.05,2)
a_max <- round(max(post_data$a1,post_data$a2)+0.05,2)
a_den <- c(density(post_data$a1)$y,density(post_data$a2)$y)
a_den_max <- round(max(a_den)+0.5,1)
# b_max <- round(max(post_data$beta.1.,post_data$beta.2.)+0.5)
b_max <- max(post_data$b1, post_data$b2)
b_den <- c(density(post_data$b1)$y,density(post_data$b1)$y)
b_den_max <- round(max(b_den)+0.0005,4)
col = 1
name <- names(post_data)[col]
ldata <- data.frame(simulation = post_data[,col])
p1 <- ggplot(ldata,aes(x=simulation, y=..density..))+
  theme(plot.title = element_text(hjust = 0.5,size = 20),
        text = element_text(size=20))+ylab("Density")+xlab(exp_name[col])+
  scale_y_continuous(breaks = c(0,a_den_max/2,a_den_max),limits = c(0,a_den_max))+
  scale_x_continuous(breaks = c(a_min,(a_min+a_max)/2,a_max),limits = c(a_min,a_max))+
  geom_histogram(fill="white", colour = "black",bins=60)+
  geom_density(alpha=.2)

col = 2
name <- names(post_data)[col]
ldata <- data.frame(simulation = post_data[,col])
p2 <- ggplot(ldata,aes(x=simulation, y=..density..))+
  theme(plot.title = element_text(hjust = 0.5,size = 20),
        text = element_text(size=20))+ylab("Density")+xlab(exp_name[col])+
  scale_y_continuous(breaks = c(0,a_den_max/2,a_den_max),limits = c(0,a_den_max))+
  scale_x_continuous(breaks = c(a_min,(a_min+a_max)/2,a_max),limits = c(a_min,a_max))+
  geom_histogram(fill="white", colour = "black",bins=60)+
  geom_density(alpha=.2, fill="#FF6666")
b_max=80
col = 3
name <- names(post_data)[col]
ldata <- data.frame(simulation = post_data[,col])
p3 <- ggplot(ldata,aes(x=simulation, y=..density..))+
  theme(plot.title = element_text(hjust = 0.5,size = 20),
        text = element_text(size=20))+ylab("Density")+xlab(exp_name[col])+
  scale_y_continuous(breaks = c(0,b_den_max/2,b_den_max),limits = c(0,b_den_max))+
  scale_x_continuous(breaks = c(0,b_max/2,b_max),limits = c(0,b_max))+
  geom_histogram(fill="white", colour = "black",bins=60)+
  geom_density(alpha=.2, fill="#FF6666")

col = 4
name <- names(post_data)[col]
ldata <- data.frame(simulation = post_data[,col])
p4 <- ggplot(ldata,aes(x=simulation, y=..density..))+
  theme(plot.title = element_text(hjust = 0.5,size = 20),
        text = element_text(size=20))+ylab("Density")+xlab(exp_name[col])+
  scale_y_continuous(breaks = c(0,b_den_max/2,b_den_max),limits = c(0,b_den_max))+
  scale_x_continuous(breaks = c(0,b_max/2,b_max),limits = c(0,b_max))+
  geom_histogram(fill="white", colour = "black",bins=60)+
  geom_density(alpha=.2, fill="#FF6666")
require(gridExtra)
setEPS()
postscript("poster1.eps")
p = grid.arrange(p1,p2,p3,p4, nrow=2, ncol=2)
dev.off()


####################### Plot the posterior distribution of change point#####################
col = 5
name <- names(post_data)[col]
ldata <- data.frame(simulation = post_data[,col])
t_min <- round(min(ldata$simulation))-10
t_max <- round(max(ldata$simulation))+10
brks <- c(t_min,(t_min+t_max)/2,t_max-10)
setEPS()
postscript("poster2.eps")
ggplot(ldata,aes(x=simulation, y=..density..))+
  theme(plot.title = element_text(hjust = 0.5,size = 20),
        text = element_text(size=20))+ylab("Density")+xlab(exp_name[col])+
  scale_x_continuous(breaks = brks,
                     labels = as.Date("1850-01-01")+7*3*4*brks,limits=c(t_min,t_max))+
  geom_histogram(bins=80, fill="white", colour = "black")
dev.off()
