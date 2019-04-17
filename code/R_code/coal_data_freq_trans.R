#################################Transform the coal data frequency###############################################
setwd("C:/Users/jerry/Dropbox/paper_PLP/coal_Data")
raw_data <- read.table("data24var.csv",sep=",",header = T)
dataset_s <- raw_data$location_id
max=40
etm <- c()
for(dataset in unique(dataset_s)){
  evetm <- as.Date(raw_data[dataset_s == dataset,"eventtime"])
  et <- rep(0,40)
  et[1:length(evetm)] <- sapply(1:length(evetm), FUN = function(x) difftime(evetm[x],as.Date("1850-01-01")))
  etm <- rbind(etm,et)
}

etm[is.na(etm)] <- 0
etm=etm/7
C <- rep(max(etm) + 1, dim(etm)[1])
nt <- t(sapply(1:dim(etm)[1], FUN = function(x) sum(etm[x,] != 0)))

write.table(etm, "eventtm1.csv", row.names = F, col.names =F,sep=',')
write.table(C, "censortm1.csv", row.names = F, col.names = F,sep=',')
write.table(nt, "NT1.csv", row.names = F, col.names = F,sep=',')
