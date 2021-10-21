library(edgeR)

# generate fake counts for drug1 and drug2
counts_drug1<-data.frame(replicate(182,sample(0:1,1000,rep=TRUE)))
counts_drug2<-data.frame(replicate(182,sample(0:1,1000,rep=TRUE)))

# time, concentration and replicate
time<-c(2,4,8,24,48)
conc<-c(0, 1, 10, 100, 1000)
replicates <- c(1,2,3,4,5,6,7)

time_s <- paste0("hrs",time)
conc_s <- paste0("c",conc)
replicates_s <- paste0("r",replicates)

# build condition matrix

df<-expand.grid(time_s,conc_s,replicates_s)

# append zero time point
zero_point <- expand.grid(c("hrs0"),c("c0"), replicates_s)
df <- rbind(df,zero_point)
names(df)<-c("Time","Conc","Rep")

group <- paste(df$Time,df$Conc,sep=".")
rep_names <- paste(df$Time,df$Conc,df$Rep,sep=".")
names(counts_drug1)<-rep_names
names(counts_drug2)<-rep_names

# build design matrix
design <- model.matrix(~0+group)
colnames(design) <- levels(factor(group))
y <- DGEList(counts=counts_drug1)
y <- estimateDisp(y, design)

fit <- glmQLFit(y, design)

# Изменение экспрессии при изменении концентрации при одном и том же времени воздействия
con <- makeContrasts(
  c1vsc0 = ((hrs4.c1 - hrs2.c1)  + (hrs8.c1 - hrs4.c1) + (hrs24.c1 - hrs8.c1)  + (hrs48.c1 - hrs24.c1))/4 - ((hrs4.c0 - hrs2.c0)  + (hrs8.c0 - hrs4.c0) + (hrs24.c0 - hrs8.c0)  + (hrs48.c0 - hrs24.c0))/4,
  c10vsc1 = ((hrs4.c10 - hrs2.c10)  + (hrs8.c10 - hrs4.c10) + (hrs24.c10 - hrs8.c10)  + (hrs48.c10 - hrs24.c10))/4 - ((hrs4.c1 - hrs2.c1)  + (hrs8.c1 - hrs4.c1) + (hrs24.c1 - hrs8.c1)  + (hrs48.c1 - hrs24.c1))/4,
  c100vsc10 = ((hrs4.c100 - hrs2.c100)  + (hrs8.c100 - hrs4.c100) + (hrs24.c100 - hrs8.c100)  + (hrs48.c100 - hrs24.c100))/4 - ((hrs4.c10 - hrs2.c10)  + (hrs8.c10 - hrs4.c10) + (hrs24.c10 - hrs8.c10)  + (hrs48.c10 - hrs24.c10))/4,
  c1000vsc100 =((hrs4.c1000 - hrs2.c1000)  + (hrs8.c1000 - hrs4.c1000) + (hrs24.c1000 - hrs8.c1000)  + (hrs48.c1000 - hrs24.c1000))/4 - ((hrs4.c100 - hrs2.c100)  + (hrs8.c100 - hrs4.c100) + (hrs24.c100 - hrs8.c100)  + (hrs48.c100 - hrs24.c100))/4,
  levels=design)

qlf <- glmQLFTest(fit, contrast=con)

topTags(qlf)

# Изменение экспрессии при изменении времени независимо от концентрации
con <- makeContrasts(
  hrs0vshrs2 = (hrs2.c1 + hrs2.c10 + hrs2.c100 + hrs2.c1000)/4 - hrs0.c0,
  hrs2vshrs4 = (hrs4.c1 + hrs4.c10 + hrs4.c100 + hrs4.c1000)/4 - (hrs2.c1 + hrs2.c10 + hrs2.c100 + hrs2.c1000)/4,
  hrs4vshrs8 = (hrs8.c1 + hrs8.c10 + hrs8.c100 + hrs8.c1000)/4 - (hrs4.c1 + hrs4.c10 + hrs4.c100 + hrs4.c1000)/4,
  hrs8vshrs24 = (hrs24.c1 + hrs24.c10 + hrs24.c100 + hrs24.c1000)/4 - (hrs8.c1 + hrs8.c10 + hrs8.c100 + hrs8.c1000)/4,
  hrs24vshrs48 = (hrs48.c1 + hrs48.c10 + hrs48.c100 + hrs48.c1000)/4 - (hrs24.c1 + hrs24.c10 + hrs24.c100 + hrs24.c1000)/4,
  levels=design)

qlf <- glmQLFTest(fit, contrast=con)
topTags(qlf)
