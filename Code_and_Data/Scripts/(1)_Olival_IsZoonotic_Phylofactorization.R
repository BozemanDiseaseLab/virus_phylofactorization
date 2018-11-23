### Downloading Olival et al. dataset and scanning NCBI FTP genomes - who do we have genomes for?
library(RCurl)
library(ape)
library(phylofactor)
library(parallel)
library(dplyr)
library(ggpubr)
library(xlsx)
library(ggtree)
source('R/fisherFactor.R')
source('R/visualization_fcns.R')
parpf <- function(reps,pf,phat){
  Pvals <- matrix(NA,nrow=reps,ncol=pf$nfactors)
  for (nn in 1:reps){
    Z <- rbinom(Ntip(pf$tree),size=1,prob=phat)
    names(Z) <- pf$tree$tip.label
    pf <- fisherFactor(Z,pf$tree,nfactors=pf$nfactors)
    Pvals[nn,] <- pf$pvals
  }
  return(Pvals)
}

parf <- function(reps,pf,phat){
  Fs <- matrix(NA,nrow=reps,ncol=pf$nfactors)
  for (nn in 1:reps){
    Z <- rbinom(Ntip(pf$tree),size=1,prob=phat)
    names(Z) <- pf$tree$tip.label
    pf <- fisherFactor(Z,pf$tree,nfactors=10)
    Pvals[nn,] <- pf$pvals
  }
  return(Pvals)
}

load('Workspaces/Viral_taxonomy_tree.rdata')
### Data from Olival et al. (2017) ###
x <- getURL('https://raw.githubusercontent.com/ecohealthalliance/HP3/master/data/viruses.csv')
VOlival <- read.csv(text=x)
tree <- di2multi(tree) ## resolve polytomies
all.equal(tree$tip.label,sapply(as.list(VOlival$vVirusNameCorrected),toString))




# Phylofactorization of IsZoonotic ----------------------------------------
Z <- VOlival$IsZoonotic
names(Z) <- tree$tip.label
phat <- mean(Z)     #fraction of viruses which are zoonotic in Olival data

pf.IsZoonotic <- fisherFactor(Z,tree,nfactors=10)
pf.IsZoonotic$factors

nms <- lapply(pf.IsZoonotic$bins,FUN=function(a,tree) tree$tip.label[a],tree=tree)


################ Null Simulation ##############
ncores=7
reps=50
cl <- phyloFcluster(ncores)
clusterExport(cl,varlist = c('fisherFactor','parpf'))
reps=as.list(rep(reps,ncores))
Ps <- parLapply(cl,reps,fun=parpf,pf=pf.IsZoonotic,phat=mean(VOlival$IsZoonotic))
stopCluster(cl)
for (i in 1:length(Ps)){
  if (i==1){Pvals <- Ps[[i]]} else{Pvals <- rbind(Pvals,Ps[[i]])}
}
rm('Ps')




# Visualization -----------------------------------------------------------
ddf <- data.frame('Pvals'=c(t(Pvals)),'factor'=rep(1:10,nrow(Pvals)),'replicate'=rep(1:nrow(Pvals),each=10))
obs <- data.frame('Pvals'=pf.IsZoonotic$pvals,'factor'=1:10,'replicate'=NA)
ddf <- rbind(ddf,obs)

cols <- c('Null Simulations'='black','Observed Phylofactorization'='red')

PLOT <- ggplot(ddf,aes(factor,Pvals))+
          stat_summary(data = ddf[!is.na(ddf$replicate),],
                       geom='ribbon',
                       fun.ymin=function(x) quantile(x,0.05),
                       fun.ymax=function(x) quantile(x,1),
                       alpha=0.2,lwd=2)+
          stat_summary(data = ddf[!is.na(ddf$replicate),],
                       aes(x=factor,y=Pvals,color='Null Simulations'),geom='line',fun.y = mean,lwd=2)+
          geom_line(data=ddf[is.na(ddf$replicate),],
                    aes(y=Pvals,color='Observed Phylofactorization'),lwd=2)+
          scale_y_continuous(trans='log',breaks = 10^(c(-8:0)))+
          scale_x_continuous(breaks=1:10)+
          scale_fill_continuous(cols)+
          theme(legend.position = c(.5,.5),title = element_text(size=18))+
          ggtitle('Phylofactorization - observed vs. null')

GG <- pf.tree(pf.IsZoonotic,
              factor.map = data.frame(1:8,rep(1,8)),
              bg.color='grey',
              bg.alpha=0.5)$ggplot
  ggarrange(GG,PLOT,ncol = 2,nrow=1)
ggsave(filename='Figures/FisherFactors_of_IsZoonotic.tiff',height=8,width=16)
save(list=ls(),file='Workspaces/IsZoonotic_phylofactorization')




# Phylofactorization by stringent Zoon ------------------------------------
Z <- VOlival$IsZoonotic.stringent
names(Z) <- tree$tip.label
phat <- mean(Z)     #fraction of viruses which are zoonotic in Olival data

pf.Stringent <- fisherFactor(Z,tree,nfactors=10)

ncores=7
reps=50
cl <- phyloFcluster(ncores)
clusterExport(cl,varlist = c('fisherFactor','parpf'))
reps=as.list(rep(reps,ncores))
Ps <- parLapply(cl,reps,fun=parpf,pf=pf.Stringent,phat=phat)
stopCluster(cl)
for (i in 1:length(Ps)){
  if (i==1){Pvals.Stringent <- Ps[[i]]} else{Pvals.Stringent <- rbind(Pvals.Stringent,Ps[[i]])}
}
rm('Ps')


############## Visualization of significance & location of factors on phylogeny
nfactors <- pf.Stringent$nfactors
ddf <- data.frame('Pvals'=c(t(Pvals.Stringent)),
                  'factor'=rep(1:nfactors,nrow(Pvals.Stringent)),
                  'replicate'=rep(1:nrow(Pvals.Stringent),each=nfactors))
obs <- data.frame('Pvals'=pf.Stringent$pvals,'factor'=1:nfactors,'replicate'=NA)
ddf <- rbind(ddf,obs)

cols <- c('Null Simulations'='black','Observed Phylofactorization'='red')

PLOT <- ggplot(ddf,aes(factor,Pvals))+
  stat_summary(data = ddf[!is.na(ddf$replicate),],
               geom='ribbon',
               fun.ymin=function(x) quantile(x,0.05),
               fun.ymax=function(x) quantile(x,1),
               alpha=0.2,lwd=2)+
  stat_summary(data = ddf[!is.na(ddf$replicate),],
               aes(x=factor,y=Pvals,color='Null Simulations'),geom='line',fun.y = mean,lwd=2)+
  geom_line(data=ddf[is.na(ddf$replicate),],
            aes(y=Pvals,color='Observed Phylofactorization'),lwd=2)+
  scale_y_continuous(trans='log',breaks = 10^(c(-8:0)))+
  scale_x_continuous(breaks=1:10)+
  scale_fill_continuous(cols)+
  theme(legend.position = c(.5,.5),title = element_text(size=18))+
  ggtitle('Phylofactorization - observed vs. null')

GG <- pf.tree(pf.Stringent,
              bg.color='grey',
              bg.alpha=0.5)$ggplot
ggarrange(GG,PLOT,ncol = 2,nrow=1)
ggsave(filename='Figures/FisherFactors_of_IsZoonotic_Stringent.tiff',height=4,width=8)




# Summarizing results -----------------------------------------------------

Clades <- bins(pf.IsZoonotic$basis)
CladeNames <- sapply(Clades,getName,VOlival=VOlival)
probs <- sapply(Clades,FUN=function(b,Z) mean(Z[b]),Z)

output <- cbind(CladeNames,probs,sapply(Clades,length),sapply(Clades,length)*probs)
colnames(output) <- c('Name','Fraction_Zoonotic','no_viruses','no_zoonotic')
write.csv(output,file='Results/phylofactor_IsZoonotic_results.csv')


Clades <- bins(pf.Stringent$basis)
sCladeNames <- sapply(Clades,getName,VOlival=VOlival)
probs <- sapply(Clades,FUN=function(b,Z) mean(Z[b]),Z)

output <- cbind(sCladeNames,probs,sapply(Clades,length),sapply(Clades,length)*probs)
colnames(output) <- c('Name','Fraction_Zoonotic','no_viruses','no_zoonotic')
write.csv(output,file='Results/phylofactor_Stringent_results.csv')

save(list=ls(),file='Workspaces/IsZoonotic_phylofactorization')








# Appending & Writing phyloBin results ------------------------------------

########## creating phylobin variable in postprocessed databse ##############
DF <- readRDS('Olival/intermediates/postprocessed_database.rds')$viruses
all.equal(DF$vVirusNameCorrected,tree$tip.label) 

B <- bins(pf.IsZoonotic$basis)
nms <- sapply(B,getName,DF)
phylobin.IsZ <- phylobin(B) %>% nms[.]
DF$pfIsZ <- phylobin.IsZ

B <- bins(pf.Stringent$basis)
nms <- sapply(B,getName,DF)
phylobin.IsZS <- phylobin(B) %>% nms[.]
DF$pfIsZS <- phylobin.IsZS

### now we append the nonZ clades, i.e. those in bins. 
### In this case, we remove all those in bins 2-6, or the monophyletic clades of factors 1-5
nonZ <- sapply(pf.IsZoonotic$groups[1:5],FUN=function(g) g[[1]]) %>% unlist

DF$nonZ <- numeric(nrow(DF))
DF$nonZ[nonZ] <- 1


write.xlsx(DF,file='Datasets/Olival_w_phylobin.xlsx')
saveRDS(DF,file='Workspaces/Olival_w_phylobin')