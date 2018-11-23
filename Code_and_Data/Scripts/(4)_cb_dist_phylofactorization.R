############ Olival cb_dist_noHoSa_maxLn factorization
library(phylofactor)
library(xlsx)
library(stringi)
library(dplyr)
library(parallel)
library(ggplot2)
library(ggpubr)
library(viridis)
source('R/fisherFactor.R')
source('R/visualization_fcns.R')

################## Loading files from Olival et al. HP3 repository #################
x <- 'https://raw.githubusercontent.com/ecohealthalliance/HP3/master'

script <- getURL(paste(x,'/R/model_reduction.R',sep=''),ssl.verifypeer = FALSE) %>%
  parse(text=.)
eval(script)
script <- getURL(paste(x,'/R/fit_gam.R',sep=''),ssl.verifypeer = FALSE) %>%
  parse(text=.)
eval(script)
script <- getURL(paste(x,'/R/relative_contributions.R',sep=''),ssl.verifypeer = FALSE) %>%
  parse(text=.)
eval(script)
script <- getURL(paste(x,'/R/cross_validation.R',sep=''),ssl.verifypeer = FALSE) %>%
  parse(text=.)
eval(script)
script <- getURL(paste(x,'/R/logp.R',sep=''),ssl.verifypeer = FALSE) %>%
  parse(text=.)
eval(script)
rm('script')

###########################
load('Workspaces/IsZoonotic_phylofactorization')
db <- readRDS("Datasets/Olival_postprocessed_database.rds")
WilcoxTest <- function(grps,tree,Z){
  
  s1 <- Z[tree$tip.label[grps[[1]]]] %>% na.omit
  s2 <- Z[tree$tip.label[grps[[2]]]] %>% na.omit
  pval <- tryCatch(wilcox.test(s1,s2)$p.value,
                   error=function(e) 1)
  return(pval)
}
parpf.cb <- function(reps,pf){
  Pvals <- matrix(NA,nrow=reps,ncol=pf$nfactors)
  for (nn in 1:reps){
    Z <- sample(pf$Data,size=length(pf$Data),replace=T)
    names(Z) <- pf$tree$tip.label
    pf <- fisherFactor(Z,pf$tree,nfactors=pf$nfactors,TestFunction = WilcoxTest)
    Pvals[nn,] <- pf$pvals
  }
  return(Pvals)
}

all.equal(tree$tip.label,sapply(as.list(db$viruses$vVirusNameCorrected),toString))

hosts <- db$hosts
viruses <- db$viruses
associations <- db$associations
rm(db)
PD_centers = names(viruses)[stri_detect_regex(names(viruses), "^(cb|st)_.*(?<!stringent)$")]

data_set = viruses %>%
  mutate_each_(funs("logp"), vars=PD_centers)
names(data_set)[names(data_set) %in% PD_centers] <- paste0(PD_centers, "Ln")

Z <- data_set$cb_dist_noHoSa_maxLn
names(Z) <- data_set$vVirusNameCorrected




pf.cb_dist_noHoSa_maxLn <- fisherFactor(Z,tree,nfactors=9,TestFunction = WilcoxTest)


# Null Simulation ---------------------------------------------------------

ncores=7
reps=50
cl <- phyloFcluster(ncores)
clusterExport(cl,varlist = c('fisherFactor','parpf','WilcoxTest'))
reps=as.list(rep(reps,ncores))
Ps <- parLapply(cl,reps,fun=parpf.cb,pf=pf.cb_dist_noHoSa_maxLn)
stopCluster(cl)
rm('cl')
for (i in 1:length(Ps)){
  if (i==1){Pvals.cb <- Ps[[i]]} else{Pvals.cb <- rbind(Pvals.cb,Ps[[i]])}
}
rm('Ps')


# Visualization -----------------------------------------------------------

nfactors <- pf.cb_dist_noHoSa_maxLn$nfactors
ddf <- data.frame('Pvals'=c(t(Pvals.cb)),
                  'factor'=rep(1:nfactors,nrow(Pvals.cb)),
                  'replicate'=rep(1:nrow(Pvals.cb),each=nfactors))
obs <- data.frame('Pvals'=pf.cb_dist_noHoSa_maxLn$pvals,
                  'factor'=1:nfactors,
                  'replicate'=NA)
ddf <- rbind(ddf,obs)

cols <- c('Null Simulations'='black','Observed Phylofactorization'='red')

PLOT <- ggplot(ddf,aes(x=factor,y=Pvals))+
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

GG <- pf.tree(pf.cb_dist_noHoSa_maxLn,
              bg.color='grey',
              bg.alpha=0.5)$ggplot
ggarrange(GG,PLOT,ncol = 2,nrow=1)
ggsave(filename='Figures/Wilcox_factorization_cb_dist.tiff',height=8,width=16)




# cb_dist by bin ----------------------------------------------

B <- bins(pf.cb_dist_noHoSa_maxLn$basis)
B_means <- sapply(B,FUN=function(b,z) mean(z[b],na.rm=T),z=Z)
nms <- sapply(B,getName,VOlival=data_set)
nms[5] <- "Order_Mononegavirales"   #from manual look - some Rhabdoviruses were listed as "Unclassified".
nms[8] <- 'Subfamily_Alphaherpesvirinae' #again, from manual look.
cbind(B_means,nms)
nm.order <- nms[order(B_means,decreasing = T)]

colfcn <- function(n)  rainbow(10)
GG <- pf.tree(pf.cb_dist_noHoSa_maxLn,Grps = B,color.fcn = colfcn,top.layer = 5,top.alpha = 0.6,bg.alpha = 1)
GTREE <- GG$ggplot+
  ggtitle('Phylogenetic Breadth Factorization')


cb.df <- data.frame('Phylogenetic_Breadth'=Z,'name'=factor(nms[phylobin(B)],levels = rev(nm.order)))

cols <- rainbow(11)[order(B_means,decreasing = F)]
alphas <- rep()
PLOT <- ggplot(cb.df,aes(x=name,y=Phylogenetic_Breadth,fill=name))+
          geom_boxplot()+
          scale_fill_manual(values=cols,guide=guide_legend(reverse=T))+
          coord_flip()+
          ggtitle('Phylogenetic Breadth')

ggarrange(GTREE,PLOT,ncol=2)
ggsave('Figures/cb_dist_Phylofactorization.tiff',height=4,width=13)


# save(list=ls(),file='Workspaces/cb_dist_and_IsZoonotic_phylofactorization')




# Appending Olival dataset with cb_dist phylobin --------------------------
DF <- readRDS('Workspaces/Olival_w_phylobin')

B <- bins(pf.cb_dist_noHoSa_maxLn$basis)
nms <- sapply(B,getName,DF)
nms[5]  <- "Order_Mononegavirales"
nms[8] <- "Subfamily_Alphaherpesvirinae"
phylobin.cbdist <- phylobin(B) %>% nms[.]
DF$pfCBdist <- factor(phylobin.cbdist)

# write.xlsx(DF,'Datasets/Olival_w_phylobin_final.xlsx',sheetName = 'Viruses')
# saveRDS(DF,file='Workspaces/Olival_w_phylobins_final')
