library(mgcv)
library(stringi)
library(dplyr)
library(purrr)
library(ReporteRs)
library(visreg)
library(gridExtra)
library(ggpubr)
x <- 'https://raw.githubusercontent.com/ecohealthalliance/HP3/master'
script <- getURL(paste(x,'/R/relative_contributions.R',sep=''),ssl.verifypeer = FALSE) %>%
  parse(text=.)
eval(script)load('Workspaces/phylofactor_comparisons')
rm('script')
dd <- readRDS('Workspaces/Olival_w_phylobins_final')
comp <- function(x) x[2]-x[1]
clo <- function(x) x/sum(x)
ilogit <- function(x) 1/(1+exp(-x))

Vmod <- vtraits$model[[1]]
Vmod_new <- vtraits_new$model[[1]]

summary(Vmod)
summary(Vmod_new)


############## Olival and Phylobin-removal datasets #############
D <- data_set_old
D$Envelope <- factor(D$Envelope)
D$Vector <- factor(D$Vector)
D$vCytoReplicTF <- factor(D$vCytoReplicTF)

D_new <-data_set_new
D_new$Envelope <- factor(D_new$Envelope)
D_new$Vector <- factor(D_new$Vector) 
D_new$vCytoReplicTF <- factor(D_new$vCytoReplicTF)

########### Olival gam ############
gam_old <- gam(formula=as.formula(Vmod$formula),data=D,family = binomial)
########### Olival formula, new dataset ##########
gam_compare <- gam(formula=as.formula(Vmod$formula),data=D_new,family = binomial)

tbl <- matrix(NA,nrow=3,ncol=3)
colnames(tbl) <- c('Original','PF_restricted','NonZ')
rownames(tbl) <- c('Vector','Envelope','CytoRepl')

ix <- setdiff(1:nrow(D),gam_old$na.action)
tbl['Vector','Original'] <- by(gam_old$linear.predictors, D$Vector[ix],mean) %>% comp
tbl['Envelope','Original'] <- by(gam_old$linear.predictors, D$Envelope[ix],mean) %>%  comp
tbl['CytoRepl','Original'] <- by(gam_old$linear.predictors, D$vCytoReplicTF[ix],mean) %>%  comp



ix <- setdiff(1:nrow(D_new),gam_compare$na.action)
tbl['Vector','PF_restricted'] <- by(gam_compare$linear.predictors, D_new$Vector[ix],mean) %>%  comp
tbl['Envelope','PF_restricted'] <- by(gam_compare$linear.predictors, D_new$Envelope[ix],mean) %>%  comp
tbl['CytoRepl','PF_restricted'] <- by(gam_compare$linear.predictors, D_new$vCytoReplicTF[ix],mean) %>%  comp

D$nonZ <- dd$nonZ
gNonZ <- gam(nonZ~Vector+Envelope+vCytoReplicTF+cb_dist_noHoSa_maxLn,data=D,family=binomial)
ix <- setdiff(1:nrow(D),gNonZ$na.action)
tbl['Vector','NonZ'] <- by(gNonZ$linear.predictors, D_new$Vector[ix],mean) %>%  comp
tbl['Envelope','NonZ'] <- by(gNonZ$linear.predictors, D_new$Envelope[ix],mean) %>% comp
tbl['CytoRepl','NonZ'] <- by(gNonZ$linear.predictors, D_new$vCytoReplicTF[ix],mean) %>% comp

tbl <- signif(tbl,digits=2)

# write.csv(tbl,'Results/Trait_Analysis_Comparisons_with_Non_Zoonotic_Clades.csv')



### GGplot will make a much more honest figure
ix <- ! 1:nrow(D) %in% gam_old$na.action
preds <- ilogit(gam_old$linear.predictors-mean(gam_old$linear.predictors))
DF_old <- cbind(D[ix,c('Vector','Envelope','cb_dist_noHoSa_maxLn','IsZoonotic','IsZoonotic.stringent')],
                'version'=rep('Olival',sum(ix)),'Normalized_Prediction'=preds)

ix <-  ! 1:nrow(D_new) %in% gam_compare$na.action
preds <- ilogit(gam_compare$linear.predictors-mean(gam_compare$linear.predictors))
DF_new <- cbind(D_new[ix,c('Vector','Envelope','cb_dist_noHoSa_maxLn','IsZoonotic','IsZoonotic.stringent')],
              'version'=rep('Phylobin_removal',sum(ix)),'Normalized_Prediction'=preds)

DF <- rbind(DF_old,DF_new)
DF$version <- factor(DF$version)

ggV <- ggplot(DF,aes(x=Vector,y=Normalized_Prediction,fill=version))+
        geom_boxplot()+
        ggtitle(label='Vector-Borne')+
        facet_grid(.~version)

ggE <- ggplot(DF,aes(x=Envelope,y=Normalized_Prediction,fill=version))+
        geom_boxplot()+
        ggtitle(label='Enveloped')+
        facet_grid(.~version)

ggP <- ggplot(DF,aes(x=cb_dist_noHoSa_maxLn,y=Normalized_Prediction,fill=version,color=version))+
        geom_point(size=5)+
        scale_x_continuous(name='Non-Human Phylogenetic Breadth',)+
        geom_smooth(size=3)+
        ggtitle(label='Phylogenetic Breadth')+
        theme(text = element_text(size=40),axis.text = element_text(size=40))
        
ggP
# ggsave('Figures/Phylogenetic_Breadth_vs_IsZoonotic.png',height=12,width=18)

# tiff('~/Bozeman/BZDEL/Figures/Olival phylofactorization/Envelope_Vector_Comparison_w_and_without_phylobin_removal.tiff',height=1200,width=400)       
# grid.arrange(ggV,ggE,ggP,layout.matrix=matrix(c(1,2,3),ncol=1))
ggarrange(ggV,ggE,ggP,ncol = 1,nrow=3)
ggarrange(ggV,ggE,ggP,ncol = 3,nrow=1)