
# phylobin classificaiton -------------------------------------------------
source('R/visualization_fcns.R')
load('Workspaces/cb_dist_and_IsZoonotic_phylofactorization')
library(mgcv)
library(stringi)
library(dplyr)
library(purrr)
library(ReporteRs)
library(visreg)
library(gridExtra)
library(ggpubr)
library(plyr)
library(xlsx)
library(ggplot2)
library(parallel)
library(multcompView)
library(phylofactor)
x <- 'https://raw.githubusercontent.com/ecohealthalliance/HP3/master'
script <- getURL(paste(x,'/R/relative_contributions.R',sep=''),ssl.verifypeer = FALSE) %>%
  parse(text=.)
eval(script)
load('Workspaces/phylofactor_comparisons')  ## gams from before

Vmod <- vtraits$model[[1]]
Vmod_new <- vtraits_new$model[[1]]
Clade <- pf.cb_dist_noHoSa_maxLn$basis %>% bins %>% phylobin %>% nms[.] %>% factor
vnames <- names(Z)

for (original in c(T,F)){
  db <- readRDS("Datasets/Olival_postprocessed_database.rds")
  
  # our modification
  if (!original){
    V <- read.xlsx('Datasets/Olival_w_phylobin_final.xlsx',sheetIndex = 1)
    db$viruses <- cbind(db$viruses,V[,c('pfIsZ','pfIsZS','nonZ')])
    db$viruses <- db$viruses[!as.logical(db$viruses$nonZ),]
  }
  
  db$viruses$Clade <- Clade[match(db$viruses$vVirusNameCorrected,vnames)]

  
  hosts <- db$hosts
  viruses <- db$viruses
  associations <- db$associations
  rm(db)
  PD_centers = names(viruses)[stri_detect_regex(names(viruses), "^(cb|st)_.*(?<!stringent)$")]
  
  model_family = binomial
  
  
  data_set = viruses %>%
    mutate_each_(funs("logp"), vars=PD_centers)
  names(data_set)[names(data_set) %in% PD_centers] <- paste0(PD_centers, "Ln")
  PD_centers <- paste0(PD_centers, "Ln")
  
  dummys = as.data.frame(with(viruses, model.matrix(~vFamily))[,-1])
  data_set = cbind(data_set, dummys)
  dummy_terms = paste0("s(", names(dummys), ", bs = 're')")
  names(dummy_terms) <- names(dummys)
  
  
  terms = list(
    PD     = paste0("s(", PD_centers, ", bs='tp', k=7)"),
    bias   = c("s(vPubMedCitesLn, bs='tp', k=7)", "s(vWOKcitesLn, bs='tp', k=7)"),
    strand = c("s(RNA, bs='re')", "s(SS, bs='re')", "s(vCytoReplicTF, bs='re')"),
    vector = "s(Vector, bs='re')",
    env = "s(Envelope, bs='re')",
    cld = "s(Clade,bs='re')",
    genome = "s(vGenomeAveLengthLn, bs='tp', k=7)",
    stringsAsFactors=FALSE)
  
  if (original){
    vtraits_Clade = fit_all_gams(data_set,
                           outcome_variable = "IsZoonotic",
                           model_family = binomial,
                           terms)
    data_set_old <- data_set
  } else {
    data_set_new <- data_set
    vtraits_Clade_new = fit_all_gams(data_set,
                               outcome_variable = "IsZoonotic",
                               model_family = binomial,
                               terms)
  }
}



# Visualiztion of Clades --------------------------------------------------

gam_Clade <- vtraits_Clade$model[[1]]
gam_Clade_new <- vtraits_Clade_new$model[[1]]

ix <- setdiff(1:nrow(data_set_old),gam_Clade$na.action)
ddf <- data.frame('IsZoonosis'=gam_Clade$linear.predictors,'Clade'=Clade[ix])

ord <- order(by(ddf$IsZoonosis,ddf$Clade,mean),decreasing = T)
ddf$Clade <- factor(ddf$Clade,levels = levels(ddf$Clade)[ord])



########## Tukey test ############

a <- aov(IsZoonosis~Clade, data=ddf)
tHSD <- TukeyHSD(a, ordered = FALSE, conf.level = 0.95)

generate_label_df <- function(HSD, flev){
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- HSD[[flev]][,4]
  Tukey.labels <- multcompLetters(Tukey.levels)['Letters']
  plot.labels <- names(Tukey.labels[['Letters']])
  
  # Get highest quantile for Tukey's 5 number summary and add a bit of space to buffer between    
  # upper quantile and label placement
  boxplot.df <- ddply(ddf, flev, function (x) max(fivenum(x$IsZoonosis)) + 0.2)
  
  # Create a data frame out of the factor levels and Tukey's homogenous group letters
  plot.levels <- data.frame(plot.labels, labels = Tukey.labels[['Letters']],
                            stringsAsFactors = FALSE)
  
  # Merge it with the labels
  labels.df <- merge(plot.levels, boxplot.df, by.x = 'plot.labels', by.y = flev, sort = FALSE)
  
  return(labels.df)
}

###############


labels.df <- generate_label_df(tHSD, 'Clade')
labels.df$plot.labels <- factor(labels.df$plot.labels,levels = levels(ddf$Clade))
names(labels.df)[1] <- 'Clade'

ggplot(data=ddf,aes(x=Clade,y=IsZoonosis,fill=Clade))+
  geom_boxplot()+
  geom_text(data = labels.df, aes(x = Clade, y = V1, label = labels))+
  theme(legend.text = element_text(size=15))
                                                       
# ggsave('Figures/cb_dist_Clades_IsZoonotic_by_original_model_selection.tiff',height=5,width=13)
