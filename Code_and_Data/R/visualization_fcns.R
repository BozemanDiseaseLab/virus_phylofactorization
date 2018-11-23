
### we need a function to grab the shortest unique taxonomy for each group.
tcollapse <- function(tax,tnms=NULL){
  if (is.null(tnms)){
    tnms <- names(tax)
  }
  return(paste(mapply(paste,tnms,tax,sep='_'),collapse='; '))
}

sutax <- function(IDs,Taxonomy,trim=F){
  nt <- length(IDs)
  nms <- colnames(Taxonomy)
  suts <- vector(mode='list',length=nt)
  if (!all(unlist(IDs) %in% Taxonomy[,1])){stop('Taxonomy improperly formatted. All IDs must be in first column of Taxonomy')}
  
  for (i in 1:nt){
    ingrp <- match(IDs[[i]],Taxonomy[,1])
    outgrp <- match(unlist(IDs[setdiff(1:nt,i)]),Taxonomy[,1])
    done=F
    ix=2
    while(!done){
      in_nms <- unique(Taxonomy[ingrp,ix])
      out_nms <- unique(Taxonomy[outgrp,ix])
      
      if (!any(in_nms %in% out_nms) | ix==ncol(Taxonomy)){
        done=T
        suts[[i]] <- unique(apply(Taxonomy[ingrp,2:ix,drop=F],1,tcollapse,tnms=colnames(Taxonomy[1,2:ix,drop=F])))
        if(trim){
          suts[[i]]  <- sapply(suts[[i]],FUN=function(x) strsplit(x,';')) %>% 
            sapply(.,FUN=function(x) x[length(x)]) %>% unique
        }
      } else {
        ix=ix+1
      }
    }
  }
  
  return(suts)
  
}

getName <- function(Clade,VOlival){
  V <- as.data.frame(VOlival[Clade,2:5,drop=F])
  ss=NULL
  if (nrow(V)==1){
    tx <- names(V)
    tx <- sapply(tx,FUN=function(x) gsub('v','',x))
    V <- sapply(V,toString)
    V <- mapply(paste,tx,V,sep='_')
    ss <- paste(V,collapse='; ')
  } else {
    for (i in ncol(V):1){
      if (length(unique(V[,i]))==1){
        if (!is.na(V[1,i])){
          tx <- gsub('v','',names(V[1,])[i])
          ss <- paste(tx,toString(V[1,i]),sep='_')
          break
        }
      }
    }
  }
  if (is.null(ss)){
    ss <- 'No Common Name'
  }
  return(ss)
}

phylobin <- function(B){
  bn <- numeric(sum(sapply(B,length)))
  for (i in 1:length(B)){
    bn[B[[i]]] <- i
  }
  return(bn)
}