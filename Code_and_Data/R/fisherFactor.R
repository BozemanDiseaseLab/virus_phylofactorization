#' @param Z
#' @param tree
#' @param nfactors
#' @param random
#' @param prob.fcn
#' @param lambda
#' @param null
#' @param TestFunction must input {grps,tree,Z} and output P-value (minimum-objective function wins)
#' 

fisherFactor <- function(Z,tree,nfactors,TestFunction=NULL,ncores=NULL,cluster.depends=NULL,random=F,prob.fcn=NULL,lambda=-1,null=F){

  if(is.null(TestFunction)){
    TestFunction <- function(grps,tree,Z){
      s1 <- Z[tree$tip.label[grps[[1]]]]
      s2 <- Z[tree$tip.label[grps[[2]]]]
      n1 <- sum(s1)
      n2 <- sum(s2)
      ft <- fisher.test(matrix(c(n1,length(s1)-n1,n2,length(s2)-n2),ncol=2),alternative='two.sided')
      return(ft$p.value)
    }
  }
  
  if (is.null(prob.fcn)){
    prob.fcn <- function(pvals,lambda){
      return(pvals^lambda)
    }
  }
  
  if(null){
    p <- sum(Z)/length(Z)
    Z <- rbinom(length(Z),1,prob = p)
    names(Z) <- tree$tip.label
  }
  if (!is.null(ncores)){
    cl <- parallel::makeCluster(ncores)
    if (!is.null(cluster.depends)){
      parallel::clusterExport(cl,varlist = 'cluster.depends')
      parallel::clusterEvalQ(cl,expr = {cluster.depends()})
    }
  }
  
  treeList <- list(tree)
  binList <- list(1:ape::Ntip(tree))
  Grps=getGroups(tree)
  output <- NULL
  pfs=0
  while (pfs < min(length(Z)-1,nfactors)){
    
    if (pfs>=1){
      treeList <- updateTreeList(treeList,binList,grp,tree,skip.check=T)
      binList <- updateBinList(binList,grp)
      Grps <- getNewGroups(tree,treeList,binList)
    }
    
    if (is.null(ncores)){
      pvals <- sapply(Grps,FUN=TestFunction,tree=tree,Z=Z)
    } else {
      pvals <- parallel::parSapply(cl,Grps,TestFunction,tree=tree,Z=Z)
    }
    
    if (!random){
      best.grp <- Grps[[which.min(pvals)]]
      P <- min(pvals,na.rm = T)
    } else {
      probs <- prob.fcn(pvals,lambda)
      ix <- sample(length(Grps),1,prob = probs)
      best.grp <- Grps[[ix]]
      P <- pvals[ix]
    }
    if (is.na(P)){
      stop('WTF')
      break
    }
    output$pvals <- c(output$pvals,P)
    
    grp <- getLabelledGrp(tree=tree,Groups=best.grp)
    output$groups <- c(output$groups,list(best.grp))
    
    grpInfo <- matrix(c(names(grp)),nrow=2)
    output$factors <- cbind(output$factors,grpInfo)
    
    pfs=pfs+1
  }
  
  
  if (!is.null(ncores)){
    parallel::stopCluster(cl)
    rm('cl')
  }
    
    output$factors <- t(output$factors)
  
    V <- NULL
    for (i in 1:length(output$groups)){
      V <- cbind(V,ilrvec(output$groups[[i]],Ntip(tree)))
    }
    
  output$basis <- V
  output$bins <- bins(V)
  output$tree <- tree
  output$nfactors <- pfs
  output$Data <- matrix(Z,ncol=1)
  rownames(output$Data) <- tree$tip.label
  class(output) <- 'phylofactor'
  
  return(output)
}