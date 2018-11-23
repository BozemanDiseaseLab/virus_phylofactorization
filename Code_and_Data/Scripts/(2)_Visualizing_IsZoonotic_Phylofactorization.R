library(phylofactor)
library(ggtree)
library(grid)
library(gridBase)
library(ggplot2)
source('R/visualization_fcns.R')

load('Workspaces/IsZoonotic_phylofactorization')

pf <- pf.IsZoonotic
Z <- VOlival$IsZoonotic

pp=pf.tree(pf)

mean(Z[pf$groups[[10]][[2]]])

B <- bins(pf$basis[,1:10])
B <- B[2:11]  ## remove paraphyletic bin
nms <- sapply(B,getName,VOlival=VOlival)
probs <- sapply(B,FUN=function(ix,Z) mean(Z[ix]),Z=Z) %>% signif(.,digits=2)
names(nms) <- probs


Legend <- pp$legend
Legend$names <- nms
P <- sapply(probs,FUN=function(x) paste('p=',toString(signif(x,digits = 2)),sep=''))
Legend$names <- mapply(paste,Legend$names,probs)
plot.new()
legend('topleft',legend=Legend$names,fill=Legend$colors,cex=1)



# The Plot ---------------------------------------------------------------
gtr <- pp$ggplot+geom_tippoint(size=3*Z,col='black')

plot.new()
grid.newpage()
pushViewport(viewport(layout = grid.layout(1,2)))

pushViewport(viewport(layout.pos.col = 1))
print(gtr, newpage = FALSE)
popViewport()

#Draw bsae plot
pushViewport(viewport(layout.pos.col = 2))
par(fig = gridFIG(), new = TRUE)
plot.new()
legend('topleft',legend=Legend$names,fill=Legend$colors,cex=1.8)
popViewport()
ggsave('Figures/IsZoonotic_fig.png')