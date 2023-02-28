# We calculate the morphological distance of the entire manus accross the entire phylogeny under an LDR hypothesis
tab.nodos <- phylogenies[[1]]$edge
LDR.full.dist <- list()
for (j in 1:100) {
  all.nodes <- 0
  for (y in 1:88) {
    pos.1 <- (tab.nodos[y,1])
    pos.2 <- (tab.nodos[y,2])
    sum.diff.sq <- 0
    diff.sq <- 0
    for (i in 1:20) {
      if (pos.1 >= 46) {
        char.val.1 <- data.LDR.full.list[[j]][[i]]$ace[which(pos.1 ==names(data.LDR.full.list[[j]][[i]]$ace))]
      } else {
        char.val.1 <-data.LDR.full[[pos.1, i]]
      }
      if (pos.2 >= 46) {
        char.val.2 <- data.LDR.full.list[[j]][[i]]$ace[which(pos.2 ==names(data.LDR.full.list[[j]][[i]]$ace))]
      } else {
        char.val.2 <- data.LDR.full[[pos.2, i]]
      }
      diff.sq <- ((char.val.1-char.val.2)^2)
      sum.diff.sq <- sum.diff.sq + diff.sq
    }
    dist.node <- unname(sqrt(sum.diff.sq))
    all.nodes <- all.nodes + dist.node
  }
  LDR.full.dist [[j]] <- all.nodes
}

# We calculate the morphological distance of the entire manus accross the entire phylogeny under a BDR hypothesis
BDR.full.dist <- list()
for (j in 1:100) {
  all.nodes <- 0
  for (y in 1:88) {
    pos.1 <- (tab.nodos[y,1])
    pos.2 <- (tab.nodos[y,2])
    sum.diff.sq <- 0
    diff.sq <- 0
    for (i in 1:20) {
      if (pos.1 >= 46) {
        char.val.1 <- data.BDR.full.list[[j]][[i]]$ace[which(pos.1 ==names(data.BDR.full.list[[j]][[i]]$ace))]
      } else {
        char.val.1 <-data.BDR.full[[pos.1, i]]
      }
      if (pos.2 >= 46) {
        char.val.2 <- data.BDR.full.list[[j]][[i]]$ace[which(pos.2 ==names(data.BDR.full.list[[j]][[i]]$ace))]
      } else {
        char.val.2 <- data.BDR.full[[pos.2, i]]
      }
      diff.sq <- ((char.val.1-char.val.2)^2)
      sum.diff.sq <- sum.diff.sq + diff.sq
    }
    dist.node <- unname(sqrt(sum.diff.sq))
    all.nodes <- all.nodes + dist.node
  }
  BDR.full.dist [[j]] <- all.nodes
}
#Next we plot the data as histograms
LDR.full.hist <- unlist(LDR.full.dist)
BDR.full.hist <- unlist(BDR.full.dist)
p1 <- hist(LDR.full.hist)
p2 <- hist(BDR.full.hist)
plot( p1, col=rgb(0,0,1,1/4), xlim=c(31.0,38.0))  # first histogram
plot( p2, col=rgb(1,0,0,1/4), xlim=c(31.0,38.0), add=T)  # second

#Next, we use t-test to test significant difference.
t.test.full <- t.test(LDR.full.hist, BDR.full.hist, paired = TRUE)
t.test.full
#Since p-value is smaller than 0.05, we can conclude that there is a significant difference.

#Next, we calculate Cohen's D to determine effect size.
library(lsr)
D.full <- cohensD(LDR.full.hist, BDR.full.hist)
D.full

vec.LDR.BDR <- rep(c("ADR", "BDR"), each = 100)
data.vec <- c(LDR.full.hist, BDR.full.hist)
violin.data.full <- data.frame(data.vec, vec.LDR.BDR)
violin.data.full
library(ggplot2)
p <- ggplot(violin.data.full, aes(x=vec.LDR.BDR, y=data.vec, fill = vec.LDR.BDR)) + 
  geom_violin(trim = F)
p + stat_summary(fun=mean, geom="point", shape=23, size=2)
p + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)


LDR.BDR.diff <- LDR.full.hist - BDR.full.hist
