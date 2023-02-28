# We calculate the morphological distance of the entire manus accross the entire phylogeny under an LDR hypothesis
tab.nodos <- phylogenies[[1]]$edge
LDR.met.dist <- list()
for (j in 1:100) {
  all.nodes <- 0
  for (y in 1:88) {
    pos.1 <- (tab.nodos[y,1])
    pos.2 <- (tab.nodos[y,2])
    sum.diff.sq <- 0
    diff.sq <- 0
    for (i in 1:5) {
      if (pos.1 >= 46) {
        char.val.1 <- data.LDR.met.list[[j]][[i]]$ace[which(pos.1 ==names(data.LDR.met.list[[j]][[i]]$ace))]
      } else {
        char.val.1 <-data.LDR.met[[pos.1, i]]
      }
      if (pos.2 >= 46) {
        char.val.2 <- data.LDR.met.list[[j]][[i]]$ace[which(pos.2 ==names(data.LDR.met.list[[j]][[i]]$ace))]
      } else {
        char.val.2 <- data.LDR.met[[pos.2, i]]
      }
      diff.sq <- ((char.val.1-char.val.2)^2)
      sum.diff.sq <- sum.diff.sq + diff.sq
    }
    dist.node <- unname(sqrt(sum.diff.sq))
    all.nodes <- all.nodes + dist.node
  }
  LDR.met.dist [[j]] <- all.nodes
}

# We calculate the morphological distance of the entire manus accross the entire phylogeny under a BDR hypothesis
BDR.met.dist <- list()
for (j in 1:100) {
  all.nodes <- 0
  for (y in 1:88) {
    pos.1 <- (tab.nodos[y,1])
    pos.2 <- (tab.nodos[y,2])
    sum.diff.sq <- 0
    diff.sq <- 0
    for (i in 1:5) {
      if (pos.1 >= 46) {
        char.val.1 <- data.BDR.met.list[[j]][[i]]$ace[which(pos.1 ==names(data.BDR.met.list[[j]][[i]]$ace))]
      } else {
        char.val.1 <-data.BDR.met[[pos.1, i]]
      }
      if (pos.2 >= 46) {
        char.val.2 <- data.BDR.met.list[[j]][[i]]$ace[which(pos.2 ==names(data.BDR.met.list[[j]][[i]]$ace))]
      } else {
        char.val.2 <- data.BDR.met[[pos.2, i]]
      }
      diff.sq <- ((char.val.1-char.val.2)^2)
      sum.diff.sq <- sum.diff.sq + diff.sq
    }
    dist.node <- unname(sqrt(sum.diff.sq))
    all.nodes <- all.nodes + dist.node
  }
  BDR.met.dist [[j]] <- all.nodes
}
#Next we plot the data as histograms
LDR.met.hist <- unlist(LDR.met.dist)
BDR.met.hist <- unlist(BDR.met.dist)
p1 <- hist(LDR.met.hist)
p2 <- hist(BDR.met.hist)
plot( p1, col=rgb(0,0,1,1/4), xlim=c(10,12.5))  # first histogram
plot( p2, col=rgb(1,0,0,1/4), xlim=c(10,12.5), add=T)  # second

#Next, we use t-test to test significant difference.
t.test.met <- t.test(LDR.met.hist, BDR.met.hist, paired = TRUE)
t.test.met
#Since p-value is smaller than 0.05, we can conclude that there is a significant difference.

#Next, we calculate Cohen's D to determine effect size.
library(lsr)
D.met <- cohensD(LDR.met.hist, BDR.met.hist)
D.met

vec.LDR.BDR <- rep(c("ADR", "BDR"), each = 100)
data.vec <- c(LDR.met.hist, BDR.met.hist)
violin.data.met <- data.frame(data.vec, vec.LDR.BDR)
violin.data.met
library(ggplot2)
p <- ggplot(violin.data.met, aes(x=vec.LDR.BDR, y=data.vec, fill = vec.LDR.BDR)) + 
  geom_violin(trim = F)
p + stat_summary(fun=mean, geom="point", shape=23, size=2)
p + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)


LDR.BDR.diff <- LDR.met.hist - BDR.met.hist
