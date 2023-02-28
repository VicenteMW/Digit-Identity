#ldp4 LDR
PCA.LDR.full <- prcomp(data.LDR.full)
PCA.BDR.full <- prcomp(data.BDR.full)

# Load package
library("factoextra")
# Create groups
group <- c(rep("Non-averostra", times=9), rep("Ceratosauria", times=4), rep("Tetanurae", times=32))
#group <- c(rep("Non-averostra", times=9), rep("Ceratosauria", times=3), rep("Allosauroidea", times=4), rep("Tyrannosauroidea", times=3), rep("Compsognathidae", times=1), rep("Ornithomimosauridae", times=5), rep("Maniraptora", times= 15)) 
# Plot
fviz_pca_biplot(PCA.LDR.full, repel=TRUE, pointsize=2, pointshape=21, col.var="red", arrowsize=0.6, labelsize=3, col.ind=group, palette=c("green2", "gold", "skyblue2"), addEllipses=TRUE, ellipse.type="convex")
fviz_pca_biplot(PCA.BDR.full, repel=TRUE, pointsize=2, pointshape=21, col.var="red", arrowsize=0.6, labelsize=3, col.ind=group, palette=c("green2", "gold", "skyblue2"), addEllipses=TRUE, ellipse.type="convex")

phy.PCA.1 <- phylogenies[[1]]
vec.1 <- rep(c(1), times=88)
phy.PCA.1$edge.length <- vec.1
library(geomorph)
tree.test <- read.tree("tree.test.txt")
PCA.w.phylo.LDR <- gm.prcomp(data.LDR.full, phy = phy.PCA.1)
plot(PCA.w.phylo.LDR, phylo = TRUE, cex = 0.5, xaxp=c(-2.5,1,7))

PCA.w.phylo.LDR <- gm.prcomp(data.LDR.full, phy = phy.avg)
plot(PCA.w.phylo.LDR, phylo = TRUE, cex = 0.5, xaxp=c(-2.5,1,7))


PCA.w.phylo.BDR <- gm.prcomp(data.BDR.full, phy = phy.PCA.1)
plot(PCA.w.phylo.BDR, phylo = TRUE, cex = 0.5, xaxp=c(-2.5,1,7))
