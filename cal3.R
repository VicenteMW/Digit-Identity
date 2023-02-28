#Loading packages

library(paleotree)

#Loading and formatting tree data, and interval data for taxa under study.
tree.data <- read.tree("phy.tree.1.txt")
mat <- scan('taxa.int.txt')
ints.taxa <- matrix(mat, ncol = 2, byrow = TRUE)
mat <- scan('intervals.txt')
ints <- matrix(mat, ncol = 2, byrow = TRUE)
vec.rownames <- c(tree.data$tip.label)
rownames(ints.taxa) <- vec.rownames
FAD.LAD <- list()
FAD.LAD[[1]] <- ints
FAD.LAD[[2]] <- ints.taxa


#Loading and formatting occurrence data from PBDB, from mesozoic saurischians
url <- 'https://paleobiodb.org/data1.2/occs/list.csv?base_name=Saurischia&taxon_reso=species&interval=Triassic,Cretaceous&show=ident,phylo'
saurischia <- read.csv(file = url)
names(saurischia)[15:16]<-c("early_age", "late_age")
sauris_sorted<-taxonSortPBDBocc(saurischia, 'species')
sauris_timelist<-occData2timeList(sauris_sorted)

#Calculation of interval mean length
interval<--apply(sauris_timelist[[1]],1,diff)
m_interval<-mean(interval, na.rm = T)


#Sampling, branching and extinction rate estimation
sauris_likFun<-make_durationFreqDisc(sauris_timelist)
spRes <- optim(
  parInit(sauris_likFun),
  sauris_likFun,
  lower = parLower(sauris_likFun),
  upper = parUpper(sauris_likFun),
  method = "L-BFGS-B",
  control = list(maxit = 1000000))
sProb <- spRes[[1]][2]
sRate <- sProb2sRate(sProb,int.length = m_interval)
divRate <- spRes[[1]][1]/m_interval
#extinction and branching rate = divRate, sampling rate = sRate


#Cal3 is used to generate 100 trees, all with the same topology but differing in branch lengths.
phylogenies <- bin_cal3TimePaleoPhy(
  tree.data,
  FAD.LAD,
  brRate = divRate,
  extRate = divRate,
  sampRate = sRate,
  anc.wt = 0,
  ntrees = 100,
  plot = TRUE)
