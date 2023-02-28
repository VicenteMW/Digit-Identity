library (phytools)

# Importing data for LDR, met hand and last phalanx missing, and transforming it into a dataframe
mat <- scan('data.LDR.met.txt')
data.LDR.met <- matrix(mat, ncol = 5, byrow = TRUE)
rownames(data.LDR.met) <- c(phylogenies[[1]]$tip.label)
colnames(data.LDR.met) <- c(1:5)
data.LDR.met.df <- as.data.frame(data.LDR.met)

# Then, we generate an empty list
data.LDR.met.list <- list()

# We use the following loop within a loop to perform a fastAnc analysis on the 100 phylogenetic trees and 19 continuous traits
for (j in 1:length(phylogenies)) {
  dummy.list <-list()
  for (i in 1:length(data.LDR.met.df)) {
    fanc <- fastAnc(phylogenies[[j]], x = setNames(data.LDR.met.df[[i]], rownames(data.LDR.met)), CI = TRUE)
    dummy.list[[i]] <- fanc                  
  }
  data.LDR.met.list[[j]] <- dummy.list
}

# Next, we do the same for BDR characters, met hand and last phalanx missing.
# Importing data for BDR and transforming it into a dataframe
mat <- scan('data.BDR.met.txt')
data.BDR.met <- matrix(mat, ncol = 5, byrow = TRUE)
rownames(data.BDR.met) <- c(phylogenies[[1]]$tip.label)
colnames(data.BDR.met) <- c(1:5)
data.BDR.met.df <- as.data.frame(data.BDR.met)

# Then, we generate an empty list
data.BDR.met.list <- list()

# We use the following loop within a loop to perform a fastAnc analysis on the 100 phylogenetic trees and 19 continuous traits
for (j in 1:length(phylogenies)) {
  dummy.list <-list()
  for (i in 1:length(data.BDR.met.df)) {
    fanc <- fastAnc(phylogenies[[j]], x = setNames(data.BDR.met.df[[i]], rownames(data.BDR.met)), CI = TRUE)
    dummy.list[[i]] <- fanc                  
  }
  data.BDR.met.list[[j]] <- dummy.list
}

