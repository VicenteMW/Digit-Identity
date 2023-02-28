library (phytools)

# Importing data for LDR, full hand and last phalanx missing, and transforming it into a dataframe
mat <- scan('data.LDR.final.txt')
data.LDR.full <- matrix(mat, ncol = 20, byrow = TRUE)
rownames(data.LDR.full) <- c(phylogenies[[1]]$tip.label)
colnames(data.LDR.full) <- c(1:20)
data.LDR.full.df <- as.data.frame(data.LDR.full)

# Then, we generate an empty list
data.LDR.full.list <- list()

# We use the following loop within a loop to perform a fastAnc analysis on the 100 phylogenetic trees and 19 continuous traits
for (j in 1:length(phylogenies)) {
  dummy.list <-list()
  for (i in 1:length(data.LDR.full.df)) {
    fanc <- fastAnc(phylogenies[[j]], x = setNames(data.LDR.full.df[[i]], rownames(data.LDR.full)), CI = TRUE)
    dummy.list[[i]] <- fanc                  
  }
  data.LDR.full.list[[j]] <- dummy.list
}

# Next, we do the same for BDR characters, full hand and last phalanx missing.
# Importing data for BDR and transforming it into a dataframe
mat <- scan('data.BDR.final.txt')
data.BDR.full <- matrix(mat, ncol = 20, byrow = TRUE)
rownames(data.BDR.full) <- c(phylogenies[[1]]$tip.label)
colnames(data.BDR.full) <- c(1:20)
data.BDR.full.df <- as.data.frame(data.BDR.full)

# Then, we generate an empty list
data.BDR.full.list <- list()

# We use the following loop within a loop to perform a fastAnc analysis on the 100 phylogenetic trees and 19 continuous traits
for (j in 1:length(phylogenies)) {
  dummy.list <-list()
  for (i in 1:length(data.BDR.full.df)) {
    fanc <- fastAnc(phylogenies[[j]], x = setNames(data.BDR.full.df[[i]], rownames(data.BDR.full)), CI = TRUE)
    dummy.list[[i]] <- fanc                  
  }
  data.BDR.full.list[[j]] <- dummy.list
}

