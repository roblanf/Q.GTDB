library(sem)
library(ggfortify)
library(devtools)
library(ggbiplot)
library(ggrepel)
library(tidyverse)
library(factoextra)

setwd("~/Documents/GitHub/Q.GTDB/Q.GTDB/PCA/")

readQ <- function(filename) {
  # read lower-diagonal matrix
  Q = readMoments(filename, diag=F)
  pi = Q[nrow(Q), 1:(ncol(Q)-1)]
  Q = Q[1:(nrow(Q)-1), 1:(ncol(Q)-1)]
  # make Q symmetric
  Q = (Q + t(Q))
  diag(Q) <- 0
  # normalise the matrix
  Q=(Q/sum(Q))*100.0
  return(rbind(Q, pi))
}

nuclear = c("WAG", "Dayhoff","JTT", "LG", "VT", "PMB", "Blosum62", "Q.pfam")
mitochondrial = c("mtREV", "mtMAM", "mtART", "mtZOA", "mtMet" , "mtVer" , "mtInv")
chloroplast = c("cpREV")
viral = c("HIVb", "HIVw", "FLU", "rtREV")
bacterial = c("Q.bacteria_class_1", "Q.bacteria_order_1", "Q.bacteria_phylum_1", "Q.GTDB_sub_1k_100", "Q.GTDB_sub_5k_100", "Q.GTDB_sub250_5k")
taxon_specific = c("Q.plant", "Q.bird", "Q.mammal", "Q.insect", "Q.yeast")

filenames = c(nuclear, mitochondrial, chloroplast, viral, bacterial, taxon_specific)
regions = c(rep("nuclear", length(nuclear)), rep("mitochondrial", length(mitochondrial)),
            rep('chloroplast', length(chloroplast)), rep('viral', length(viral)),
            rep('bacterial',length(bacterial)), rep('taxon-specific',length(taxon_specific)))

allQ = NULL
allF = NULL

for (f in filenames) {
  print(f)
  Q_pi = readQ(f)
  Q = Q_pi[1:(nrow(Q_pi)-1),]
  pi = Q_pi[nrow(Q_pi),]
  
  n <- nrow(Q)
  
  allQ = rbind(allQ, Q[lower.tri(Q)])
  allF = rbind(allF, pi)
}

rownames(allQ) <- filenames
rownames(allF) <- filenames

pcaQ <- prcomp(allQ, scale = TRUE)
pca_summary <- summary(pcaQ)
pc1_var <- pca_summary$importance[2, 1]  # Proportion of variance explained by PC1
pc2_var <- pca_summary$importance[2, 2]  # Proportion of variance explained by PC2
pc1_label <- sprintf("PC1 (%.2f%% Variance Explained)", pc1_var * 100)
pc2_label <- sprintf("PC2 (%.2f%% Variance Explained)", pc2_var * 100)

ggplot(pcaQ$x, aes(x = PC1, y = PC2, colour = as.factor(regions))) + 
  geom_point() +
  geom_text_repel(aes(label = rownames(pcaQ$x))) +
  xlab(pc1_label) +
  ylab(pc2_label)



pcaF <- prcomp(allF, scale = TRUE)
pca_summary <- summary(pcaF)
pc1_var <- pca_summary$importance[2, 1]  # Proportion of variance explained by PC1
pc2_var <- pca_summary$importance[2, 2]  # Proportion of variance explained by PC2
pc1_label <- sprintf("PC1 (%.2f%% Variance Explained)", pc1_var * 100)
pc2_label <- sprintf("PC2 (%.2f%% Variance Explained)", pc2_var * 100)

ggplot(pcaF$x, aes(x = PC1, y = PC2, colour = as.factor(regions))) + 
  geom_point() +
  geom_text_repel(aes(label = rownames(pcaF$x))) +
  xlab(pc1_label) +
  ylab(pc2_label)
