
# 1. PACKAGES

library(ape)
library(phylobase)
library(data.table)
library(ggplot2)

# 2. INPUT DATA (CSV)

species_data <- read.csv("species_data.csv", sep = ",", stringsAsFactors = FALSE)

# Keep only the first two columns
species_data <- species_data[, 1:2]

# Rename colums
colnames(species_data) <- c("species","RL")

str(species_data)
# Remove leading and trailing whitespace
species_data$species <- trimws(species_data$species)

# 3. LOAD BEAST TREES

trees <- read.nexus("cantha.trees")

# Check the number of trees
length(trees)

# 4. REMOVE BURN-IN

burnin <- 0.1

start <- floor(length(trees)*burnin) + 1
trees <- trees[start:length(trees)]

# 5. SAMPLE TREES

set.seed(123)
trees_sample <- sample(trees, 1000)

# Clean tree tip labels
trees_sample <- lapply(trees_sample, function(tr){
  tr$tip.label <- trimws(tr$tip.label)
  return(tr)
})

# 6. CHECK SPECIES MATCHING

tree_test <- trees_sample[[1]]

length(intersect(tree_test$tip.label, species_data$species))

# 7. GENERATE THE GE2 DISTRIBUTION

pext.vals <- data.frame(
  rl.cat = rev(c("CR","EN","VU","NT","LC")),
  pext = rev(c(0.97, 0.97/2, 0.97/4, 0.97/8, 0.97/16))
)

GE.2.calc <- function(pext){
  iucn <- sample(1:5, size=1000000, replace=TRUE)
  data <- data.frame(species=c(1:1000000), pext=pext$pext[iucn])
  data <- data[order(data$pext),]
  data$rank <- seq_len(nrow(data))
  rank <- c(0, with(data, tapply(rank, pext, median)))
  pext.tap <- c(0, pext$pext)
  rank.sq <- rank^2; rank.cub <- rank^3; rank.qu <- rank^4
  model <- lm(pext.tap ~ rank + rank.sq + rank.cub + rank.qu)
  data$rank.sq <- data$rank^2
  data$rank.cub <- data$rank^3
  data$rank.qu <- data$rank^4
  data$rank.pext <- predict(model, data)
  data$rank.pext[data$rank.pext <= 0] <- 0.0001
  data$rank.pext[data$rank.pext >= 1] <- 0.9999
  
  pext.LC <- data.frame(RL.cat="LC", pext=data$rank.pext[data$pext==pext.tap[2]])
  pext.NT <- data.frame(RL.cat="NT", pext=data$rank.pext[data$pext==pext.tap[3]])
  pext.VU <- data.frame(RL.cat="VU", pext=data$rank.pext[data$pext==pext.tap[4]])
  pext.EN <- data.frame(RL.cat="EN", pext=data$rank.pext[data$pext==pext.tap[5]])
  pext.CR <- data.frame(RL.cat="CR", pext=data$rank.pext[data$pext==pext.tap[6]])
  
  pext.obj <- rbind(pext.CR,pext.EN,pext.VU,pext.NT,pext.LC)
  
  pext.obj.sample <- NULL
  
  for(i in unique(pext.obj$RL.cat)){
    a <- pext.obj$pext[pext.obj$RL.cat == i]
    
    if(median(a) < pext$pext[pext$rl.cat == i]){
      while(median(a) < pext$pext[pext$rl.cat == i]){
        a <- a[-sample(1:length(a),50,prob=rev(a))]
      }
    } else {
      while(median(a) > pext$pext[pext$rl.cat == i]){
        a <- a[-sample(1:length(a),50,prob=a)]
      }
    }
    
    pext.obj.sample <- rbind(pext.obj.sample, data.frame(RL.cat=i, pext=a))
  }
  
  return(pext.obj.sample)
}

pext_dist <- GE.2.calc(pext.vals)


# 8. EDGE2 AUXILIARY FUNCTIONS

into_order <- function(tree, pext){
  new_pext <- pext[match(tree$tip.label, pext$species),]
  return(new_pext)
}

reorder_tree <- function(tree, ordering){
  tree@edge.length <- tree@edge.length[ordering]
  tree@edge <- tree@edge[ordering,]
  return(tree)
}

# 9. EDGE2 FUNCTION

EDGE2_mod <- function(tree, pext){
  
  require(phylobase)
  require(data.table)
  
  names(pext) <- c("species","pext")
  
  N_species <- length(tree$tip.label)
  N_nodes <- tree$Nnode
  N_tot <- N_species + N_nodes
  
  # Ensure extinction probabilities are given in same order as tree$tip.label.
  if (!identical(tree$tip.label, pext$species)){
    pext <- into_order(tree, pext)
  }
  
  if(!class(tree) == "phylo"){
    tree <- as(tree, "phylo")
  }
  
  tree_dat <- data.frame(Species = as.character(tree$tip.label),
                         TBL = NA, 
                         pext = pext$pext, ED = NA, EDGE = NA)
  ePD.dat <- data.frame(PD = sum(tree$edge.length),ePDloss = NA)
  
  tree <- as(tree, "phylo4")
  root <- rootNode(tree)
  nodes <- c(root, descendants(tree, root, "all"))
  
  # Reorder tree components more instinctively, such that nodes are easier to find
  ord <- order(nodes)
  tree <- reorder_tree(tree, ord)
  nodes <- nodes[ord]
  
  tree_dat$TBL <- tree@edge.length[1:N_species]
  
  node_data <- data.frame(Node = 1:N_tot, Pext = rep(1, N_tot), Edge_Sum = NA)
  node_data[1:N_species, 2] <- pext[,2]
  
  # Assign the product of its descendant tips to each node
  for (i in c(1:length(tree@label), N_tot:(root+1))){       # for each node, beginning with tips
    anc <- tree@edge[i,1]                                   # find ancestor of node
    node_data[anc, 2] <- node_data[anc, 2]*node_data[i,2]   # muliply ancestor value by node "pext" 
  }
  
  # Multiply each edge.length by each pext caluclated above
  for(i in 1:length(nodes)){
    tree@edge.length[i] <- tree@edge.length[i]*node_data[i,2]
  }
  save(tree, file ="tree.rda")
  
  if (is.na(tree@edge.length[root])){
    tree@edge.length[root] <- 0
  }
  node_data$Edge_Sum[root] <- tree@edge.length[root]
  
  # For each internal node, summate ancesteral edgelengths
  for (i in (root+1):N_tot){
    ans <- tree@edge[i,1]
    node_data$Edge_Sum[i] <- node_data$Edge_Sum[ans] + tree@edge.length[i]
  }
  
  # For each tip, summate ancesteral edgelengths to find EDGE2 score
  for (i in 1:N_species){
    ans <- tree@edge[i,1]
    tree_dat$EDGE[i] <- node_data$Edge_Sum[ans] + tree@edge.length[i]
  }  
  
  tree_dat$ED <- tree_dat$EDGE / tree_dat$pext
  # Reorder tree
  tree <- reorder_tree(tree, order(ord))
  
  tree <- as(tree, "phylo")
  ePD.dat$ePDloss <- sum(tree$edge.length)
  edge.res <- list(tree_dat,tree,ePD.dat)
  return(edge.res)
}

# 10. MAIN ANALYSIS LOOP

n_iter <- 1000
results <- list()

for(i in 1:n_iter){
  
  cat("Iteracao:", i, "\n")
  
  species_iter <- species_data
  tree_i <- trees_sample[[i]]
  
  # Match species to tree tip labels
  species_iter <- species_iter[match(tree_i$tip.label, species_iter$species),]
  
  # Check for matching errors
  if(any(is.na(species_iter$species))){
    print(paste("Matching error in iteration", i))
    next
  }
  
  # Identify Data Deficient (DD) species
  dd_index <- which(species_iter$RL == "DD")
  
  # Generate extinction probabilities (pext)
  species_iter$pext <- sapply(species_iter$RL, function(cat){
    if(cat != "DD"){
      sample(pext_dist$pext[pext_dist$RL.cat == cat], 1)
    } else {
      NA
    }
  })
  
  # Assign random extinction probabilities to DD species
  if(length(dd_index) > 0){
    species_iter$pext[dd_index] <- sample(
      pext_dist$pext,
      size = length(dd_index),
      replace = TRUE
    )
  }
  
  # Prepare EDGE2 input
  pext_input <- data.frame(
    species = species_iter$species,
    pext = species_iter$pext
  )
  
  # Run EDGE2
  res <- tryCatch({
    EDGE2_mod(tree_i, pext_input)
  }, error = function(e){
    print(paste("Error in iteration", i))
    return(NULL)
  })
  
  results[[i]] <- res
}

length(results)

# 11. CREATE THE FINAL EDGE MATRIX

results <- results[!sapply(results, is.null)]

EDGE_mat <- do.call(cbind, lapply(results, function(x) x[[1]]$EDGE))

rownames(EDGE_mat) <- results[[1]][[1]]$Species

# 12. CREATE THE PEXT (GE2) MATRIX

# Combine extinction probability values from all simulations
pext_mat <- do.call(cbind, lapply(results, function(x) x[[1]]$pext))

# 13. CALCULATE EVOLUTIONARY DISTINCTIVENESS (ED)

# ED is obtained by dividing EDGE by extinction probability (pext)
ED_mat <- EDGE_mat / pext_mat

# 14. SUMMARIZE RESULTS (MEDIANS)

# Calculate the median across simulations for each species
EDGE_mean   <- rowMeans(EDGE_mat, na.rm = TRUE)
EDGE_median <- apply(EDGE_mat, 1, median, na.rm = TRUE)
ED_median   <- apply(ED_mat, 1, median, na.rm = TRUE)
GE2_median  <- apply(pext_mat, 1, median, na.rm = TRUE)

# 15. CREATE THE FINAL RESULTS TABLE

# Combine all metrics into a single table
final_full <- data.frame(
  Species = rownames(EDGE_mat),
  EDGE = EDGE_median,
  ED = ED_median,
  GE2 = GE2_median
)

# 16. SORT SPECIES BY EDGE PRIORITY

# Sort species from highest to lowest EDGE score
final_full <- final_full[order(final_full$EDGE, decreasing = TRUE),]

# 17. CREATE RANKED RESULTS

EDGE_rank <- sort(EDGE_median, decreasing = TRUE)

final_results <- data.frame(
  Species = names(EDGE_rank),
  EDGE_median = EDGE_rank
)

head(final_results)

# 18. CALCULATE THE GLOBAL MEDIAN

# Global median EDGE score

global_median <- median(EDGE_median, na.rm = TRUE)

# 19. CALCULATE THE FREQUENCY ABOVE THE GLOBAL MEDIAN

# Calculate the proportion of simulations in which each species scored above the global median EDGE value
freq_above <- apply(EDGE_mat, 1, function(x){
  mean(x > global_median, na.rm = TRUE)
})

# 20. ADD THE FREQUENCY TO THE FINAL TABLE

final_full$Frequency <- freq_above

# 21. CREATE PRIORITY LISTS

# Main EDGE list (threatened species)
EDGE_list <- final_full[
  final_full$Frequency >= 0.95 &
    species_data$RL[match(final_full$Species, species_data$species)] %in% c("VU","EN","CR"),
]

# Research list (Data Deficient species)
Research_list <- final_full[
  final_full$Frequency >= 0.95 &
    species_data$RL[match(final_full$Species, species_data$species)] == "DD",
]

# Watch list (Least Concern and Near Threatened species)
Watch_list <- final_full[
  final_full$Frequency >= 0.95 &
    species_data$RL[match(final_full$Species, species_data$species)] %in% c("LC","NT"),
]

# Borderline list
Borderline_list <- final_full[
  final_full$Frequency >= 0.80 &
    final_full$Frequency < 0.95 &
    species_data$RL[match(final_full$Species, species_data$species)] %in% c("VU","EN","CR"),
]

# 22. SAVE RESULTS

# Save the complete results table

write.csv(final_full, "EDGE_full_com_frequency.csv", row.names = FALSE)

# Save priority lists

write.csv(EDGE_list, "EDGE_list.csv", row.names = FALSE)
write.csv(Research_list, "Research_list.csv", row.names = FALSE)
write.csv(Watch_list, "Watch_list.csv", row.names = FALSE)
write.csv(Borderline_list, "Borderline_list.csv", row.names = FALSE)

# ave ranked results

write.csv(final_results, "EDGE_final_results.csv", row.names = FALSE)

# Save the complete R workspace

save.image("EDGE_completo.Rdata")

# PLOTS

# 23. LOAD PACKAGES

library(ggplot2)
library(ape)

# 24. LOAD DATA

dados <- final_full
tree <- read.tree("cantha_median_clean.tree")

# 25. CLEAN SPECIES NAMES

dados$Species <- trimws(dados$Species)
species_data$species <- trimws(species_data$species)

tree$tip.label <- trimws(tree$tip.label)
tree$tip.label <- gsub("'", "", tree$tip.label)

# 26. CHECK DATA (SHOULD BE 108)

cat("Tree tips:", length(tree$tip.label), "\n")
cat("Matching species:", length(intersect(tree$tip.label, dados$Species)), "\n")

# 27. REMOVE INVALID TIPS (IF NECESSARY)

tips_invalidos <- setdiff(tree$tip.label, dados$Species)

if(length(tips_invalidos) > 0){
  tree <- drop.tip(tree, tips_invalidos)
}

# 28. MATCH DATA TO THE TREE

dados_plot <- dados[match(tree$tip.label, dados$Species),]

# 29. DD IUCN RED LIST CATEGORIES

dados_plot$RL <- species_data$RL[
  match(dados_plot$Species, species_data$species)
]

# 30. DEFINE PHYLOGENETIC ORDER

dados_plot$Species <- factor(
  dados_plot$Species,
  levels = rev(tree$tip.label)  # mantém ordem da árvore
)

# 31. DEFINE IUCN CATEGORY ORDER

dados_plot$RL <- factor(
  dados_plot$RL,
  levels = c("DD","LC","NT","VU","EN","CR") 
)

# 32. DEFINE IUCN COLORS

cores_iucn <- c(
  "DD" = "#D1D1C6",
  "LC" = "#60C659",
  "NT" = "#CCE226",
  "VU" = "#F9E814",
  "EN" = "#FC7F3F",
  "CR" = "#D81E05"
)


# 33. CREATE THE FINAL PLOT

p <- ggplot(dados_plot, aes(x = EDGE, y = Species, fill = RL)) +
  geom_col(width = 0.7) +
  
  scale_fill_manual(
    values = cores_iucn,
    drop = FALSE,  #  Keep all IUCN categories in the legend, even if absent from the data
    name = "IUCN",
    labels = c(
      "DD" = "Data Deficient (DD)",
      "LC" = "Least Concern (LC)",
      "NT" = "Near Threatened (NT)",
      "VU" = "Vulnerable (VU)",
      "EN" = "Endangered (EN)",
      "CR" = "Critically Endangered (CR)"
    )
  ) +
  
  labs(x = "EDGE2 score", y = "") +
  
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 6),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

# 34. DISPLAY THE PLOT

print(p)

# 35. SAVE THE PLOT

ggsave("EDGE_plot_final.png", p, width = 8, height = 12, dpi = 300)
ggsave("EDGE_plot_final.pdf", p, width = 8, height = 12)
