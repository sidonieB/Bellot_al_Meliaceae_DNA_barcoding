library(ape)
library(phytools)
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(VennDiagram)


##################################
# Script used to analyse gene trees and classify the genes depending on recovery across samples and informativeness
##################################

setwd("rooted")    # directory with the ROOTED gene trees
input_sp <- "list.txt"       # list of genera, species names and number of individuals per species
pattern_trees <- "Rooted.tre"                         # whatever pattern is common to all your ROOTED tree files AND ONLY TO THEM

sp_num <- read.table(input_sp, header=T, sep="\t", stringsAsFactors = T)

sp_num_sub <- split.data.frame(sp_num, sp_num$genus)

trees <- list.files(".", pattern = pattern_trees, full.names = T) 

#the loop below can take forever if many combinations exist among which to find the largest


for (G in 1:length(sp_num_sub)){
  
  genus <- as.character(sp_num_sub[[G]][1,1])
  
  res2 <- as.data.frame(matrix(data=NA, nrow=0, ncol=8))
  names(res2) <- c("genus", "gene", "species", "intraspecies_recovery_number", "intraspecies_recovery_percent", "individuals_in_largest_clade_number", "individuals_in_largest_clade_percent", "largest_clade_BP")
  
 for (z in 1:length(trees)){
   
   gene <- gsub("_r_CI85_T_g.fasta.contree_PxrrRooted.tre", "", trees[z]) # modify this depending on your tree names so that only keep simple gene names
   gene <- gsub("./", "", gene)
   
   t <- read.tree(trees[[z]])
   
   res <- as.data.frame(matrix(data=NA, nrow=(length(sp_num_sub[[G]][,1]))+2, ncol=8))
   names(res) <- c("genus", "gene", "species", "intraspecies_recovery_number", "intraspecies_recovery_percent", "individuals_in_largest_clade_number", "individuals_in_largest_clade_percent", "largest_clade_BP")
   
   
  for (s in 1:length(sp_num_sub[[G]][,1])) {
    
    sp <-  as.character(sp_num_sub[[G]][s,2])
    indiv_num <- sp_num_sub[[G]][s,3]
    indiv_pres_num <- length(grep(sp, t$tip.label))
    indiv_pres_per <- (indiv_pres_num*100)/indiv_num
    
    if (indiv_pres_num > 1) {
    
    present_indiv <- t$tip.label[grep(sp, t$tip.label)]
   
    if (is.monophyletic(t, present_indiv)) {
 
      length_final_indiv <- length(present_indiv)
      length_final_indiv_per <- 100
      largest_clade_node <- getMRCA(t,present_indiv)
      node_index <- largest_clade_node - length(t$tip.label)
      largest_clade_BP <- t$node.label[node_index]   

    
    } else {
      
      length_final_indiv <- 0
      indiv_combs <- do.call("c", lapply(seq_along(present_indiv), function(i) combn(present_indiv, i, FUN = list)))
      for (c in 1:length(indiv_combs)) {
        if (is.monophyletic(t, indiv_combs[[c]])) {
          if (length(indiv_combs[[c]]) > length_final_indiv) {
            final_indiv <- indiv_combs[[c]]
            length_final_indiv <- length(indiv_combs[[c]])
          }
        }
      }
      
      length_final_indiv_per <- (length_final_indiv * 100)/ length(present_indiv)
      
      if (length_final_indiv > 1) {
        largest_clade_node <- getMRCA(t,final_indiv)
        node_index <- largest_clade_node - length(t$tip.label)
        largest_clade_BP <- t$node.label[node_index]
      } else {
        largest_clade_BP <- NA
      }     
      
    }
    
    } else {
      length_final_indiv <- NA
      length_final_indiv_per <- NA
      largest_clade_BP <- NA
    }
    
    res[s,] <- c(genus, gene, sp, indiv_pres_num, indiv_pres_per, length_final_indiv, length_final_indiv_per, largest_clade_BP)
    
    
  }
   
  
  res[(length(sp_num_sub[[G]][,1]))+1,1:3] <- c(genus, gene, "Average")
  res[(length(sp_num_sub[[G]][,1]))+2,1:3] <- c(genus, gene, "Median")
  for (x in 4:8) {
    res[(length(sp_num_sub[[G]][,1]))+1,x] <- mean(as.numeric(res[1:(length(sp_num_sub[[G]][,1])),x]), na.rm = T)
    res[(length(sp_num_sub[[G]][,1]))+2,x] <- median(as.numeric(res[1:(length(sp_num_sub[[G]][,1])),x]), na.rm = T)
  }

res2 <- rbind(res2, res)
  
 }

# export table for this genus
  write.table(res2, paste0(genus, "_barcoding_potential_assessment.txt"))

}



### Plots


Kg <- read.table("Khaya_barcoding_potential_assessment.txt", header=T, sep=" ")
Eg <- read.table("Entandrophragma_barcoding_potential_assessment.txt", header=T, sep=" ")
Lg <- read.table("Lovoa_barcoding_potential_assessment.txt", header=T, sep=" ")
Sg <- read.table("Swietenia_barcoding_potential_assessment.txt", header=T, sep=" ")

head(Kg)

Kg2 <- subset(Kg, Kg$species != "Average")
Kg3 <- subset(Kg2, Kg2$species != "Median")
Eg2 <- subset(Eg, Eg$species != "Average")
Eg3 <- subset(Eg2, Eg2$species != "Median")
Sg2 <- subset(Sg, Sg$species != "Average")
Sg3 <- subset(Sg2, Sg2$species != "Median")
Lg2 <- subset(Lg, Lg$species != "Average")
Lg3 <- subset(Lg2, Lg2$species != "Median")

All_ass <- rbind(Kg3, Eg3, Lg3, Sg3)

genes <- unique(All_ass$gene)
genes_tab <- as.data.frame(genes)
genes_tab$intraspecies_recovery_percent_Med <- rep(NA, length(genes_tab$genes))
genes_tab$intraspecies_recovery_percent_Av <- rep(NA, length(genes_tab$genes))
genes_tab$num_monophyletic_species <- rep(NA, length(genes_tab$genes))
genes_tab$num_monophyletic_species_relaxed80 <- rep(NA, length(genes_tab$genes))
genes_tab$species_monophylousness_Av <- rep(NA, length(genes_tab$genes))
genes_tab$species_monophylousness_Med <- rep(NA, length(genes_tab$genes))

for (g in 1:length(genes_tab$genes)) {
  
  gt <- subset(All_ass, All_ass$gene == genes_tab$genes[g])
  
  genes_tab$intraspecies_recovery_percent_Med[g] <- median(gt$intraspecies_recovery_percent, na.rm=T)
  genes_tab$intraspecies_recovery_percent_Av[g] <- mean(gt$intraspecies_recovery_percent, na.rm=T)
  genes_tab$num_monophyletic_species[g] <- length(which(gt$individuals_in_largest_clade_percent > 100))
  genes_tab$num_monophyletic_species_relaxed80[g] <- length(which(gt$individuals_in_largest_clade_percent > 80))
  genes_tab$species_monophylousness_Med[g] <- median(gt$individuals_in_largest_clade_percent, na.rm=T)
  genes_tab$species_monophylousness_Av[g] <- mean(gt$individuals_in_largest_clade_percent, na.rm=T)
}

# do some gene name cleaning for easier display, modify depending on gene names
genes_tab$region <- gsub("_.*", "", genes_tab$genes)
genes_tab$region2 <- gsub("g.*", "", genes_tab$region)
genes_tab$genes2 <- gsub("ANGIO353", "", genes_tab$genes)
genes_tab$genes3 <- gsub("_ING", "", genes_tab$genes2)
genes_tab$genes4 <- gsub("_IPIM", "", genes_tab$genes3)
genes_tab$genes5 <- gsub("PLANCS_", "", genes_tab$genes4)
genes_tab$genes6 <- gsub("PLACDS_", "", genes_tab$genes5)
genes_tab$genes7 <- gsub("PLARRN_", "", genes_tab$genes6)

# Plot all genes coloring by nuclear Angiosperms353/plastid/ITS
p1 <- ggplot(genes_tab, aes(x=species_monophylousness_Av, y=intraspecies_recovery_percent_Av, color=region2)) +
  geom_point()+
  theme_classic()+
  theme(legend.position = "bottom")

# plot best genes with names
p2 <- ggplot(genes_tab, aes(x=species_monophylousness_Av, y=intraspecies_recovery_percent_Av, color=region2)) +
  geom_point() + xlim(70,100) + ylim(70,100) +
  geom_text_repel(label=genes_tab$genes7, size=2, max.overlaps = 100)+
  theme_classic()+
  theme(legend.position = "bottom")
  








