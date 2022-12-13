#######load packages/working directory

setwd("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Body_size_evolution/Data/PIC_script")

library(tidyverse)
library(ggpubr)
library(ape)
library(ggtree)
library(ggtreeExtra)
library(treeio)
library(grid)
library(gridExtra)
library(picante)

#####load morphological data

# read in data
morpho = read_csv("20221206_phy_info_dimorphism_add_gonopodia.csv")

# add sword index
morpho = morpho %>%
  mutate(swd_index = swd_length/std_length)

# create dataframe for males
m_morpho = morpho %>%
  filter(sex == "M") %>%
  filter(!species %in% c("Xkallmani"))

#####load modified phylogeny

swdtree = read.nexus("tobir3concat_pruned_v4_revised_August2022_branchlengths.tre")
swdtree = root(swdtree, outgroup = "Psjonesii")
species = swdtree$tip.label

#remove duplicate ids/missing species
species = subset(species, species != "Xbirref" &
                   species != "Priapella")
pruned_swdtree = drop.tip(swdtree,swdtree$tip.label[-match(species, swdtree$tip.label)])

########################
###test assumptions for traits of interest re Garland et al 1992
###following workflow here: http://phylo.wikidot.com/phylogenetically-independent-contrasts
########################

diagnostics(m_morpho$lower_sword_edge_width,pruned_swdtree) 
diagnostics(m_morpho$swd_index,pruned_swdtree)

#transformation for branch lengths to improve diagnostics
pruned_swdtree_transform<-pruned_swdtree
pruned_swdtree_transform$edge.length<-pruned_swdtree$edge.length^2

#recheck diagnostics
diagnostics(m_morpho$lower_sword_edge_width,pruned_swdtree_transform)

######################
####Run PIC analysis for traits of interest

###set traits for analysis
v_m_morpho_trait1 <- with(m_morpho, setNames(m_morpho$lower_sword_edge_width, m_morpho$species))
v_m_morpho_trait2 <- with(m_morpho, setNames(m_morpho$std_length, m_morpho$species))

tree_trait1<-pruned_swdtree_transform
tree_trait2<-pruned_swdtree

#run PIC
picres1 = pic(v_m_morpho_trait1, tree_trait1)
picres2 = pic(v_m_morpho_trait2, tree_trait2)
summary(lm(picres1 ~ 0 + picres2))

####plot results
rl<-lm(picres1 ~ 0 + picres2)
plot(picres2,picres1,pch=20,cex=3,cex.axis=1.2,ylab="PIC lower sword edge",xlab="PIC body size",col=rgb(50/255,96/255,168/255,alpha=0.4),cex.lab=1.2)
abline(rl,col="darkgray",lty=2,lwd=2)

