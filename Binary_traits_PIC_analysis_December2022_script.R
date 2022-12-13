#######load packages/working directory

setwd("~/Swordtail Dropbox/Schumer_lab_resources/Project_files/Body_size_evolution/Data/PIC_script")

# R packages
library(phytools)
library(geiger)
library(ape)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggtree)
library(ggtreeExtra)
library(treeio)
library(grid)
library(gridExtra)
##devtools::install_github("lamho86/phylolm")
library(phylolm)

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
#####Running phylolm
########################

#convert sword index to binary variable
sword_index_threshold<- 0.01

m_morpho$sword_binary<-rep(0,length(m_morpho$std_length))
tmpindex<-m_morpho$swd_length/m_morpho$std_length
m_morpho$sword_binary[tmpindex < sword_index_threshold] <- 0
m_morpho$sword_binary[tmpindex > sword_index_threshold] <- 1

#use Garland & Ives pglm 
#set traits for analysis
v_m_morpho_trait1 <- with(m_morpho, setNames(m_morpho$sword_binary, m_morpho$species))  
v_m_morpho_trait2 <- with(m_morpho, setNames(m_morpho$std_length, m_morpho$species))  

dat=data.frame(y=v_m_morpho_trait1,x=v_m_morpho_trait2) 
summary(phyloglm(y~x,dat,pruned_swdtree,method="logistic_IG10") )

##############
##Also run analysis drop Pjonesii
##############

species = subset(species, species != "Xbirref" &
                   species != "Priapella" & species != "Psjonesii")
pruned_swdtree2 = drop.tip(swdtree,swdtree$tip.label[-match(species, swdtree$tip.label)])

# create dataframe for males
morpho2 = morpho %>%
  filter(sex == "M") %>%
  filter(!species %in% c("Xkallmani")) %>%
  filter(!species %in% c("Psjonesii"))
  
 #analysis for gonopodial traits
 
v_m_morpho_trait1 <- with(morpho2, setNames(morpho2$spine_angle_ray3, morpho2$species)) #significant
v_m_morpho_trait2 <- with(morpho2, setNames(morpho2$std_length_literature, morpho2$species))  

dat=data.frame(y=v_m_morpho_trait1,x=v_m_morpho_trait2) 
summary(phyloglm(y~x,dat,pruned_swdtree2,method="logistic_IG10") )
