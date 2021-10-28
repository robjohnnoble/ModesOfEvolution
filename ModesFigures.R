library(ape)
library(readr)
library(dplyr)
library(ggmuller)
library(dils)
library(Rgraphviz)
library(ggplot2)
library(gridExtra)
library(demonanalysis)
library(ggrepel)
library(data.table)

move_down_mod <- function(edges, parent) {
  edges <- filter(edges, Parent != 0 | Identity != 0)
  if (!(parent %in% edges$Identity) & !(parent %in% edges$Parent)) 
    stop("Invalid parent.")
  daughters <- filter(edges, Parent == parent)$Identity
  if (length(daughters) == 0) 
    return(parent)
  if (is.factor(daughters)) 
    daughters <- levels(daughters)[daughters]
  return(sort(daughters))
}

get_ITH_from_tree <- function(tree) {
  node <- 0
  clonal <- 0
  for(i in 1:100) {
    node <- move_down_mod(tree, node)
    clonal <- clonal + 1
    if(length(node) > 1) return((dim(tree)[1] - clonal) / clonal)
  }
  return(0)
}

metrics_extended <- function(data) {
  return(c(metrics(data), ITH = get_ITH_from_tree(data)))
}

cutoff_tree <- function(tree_df, cutoff) {
  pop_cutoff <- cutoff * sum(tree_df$Population)
  for(node in tree_df$Identity) {
    parent_node <- unlist(tree_df[which(tree_df$Identity == node), "Parent"])
    node_desc <- unlist(tree_df[which(tree_df$Identity == node), "Descendants"])
    parent_desc <- unlist(tree_df[which(tree_df$Identity == parent_node), "Descendants"])
    parent_pop <- unlist(tree_df[which(tree_df$Identity == parent_node), "Population"])
    if(node_desc < pop_cutoff && parent_desc >= pop_cutoff)
      tree_df[which(tree_df$Identity == parent_node), "Population"] <- parent_pop + node_desc
  }
  tree_df <- filter(tree_df, Descendants >= pop_cutoff) %>% 
    select(Parent, Identity, Population, Descendants)
  return(tree_df)
}

######

read_one_set <- function(folder, pattern, dataset, minimal) {
  lapply(list.files(folder, pattern = pattern), 
         function(x) {
           read.csv(file = paste0(folder, "/", x), stringsAsFactors=FALSE) %>% 
             mutate(frequency = Population / sum(Population), 
                    tumour = gsub("\\..*", "", x),
                    dataset = dataset,
                    minimal = minimal)
         }
  ) %>% 
    rbindlist
}

###### Read and process real tumour tree data:

r1 <- read_one_set("RealTumourTreesData/TRACERx_Trees", "K....csv", "kidney", 0)
r2 <- read_one_set("RealTumourTreesData/TRACERx_Trees/Minimal", "K....csv", "kidney", 1)
r3 <- read_one_set("RealTumourTreesData/TRACERx_Trees", "CRUK", "lung", 0)
r4 <- read_one_set("RealTumourTreesData/TRACERx_Trees/Minimal", "CRUK", "lung", 1)
r5 <- read_one_set("RealTumourTreesData/YatesEtAl_Trees", "PD", "breast", 0)
r6 <- read_one_set("RealTumourTreesData/YatesEtAl_Trees/Minimal", "PD", "breast", 1)
real_trees <- rbind(r1, r2, r3, r4, r5, r6)

K153 <- read.csv("RealTumourTreesData/TRACERx_Trees/K153.csv", stringsAsFactors=FALSE)
K255 <- read.csv("RealTumourTreesData/TRACERx_Trees/K255.csv", stringsAsFactors=FALSE)
K448 <- read.csv("RealTumourTreesData/TRACERx_Trees/K448.csv", stringsAsFactors=FALSE)
K252 <- read.csv("RealTumourTreesData/TRACERx_Trees/K252.csv", stringsAsFactors=FALSE)
K136 <- read.csv("RealTumourTreesData/TRACERx_Trees/K136.csv", stringsAsFactors=FALSE)
K153min <- read.csv("RealTumourTreesData/TRACERx_Trees/Minimal/K153.csv", stringsAsFactors=FALSE)
K255min <- read.csv("RealTumourTreesData/TRACERx_Trees/Minimal/K255.csv", stringsAsFactors=FALSE)
K448min <- read.csv("RealTumourTreesData/TRACERx_Trees/Minimal/K448.csv", stringsAsFactors=FALSE)
K252min <- read.csv("RealTumourTreesData/TRACERx_Trees/Minimal/K252.csv", stringsAsFactors=FALSE)
K136min <- read.csv("RealTumourTreesData/TRACERx_Trees/Minimal/K136.csv", stringsAsFactors=FALSE)

CRUK0029 <- read.csv("RealTumourTreesData/TRACERx_Trees/CRUK0029.csv", stringsAsFactors=FALSE)
CRUK0062 <- read.csv("RealTumourTreesData/TRACERx_Trees/CRUK0062.csv", stringsAsFactors=FALSE)
CRUK0065 <- read.csv("RealTumourTreesData/TRACERx_Trees/CRUK0065.csv", stringsAsFactors=FALSE)
CRUK0071 <- read.csv("RealTumourTreesData/TRACERx_Trees/CRUK0071.csv", stringsAsFactors=FALSE)
CRUK0096 <- read.csv("RealTumourTreesData/TRACERx_Trees/CRUK0096.csv", stringsAsFactors=FALSE)
CRUK0029min <- read.csv("RealTumourTreesData/TRACERx_Trees/Minimal/CRUK0029.csv", stringsAsFactors=FALSE)
CRUK0062min <- read.csv("RealTumourTreesData/TRACERx_Trees/Minimal/CRUK0062.csv", stringsAsFactors=FALSE)
CRUK0065min <- read.csv("RealTumourTreesData/TRACERx_Trees/Minimal/CRUK0065.csv", stringsAsFactors=FALSE)
CRUK0071min <- read.csv("RealTumourTreesData/TRACERx_Trees/Minimal/CRUK0071.csv", stringsAsFactors=FALSE)
CRUK0096min <- read.csv("RealTumourTreesData/TRACERx_Trees/Minimal/CRUK0096.csv", stringsAsFactors=FALSE)

PD9849 <- read.csv("RealTumourTreesData/YatesEtAl_Trees/PD9849.csv", stringsAsFactors=FALSE)
PD9694 <- read.csv("RealTumourTreesData/YatesEtAl_Trees/PD9694.csv", stringsAsFactors=FALSE)
PD9852 <- read.csv("RealTumourTreesData/YatesEtAl_Trees/PD9852.csv", stringsAsFactors=FALSE)
PD9849min <- read.csv("RealTumourTreesData/YatesEtAl_Trees/Minimal/PD9849.csv", stringsAsFactors=FALSE)
PD9694min <- read.csv("RealTumourTreesData/YatesEtAl_Trees/Minimal/PD9694.csv", stringsAsFactors=FALSE)
PD9852min <- read.csv("RealTumourTreesData/YatesEtAl_Trees/Minimal/PD9852.csv", stringsAsFactors=FALSE)

AML02 <- read.csv("RealTumourTreesData/AML_Trees/AML-02-001.csv", stringsAsFactors=FALSE)
AML05 <- read.csv("RealTumourTreesData/AML_Trees/AML-05-001.csv", stringsAsFactors=FALSE)
AML16 <- read.csv("RealTumourTreesData/AML_Trees/AML-16-001.csv", stringsAsFactors=FALSE)
AML33 <- read.csv("RealTumourTreesData/AML_Trees/AML-33-001.csv", stringsAsFactors=FALSE)
AML35 <- read.csv("RealTumourTreesData/AML_Trees/AML-35-001.csv", stringsAsFactors=FALSE)
AML55 <- read.csv("RealTumourTreesData/AML_Trees/AML-55-001.csv", stringsAsFactors=FALSE)
AML73 <- read.csv("RealTumourTreesData/AML_Trees/AML-73-001.csv", stringsAsFactors=FALSE)
AML77 <- read.csv("RealTumourTreesData/AML_Trees/AML-77-001.csv", stringsAsFactors=FALSE)

UMM059 <- read.csv("RealTumourTreesData/DuranteEtAl_Trees/UMM059.csv", stringsAsFactors=FALSE)
UMM061 <- read.csv("RealTumourTreesData/DuranteEtAl_Trees/UMM061.csv", stringsAsFactors=FALSE)
UMM062 <- read.csv("RealTumourTreesData/DuranteEtAl_Trees/UMM062.csv", stringsAsFactors=FALSE)
UMM063 <- read.csv("RealTumourTreesData/DuranteEtAl_Trees/UMM063.csv", stringsAsFactors=FALSE)
UMM064 <- read.csv("RealTumourTreesData/DuranteEtAl_Trees/UMM064.csv", stringsAsFactors=FALSE)
UMM065 <- read.csv("RealTumourTreesData/DuranteEtAl_Trees/UMM065.csv", stringsAsFactors=FALSE)
UMM066 <- read.csv("RealTumourTreesData/DuranteEtAl_Trees/UMM066.csv", stringsAsFactors=FALSE)
UMM069 <- read.csv("RealTumourTreesData/DuranteEtAl_Trees/UMM069.csv", stringsAsFactors=FALSE)

MED001 <- read.csv("RealTumourTreesData/ZhangEtAl_Trees/MED001.csv", stringsAsFactors=FALSE)
MED012 <- read.csv("RealTumourTreesData/ZhangEtAl_Trees/MED012.csv", stringsAsFactors=FALSE)
MED023 <- read.csv("RealTumourTreesData/ZhangEtAl_Trees/MED023.csv", stringsAsFactors=FALSE)
MED024 <- read.csv("RealTumourTreesData/ZhangEtAl_Trees/MED024.csv", stringsAsFactors=FALSE)
MED027 <- read.csv("RealTumourTreesData/ZhangEtAl_Trees/MED027.csv", stringsAsFactors=FALSE)
MED034 <- read.csv("RealTumourTreesData/ZhangEtAl_Trees/MED034.csv", stringsAsFactors=FALSE)

TN1 <- read.csv("RealTumourTreesData/MinussiEtAl_Trees/TN1.csv", stringsAsFactors=FALSE)
TN2 <- read.csv("RealTumourTreesData/MinussiEtAl_Trees/TN2.csv", stringsAsFactors=FALSE)
TN3 <- read.csv("RealTumourTreesData/MinussiEtAl_Trees/TN3.csv", stringsAsFactors=FALSE)
TN4 <- read.csv("RealTumourTreesData/MinussiEtAl_Trees/TN4.csv", stringsAsFactors=FALSE)
TN5 <- read.csv("RealTumourTreesData/MinussiEtAl_Trees/TN5.csv", stringsAsFactors=FALSE)
TN6 <- read.csv("RealTumourTreesData/MinussiEtAl_Trees/TN6.csv", stringsAsFactors=FALSE)
TN7 <- read.csv("RealTumourTreesData/MinussiEtAl_Trees/TN7.csv", stringsAsFactors=FALSE)
TN8 <- read.csv("RealTumourTreesData/MinussiEtAl_Trees/TN8.csv", stringsAsFactors=FALSE)

real_points <- as.data.frame(rbind(metrics_extended(K153),
                                   metrics_extended(K255),
                                   metrics_extended(K448),
                                   metrics_extended(K252),
                                   metrics_extended(K136),
                                   metrics_extended(K153min),
                                   metrics_extended(K255min),
                                   metrics_extended(K448min),
                                   metrics_extended(K252min),
                                   metrics_extended(K136min),
                                   metrics_extended(CRUK0029),
                                   metrics_extended(CRUK0062),
                                   metrics_extended(CRUK0065),
                                   metrics_extended(CRUK0071),
                                   metrics_extended(CRUK0096),
                                   metrics_extended(CRUK0029min),
                                   metrics_extended(CRUK0062min),
                                   metrics_extended(CRUK0065min),
                                   metrics_extended(CRUK0071min),
                                   metrics_extended(CRUK0096min),
                                   metrics_extended(PD9852),
                                   metrics_extended(PD9849),
                                   metrics_extended(PD9694),
                                   metrics_extended(PD9852min),
                                   metrics_extended(PD9849min),
                                   metrics_extended(PD9694min),
                                   metrics_extended(AML02),
                                   metrics_extended(AML05),
                                   metrics_extended(AML16),
                                   metrics_extended(AML33),
                                   metrics_extended(AML35),
                                   metrics_extended(AML55),
                                   metrics_extended(AML73),
                                   metrics_extended(AML77),
                                   metrics_extended(UMM059),
                                   metrics_extended(UMM061),
                                   metrics_extended(UMM062),
                                   metrics_extended(UMM063),
                                   metrics_extended(UMM064),
                                   metrics_extended(UMM065),
                                   metrics_extended(UMM066),
                                   metrics_extended(UMM069),
                                   metrics_extended(MED001),
                                   metrics_extended(MED012),
                                   metrics_extended(MED023),
                                   metrics_extended(MED024),
                                   metrics_extended(MED027),
                                   metrics_extended(MED034),
                                   metrics_extended(TN1),
                                   metrics_extended(TN2),
                                   metrics_extended(TN3),
                                   metrics_extended(TN4),
                                   metrics_extended(TN5),
                                   metrics_extended(TN6),
                                   metrics_extended(TN7),
                                   metrics_extended(TN8)
))
real_points$tumour <- c(rep(c("K153", "K255", "K448", "K252", "K136"), 2), 
                        rep(c("CRUK0029", "CRUK0062", "CRUK0065", "CRUK0071", "CRUK0096"), 2), 
                        rep(c("PD9852", "PD9849", "PD9694"), 2), 
                        c("AML-02-001","AML-05-001","AML-16-001","AML-33-001","AML-35-001","AML-55-001","AML-73-001","AML-77-001"), 
                        c("UMM059","UMM061","UMM062","UMM063","UMM064","UMM065","UMM066","UMM069"), 
                        c("MED001","MED012","MED023","MED024","MED027","MED034"), 
                        c("TN1","TN2","TN3","TN4","TN5","TN6","TN7","TN8"))
real_points$tumourshort <- c(rep(c("K153", "K255", "K448", "K252", "K136"), 2), 
                             rep(c("C29", "C62", "C65", "C71", "C96"), 2), 
                             rep(c("P852", "P849", "P694"), 2), 
                             c("A02","A05","A16","A33","A35","A55","A73","A77"), 
                             c("U59","U61","U62","U63","U64","U65","U66","U69"), 
                             c("M01","M12","M23","M24","M27","M34"), 
                             c("TN1","TN2","TN3","TN4","TN5","TN6","TN7","TN8"))
real_points$dataset <- c(rep("kidney", 10), rep("lung", 10), rep("breast", 6), rep("AML", 8), rep("uveal", 8), rep("mesothelioma", 6), rep("breast_SC", 8))
real_points$minimal <- c(rep(0, 5), rep(1, 5), rep(0, 5), rep(1, 5), rep(0, 3), rep(1, 3), rep(0, 8), rep(0, 8), rep(0, 6), rep(0, 8))

##### plot real tumour trees:

# Figure S10:
par(mfrow = c(3,5))
plot_tree(K153min) # stage IV
plot_tree(K255min) # stage III
plot_tree(K448min) # stage IV
plot_tree(K252min) # stage III
plot_tree(K136min) # stage II
plot_tree(CRUK0029min)
plot_tree(CRUK0062min)
plot_tree(CRUK0065min)
plot_tree(CRUK0071min)
plot_tree(CRUK0096min)
plot_tree(PD9852min)
plot_tree(PD9849min)
plot_tree(PD9694min)

# Figure S11:
par(mfrow = c(3,5))
plot_tree(K153) # stage IV
plot_tree(K255) # stage III
plot_tree(K448) # stage IV
plot_tree(K252) # stage III
plot_tree(K136) # stage II
plot_tree(CRUK0029)
plot_tree(CRUK0062)
plot_tree(CRUK0065)
plot_tree(CRUK0071)
plot_tree(CRUK0096)
plot_tree(PD9852)
plot_tree(PD9849)
plot_tree(PD9694)
plot.new()
plot.new()

######## read simulated driver phylogenies:

caseA <- read_table2("ModelOutput/NonspatialModel/output_driver_genotype_properties.dat")
caseB <- read_table2("ModelOutput/FissionModelK8192/output_driver_genotype_properties.dat")
caseC <- read_table2("ModelOutput/InvasiveK512/output_driver_genotype_properties.dat")
caseD <- read_table2("ModelOutput/EdenModel/output_driver_genotype_properties.dat")

caseC0 <- read_table2("ModelOutput/InvasiveK512_extraseeds/seed_0/output_driver_genotype_properties.dat")
caseC1 <- read_table2("ModelOutput/InvasiveK512_extraseeds/seed_1/output_driver_genotype_properties.dat")
caseC2 <- read_table2("ModelOutput/InvasiveK512_extraseeds/seed_2/output_driver_genotype_properties.dat")
caseC3 <- read_table2("ModelOutput/InvasiveK512_extraseeds/seed_3/output_driver_genotype_properties.dat")
caseC4 <- read_table2("ModelOutput/InvasiveK512_extraseeds/seed_4/output_driver_genotype_properties.dat")

case_supp5 <- read_table2("ModelOutput/supplementary/K8192_fission_migrationevolves_seed18_output_driver_genotype_properties.dat")
case_supp2 <- read_table2("ModelOutput/supplementary/K512_invasion_seed51_output_driver_genotype_properties.dat")
case_supp1 <- read_table2("ModelOutput/supplementary/K1_voter_seed91_output_driver_genotype_properties.dat")
case_supp3 <- read_table2("ModelOutput/supplementary/K2048_invasion_seed21_output_driver_genotype_properties.dat")
case_supp4 <- read_table2("ModelOutput/supplementary/K2048_fission_seed64_output_driver_genotype_properties.dat")

# filter for only the essential columns:
edgesA <- filter(caseA, Descendants >= 1e6 / 100) %>% dplyr::select(Parent, Identity, Population)
edgesB <- filter(caseB, Descendants >= 1e6 / 100) %>% dplyr::select(Parent, Identity, Population)
edgesC <- filter(caseC, Descendants >= 1e6 / 100) %>% dplyr::select(Parent, Identity, Population)
edgesD <- filter(caseD, Descendants >= 1e6 / 100) %>% dplyr::select(Parent, Identity, Population)

edges_supp5 <- filter(case_supp5, Descendants >= 1e6 / 100) %>% dplyr::select(Parent, Identity, Population)
edges_supp2 <- filter(case_supp2, Descendants >= 1e6 / 100) %>% dplyr::select(Parent, Identity, Population)
edges_supp1 <- filter(case_supp1, Descendants >= 1e6 / 100) %>% dplyr::select(Parent, Identity, Population)
edges_supp3 <- filter(case_supp3, Descendants >= 1e6 / 100) %>% dplyr::select(Parent, Identity, Population)
edges_supp4 <- filter(case_supp4, Descendants >= 1e6 / 100) %>% dplyr::select(Parent, Identity, Population)

edgesC0_cutoff05pc <- cutoff_tree(caseC0, 0.005) # filter(caseC0, Descendants >= 1e6 / 200) %>% dplyr::select(Parent, Identity, Population)
edgesC1_cutoff05pc <- cutoff_tree(caseC1, 0.005) # filter(caseC1, Descendants >= 1e6 / 200) %>% dplyr::select(Parent, Identity, Population)
edgesC2_cutoff05pc <- cutoff_tree(caseC2, 0.005) # filter(caseC2, Descendants >= 1e6 / 200) %>% dplyr::select(Parent, Identity, Population)
edgesC3_cutoff05pc <- cutoff_tree(caseC3, 0.005) # filter(caseC3, Descendants >= 1e6 / 200) %>% dplyr::select(Parent, Identity, Population)
edgesC4_cutoff05pc <- cutoff_tree(caseC4, 0.005) # filter(caseC4, Descendants >= 1e6 / 200) %>% dplyr::select(Parent, Identity, Population)
edgesC0_cutoff1pc <- cutoff_tree(caseC0, 0.01) # filter(caseC0, Descendants >= 1e6 / 100) %>% dplyr::select(Parent, Identity, Population)
edgesC1_cutoff1pc <- cutoff_tree(caseC1, 0.01) # filter(caseC1, Descendants >= 1e6 / 100) %>% dplyr::select(Parent, Identity, Population)
edgesC2_cutoff1pc <- cutoff_tree(caseC2, 0.01) # filter(caseC2, Descendants >= 1e6 / 100) %>% dplyr::select(Parent, Identity, Population)
edgesC3_cutoff1pc <- cutoff_tree(caseC3, 0.01) # filter(caseC3, Descendants >= 1e6 / 100) %>% dplyr::select(Parent, Identity, Population)
edgesC4_cutoff1pc <- cutoff_tree(caseC4, 0.01) # filter(caseC4, Descendants >= 1e6 / 100) %>% dplyr::select(Parent, Identity, Population)
edgesC0_cutoff2pc <- cutoff_tree(caseC0, 0.02) # filter(caseC0, Descendants >= 1e6 / 50) %>% dplyr::select(Parent, Identity, Population)
edgesC1_cutoff2pc <- cutoff_tree(caseC1, 0.02) # filter(caseC1, Descendants >= 1e6 / 50) %>% dplyr::select(Parent, Identity, Population)
edgesC2_cutoff2pc <- cutoff_tree(caseC2, 0.02) # filter(caseC2, Descendants >= 1e6 / 50) %>% dplyr::select(Parent, Identity, Population)
edgesC3_cutoff2pc <- cutoff_tree(caseC3, 0.02) # filter(caseC3, Descendants >= 1e6 / 50) %>% dplyr::select(Parent, Identity, Population)
edgesC4_cutoff2pc <- cutoff_tree(caseC4, 0.02) # filter(caseC4, Descendants >= 1e6 / 50) %>% dplyr::select(Parent, Identity, Population)
edgesC0_cutoff5pc <- cutoff_tree(caseC0, 0.05) # filter(caseC0, Descendants >= 1e6 / 20) %>% dplyr::select(Parent, Identity, Population)
edgesC1_cutoff5pc <- cutoff_tree(caseC1, 0.05) # filter(caseC1, Descendants >= 1e6 / 20) %>% dplyr::select(Parent, Identity, Population)
edgesC2_cutoff5pc <- cutoff_tree(caseC2, 0.05) # filter(caseC2, Descendants >= 1e6 / 20) %>% dplyr::select(Parent, Identity, Population)
edgesC3_cutoff5pc <- cutoff_tree(caseC3, 0.05) # filter(caseC3, Descendants >= 1e6 / 20) %>% dplyr::select(Parent, Identity, Population)
edgesC4_cutoff5pc <- cutoff_tree(caseC4, 0.05) # filter(caseC4, Descendants >= 1e6 / 20) %>% dplyr::select(Parent, Identity, Population)
edgesC0_cutoff10pc <- cutoff_tree(caseC0, 0.1) # filter(caseC0, Descendants >= 1e6 / 10) %>% dplyr::select(Parent, Identity, Population)
edgesC1_cutoff10pc <- cutoff_tree(caseC1, 0.1) # filter(caseC1, Descendants >= 1e6 / 10) %>% dplyr::select(Parent, Identity, Population)
edgesC2_cutoff10pc <- cutoff_tree(caseC2, 0.1) # filter(caseC2, Descendants >= 1e6 / 10) %>% dplyr::select(Parent, Identity, Population)
edgesC3_cutoff10pc <- cutoff_tree(caseC3, 0.1) # filter(caseC3, Descendants >= 1e6 / 10) %>% dplyr::select(Parent, Identity, Population)
edgesC4_cutoff10pc <- cutoff_tree(caseC4, 0.1) # filter(caseC4, Descendants >= 1e6 / 10) %>% dplyr::select(Parent, Identity, Population)

####### plot simulated trees:

# For Figure 2:
plot_tree(edgesA) # first row
plot_tree(edgesB) # second row
plot_tree(edgesC) # third row
plot_tree(edgesD) # fourth row

# For Extended Data Figure 2:
plot_tree(edges_supp5) # first row
plot_tree(edges_supp2) # second row
plot_tree(edges_supp1) # third row
plot_tree(edges_supp3) # fourth row
plot_tree(edges_supp4) # fifth row

# Extended data figure 7:
par(mfrow = c(5,5))
plot_tree(edgesC0_cutoff05pc)
plot_tree(edgesC1_cutoff05pc)
plot_tree(edgesC2_cutoff05pc)
plot_tree(edgesC3_cutoff05pc)
plot_tree(edgesC4_cutoff05pc)
plot_tree(edgesC0_cutoff1pc)
plot_tree(edgesC1_cutoff1pc)
plot_tree(edgesC2_cutoff1pc)
plot_tree(edgesC3_cutoff1pc)
plot_tree(edgesC4_cutoff1pc)
plot_tree(edgesC0_cutoff2pc)
plot_tree(edgesC1_cutoff2pc)
plot_tree(edgesC2_cutoff2pc)
plot_tree(edgesC3_cutoff2pc)
plot_tree(edgesC4_cutoff2pc)
plot_tree(edgesC0_cutoff5pc)
plot_tree(edgesC1_cutoff5pc)
plot_tree(edgesC2_cutoff5pc)
plot_tree(edgesC3_cutoff5pc)
plot_tree(edgesC4_cutoff5pc)
plot_tree(edgesC0_cutoff10pc)
plot_tree(edgesC1_cutoff10pc)
plot_tree(edgesC2_cutoff10pc)
plot_tree(edgesC3_cutoff10pc)
plot_tree(edgesC4_cutoff10pc)

########### read and process simulation data:

dataForMetricPlots<-read_csv("dataForMetricPlots.csv", guess_max = 1E4)
combined_cases <- read.csv("DivMutation_Allcombined_cases.csv",sep="")

InvasiveGlandular_K512_EdgeOnly0 <- subset(dataForMetricPlots, 
                                           dataForMetricPlots$K==512 & 
                                             dataForMetricPlots$s_driver_birth==0.1 & 
                                             dataForMetricPlots$mu_driver_birth==10^-5) %>% 
  mutate(case = "caseC") %>% 
  select(Drivers, DriverDiversity, case)
colnames(InvasiveGlandular_K512_EdgeOnly0)[which(colnames(InvasiveGlandular_K512_EdgeOnly0)=="DriverDiversity")]<-"Diversity"

combined_cases$case <- as.character(combined_cases$case)

combined_cases <- combined_cases %>% mutate(in_chart = ifelse((migration_type == 2 & migration_edge_only == 0 & K > 1 & s_driver_birth > 0) |
                                                                (migration_type == 0 & normal_birth_rate == 0.9 & s_driver_birth > 0), 0.5, 0))

combined_cases <- combined_cases %>% mutate(in_chart = ifelse((K == 1 & migration_type == 2 & migration_edge_only == 1)
                                                              | (K == 2^20 & migration_type == 0 & migration_edge_only == 0 & baseline_death_rate == 0.98)
                                                              | (K == 2^13 & migration_type == 2 & migration_edge_only == 0 & s_driver_migration == 0), 1, in_chart))
combined_cases <- combined_cases %>% mutate(in_chart = case_when((s_driver_migration == 0.5 & s_driver_birth > 0 & migration_edge_only == 0) ~ 11, 
                                                                 (K == 512 & migration_type == 0 & migration_edge_only == 0 & s_driver_birth > 0) ~ 12,
                                                                 (K == 1 & migration_type == 0 & migration_edge_only == 0 & s_driver_birth > 0) ~ 13,
                                                                 (K == 2048 & migration_type == 0 & migration_edge_only == 1 & s_driver_birth > 0) ~ 14,
                                                                 (K == 2048 & migration_type == 2 & migration_edge_only == 0 & s_driver_birth > 0) ~ 15,
                                                                 TRUE ~ in_chart))

combined_cases <- combined_cases %>% mutate(case = case_when(s_driver_birth == 0 ~ "neutral", 
                                                             case == "caseC" & migration_edge_only == 0 ~ "caseE",
                                                             TRUE ~ case))
chart_df_1 <- filter(combined_cases, in_chart == 1)

DataMetricPlot_InvasiveGlandularK512<-rbind(select(chart_df_1,Drivers,Diversity, case ),
                                            select(InvasiveGlandular_K512_EdgeOnly0,Drivers,Diversity, case ))

###########

topx <- 25
topy <- 400
get_n_and_D <- function(ff) {
  ff <- sort(ff / sum(ff), decreasing = TRUE)
  n <- sum(ff*(1:length(ff)))
  D <- 1/sum(ff^2)
  return(c(n = n, D = D))
}
intermediate_curve <- function(n) {
  if(n == 1) return(1)
  q <- ceiling(3 * n - 2)
  return(2/((q-1)*q) * (2*q-3*n+1+3*(2*n-q-1)*(1:q)/(q+1)))
}
sweeps_func <- function(n) 1 / ((n %% 1 - 1)^2 + (n %% 1)^2)
curve_df <- data.frame(linex = seq(1, topx, length = 19000))
curve_df$maxy <- ifelse(curve_df$linex < 2, 1 / (2 - curve_df$linex)^2, NA)
curve_df$sweeps <- sapply(curve_df$linex, sweeps_func)
curve_df$int_curve <- sapply(sapply(curve_df$linex, intermediate_curve), get_n_and_D)["D",]

# Figure 3b:
minim <- 0
ggplot() +
  geom_ribbon(aes(x = linex, ymin = maxy, ymax = topy), curve_df, lty = 1, fill = "grey90") +
  scale_y_log10(limits = c(1, topy), name = expression(paste("Clonal diversity, ", italic("D"))), expand = c(0, 0)) +
  scale_x_log10(limits = c(1, topx), name = expression(paste("Mean driver mutations per cell, ", italic("n"))), expand = c(0, 0)) +
  theme_classic() +
  geom_line(aes(x = linex, y = sweeps), curve_df, lty = 1, color = "pink") +
  geom_line(aes(x = linex, y = int_curve), curve_df, lty = 1, color = "skyblue") +
  geom_line(aes(x = linex, y = maxy), curve_df, lty = 1, color = "grey") +
  geom_point(aes(x = n, y = D, group = dataset), filter(real_points, minimal == minim, dataset == "kidney")) +
  geom_point(aes(x = n, y = D, group = dataset), filter(real_points, minimal == minim, dataset == "lung")) +#, shape = 15) +
  geom_point(aes(x = n, y = D, group = dataset), filter(real_points, minimal == minim, dataset == "breast")) +#, shape = 17) + 
  geom_point(aes(x = n, y = D, group = dataset), filter(real_points, minimal == minim, dataset == "AML"), shape = 15, col = "purple") +
  scale_color_manual(values = c("red", "yellow4", "dodgerblue", "chocolate", "cyan"), name = "deme\ncarrying\ncapacity") +
  theme(legend.position="none") + 
  geom_point(aes(x = n, y = D, group = dataset), filter(real_points, minimal == minim, dataset == "uveal"), col = "black", shape =1) +
  geom_point(aes(x = n, y = D, group = dataset), filter(real_points, minimal == minim, dataset == "mesothelioma"), col = "black", shape =1) +
  geom_point(aes(x = n, y = D, group = dataset), filter(real_points, minimal == minim, dataset == "breast_SC"), col = "black", shape =1) +
  geom_text_repel(aes(x = n, y = D, label = tumourshort), filter(real_points, minimal == minim, dataset != "AML", dataset != "thyroid"), 
                  size = 2.5/2, box.padding = 0.15/1.5) +
  geom_text_repel(aes(x = n, y = D, label = tumourshort), filter(real_points, minimal == minim, dataset == "AML"), 
                  size = 2.5/2, box.padding = 0.15/1.5, col = "purple")

# Figure 3c:
ggplot() +
  geom_ribbon(aes(x = linex, ymin = maxy, ymax = topy), curve_df, lty = 1, fill = "grey90") + 
  geom_line(aes(x = linex, y = sweeps), curve_df, lty = 1, color = "grey") +
  geom_line(aes(x = linex, y = toplinear_yannick), curve_df, lty = 1, color = "grey") +
  geom_line(aes(x = linex, y = maxy), curve_df, lty = 1, color = "grey") +
  geom_point(aes(x=Drivers+1, y=Diversity, color=as.factor(case), shape =as.factor(case)), size=2 , subset(DataMetricPlot_InvasiveGlandularK512, DataMetricPlot_InvasiveGlandularK512$case=="caseD_new" ) ) +
  geom_point(aes(x=Drivers+1, y=Diversity, color=as.factor(case), shape =as.factor(case)), size=2 , subset(DataMetricPlot_InvasiveGlandularK512, ! DataMetricPlot_InvasiveGlandularK512$case=="caseD_new" ) ) +
  scale_y_log10(limits = c(1, topy), name = "Clonal diversity", expand = c(0, 0)) +
  scale_x_log10(limits = c(1, topx), name = "Mean driver mutations per cell", expand = c(0, 0)) +
  theme_classic(base_size = 18)+
  theme(aspect.ratio=1)+
  guides(fill=FALSE,shape=FALSE,color= FALSE)+
  scale_color_manual(values = c( "caseA"="red", "caseB"="yellow4", "caseC"="dodgerblue", 
                                 "caseD_new"="chocolate", "neutral"="cyan"), name = "Oncoevotypes", 
                     labels=c("caseA"="non-spatial", "caseB"="gland fission", "caseC"="invasive glandular", 
                              "caseD_new"="boundary growth", "neutral"="neutral")) +
  scale_shape_manual(values=c(0, 1, 2,3,4), name = "Oncoevotypes", labels=c("non-spatial", "gland fission", 
                                                                            "invasive glandular", "boundary growth", "neutral"))

########

# Extended data figure 10i-l:
par(mar = c(4,4,1,1))
par(mfrow = c(1, 4))
plot_counts("ModelOutput/NonspatialModel/output_allele_counts.dat", ylim = c(0, 100), xlab = "Frequency")
plot_counts("ModelOutput/FissionModelK8192/output_allele_counts.dat", ylim = c(0, 100), xlab = "Frequency")
plot_counts("ModelOutput/InvasiveK512/output_allele_counts.dat", ylim = c(0, 100), xlab = "Frequency")
plot_counts("ModelOutput/EdenModel/output_allele_counts.dat", ylim = c(0, 100), xlab = "Frequency")

########

image_df_from_grid_file_alt <- function (file, trim = -1, as_matrix = FALSE) 
{
  if (!file.exists(file)) {
    warning(paste0(file, " not found"))
    return(NA)
  }
  res <- read_delim(file, "\t", escape_double = FALSE, col_names = FALSE, 
                    trim_ws = TRUE)
  res <- data.matrix(res)
  res <- res[, -ncol(res)]
  if (length(res) > 1) {
    if (trim < 0) {
      res <- res[, colSums(is.na(res)) < nrow(res)]
      res <- res[rowSums(is.na(res)) < ncol(res), ]
    }
    else res <- res[(trim + 1):(nrow(res) - trim), (trim + 
                                                      1):(ncol(res) - trim)]
  } else {
    return(data.frame(x = 1, y = 1, z = 1))
  }
  if (as_matrix) 
    return(res)
  df <- expand.grid(x = 1:ncol(as.data.frame(res)), y = 1:nrow(as.data.frame(res)))
  df$z <- as.vector(t(res))
  return(df)
}

plot_figure2_alt <- function (path, output_filename = NA, file_type = "png", output_dir = NA, 
                              trim = -1, cutoff = 0) 
{
  if (substr(path, nchar(path), nchar(path)) != "/") 
    path <- paste0(path, "/")
  if (!is.na(output_dir)) 
    if (substr(output_dir, nchar(output_dir), nchar(output_dir)) != 
        "/") 
      output_dir <- paste0(output_dir, "/")
  
  image_df <- image_df_from_grid_file_alt(paste0(path, "/output_driversgrid.dat"), trim)
  output_genotype_properties <- read_delim_special(paste0(path, "/output_genotype_properties.dat"))
  driver_phylo <- read_delim_special(paste0(path, "/driver_phylo.dat"))
  
  t1 <- unique(driver_phylo$DriverIdentity)
  t2 <- unique(output_genotype_properties$DriverIdentity)
  t3 <- as.numeric(unique(image_df$z))
  drivers_to_use <- 0:max(t1, t2, t3, na.rm = TRUE)
  
  long_palette <- c("#8A7C64", "#599861", "#89C5DA", "#DA5724", "#74D944", "#CE50CA", 
                    "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", 
                    "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
                    "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", 
                    "#5E738F", "#D1A33D")
  long_palette <- rep(long_palette, ceiling(length(drivers_to_use)/length(long_palette)))
  dd <- drivers_to_use
  dd.col <- long_palette
  names(dd.col) <- dd
  
  muller_df <- muller_df_from_file(paste0(path, "/driver_phylo.dat"), cutoff = cutoff)
  if (class(muller_df) != "data.frame") return(NA)
  muller_df_filtered <- filter(muller_df, DriverIdentity %in% drivers_to_use)
  g1 <- Muller_plot(muller_df_filtered, colour_by = "Identity", palette = dd.col) +
    theme_void() + theme(legend.position = "none")
  
  image_df[which(image_df$z > 0), "z"] <- as.character(image_df[which(image_df$z > 0), "z"])
  g2 <- grid_plot(image_df, palette = dd.col, discrete = TRUE)
  
  output_genotype_properties_filtered <- filter(output_genotype_properties, DriverIdentity %in% drivers_to_use)
  output_genotype_properties_filtered <- filter(output_genotype_properties_filtered, Descendants > 100)
  g3 <- plot_allelecount_vs_origintime(output_genotype_properties_filtered,
                                       colour_by = "DriverIdentity",
                                       palette = dd.col,
                                       discrete = TRUE, print_plot = FALSE)  +
    theme_classic() + theme(axis.text.y=element_blank(),
                            axis.text.x = element_blank(),
                            axis.ticks=element_blank(),
                            axis.title.x=element_blank(),
                            axis.title.y=element_blank(),
                            legend.position="none")
  
  if (!is.na(output_filename)) 
    print(paste0("Created all plots for file ", output_filename), 
          quote = FALSE)
  if (!is.na(output_filename) & !is.na(output_dir)) {
    if (file_type == "png") 
      png(paste0(output_dir, output_filename, ".png"), 
          width = 24, height = 6, units = "cm", res=500)
    else pdf(paste0(output_dir, output_filename, ".pdf"), 
             width = 20, height = 7)
  }
  
  lay <- rbind(c(1,1,2,3))
  print(grid.arrange(g1, g2, g3, layout_matrix = lay))
  if (!is.na(output_filename) & !is.na(output_dir)) 
    dev.off()
  if (!is.na(output_filename)) 
    print("Saved the plot", quote = FALSE)
}

# Figure 2 and Extended data figure 10e-h:
# (Eden model omitted here due to prohibitively large file size)
plot_figure2_alt("ModelOutput/NonspatialModel")
plot_figure2_alt("ModelOutput/FissionModelK8192")
plot_figure2_alt("ModelOutput/InvasiveK512")

##########

plot_all_images_alt <- function(path, output_filename = NA, file_type = "png", output_dir = NA, trim = -1, 
                                cutoff = 0, min_birth_rate = NA, max_birth_rate = NA) {
  if(substr(path, nchar(path), nchar(path)) != "/") path <- paste0(path, "/")
  if(!is.na(output_dir)) if(substr(output_dir, nchar(output_dir), nchar(output_dir)) != "/") output_dir <- paste0(output_dir, "/")
  
  Muller_df <- muller_df_from_file(paste0(path, "driver_phylo.dat"), cutoff = cutoff)
  if(class(Muller_df) != "data.frame") return(NA)
  
  b_grid <- image_df_from_grid_file(paste0(path, "output_birthratesgrid.dat"), trim)
  
  if(is.na(min_birth_rate)) min_birth_rate <- min(c(b_grid$z, Muller_df$BirthRate), na.rm = TRUE)
  if(is.na(max_birth_rate)) max_birth_rate <- max(c(b_grid$z, Muller_df$BirthRate), na.rm = TRUE)
  h3 <- Muller_plot(Muller_df, colour_by = "BirthRate", add_legend = FALSE) + 
    scale_fill_distiller(palette = "RdBu", direction = -1, 
                         limits = c(min_birth_rate, max_birth_rate)) + 
    theme(line = element_blank(), rect = element_blank()) + 
    scale_x_continuous(breaks = c(0, round(max(Muller_df$Generation)))) + 
    scale_y_continuous(breaks = NULL)
  
  g2 <- grid_plot(b_grid, add_legend = TRUE, legend_title = "Mean cell\nproliferation rate   ") + 
    scale_fill_distiller(name = "Mean cell\nproliferation rate   ", palette ="RdBu", 
                         direction = -1, na.value="white", 
                         limits = c(min_birth_rate, max_birth_rate))
  
  image_df <- image_df_from_grid_file(paste0(path, "output_passengersgrid.dat"), trim)
  g3 <- grid_plot(image_df, add_legend = TRUE, legend_title = "Mean passenger\nmutations per cell")
  
  if(!is.na(output_filename)) print(paste0("Created all plots for file ", output_filename), quote = FALSE)
  
  if(!is.na(output_filename) & !is.na(output_dir)) {
    if(file_type == "png") png(paste0(output_dir,output_filename,".png"), width = 1000, height = 180, res = 100)
    else pdf(paste0(output_dir,output_filename,".pdf"), width = 10, height = 1.8)
  }
  lay <- rbind(c(1,2,3))
  print(grid.arrange(h3, g2, g3, layout_matrix = lay))
  if(!is.na(output_filename) & !is.na(output_dir)) dev.off()
  
  if(!is.na(output_filename)) print("Saved the plot", quote = FALSE)
}

# Extended data figure 1:
plot_all_images_alt("/Users/rnoble/Documents/MontpellierDocuments/Data/testdata/caseB_seed_76")
plot_all_images_alt("/Users/rnoble/Documents/MontpellierDocuments/Data/testdata/caseC_edgeonly0_seed51")
plot_all_images_alt("/Users/rnoble/Documents/MontpellierDocuments/Data/testdata/caseD_seed_10")



