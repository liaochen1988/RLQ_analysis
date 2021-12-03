library(ade4)
library(adegraphics)
library(dplyr)

###########################
# Step 1: Read R,L,Q tables
###########################

# R: sites by environmental variables
# in this application, sites = Pseudomonas strains, environmental variables = phenotypes, values are phenotypic measurements
# note that phenotypes can be both categorical (factor) or numeric (double)
tableR <- read.csv(file = 'input/tableR_strain_by_phenotype.csv', row.names = 'StrainID', colClasses=c('factor','factor','factor',rep('double',121)))

# L: sites by species
# in this application, sites = Pseudomonas strains, species = genes, values are normalized gene expression (RNAseq)
tableL <- read.csv(file = 'input/tableL_strain_by_gene.csv', row.names = 'StrainID')

# Q: species by traits
# in this application, species = genes, traits = KEGG pathways, values indicates whether a gene belongs to a pathway
tableQ <- read.csv(file = 'input/tableQ_gene_by_pathway.csv', row.names = 'GeneID', colClasses=c(rep('factor',1)))

#################################################
# Step 2: Perform separate analyses of each table
#################################################

# Correspondence analysis 
coa1 <- dudi.coa(tableL, scannf = FALSE, nf = 2)

# Note seperate analyses of traits (Q) and enviornmental variables (R) should be weighted by the sites and species weights derived form the correspondence analysis of L.

# HillSmith
dudiphe <- dudi.hillsmith(tableR, scannf = FALSE, nf = 2, row.w = coa1$lw)

# Multiple correspondence analysis
# Principle component analysis is also fine (dudi.pca)
# Use hillsmith instead if mixed types of variables
dudigfunc <- dudi.acm(tableQ, scannf = FALSE, nf = 2, row.w = coa1$cw)

##############################
# Step 3: Perform RLQ analysis
##############################

rlq1<- rlq(dudiphe, coa1, dudigfunc, scannf = FALSE, nf = 2)

##############################
# Step 4: Write to files
##############################

write.table(rlq1$lQ, 'output/lQ_scores_of_genes.csv', sep=',')
write.table(rlq1$l1, 'output/l1_loadings_of_phenotypes.csv', sep=',')
write.table(rlq1$c1, 'output/c1_loadings_of_pathways.csv', sep=',')
write.table(rlq1$lR, 'output/lR_scores_of_strains.csv', sep=',')
write.table(rlq1$eig, 'output/eigenvalues.csv', sep=',')
write.table(rlq1$aQ, 'output/aQ_axes_of_pathways_on_RLQ_axes.csv', sep=',')
write.table(rlq1$aR, 'output/aR_axes_of_phenotypes_on_RLQ_axes.csv', sep=',')
