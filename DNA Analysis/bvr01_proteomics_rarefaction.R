
#clearing all objects from workspace
rm(list=ls())

##SET DIRECTORY
#
setwd("C:/Users/Elliot/Desktop/Taxonomy Analysis/")

library('ggplot2')
library("vegan")
library(RColorBrewer)
library(phyloseq)
library("plyr")

scaleyellowred <- colorRampPalette(c("lightyellow", "red"), space = "rgb")(100)

samplet <- read.csv("bvr01_proteomics_samplesheet.csv", header=TRUE, row.names =1, sep=",")

################################
#16S DATA
#################################

otut16 <- read.csv("bvr01_proteomics_abs_16s.csv", header=TRUE,row.names=1, sep=",")
taxat16 <- read.csv("bvr01_proteomics_taxatable_16s.csv", header=TRUE, row.names =1, sep=",")

OTU = otu_table(as.matrix(otut16), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(taxat16))
SAMPLE <- sample_data(samplet)

physeq = phyloseq(OTU, TAX)
physeq
physeq1 = merge_phyloseq(physeq, SAMPLE)

##########################
#Pruning taxa > 0
###########################
physeq2 <- prune_taxa( taxa_sums(physeq1) > 0, physeq1)
physeq2

######################

rank_name <- 'Species'
phy.organism <- phyloseq::tax_glom(physeq2, rank_name)

#Creating relative abundance
phy.organism.rel <- transform_sample_counts(phy.organism, function(OTU) OTU/sum(OTU))

#################
# Alpha diversity plots
#################
p = plot_richness(phy.organism, x="subject_id", color = "is_case", measures=c("Chao1", "Shannon"))
p + geom_point(size=6, alpha=0.2)+ geom_text(mapping = aes(label = subject_id), size = 2, color="black")

##################
# Alpha diversity exports
##################

phy.organism.chao1 = estimate_richness(phy.organism, measures=c("Chao1"))
phy.organism.shannon = estimate_richness(phy.organism, measures=c("Shannon"))

phy.organism.chao1$id <- row.names(phy.organism.chao1)
phy.organism.shannon$id <- row.names(phy.organism.shannon)

phy.organism.shannon$id2 = gsub("-", ".", phy.organism.shannon$id)
phy.organism.chao1$id2 = phy.organism.chao1$id

phy.organism.alphad.merge <- merge(phy.organism.chao1,phy.organism.shannon, by = "id2")

write.csv(phy.organism.alphad.merge, file = "phy_organism_alphad_merge_16s.csv")

##############################################

##################
#Beta diversity - MDS with Bray distance
##################

phy.subset.dist <- phy.organism.rel
braydist <- distance(phy.subset.dist, method="bray")
braymds  <- ordinate(phy.subset.dist, "MDS", distance=braydist)
pdist1 <- NULL
pdist1 <- plot_ordination(phy.subset.dist, braymds, color="is_case")
pdist1 <- pdist1 + ggtitle(paste("MDS using bray distance"))
pdist1 + geom_point(size = 7)+ geom_text(mapping = aes(label = method_manuscript_id), size = 2, color="black")

##################
#Beta diversity - Multiple methods
##################

#####
# Relative abundance
#####
phy.organism.ord.rel <- ordinate(phy.organism.rel, "NMDS", "euclidean")
p4.rel = plot_ordination(phy.organism.rel, phy.organism.ord.rel, type="samples", color="is_case", shape="vag_fluid")
p4.rel + geom_point(size=7) + geom_text(mapping = aes(label = subject_id), size = 2, color="black")+ ggtitle("samples")
######

dist = "bray"
ord_meths = c("DCA", "CCA", "RDA", "NMDS", "MDS", "PCoA")
plist = llply(as.list(ord_meths), function(i, phy.organism, dist){
  ordi = ordinate(phy.organism, method=i, distance=dist)
  plot_ordination(phy.organism, ordi, "samples", color="is_case")
}, phy.organism, dist)

names(plist) <- ord_meths

pdataframe = ldply(plist, function(x){
  df = x$data[, 1:2]
  colnames(df) = c("Axis_1", "Axis_2")
  return(cbind(df, x$data))
})
names(pdataframe)[1] = "method"

p5 = ggplot(pdataframe, aes(Axis_1, Axis_2, color=is_case, shape=vag_fluid, fill="is_case", label=method_manuscript_id))
p5 = p5 + geom_point(size=2) +geom_text(hjust=-.3, vjust=-.8, size=1.5)
p5 = p5 + facet_wrap(~method, scales="free")
p5 = p5 + scale_fill_brewer(type="qual", palette="Set1")
p5 = p5 + scale_colour_brewer(type="qual", palette="Set1")
p5

######
library(vegan)
plot_bar(phy.organism.rel, fill="Species")+
  theme(legend.position = "none")+ ggtitle(paste("Plot of Relative Abundances, Not Rarefied"))

library(vegan)
plot_bar(phy.organism.rel, fill="Species")+
  facet_wrap(~amsel, scales="free_x", nrow=1)+
  theme(legend.position = "none")+ ggtitle(paste("Plot of Relative Abundances, Not Rarefied"))


phy.organism.rarefied <- rarefy_even_depth(phy.organism, rngseed=1, sample.size=0.9*min(sample_sums(phy.organism)), replace=F)
plot_bar(phy.organism.rarefied, fill="Species")+
  theme(legend.position = "none")+ ggtitle(paste("Plot of Rarefied Abundances, set to 11628 representing 90% of the min read depth"))

phy.organism.rarefied <- rarefy_even_depth(phy.organism, rngseed=1, sample.size=0.9*min(sample_sums(phy.organism)), replace=F)
plot_bar(phy.organism.rarefied, fill="Species")+
  facet_wrap(~amsel, scales="free_x", nrow=1)+
  theme(legend.position = "none")+ ggtitle(paste("Plot of Rarefied Abundances, set to 11628 representing 90% of the min read depth"))


raredf <- as.data.frame(otu_table(phy.organism))

raretab <- as.matrix(t(raredf))

rarecurve(raretab, step=50, cex=0.5, label = FALSE)
ordilabel(cbind(rowSums(raretab), specnumber(raretab)), labels=samplet$method_manuscript_id)



################################
#SHOTGUN SEQUENCING DATA
#################################

rm(list=ls())

library('ggplot2')
library("vegan")
library(RColorBrewer)
library(phyloseq)
library("plyr")

scaleyellowred <- colorRampPalette(c("lightyellow", "red"), space = "rgb")(100)

samplet <- read.csv("bvr01_proteomics_samplesheet.csv", header=TRUE, row.names =1, sep=",")
otut_shotgun <- read.csv("bvr01_proteomics_abs_shotgun.csv", header=TRUE,row.names=1, sep=",")
taxat_shotgun <- read.csv("bvr01_proteomics_taxatable_shotgun.csv", header=TRUE, row.names =1, sep=",")

OTU = otu_table(as.matrix(otut_shotgun), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(taxat_shotgun))
SAMPLE <- sample_data(samplet)

physeq = phyloseq(OTU, TAX)
physeq
physeq1 = merge_phyloseq(physeq, SAMPLE)

##########################
#Pruning taxa > 0
###########################
physeq2 <- prune_taxa( taxa_sums(physeq1) > 0, physeq1)
physeq2

######################

rank_name <- 'Species'
phy.organism <- phyloseq::tax_glom(physeq2, rank_name)

#Creating relative abundance
phy.organism.rel <- transform_sample_counts(phy.organism, function(OTU) OTU/sum(OTU))

#################
# Alpha diversity plots
#################
p = plot_richness(phy.organism, x="subject_id", color = "is_case", measures=c("Chao1", "Shannon"))
p + geom_point(size=6, alpha=0.2)+ geom_text(mapping = aes(label = subject_id), size = 2, color="black")

##################
# Alpha diversity exports
##################

phy.organism.chao1 = estimate_richness(phy.organism, measures=c("Chao1"))
phy.organism.shannon = estimate_richness(phy.organism, measures=c("Shannon"))

phy.organism.chao1$id <- row.names(phy.organism.chao1)
phy.organism.shannon$id <- row.names(phy.organism.shannon)

phy.organism.shannon$id2 = gsub("-", ".", phy.organism.shannon$id)
phy.organism.chao1$id2 = phy.organism.chao1$id

phy.organism.alphad.merge <- merge(phy.organism.chao1,phy.organism.shannon, by = "id2")

write.csv(phy.organism.alphad.merge, file = "phy_organism_alphad_merge_shotgun.csv")

##############################################

##################
#Beta diversity - MDS with Bray distance
##################

phy.subset.dist <- phy.organism.rel
braydist <- distance(phy.subset.dist, method="bray")
braymds  <- ordinate(phy.subset.dist, "MDS", distance=braydist)
pdist1 <- NULL
pdist1 <- plot_ordination(phy.subset.dist, braymds, color="is_case")
pdist1 <- pdist1 + ggtitle(paste("MDS using bray distance"))
pdist1 + geom_point(size = 7)+ geom_text(mapping = aes(label = method_manuscript_id), size = 2, color="black")

##################
#Beta diversity - Multiple methods
##################

#####
# Relative abundance
#####
phy.organism.ord.rel <- ordinate(phy.organism.rel, "NMDS", "euclidean")
p4.rel = plot_ordination(phy.organism.rel, phy.organism.ord.rel, type="samples", color="is_case", shape="vag_fluid")
p4.rel + geom_point(size=7) + geom_text(mapping = aes(label = subject_id), size = 2, color="black")+ ggtitle("samples")
######

dist = "bray"
ord_meths = c("DCA", "CCA", "RDA", "NMDS", "MDS", "PCoA")
plist = llply(as.list(ord_meths), function(i, phy.organism, dist){
  ordi = ordinate(phy.organism, method=i, distance=dist)
  plot_ordination(phy.organism, ordi, "samples", color="is_case")
}, phy.organism, dist)

names(plist) <- ord_meths

pdataframe = ldply(plist, function(x){
  df = x$data[, 1:2]
  colnames(df) = c("Axis_1", "Axis_2")
  return(cbind(df, x$data))
})
names(pdataframe)[1] = "method"

p5 = ggplot(pdataframe, aes(Axis_1, Axis_2, color=is_case, shape=vag_fluid, fill="is_case", label=method_manuscript_id))
p5 = p5 + geom_point(size=2) +geom_text(hjust=-.3, vjust=-.8, size=1.5)
p5 = p5 + facet_wrap(~method, scales="free")
p5 = p5 + scale_fill_brewer(type="qual", palette="Set1")
p5 = p5 + scale_colour_brewer(type="qual", palette="Set1")
p5

######
library(vegan)
plot_bar(phy.organism.rel, fill="Species")+
  theme(legend.position = "none")+ ggtitle(paste("Plot of Relative Abundances, Not Rarefied"))

library(vegan)
plot_bar(phy.organism.rel, fill="Species")+
  facet_wrap(~amsel, scales="free_x", nrow=1)+
  theme(legend.position = "none")+ ggtitle(paste("Plot of Relative Abundances, Not Rarefied"))


phy.organism.rarefied <- rarefy_even_depth(phy.organism, rngseed=1, sample.size=0.9*min(sample_sums(phy.organism)), replace=F)
plot_bar(phy.organism.rarefied, fill="Species")+
  theme(legend.position = "none")+ ggtitle(paste("Plot of Rarefied Abundances, set to 11628 representing 90% of the min read depth"))

phy.organism.rarefied <- rarefy_even_depth(phy.organism, rngseed=1, sample.size=0.9*min(sample_sums(phy.organism)), replace=F)
plot_bar(phy.organism.rarefied, fill="Species")+
  facet_wrap(~amsel, scales="free_x", nrow=1)+
  theme(legend.position = "none")+ ggtitle(paste("Plot of Rarefied Abundances, set to 11628 representing 90% of the min read depth"))


raredf <- as.data.frame(otu_table(phy.organism))

raretab <- as.matrix(t(raredf))

rarecurve(raretab, step=1000, cex=0.5, label = FALSE)
ordilabel(cbind(rowSums(raretab), specnumber(raretab)), labels=samplet$method_manuscript_id)
