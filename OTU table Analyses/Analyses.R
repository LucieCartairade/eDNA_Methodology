#### Environment #### 

setwd("C:/Users/lcartair/Documents/eDNA_FromWaterToBiodiversity/2209_Papetoai_Manip_PorositeVolume/")
Datas_path <- "C:/Users/lcartair/Documents/eDNA_FromWaterToBiodiversity/2209_Papetoai_Manip_PorositeVolume/Res_Decona/"
Images_path <- "C:/Users/lcartair/Documents/eDNA_FromWaterToBiodiversity/2209_Papetoai_Manip_PorositeVolume/ResultsR/"

source("../Functions.R")

#library always needed : 
library(ggplot2)
theme_set(theme_bw())
library(dplyr)

#### Reading Files #### 

##### Datas #####
metadatas <- read.csv(file = "Datas_Pap.csv", header = T, sep = ";", dec = ",", na.strings = "NA", fileEncoding = "ISO-8859-1")
metadatas <- metadatas[which(metadatas$Barcod != ""),]
metadatas$Barcod <- stringr::str_sub(metadatas$Barcod,3,4)
#metadatas$Sample.ID <- row.names(metadatas)
# match sample names
metadatas$Run_Barcod <- paste0(metadatas$Run.name,"_barcode",metadatas$Barcod,"_concatenated")

metadatas <- metadatas[metadatas$Sample.Type == "Porosity" | metadatas$Sample.Type == "Control",]

###############################################################################b
#####                    Results From Decona                               #####
###############################################################################b

Res <- read.table(file='Res_Decona/BLAST_out_reclustered_summary_tax_seq_counts_250128.txt',sep ="\t", header = T, na.string = "")
names(Res)[14:dim(Res)[2]] <- gsub("\\.", "-", names(Res)[14:dim(Res)[2]])

# Retirer la contamination Homo Sapiens 
Res <- Res[-which(Res$tax.id == "9606"),]

# Normaliser
colnames(Res)
colnames(Res[,c(14:dim(Res)[2])])
Res_N <- t(t(Res[,c(14:dim(Res)[2])])/colSums(Res[,c(14:dim(Res)[2])]) * 100)
Res_N <- cbind(Res[,c(1:13)], Res_N)


# Melting function 
Res_melt <- Melting(Res_N, c("Run_Barcod","Sample.ID","Sample.Type","Filtration.volume", "Size.fraction","Replica"))
Res_melt <- Melting(Res, c("Run_Barcod","Sample.ID","Sample.Type","Filtration.volume", "Size.fraction","Replica"))

unique(Res_melt[is.na(Res_melt$Sample.ID),"Run_Barcod"])
Res_melt <- Res_melt[!is.na(Res_melt$Sample.ID),]
unique(Res_melt[is.na(Res_melt$Nb.reads),])

# Look at reads an run quality
dim(Res_melt)
# Sur des données normalisées : 
dim(Res_melt[Res_melt$Nb.reads <1,]) # 
dim(Res_melt[Res_melt$Nb.reads <1,]) / dim(Res_melt) # 74 % des OTUs représentent moins de 1% des reads de chaque échantillon 
#Res_melt[Res_melt$Nb.reads >= 1,]

sum(Res_melt$Nb.reads) #599714
sum(subset(Res_melt,Res_melt$Family=="unknown")$Nb.reads) #33520
(sum(subset(Res_melt,Res_melt$Family=="unknown")$Nb.reads)/sum(Res_melt$Nb.reads))*100 # 5.59% de unknown

length(unique(Res_melt$Taxon)) # 171
mean(Res_melt$X.ID, na.rm = T)
mean(Res_melt$bit.score, na.rm = T)

sum(Res_melt[ Res_melt$X.ID > 97, "Nb.reads"], na.rm = T) / sum(Res_melt$Nb.reads)# 79% des reads ont une assignation taxonomique à plus de 97% d'identité
Res_melt[ Res_melt$X.ID > 97,]

#Metric(Res)
#Violinplot(Res_melt)

##### Supprimer les assignations dont la longeur d'alignement est inf?rieure ? 160 #### 
pdf(paste0(Images_path,"Before-After_filter-on-align-length.pdf"))
Boxplot(Res_melt)
Res_melt[which(Res_melt$alignment.length < 160),c("Family", "Genus", "Species", "Taxon")] <- "unknown"
Boxplot(Res_melt)
dev.off()
dev.off()

# Look at reads an run quality
sum(Res_melt$Nb.reads) #16294993
sum(subset(Res_melt,Res_melt$Family=="unknown")$Nb.reads) #469846
(sum(subset(Res_melt,Res_melt$Family=="unknown")$Nb.reads)/sum(Res_melt$Nb.reads))*100 # 2.88% de unknown




###############################################################################b
####                                Tab                                    #####
###############################################################################b

Res_melt$Taxon <- paste(Res_melt$Family, Res_melt$Genus,Res_melt$Species, sep = "_")

Tab_OTU <- reshape2::acast(Res_melt, value.var = "Nb.reads",Res_melt$clusters.id~Res_melt$Sample.ID)
Tab_Taxon <- reshape2::acast(Res_melt, value.var = "Nb.reads", Res_melt$Taxon~Res_melt$Sample.ID, fun.aggregate = sum)
colnames(Tab_Taxon)

sum(Tab_OTU,na.rm = T)
sum(Tab_Taxon,na.rm = T)

Tab_01 <- ifelse(Tab_Taxon == 0,FALSE,TRUE)
colSums(Tab_01,na.rm = T)
Tab_01 <- ifelse(Tab_OTU == 0,FALSE,TRUE)
colSums(Tab_01,na.rm = T)



###############################################################################b
####                      Rarefaction curves                               #####
###############################################################################b

# Il ne faut pas utiliser les données normalisées pour les courbes de rarefaction !!!!!!

Tab_raw <- Tab_Taxon
#Tab_raw <- Tab_OTU
Tab_raw[is.na(Tab_raw)] <- 0 

#total number of species at each site (row of data)
S <- vegan::specnumber(t(Tab_raw))

# Number of INDIVIDULS per site (?)
raremax <- min(rowSums(t(Tab_raw), na.rm = T)) 

# rarefy, w/ raremax as input (?)
Srare <- vegan::rarefy(t(Tab_raw), raremax)
Tab_raw

#Plot rarefaction results
par(mfrow = c(1,2))
plot(S, Srare, xlab = "Observed No. of Species", 
     ylab = "Rarefied No. of Species",
     main = " plot(rarefy(Tab, raremax))")
abline(0, 1)
#vegan::rarecurve(t(Tab_raw), step = 20, 
#                 sample = raremax, 
#                 col = "blue", 
#                 cex = 0.6,
#                 main = "rarecurve()")
dev.off()


pdf(paste0(Images_path,"rarefaction.pdf"), width = 9, height = 6)
vegan::rarecurve(t(Tab_raw), step = 20, sample = raremax, col = "blue", cex = 0.6,
                 main = "rarecurve() on subset of data")
dev.off()

###############################################################################b
####               Correlations between replicates                        ######
###############################################################################b

#Pearson Correlation (linear correlation)
corP <- cor(as.matrix(Tab_Taxon), method = "pearson")
corP 
corrplot::corrplot(corP, type="upper", tl.col="black", tl.srt=45)
colMeans(corP)

#Spearman Correlation (non linear correlation)
corS <- cor(as.matrix(Tab_Taxon), method = "spearman")
corS
corrplot::corrplot(corS, type="upper", tl.col="black", tl.srt=45)

#### Hierarchical Ascending Classification ####
HA(Tab_Taxon, file_path = paste0(Images_path,"/Dendrogramme.svg"), line = 5.7)

###############################################################################b
####                           Heat map                                   ######
###############################################################################b

HeatMap <- function()
{
  Res_melt$Taxon <- factor(Res_melt$Taxon, levels = rev(unique(c("unknown_unknown_unknown", as.vector(Res_melt$Taxon[order(Res_melt$Nb.reads, decreasing = T)])))))
  ggplot(as.data.frame(Res_melt)) +
    geom_tile(aes(x=Sample.ID , y= Taxon, fill = log(Nb.reads))) +
    xlab("20L Replicates") + 
    facet_grid(~Size.fraction, scales = "free", space = "free") +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    scale_fill_gradient(low = "#55aff4", high = "#203850", na.value = NA)
  #write.table(Tab, "Tab.csv", row.names = T, sep = ";") #Pour Gilles
}

HeatMap()
ggsave(path = Images_path, filename = "HeatMap.svg", width = 15, height = 25) #, device='tiff', dpi=700)

###############################################################################b
####                            Physolseq                                  #####
###############################################################################b

##### Phyloseq on OTUs #### 

colnames(Res)
colnames(Res[,c(20:dim(Res)[2])])
metadatas[,c("Sample.ID","Run_Barcod")]
Res[is.na(Res$X.ID),"alignment.length"] <- 0 
Res[Res$alignment.length < 160, c("Family","Genus","Species")] <- "unknown"

physeq <- Phyloseq(Res, rarefied = F, sampleDBool = T)

# Quality check
p <- phyloseq::plot_bar(physeq, fill = "Family", x = "Sample") + scale_color_manual(na.value = "grey50") + scale_fill_manual(values= Palette(physeq,"Family"))
order <- rev(unique(c("unknown", as.vector(p$data$Family[order(p$data$Abundance, decreasing = T)]))))
p$data$Family <- factor(p$data$Family, levels = order)
#p$data$Sample <- factor(p$data$Sample, levels = colnames(Tab))
p$data$Sample <- paste(p$data$Sample.Type, p$data$Replica)
p$data[which(p$data$OTUs == "unknown unknown"),c("OTUs", "Family", "Genus","Species")] <- NA
p + scale_x_discrete(guide = guide_axis(angle = 45))
ggsave(path = Images_path, filename = "Barplot_Phyloseq_Samples.svg", width = 15, height = 8)

# Rarefy
physeq = phyloseq::rarefy_even_depth(physeq, rngseed=1, sample.size=0.9*min(phyloseq::sample_sums(physeq)), replace=F)

# Plot Phyloseq
p <- phyloseq::plot_bar(physeq, fill = "Family", x = "Replica") + scale_color_manual(na.value = "grey50") + scale_fill_manual(values= Palette(physeq,"Family"))
p$data$Family <- factor(p$data$Family, levels = order) # orderring by family
p$data[which(p$data$OTUs == "unknown unknown"),c("OTU", "Family", "Genus","Species")] <- NA # set unknown to na.value
p <- p + facet_grid(~Sampling.Site, scales = "free", space = "free_x")
p <- p + xlab("20L Replicates") + guides(fill=guide_legend(ncol=2))
p <- p + theme(text=element_text(size = 20)) 
print(p)
ggsave(path = Images_path, filename = "Barplot_Phyloseq_Rarefied_OTUs.svg", width = 20, height = 10)

#p <- plot_composition(physeq, "tax_level_selection", "tax_name_selection", "tax_level_aggregation", numberOfTaxa = NB_OTU, fill = "tax_color")
microbiome::plot_composition(physeq, "Family", "Pomacentridae", "Genus", numberOfTaxa = 5, fill = "Genus")
phyloseq::tax_table(physeq)

##### More abondant taxa barplot #### 

unique(phyloseq::tax_table(physeq))
top_nested <- fantaxtic::nested_top_taxa(physeq,
                                         top_tax_level = "Family",
                                         nested_tax_level = "OTUs",
                                         n_top_taxa = 6, 
                                         n_nested_taxa = 3, 
                                         include_na_taxa = T)
# Enlever le groupe "Other"
phyloseq::tax_table(top_nested$ps_obj) <- phyloseq::tax_table(top_nested$ps_obj)[-which(phyloseq::tax_table(top_nested$ps_obj)[,"Family"] == "Other"),]
phyloseq::tax_table(top_nested$ps_obj) <- phyloseq::tax_table(top_nested$ps_obj)[-which(stringr::str_split(phyloseq::tax_table(top_nested$ps_obj)[,"Genus"], " " ,simplify = T)[,1] == "Other"),]

plot_nested_bar_Lucie(ps_obj = top_nested$ps_obj, top_level = "Family", nested_level = "OTUs", x_value= "Replica",
                      palette = c(unknown = "gray50"),
                      #na_taxon_label = "unknown unknown",
                      merged_clr = "black",
                      legend_title = "Species") +
  facet_grid(~Sampling.Site, scales = "free_x", space = "free_x")
ggsave(path = Images_path, filename = "Barplot_Phyloseq_Nested_top6_without_Other.svg", width = 15, height = 8)

top_taxa <- top_nested$top_taxa


###############################################################################b
####                           Phyloseq Tree                                 #####
###############################################################################b
#https://joey711.github.io/phyloseq/merge.html

t <- Tree(ps = physeq)
t <- Tree(ps = physeq, color = "Loc", shape = "NM")#, label = "taxa_names")

library(ape)
random_tree = ape::rtree(ntaxa(physeq), rooted = T, tip.label = taxa_names(physeq))

ps_tree = merge_phyloseq(physeq, random_tree)
sample_data(ps_tree)
plot_tree(ps_tree, color ="Habitat", label.tips= "taxa_names", ladderize="left")#, text.size = 5)



print(t)
rank_names(physeq)
tax_table(physeq)


ggsave(t, path = Images_path, filename = "Tree.svg", width = 25, height = 45)


###############################################################################b
####                          Alpha diversity                                #####
###############################################################################b

a <- phyloseq::plot_richness(physeq, x = "Sampling.Site", color = "Sampling.Site",  measures=c("Observed","Shannon")) + 
  geom_boxplot(alpha = 0.2)
#a <- a + geom_boxplot(aes(fill = Sampling.Site))
#a$data$Rep <- factor(a$data$Rep, levels = c(1,3,5,10))
print(a)
ggsave(path = Images_path, filename = "AlphaDiv.svg", width = 7, height = 5)

rich = phyloseq::estimate_richness(physeq, measures = c("Observed", "Chao1", "Shannon", "InvSimpson"))

data <- cbind(phyloseq::sample_data(physeq), rich)

physeq_anova <- aov(Observed ~ Replica, data) 
physeq_anova <- aov(Observed ~ Sampling.Site, subset(data,Sampling.Site != "Samoa"))
summary(physeq_anova) # il y a un effet siginificatif du nombre de rep sur le nombre d'OTUs observé 
physeq_anova <- aov(Shannon ~ Rep, data) 
summary(physeq_anova) # il n'y a pas d'effet sur l'uniformité 
physeq_anova <- aov(InvSimpson ~ Rep, data) 
summary(physeq_anova) # ni sur la probalité de tirer deux fois la même espèce. 



###############################################################################b
####                                PCA                                      #####
###############################################################################b


cor <- cor(Res[,c(20:25,28:39)], method = "pearson")
cor <- cor(Tab_Taxon, method = "pearson")
cor
corrplot::corrplot(cor)

ord_meths = c("DCA", "CCA", "RDA", "DPCoA", "NMDS", "MDS", "PCoA")

dist.bc <- phyloseq::distance(physeq, method = "bray")
ord <- phyloseq::ordinate(physeq, method = "MDS", distance = dist.bc)


#p <- phyloseq.extended::plot_dist_as_heatmap(dist.bc, title = "Bray-Curtis")
#plot(p)


p1 = phyloseq::plot_ordination(physeq, ord, type="taxa", color="Family", title="taxa")
p1
p1 + facet_wrap(~Family, 3)

p2 = phyloseq::plot_ordination(physeq, ord, type="samples", color="Sampling.Site") 
p2 + geom_polygon(aes(fill=Sampling.Site)) + geom_point() #+ ggtitle("samples")

ggsave(path = Images_path, filename = "BrayCurtis_polygon.svg", width = 6, height = 5)



#install.packages(c("FactoMineR", "factoextra"))

#library("FactoMineR")
#library("factoextra")

res.PCA <- FactoMineR::PCA(Tab_Taxon, scale.unit = T)
FactoMineR::plot.PCA(res.PCA, choix = "var")
FactoMineR::plot.PCA(res.PCA, choix = "ind")



###############################################################################b
####                            Upset Taxon                                  #####
###############################################################################b
# https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html
# https://jokergoo.github.io/ComplexHeatmap/reference/make_comb_mat.html


Tab_Taxon_Site[paste(top_taxa$Family, top_taxa$Genus, top_taxa$Species, sep = "_"),]

Tab_melt_Upset <- Res_melt
#Tab_melt_Upset <- subset(Res_melt, Taxon == paste(top_taxa$Family, top_taxa$Genus, top_taxa$Species, sep = "_")) # c'est nul
unique(Tab_melt_Upset$Sampling.Site)

lt <- list( Motu_Ahi = subset(Tab_melt_Upset,Sampling.Site == "MotuAhi")$Taxon,
            Afareaitu = subset(Tab_melt_Upset,Sampling.Site == "Afareaitu")$Taxon,
            Pihaena = subset(Tab_melt_Upset,Sampling.Site == "Pihaena")$Taxon,
            Entre_2_baies = subset(Tab_melt_Upset,Sampling.Site == "Entre2Baies")$Taxon,
            Taotaha = subset(Tab_melt_Upset,Sampling.Site == "Taotaha")$Taxon,
            Haapiti = subset(Tab_melt_Upset,Sampling.Site == "Haapiti")$Taxon)

ComplexHeatmap::list_to_matrix(lt)
m <- ComplexHeatmap::make_comb_mat(lt)
#ComplexHeatmap::comb_name(m)
#t(m)
m2 <- m
#m2 <- m[ComplexHeatmap::comb_size(m) > 6]
upsetPlot <- ComplexHeatmap::UpSet(m2,comb_order = order(ComplexHeatmap::comb_size(m2), decreasing = T), 
                                   top_annotation = ComplexHeatmap::upset_top_annotation(m2, add_numbers = TRUE),
                                   right_annotation = ComplexHeatmap::upset_right_annotation(m2, add_numbers = TRUE))
upsetPlot
ggsave(path = Images_path, ggplotify::as.ggplot(upsetPlot), file = "UpSetR_Taxon.svg", width = 9, height = 3) 


mat <- ComplexHeatmap::list_to_matrix(lt)
#Cosmopolites
mat[which(mat[,"Hiva_Oa"] == 1 & mat[,"Samoa"] == 1 & mat[,"Tonga"] == 1  & mat[,"Rarotonga"] == 1 & mat[,"Tikehau"] == 1 & mat[,"Christmas_Island"] == 1), ]
# Cosmopolites mais pas Samoa
mat[which(mat[,"Hiva_Oa"] == 1 & mat[,"Samoa"] == 0 & mat[,"Tonga"] == 1  & mat[,"Rarotonga"] == 1 & mat[,"Tikehau"] == 1 & mat[,"Christmas_Island"] == 1), ]
# Specific to Tonga 
mat[which(mat[,"Hiva_Oa"] == 0 & mat[,"Samoa"] == 0 & mat[,"Tonga"] == 1  & mat[,"Rarotonga"] == 0 & mat[,"Tikehau"] == 0 & mat[,"Christmas_Island"] == 0), ]
# Specific to Rarotonga 
mat[which(mat[,"Hiva_Oa"] == 0 & mat[,"Samoa"] == 0 & mat[,"Tonga"] == 0  & mat[,"Rarotonga"] == 1 & mat[,"Tikehau"] == 0 & mat[,"Christmas_Island"] == 0), ]
# Specific to Christmas Island
mat[which(mat[,"Hiva_Oa"] == 0 & mat[,"Samoa"] == 0 & mat[,"Tonga"] == 0  & mat[,"Rarotonga"] == 0 & mat[,"Tikehau"] == 0 & mat[,"Christmas_Island"] == 1), ]



###############################################################################b
####                           Upset by OTUs                                 #####
###############################################################################b

lt <- list( Hiva_Oa = subset(Res_melt,Sampling.Site == "Hiva Oa")$clusters.id,
            Christmas_Island = subset(Res_melt,Sampling.Site == "Christmas Island")$clusters.id,
            Tikehau = subset(Res_melt,Sampling.Site == "Tikehau")$clusters.id,
            Rarotonga = subset(Res_melt,Sampling.Site == "Rarotonga")$clusters.id,
            Tonga = subset(Res_melt,Sampling.Site == "Tonga")$clusters.id,
            Samoa = subset(Res_melt,Sampling.Site == "Samoa")$clusters.id)

ComplexHeatmap::list_to_matrix(lt)
m <- ComplexHeatmap::make_comb_mat(lt)
#ComplexHeatmap::comb_name(m)
#t(m)
m2 <- m[ComplexHeatmap::comb_size(m) > 6]
upsetPlot <- ComplexHeatmap::UpSet(m2, comb_order = order(ComplexHeatmap::comb_size(m2), decreasing = T), 
                                   top_annotation = ComplexHeatmap::upset_top_annotation(m2, add_numbers = TRUE),
                                   right_annotation = ComplexHeatmap::upset_right_annotation(m2, add_numbers = TRUE))
upsetPlot
ggsave(path = Images_path, ggplotify::as.ggplot(upsetPlot), file = "UpSetR_OTUs.svg", width = 9, height = 4) 

###############################################################################b
####                            Venn Euler                                 #####
###############################################################################b

# Considérer le taxon comme présent dans le site quand il comporte plus de 100 reads
Tab_Euler <- ifelse(Tab_Taxon_Site < 100,FALSE,TRUE)
Tab_Euler <- ifelse(Tab_Taxon_Site == 0,FALSE,TRUE)

Tab_Taxon_Site[paste(top_taxa$Family, top_taxa$Genus, top_taxa$Species, sep = "_"),]

set.seed(10)
pdf(file = paste0(Images_path,"Euler_plot.pdf"))
#plot(eulerr::euler(ifelse(is.na(Tab_OTU),FALSE,TRUE), shape = "ellipse"), quantities = TRUE) # ne tourne pas 
plot(eulerr::euler(Tab_Euler, shape = "ellipse"), quantities = TRUE)
plot(eulerr::venn(Tab_Euler[,1:5]))
plot(eulerr::venn(Tab_Euler[,2:6]))
dev.off()
