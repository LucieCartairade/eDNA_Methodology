# Commons scripts to all figures
## Loading the working environment
```r
source("Functions.R")

library(ggplot2)
theme_set(theme_bw())
library(dplyr)
```
## Reading metadatas file 
```r
metadatas <- read.csv(file = "Datas.csv", header = T, sep = ";", dec = ",", na.strings = "NA", fileEncoding = "ISO-8859-1")
metadatas <- metadatas[which(!is.na(metadatas$Barcod)),]
metadatas$Barcod <- sprintf("%02d",metadatas$Barcod)
# match sample names with OTU table
metadatas$Run_Barcod <- paste0(metadatas$Run.name,"_barcode",metadatas$Barcod,"_concatenated")
```

## Reading OTU table file 
```r
Res <- read.table(file='BLAST_out_reclustered_summary_tax_seq_counts.txt',sep ="\t", header = T, na.string = "")

# Remove Homo sapiens reads 
Res <- Res[-which(Res$tax.id == "9606"),]

Res_melt <- Melting(Res, metadatas_selected_col = c("Run_Barcod","Sample.ID","Replica"))

100 * sum(subset(Res_melt,Res_melt$Family=="unknown")$Nb.reads)/sum(Res_melt$Nb.reads) # 2.29% of unknown

Tax_melt <- Res_melt %>%
  group_by(tax.id, Sample.ID, Replica, Family) %>%
  summarise(
    Nb.reads_sum = sum(Nb.reads),
    X.ID_mean = weighted.mean(X.ID, Nb.reads),
    alignment.length_mean = weighted.mean(alignment.length, Nb.reads),
    mismatches_mean = weighted.mean(mismatches, Nb.reads),
    gap.opens_mean = weighted.mean(gap.opens, Nb.reads),
    evalue_mean = weighted.mean(evalue, Nb.reads),
    bit.score_mean = weighted.mean(bit.score, Nb.reads),
    qcovs_mean = weighted.mean(qcovs, Nb.reads),
    sequence_max = sequence[which.max(Nb.reads)][1]
  ) %>%
  group_by(Sample.ID) %>%
  mutate(
    relative_biomass = 100 * Nb.reads_sum / sum(Nb.reads_sum)
  ) %>%
  ungroup()

Tax_melt[which(Tax_melt$bit.score_mean < 200),"Taxon"] <- "unknown unknown"
Tax_melt[which(Tax_melt$bit.score_mean < 200),"Family"] <- "unknown"
Tax_melt[which(Tax_melt$bit.score_mean < 200),"X.ID_mean"] <- NA

100 * sum(subset(Tax_melt,Tax_melt$Family=="unknown")$Nb.reads_sum)/sum(Tax_melt$Nb.reads_sum) # 3.77% of unknown

Tax_table <- reshape2::acast(Tax_melt, value.var = "Nb.reads_sum", Taxon~Sample.ID, fill = 0, fun.aggregate = sum)
```

## Creating a phyloseq object
```r
# OTUs object
OTUs <- phyloseq::otu_table(as.data.frame(Tax_table), taxa_are_rows = T)

# TAX object
TAX <- data.frame(unique(Tax_melt[c("Family", "Taxon")]), row.names = unique(Tax_melt$Taxon))
names(TAX) <-  c("Family", "Species")
TAX <- phyloseq::tax_table(as.matrix(TAX))

# SAMPLE object
sample <- data.frame(metadatas, row.names = metadatas$Sample.ID)
SAMPLE <- phyloseq::sample_data(sample)

physeq <- phyloseq::phyloseq(OTUs, TAX, SAMPLE)

physeq <- phyloseq::subset_samples(physeq, Sample.Type != "Control")

# Rarefy
physeq <- phyloseq::rarefy_even_depth(physeq, rngseed=1, sample.size=0.9*min(phyloseq::sample_sums(physeq)), replace=F)
```

# Figure 1 : Porosity - Alpha Diversity 
![image](Figures/Figure1.png)
```r
rich <- phyloseq::estimate_richness(physeq, measures = c("Observed", "Chao1", "Shannon", "InvSimpson"))
data <- merge(as.data.frame(phyloseq::sample_data(physeq)), rich, by.x = "row.names", by.y = "row.names")
data$Group <- ifelse(data$Size.fraction %in% c("0.2-0.8", "0.2-1.2", "0.2-3"), "Group1", "Group2")
data[,c("Size.fraction","Group","Observed","Shannon")]

# Graphique A : Observed
p1 <- ggplot(subset(data, Sample.Type != "Control"), aes(x = Size.fraction, y = Observed, fill = Group)) +
  stat_summary(fun = mean, geom = "bar", color = "black", width = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  labs(y = "Observed richness", x = "Size fraction", title = "Observed") + 
  scale_fill_manual(values = c("Group1" = "gray80", "Group2" = "gray60")) +
  scale_x_discrete(guide = guide_axis(angle = 45)) + theme(legend.position = "none")

# Graphique B : Shannon
p2 <- ggplot(subset(data, Sample.Type != "Control"), aes(x = Size.fraction, y = Shannon, fill = Group)) +
  stat_summary(fun = mean, geom = "bar", color = "black", width = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  labs(y = "Shannon index", x = "Size fraction", title = "Shannon") + 
  scale_fill_manual(values = c("Group1" = "gray80", "Group2" = "gray60")) +
  scale_x_discrete(guide = guide_axis(angle = 45)) + theme(legend.position = "none")

(p1 | p2) + plot_annotation(tag_levels = 'A')
ggsave(path = Images_path, filename = "AlphaDiv.svg", width = 6, height = 4)
```

# Figure 2 : Barplot - Robot vs Tripod
![image](Figures/Figure2.png)
```r
top_nested <- fantaxtic::nested_top_taxa(physeq, top_tax_level = "Family", nested_tax_level = "Species", n_top_taxa = 7, n_nested_taxa = 8, include_na_taxa = T)

plot_nested_bar_Lucie(ps_obj = top_nested$ps_obj, top_level = "Family", nested_level = "Species", x_value= "Rep",
                      palette = c(unknown = "gray50"), merged_clr = "black", legend_title = "Species") +
  facet_grid(~Method, scales = "free_x", space = "free_x") + 
  xlab("20L Replica")

ggsave(path = Images_path, filename = "Barplot_Phyloseq_Nested_top8.svg", width = 12, height = 8)
```

# Figure 3 : Euler plot - Robot vs Tripod
![image](Figures/Figure3.png)
```r
Tab_Euler <- Tax_melt_wV
Tab_Euler <- aggregate(as.numeric(Tab_Euler$relative_abundance), by= list(Method = Tab_Euler$Method, Species = Tab_Euler$Taxon), mean)
Tab_Euler <- acast(Tab_Euler, value.var = "x", Tab_Euler$Species~Tab_Euler$Method, fill = 0)
Tab_Euler <- as.data.frame(Tab_Euler) %>% 
  select("Robot", "Tripode", "Visual Census")
colnames(Tab_Euler) <- c("Robot", "Tripod", "Visual Census")

Tab_Euler_final <- ifelse(Tab_Euler == 0,FALSE,TRUE)

pdf(file = paste0(Images_path,"Euler_plot.pdf"), width = 5, height = 5)
set.seed(19980821)
plot(eulerr::euler(Tab_Euler_final, shape = "ellipse"), fills = c("#4A90E2", "#F5A623", "#50E3C2"), quantities = TRUE, alpha = 0.5)
dev.off()
```

# Figure 4 : Barplot for each taxon - Robot vs Tripode
![image](Figures/Figure4.png)
```r
df <- subset(Tax_melt_wV, Replicas.Nb != "Control" & Replicas.Nb != 0) %>%
  group_by(Sample.ID) %>%
  mutate(total_abundance = sum(Nb.reads_sum)) %>%  # Compute total abundance for each sample
  ungroup() %>%
  mutate(relative_abundance = (Nb.reads_sum / total_abundance) * 100) %>%  # Compute relative abundance
  # Complete the data for each combination of Taxon and Method, adding zeros where data is missing
  tidyr::complete(Taxon, Method, fill = list(relative_abundance = 0))  %>% # Ensure missing combinations have relative_abundance set to 0
  # Sort by relative_abundance within each method
  arrange(Method, relative_abundance) %>%
  # Force the order of Taxa based on relative abundance in each method
  mutate(Taxon = factor(Taxon, levels = unique(Taxon)))  %>% # Update factor levels for Taxon
  mutate(Method = recode(Method, "Tripode" = "Tripod"))

# Compute the maximum relative abundance per taxon
taxon_order <- df %>%
  group_by(Taxon) %>%
  summarise(max_abundance = max(relative_abundance, na.rm = TRUE)) %>%
  arrange(max_abundance) %>%
  pull(Taxon)

# Update Taxon factor levels based on maximum abundance
df_ordered <- df %>%
  mutate(Taxon = factor(Taxon, levels = taxon_order))

sum(df$relative_abundance, na.rm = T)

ggplot(df_ordered, aes(x = Taxon, y = relative_abundance, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge",  width = 0.5) +
  scale_fill_manual(values = c("Robot" = "#4A90E2", "Tripod" = "#F5A623", "Visual Census" = "#50E3C2")) +
  labs(
    title = "Average Relative Abundance of Species by Method",
    x = "Species",
    y = "Average Relative Abundance",
    fill = "Method"
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_flip()  # Flip the axes for better readability

ggsave(path = Images_path, filename = "RelativeAbundance_wVC.svg", width = 12, height = 6)
```

# Figure 5 : Alpha Diversity - Volume
![image](Figures/Figure5.png)
```r
rich = phyloseq::estimate_richness(physeq_wout_Ctrl, measures = c("Observed", "Chao1", "Shannon", "InvSimpson"))
data <- merge(as.data.frame(phyloseq::sample_data(physeq_wout_Ctrl)), rich, by.x = "row.names", by.y = "row.names")
data[,c("Filtration.volume","Observed","Shannon")]

data <- cbind(data,Nb.reads = t(t(colSums(phyloseq::otu_table(physeq_wout_Ctrl)))))

# Graphique A : Observed
p1 <- ggplot(data, aes(x = factor(Filtration.volume), y = Observed, fill = factor(Filtration.volume))) +
  stat_summary(fun = mean, geom = "bar", color = "black", width = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  #geom_jitter(width = 0.1, size = 2, alpha = 0.6) +
  labs(y = "Observed richness", x = "Filtration volume (L)", title = "Observed") + 
  scale_fill_manual(values = c("gray80", "gray70", "gray60")) + theme(legend.position = "none")

# Graphique B : Shannon
p2 <- ggplot(data, aes(x = factor(Filtration.volume), y = Shannon, fill = factor(Filtration.volume))) +
  stat_summary(fun = mean, geom = "bar", color = "black", width = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  #geom_jitter(width = 0.1, size = 2, alpha = 0.6) +
  labs(y = "Shannon index", x = "Filtration volume (L)", title = "Shannon") + 
  scale_fill_manual(values = c("gray80", "gray70", "gray60"))+ theme(legend.position = "none")
p2

# Combine les deux graphiques
(p1 | p2 ) + plot_annotation(tag_levels = 'A')
ggsave(path = Images_path, filename = "AlphaDiv.svg", width = 4, height = 4)
```

# Figure 6 : Accumulation curves - Sampling Replicates
![image](Figures/Figure6.png)
```r
sp <- vegan::specaccum(t(Tax_table), method = "collector")

pdf(paste0(Images_path,"spaceaccum.pdf"),width = 5, height = 5)
plot(sp, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue", 
    main = "Species Accumulation curve", xlab = "Replicas", ylab = "Number of species")
dev.off()
```

# Figure 7 : PCR Replicates
![image](Figures/Figure7.png)
```r
Nb <- data.frame(Nb_OTUs = colSums(Tax_table != 0 ), 
                Origin = stringr::str_split(colnames(Tax_table), "_", simplify= T)[,1],
                Rep = as.numeric(stringr::str_split(colnames(Tax_table), "_", simplify= T)[,2]))

summary_data <- Nb %>%
  group_by(Rep) %>%
  summarise(
    mean_OTUs = mean(Nb_OTUs),
    sd_OTUs = sd(Nb_OTUs)
  )

p1 <- ggplot(summary_data, aes(x = Rep, y = mean_OTUs)) +
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(aes(ymin = mean_OTUs - sd_OTUs, ymax = mean_OTUs + sd_OTUs), width = 0.3) +
  scale_x_continuous(breaks = c(1, 3, 5, 10)) +
  labs(
    x = "Number of PCR replicates",
    y = "Mean number of OTUs (± SD)",
    title = "Observed"
  )

shannon_df <- Tax_melt %>%
  group_by(Sample.ID) %>%
  mutate(p = Nb.reads_sum / sum(Nb.reads_sum)) %>%
  summarise(Shannon = -sum(p * log(p))) %>%
  tidyr::separate(Sample.ID, into = c("Origin", "Rep"), sep = "_") %>%
  mutate(
    Origin = as.factor(Origin),
    Rep = as.numeric(Rep)
  )

summary_data <- shannon_df %>%
  group_by(Rep) %>%
  summarise(
    mean_OTUs = mean(Shannon),
    sd_OTUs = sd(Shannon)
  )

p2 <- ggplot(summary_data, aes(x = Rep, y = mean_OTUs)) +
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(aes(ymin = mean_OTUs - sd_OTUs, ymax = mean_OTUs + sd_OTUs), width = 0.3) +
  scale_x_continuous(breaks = c(1, 3, 5, 10)) +
  labs(
    x = "Number of PCR replicates",
    y = "Mean Shannon index (± SD)",
    title = "Shannon index"
  )

(p1 | p2) + plot_annotation(tag_levels = 'A')
ggsave(path = Images_path, filename = "PCR_replicates.svg", width = 8, height = 4)
```

# Figure 8 : Sequencing depth
![image](Figures/Figure8.png)
```r
Tab_raw <- Tax_table
Tab_raw[is.na(Tab_raw)] <- 0 

#total number of species at each site (row of data)
S <- vegan::specnumber(t(Tab_raw))

# Number of Taxon per sample
raremax <- min(rowSums(t(Tab_raw), na.rm = T)) 

# rarefy, w/ raremax as input
Srare <- vegan::rarefy(t(Tab_raw), raremax)

#Plot rarefaction results
par(mfrow = c(1,2))
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
vegan::rarecurve(t(Tab_raw), step = 20,  sample = raremax,  col = "blue",  cex = 0.6)
dev.off()
```

# Figure 9 : Distance matrix - Tiahura 
![image](Figures/Figure9.png)
```r
dist.jc. <- betapart::beta.pair(t(ifelse(Tax_table != 0, 1, 0)), index.family="jaccard")
dist.bc. <- vegan::vegdist(t(Tax_table), method = "bray")

pheatmap::pheatmap(as.matrix(dist.jc$beta.jac), cluster_rows = F, cluster_cols = F, cellwidth = 10, cellheight = 10, legend = TRUE, main = "Jaccard")
pheatmap::pheatmap(as.matrix(dist.bc), cluster_rows = F, cluster_cols = F, cellwidth = 10, cellheight = 10, legend = TRUE, main = "BrayCurtis")
```

# Figure 10 : PCoA - Tiahura
![image](Figures/Figure10.png)
```r
#Performing PCOA
pcoa <- cmdscale(dist.jc.both$beta.jac, eig = T, add = T)
position <- pcoa$points[,c(1,2)]
colnames(position) <- c("Dim1","Dim2")

# Naming axis
percent_explained <- round(100 * pcoa$eig / sum(pcoa$eig), digits = 1)
percent_explained[1:2]
labs <- c(glue::glue("PCo 1 ({percent_explained[1]}%)"), glue::glue("PCo 2 ({percent_explained[2]}%)"))

#Merging sample_data
pcoa_data <- merge(as.data.frame(position), as.data.frame(Tax_melt_wVC), by.x = 0, by.y = "Sample.ID")

pcoa_data %>% ggplot(aes_string(x = "Dim1", y = "Dim2", shape = "Sample.Type")) + 
  geom_point(aes_string(color = "Habitat"),size = 3)  + 
  labs(x = labs[1], y = labs[2]) +
  facet_wrap(~paste("Month",Month)) + 
  stat_ellipse(aes_string(fill = "Sample.Type"),geom = "polygon", type = "norm", level = 0.9, alpha = 0.2) + 
  scale_color_brewer(palette = "Paired")
```

# Figure 11 : Barplot Activity - Along24h
![image](Figures/Figure11.png)
```r
p <- phyloseq::plot_bar(physeq, fill = "Activity", x = "Replica") +
  scale_color_manual(na.value = "grey50") + 
  scale_fill_manual(values= Palette(physeq, "Activity"))
p$data[,"Activity"] <- factor(p$data[,"Activity"], levels = c("nocturnal", "both", "diurnal"))
p + xlab("20L Replicates") +
  theme(text=element_text(size = 20)) + scale_x_discrete(guide = guide_axis(angle = 0)) + geom_col(color = "black", size = 0.05) + 
  ggh4x::facet_nested(~ Sampling.Time + paste("Day",Sampling.Day), scales = "free", space = "free_x")

ggsave(path = Images_path, "Barplot_Activity.svg", width = 15, height = 9)
```

# Figure 12 : Nocturnal activity ratio
![image](Figures/Figure12.png)
```r
data <- Tax_melt[,c("Sample.ID","Replica","Family","Taxon","Sampling.Time","Sampling.Day","Nb.reads_sum","Activity")]

ratio <- data %>%
  group_by(Sample.ID, Sampling.Time) %>%
  summarise(
    reads_nocturnal = sum(ifelse(Activity == "nocturnal", Nb.reads_sum, 0), na.rm = TRUE),
    reads_nocturnal_both = sum(ifelse(Activity == "nocturnal" | Activity == "both" , Nb.reads_sum, 0), na.rm = TRUE),
    reads_total = sum(Nb.reads_sum, na.rm = TRUE),
    ratio_nocturnal = reads_nocturnal / reads_total, 
    ratio_nocturnal_both = reads_nocturnal_both / reads_total
  )

shapiro.test(ratio$ratio_nocturnal)
kruskal.test(ratio_nocturnal ~ Sampling.Time, data = ratio)
FSA::dunnTest(ratio_nocturnal ~ Sampling.Time, data = ratio, method = "bonferroni")

p <- ggplot(ratio, aes(x = Sampling.Time, y = ratio_nocturnal)) +
  geom_boxplot(width = 0.6) +
  stat_summary(fun = mean, geom = "point", shape = 18, color = "red") +
  labs(y = "Nocturnal species reads ratio", x = "Sampling time") +
  theme(legend.position = "none")

stat <- data.frame( group1 = c("00:30", "00:30"), group2 = c("12:30", "18:30"), y.position = c(0.6, 0.65),  p.adj = c(0.00065, 0.00019), p.adj.signif = c("***", "***"))
p + ggpubr::stat_pvalue_manual(data = stat, label = "p.adj.signif", tip.length = 0.01)
```