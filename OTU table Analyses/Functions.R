#### Managing dataframe #### 

# Melting function 
Melting <- function(Res, metadatas_selected_col)
{

  Res_melt <- reshape2::melt(Res[,c(1,14:dim(Res)[2])], id = "clusters.id", variable.name = "Run_Barcod", value.name = "Nb.reads")
  Res_melt <- Res_melt[Res_melt$Nb.reads != 0,]
  Res_melt <- merge(Res_melt,Res[,c(1:13)], by = "clusters.id", all = T)
  Res_melt <- merge(Res_melt, metadatas[,metadatas_selected_col], by = "Run_Barcod", all = T)
  # Supprimer les échantillons pour lesquels on a aucun résultat
  Res_melt <- Res_melt[!is.na(Res_melt$clusters.id),]
  # Creating Taxon column
  Res_melt[which(is.na(Res_melt$X.ID)),c("Family","Genus","Species")] <- "unknown"
  Res_melt[which(is.na(Res_melt$Family)),c("Family","Genus","Species")] <- "unknown"
  Res_melt$Taxon <- ifelse(is.na(Res_melt$Genus), Res_melt$Family,paste(Res_melt$Genus, Res_melt$Species))
  unique(Res_melt$Taxon)
  return(Res_melt)
}

Boxplot <- function(Res_melt)
{
  par(mfrow = c(1,2))
  boxplot(Res_melt$alignment.length, ylim = c(0,200), 
          pch = 20, cex = 0.7, xlab = "with unknown", ylab = "alignement length")
  boxplot(subset(Res_melt,Taxon != "unknown")$alignment.length, ylim = c(0,200), 
          pch = 20, cex = 0.7, xlab = "without unknown", ylab = "alignement length")
  
  par(fig=c(0,0.8,0,0.8))
  plot(subset(Res_melt,Taxon != "unknown")$X.ID, xlab = "% identity",
       subset(Res_melt,Taxon != "unknown")$bit.score, ylab = "bit score",
       ylim = c(50,350), pch = 20, cex = 0.7)
  par(fig=c(0,0.8,0.55,1), new = TRUE)
  boxplot(subset(Res_melt,Taxon != "unknown")$X.ID,
          horizontal=TRUE, axes = FALSE, pch = 20, cex = 0.7)
  par(fig=c(0.65,1,0,0.8), new = T)
  boxplot(subset(Res_melt,Taxon != "unknown")$bit.score, 
          axes=FALSE, pch = 20, cex = 0.7, ylim = c(50,350))
}

#### Traitments ####

My_normalize <- function(Tab){
  return(t(t(Tab)/colSums(Tab) *100))
}

normalize_to_100 <- function(x) {
  x / sum(x) * 100
}


Metrics2Groups <- function(Res,Group,Char_g1,Char_g2){
  library(tidyr)
  library(ggpubr)
  library(VennDiagram)
  
  cat("\t\t\t",Char_g1,"\t\t\t",Char_g2,"\n",
      "Mean % ID \t\t", mean(subset(Res,Group==Char_g1)$X.ID,na.rm = T), "\t\t", mean(subset(Res,Group==Char_g2)$X.ID,na.rm = T),"\n",
      "Mean evalue\t\t", mean(subset(Res,Group==Char_g1)$evalue,na.rm = T), "\t\t", mean(subset(Res,Group==Char_g2)$evalue,na.rm = T),"\n"
      )
  
  ID <- ggplot(Res, aes(x=Group, y=X.ID )) + 
    geom_violin() + 
    stat_summary(fun = "mean", geom = "crossbar", width = 0.5, colour = "red")
  Align <- ggplot(Res, aes(x=Group, y=alignment.length )) + 
    geom_violin()+ 
    stat_summary(fun = "mean", geom = "crossbar", width = 0.5, colour = "red")
  mis <- ggplot(Res, aes(x=Group, y=mismatches )) + 
    geom_violin()+ 
    stat_summary(fun = "mean", geom = "crossbar", width = 0.5, colour = "red")
  eval<- ggplot(Res, aes(x=Group, y=evalue)) + 
    geom_violin()+ 
    stat_summary(fun = "mean", geom = "crossbar", width = 0.5, colour = "red")
  
  ggarrange(ID, Align, mis, eval, 
            labels = c("A", "B", "C","D"),
            ncol = 2, nrow = 2)
  ggsave(path = Images_path, filename = "Metrics_ViolinPlot_4Metrics_byGroup.svg", width = 7, height = 7)

  
  ggplot(Res, aes(x=Group, y=Nb_reads)) + 
    geom_violin()
  ggsave(path = Images_path, filename = "Metrics_ViolinPlot_Nb_reads_byGroup.svg", width = 20, height = 7) 
  
  Sum_Tot_1 <- sum(Res$Nb_reads[which(Group==Char_g1)])
  Sum_Tot_2 <- sum(Res$Nb_reads[which(Group==Char_g2)])
  Sum_UnK_1 <- sum(subset(Res,Group == Char_g1 & Specie == "unknown")$Nb_reads)
  Sum_UnK_2 <- sum(subset(Res,Group == Char_g2 & Specie == "unknown")$Nb_reads)
  Nb_Clusters_1 <- dim(subset(Res,Group == Char_g1))[1]
  Nb_Clusters_2 <- dim(subset(Res,Group == Char_g2))[1]
  Nb_Clusters_F_G_1 <- dim(subset(Res,Group == Char_g1 & Specie == ""))[1]
  Nb_Clusters_F_G_2 <- dim(subset(Res,Group == Char_g2 & Specie == ""))[1]
  Nb_Clusters_F__1 <- dim(subset(Res,Group == Char_g1 & Genus == ""))[1]
  Nb_Clusters_F__2 <- dim(subset(Res,Group == Char_g2 & Genus == ""))[1] 
  Nb_Clusters_unK_1 <- dim(subset(Res,Group == Char_g1 & Specie == "unknown"))[1]
  Nb_Clusters_unK_2 <- dim(subset(Res,Group == Char_g2 & Specie == "unknown"))[1]
  
  Res_Agg <- aggregate(as.numeric(Res$Nb_reads), by= list(Group = Group, Taxon =  paste(Res$Family, Res$Genus, Res$Specie, sep = "_")), sum)
  
  
  cat(" Sum Nb reads Tot \t", Sum_Tot_1, "\t\t", Sum_Tot_2, "\n", 
      "Sum Nb reads unknown \t",Sum_UnK_1, "(",sprintf("%.2f",Sum_UnK_1*100/Sum_Tot_1),"% )\t", Sum_UnK_2, "(",sprintf("%.2f",Sum_UnK_2*100/Sum_Tot_2),"% )\n",
      "Nb unique(Taxons) \t", length(Res_Agg$Taxon[which(Res_Agg$Group == Char_g1)]) , "\t\t\t", length(Res_Agg$Taxon[which(Res_Agg$Group == Char_g2)]) ,"\n",
      "Nb unique(Family) \t", length(unique(Res$Family[which(Group==Char_g1)])), "\t\t\t", length(unique(Res$Family[which(Group==Char_g2)])), "\n",
      "Nb unique(Genus) \t", length(unique(Res$Genus[which(Group==Char_g1)])), "\t\t\t", length(unique(Res$Genus[which(Group==Char_g2)])), "\n",
      "Nb unique(Species) \t", length(unique(Res$Specie[which(Group==Char_g1)])), "\t\t\t", length(unique(Res$Specie[which(Group==Char_g2)])), "\n",
      "Nb clusters \t\t", Nb_Clusters_1,"\t\t\t", Nb_Clusters_2,"\n",
      "Nb clusters F_G_ \t", Nb_Clusters_F_G_1,"(",sprintf("%.2f",Nb_Clusters_F_G_1*100/Nb_Clusters_1),"% )\t\t", Nb_Clusters_F_G_2,"(",sprintf("%.2f",Nb_Clusters_F_G_2*100/Nb_Clusters_2),"% )\n",
      "Nb Clusters F__ \t", Nb_Clusters_F__1,"(",sprintf("%.2f",Nb_Clusters_F__1*100/Nb_Clusters_1),"% )\t\t", Nb_Clusters_F__2,"(",sprintf("%.2f",Nb_Clusters_F__2*100/Nb_Clusters_2),"% )\n",
      "Nb Clusters \"unknown\" \t", Nb_Clusters_unK_1,"(",sprintf("%.2f",Nb_Clusters_unK_1*100/Nb_Clusters_1),"% )\t", Nb_Clusters_unK_2,"(",sprintf("%.2f",Nb_Clusters_unK_2*100/Nb_Clusters_2),"% )\n"
  )
  
  #### Venn Diagram for compare OTUs presence from two groups ####
  
  x_sp <- list(
    Char_g1 = as.vector(Res_Agg[which(Res_Agg$Group == Char_g1),"Taxon"]),
    Char_g2 = as.vector(Res_Agg[which(Res_Agg$Group == Char_g2),"Taxon"])
  )
  ggVennDiagram(x_sp, label_alpha = 0, category.names = c(Char_g1,Char_g2)) +
    scale_fill_gradient(low="white",high = "#0073C2FF")#, label = "count", show_intersect = TRUE,fill = "white") #+ ggplot2::scale_fill_gradient(low="white",high = "blue")
  

  venn.diagram(x_sp, 
               filename="ResultsR/Metrics_VennDiagram.png", output=F,
               lwd = 1,cex = 1)
  cat("\nOTUs list present in ", Char_g1, " but not in ", Char_g2, " :\n",
    setdiff(x_sp[[1]],x_sp[[2]]), "\n",
    "\nOTUs list present in ", Char_g2, " but not in ", Char_g1, " :\n",
    setdiff(x_sp[[2]],x_sp[[1]]),"\n")
  
  #### Number of reads distribution by OTUs (separated in two graph because ther's often too muxh of OTUs for one)  #### 
  
  median <- median(Res_Agg$x)
  
  SupMedian <- ggplot(data=subset(Res_Agg,x>median), aes(x=Taxon, y=x, fill = Group)) +
    geom_bar(stat="identity", position=position_dodge()) +
    scale_x_discrete(guide = guide_axis(angle = 45))
  
  InfMedian <- ggplot(data=subset(Res_Agg,x<median), aes(x=Taxon, y=x, fill = Group)) +
    geom_bar(stat="identity", position=position_dodge()) +
    scale_x_discrete(guide = guide_axis(angle = 45))
  
  ggarrange(SupMedian, InfMedian, 
            labels = c("> to median", "< to median"),
            ncol = 1, nrow = 2, common.legend = TRUE, legend = "right")
  ggsave(path = Images_path, filename = "Metrics_Barplot_NbReads_byOTUs.svg", width = 20, height = 10)
  
  #### % ID distribution by OTUs (separated in two graph because ther's often too muxh of OTUs for one) #### 
  
  ggplot(subset(Res,Group == Char_g2), aes(x=Family, y=X.ID )) + 
    geom_violin() + 
    stat_summary(fun = "mean", geom = "crossbar", width = 0.5, colour = "red")+
    scale_x_discrete(guide = guide_axis(angle = 45)) 
  ggsave(path = Images_path, filename = "Metrics_ViolinPlot_XID_byFamily.svg", width = 20, height = 7)
  
  
  Res_Agg <- as.data.frame(aggregate(as.numeric(Res$evalue), by= list(Group = Group, Taxon =  paste(Res$Family, Res$Genus, Res$Specie, sep = "_"), Family = Res$Family), mean))
  Res_Agg1 <- as.data.frame(aggregate(as.numeric(Res$evalue), by= list(Family = Res$Family), mean))
  Res_Agg2 <- as.data.frame(aggregate(as.numeric(Res$evalue), by= list(Group = Group, Family = Res$Family), mean))
  Res_Agg2 <- as.data.frame(aggregate(as.numeric(Res_Agg2$x), by= list(Family = Res_Agg2$Family), median))
  
  
  #Res_Agg <- Res_Agg2[order(Res_Agg2$x, decreasing = T),]
  Res_Agg$Family <- factor(Res_Agg$Family, levels = Res_Agg2[order(Res_Agg2$x, decreasing = F),"Family"])
  
  #Res_Agg$Taxon <- factor(Res_Agg$Taxon, levels = unique(Res_Agg$Taxon))
  
  levels(Res_Agg$Taxon)
  levels(Res_Agg$Family)

  median_list <- as.vector(levels(Res_Agg$Family)[1:as.numeric(length(levels(Res_Agg$Family))%/%2)])
  
  FirstList <- ggplot(Res_Agg, aes(x=Family, y=x, fill = Group)) + 
    geom_violin() + 
    scale_x_discrete(limits = c(median_list), 
                     guide = guide_axis(angle = 45)) +
    stat_summary(fun = "mean", geom = "crossbar", width = 0.5, colour = "red")
    
  median_list <- levels(Res_Agg$Family)[(as.numeric(length(levels(Res_Agg$Family))%/%2)+1):length(levels(Res_Agg$Family))]
    
  SecondList <- ggplot(Res_Agg, aes(x=Family, y=x, fill = Group)) + 
    geom_violin() + 
    scale_x_discrete(limits = c(median_list),
                     guide = guide_axis(angle = 45)) + 
    stat_summary(fun = "mean", geom = "crossbar", width = 0.5, colour = "red")
  
  ggarrange(FirstList, SecondList, 
            #labels = c("A", "B"),
            ncol = 1, nrow = 2, common.legend = TRUE, legend = "right")
  ggsave(path = Images_path, filename = "Metrics_ViolinPlot_evalue_byFamily.svg", width = 20, height = 10)
  
  
}

Metrics <- function(Res){
  Images_path <- paste0(Images_path, "/Metrics")
  library(tidyr)
  library(ggpubr)
  library(VennDiagram)
  library(ggVennDiagram)
  
  for (g in unique(Res$Group)) cat("\t\t\t",g) 
  cat("\nMean % ID\t\t")
  for (g in unique(Res$Group)) cat(mean(subset(Res,Group==g)$X.ID,na.rm = T),"\t\t")
  cat("\nMean evalue\t\t")
  for (g in unique(Res$Group)) cat(mean(subset(Res,Group==g)$evalue,na.rm = T),"\t\t")
  cat("\nMean bit score\t\t")
  for (g in unique(Res$Group)) cat(mean(subset(Res,Group==g)$bit.score,na.rm = T),"\t\t")

  
  ID <- ggplot(Res, aes(x=Group, y=X.ID )) + 
    geom_boxplot(na.rm = T) + 
    stat_summary(fun = "mean", geom = "point", colour = "red", na.rm = T) + labs(x="")
  Align <- ggplot(Res, aes(x=Group, y=alignment.length )) + 
    geom_boxplot(na.rm = T)+ 
    stat_summary(fun = "mean", geom = "point", colour = "red", na.rm = T) + labs(x="")
  mis <- ggplot(Res, aes(x=Group, y=mismatches )) + 
    geom_boxplot(na.rm = T)+ 
    stat_summary(fun = "mean", geom = "point", colour = "red", na.rm = T) + labs(x="")
  evalue <- ggplot(Res, aes(x=Group, y=evalue)) + 
    geom_boxplot(na.rm = T)+ 
    stat_summary(fun = "mean", geom = "point", colour = "red", na.rm = T) + labs(x="")
  log_evalue <- ggplot(Res, aes(x=Group, y=log(evalue))) + 
    geom_boxplot(na.rm = T)+ 
    stat_summary(fun = "mean", geom = "point", colour = "red", na.rm = T) + labs(x="")
  bit.score <- ggplot(Res, aes(x=Group, y=bit.score)) + 
    geom_boxplot(na.rm = T) + 
    stat_summary(fun = "mean", geom = "point", colour = "red", na.rm = T) + labs(x="")
  
  ggarrange(ID, Align, mis, evalue, log_evalue, bit.score,
            labels = c("A", "B", "C","D", "E","F"),
            ncol = 3, nrow = 2) 
  ggsave(path = Images_path, filename = "Boxplots_Metrics_byGroup.svg", width = 12, height = 8)
  
  table(Res$Family)

  
  
  ggplot(Res, aes(x=Group, y=log(Nb_reads))) + 
    geom_boxplot(na.rm=T)
  ggsave(path = Images_path, filename = "Boxplot_Nb_reads_byGroup.svg", width = 14, height = 7) 
  
  cat("\nSum Nb reads Tot \t")
  for (g in unique(Res$Group)) {Sum_Tot <- sum(Res$Nb_reads[which(Res$Group==g)], na.rm=T); cat(Sum_Tot,"\t\t")}
  cat("\nSum Nb reads unknown \t")
  for (g in unique(Res$Group)) 
  {
    Sum_UnK <- sum(subset(Res,Group == g & Specie == "unknown")$Nb_reads, na.rm=T)
    Sum_Tot <- sum(Res$Nb_reads[which(Res$Group==g)], na.rm=T)
    cat(Sum_UnK,"(",sprintf("%.2f",Sum_UnK*100/Sum_Tot),"% )\t\t")
  }
  cat("\nNb unique(Taxons) \t")
  for (g in unique(Res$Group)) cat(length(unique(Res$Taxon[which(Res$Group == g)])) , "\t\t\t")
  cat("\nNb unique(Family) \t")
  for (g in unique(Res$Group)) cat(length(unique(Res$Family[which(Res$Group == g)])) , "\t\t\t")
  cat("\nNb unique(Genus) \t")
  for (g in unique(Res$Group)) cat(length(unique(paste(Res$Family[which(Res$Group == g)],Res$Genus[which(Res$Group == g)]))) , "\t\t\t")
  cat("\nNb unique(Species) \t")
  for (g in unique(Res$Group)) cat(length(unique(Res$Specie[which(Res$Group == g)])) , "\t\t\t")
  cat("\nNb clusters \t\t")
  for (g in unique(Res$Group)) cat(dim(subset(Res,Group == g))[1] , "\t\t\t")
  cat("\nNb clusters F_G_ \t")
  for (g in unique(Res$Group)) {
    Nb_Clusters <- dim(subset(Res,Group == g))[1]
    Nb_Clusters_F_G <- dim(subset(Res,Group == g & is.na(Res$Specie)))[1]
    cat( Nb_Clusters_F_G,"(",sprintf("%.2f",Nb_Clusters_F_G*100/Nb_Clusters),"% )\t\t")
  }
  cat("\nNb clusters F__ \t")
  for (g in unique(Res$Group)) {
    Nb_Clusters <- dim(subset(Res,Group == g))[1]
    Nb_Clusters_F__ <- dim(subset(Res,Group == g & is.na(Res$Genus)))[1]
    cat( Nb_Clusters_F__,"(",sprintf("%.2f",Nb_Clusters_F__*100/Nb_Clusters),"% )\t\t")
  }
  cat("\nNb clusters \"unknown\" \t")
  for (g in unique(Res$Group)) {
    Nb_Clusters <- dim(subset(Res,Group == g))[1]
    Nb_Clusters_unK <- dim(subset(Res,Group == g & is.na(Res$Specie)))[1]
    cat( Nb_Clusters_unK,"(",sprintf("%.2f",Nb_Clusters_unK*100/Nb_Clusters),"% )\t\t")
  }
  cat("\n")

  #### Number of reads distribution by OTUs (separated in two graph because ther's often too much of OTUs for one)  #### 
  Res_Agg <- aggregate(as.numeric(Res$Nb_reads), by= list(Group = Res$Group, Family = Res$Family, Genus = Res$Genus, Specie = Res$Specie), sum)
  Res_Agg$Taxon <- ifelse(is.na(Res_Agg$Genus),Res_Agg$Family,paste(Res_Agg$Genus, Res_Agg$Specie))
  #print("coucouc2")
  Res_Agg$Taxon <- factor(Res_Agg$Taxon, levels = unique(Res_Agg[order(Res_Agg$x, decreasing = T),"Taxon"]))
  #print("coucou3")
  median <- median(Res_Agg$x)
  
  SupMedian <- ggplot(data=subset(Res_Agg,x>median), aes(x=Taxon, y=x, fill = Group)) +
    geom_bar(stat="identity", position=position_dodge()) +
    scale_x_discrete(guide = guide_axis(angle = 45)) + 
    ylab("Nb_reads")
  
  InfMedian <- ggplot(data=subset(Res_Agg,x<median), aes(x=Taxon, y=x, fill = Group)) +
    geom_bar(stat="identity", position=position_dodge()) +
    scale_x_discrete(guide = guide_axis(angle = 45)) + 
    ylab("Nb_reads")
  
  ggarrange(SupMedian, InfMedian, 
            labels = c(paste0("> to ",median), paste0("< to ",median)),
            ncol = 1, nrow = 2, common.legend = TRUE, legend = "right")
  ggsave(path = Images_path, filename = "Barplot_NbReads_byOTUs.svg", width = 20, height = 10)
  
  #### Metrics Distribution by OTUs (separated in two graph because ther's often too much of OTUs for one) #### 
  
  for (m in c("evalue","X.ID","alignment.length","bit.score"))
  {
    for (g in unique(Res$Group))
    {
      #cat("coucou",m,g,"\n")
      Res_temp <- subset(Res, Group == g)
      Res_Agg <- aggregate(as.numeric(Res_temp[,m]), by= list(Taxon = Res_temp$Taxon ), mean)
      #Res_Agg$Taxon <- ifelse(is.na(Res_Agg$Genus),Res_Agg$Family,paste(Res_Agg$Genus, Res_Agg$Specie))
      #cat("coucou2",m,g,"\n")
      Res_Agg$Taxon <- factor(Res_Agg$Taxon, levels = unique(Res_Agg[order(Res_Agg$x, decreasing = F),"Taxon"]))
      #cat("coucou3",m,g,"\n")
      Res_temp$Taxon <- factor(Res_temp$Taxon, levels = levels(Res_Agg$Taxon))
      #cat("coucou4",m,g,"\n")
      cut <- length(levels(as.factor(Res_temp$Taxon)))

      
      SupMedian <- ggplot(Res_temp, aes(x = Taxon, y=Res_temp[,m])) +
        geom_jitter(size = 0.2, na.rm = T) + 
        stat_summary(fun = "mean", geom = "point", colour = "red", na.rm = T) +
        scale_x_discrete(guide = guide_axis(angle = 45), limit = c(levels(as.factor(Res_temp$Taxon))[1:(cut%/%2)])) +
        #scale_y_continuous(breaks = seq(1, 10, 1)) + 
        ylab(m)
      
      InfMedian <- ggplot(Res_temp, aes(x = Taxon, y=Res_temp[,m])) + 
        geom_jitter(size = 0.2, na.rm = T) +  
        stat_summary(fun = "mean", geom = "point", colour = "red", na.rm = T) +
        scale_x_discrete(guide = guide_axis(angle = 45), limit = c(levels(as.factor(Res_temp$Taxon))[((cut%/%2)+1):cut])) +
        #scale_y_continuous(breaks = seq(1, 10, 1)) + 
        ylab(m)
      
      ggarrange(SupMedian, InfMedian, ncol = 1, nrow = 2) 
      ggsave(path = Images_path, filename = paste0("DotPlot_",m,"_byFamily_Group",g,".svg"), width = 20, height = 10)
    }
  }
    
  ggplot(Res, aes(x = Taxon, y=evalue)) + 
    geom_jitter(aes(colour = Group),size = 0.2, na.rm = T) + 
    stat_summary(data = subset(Res,Group == "0.80"),fun = "mean", geom = "point", size = 1, colour = "red", na.rm = T) +
    stat_summary(data = subset(Res,Group == "0.95"),fun = "mean", geom = "point", size = 1, colour = "green", na.rm = T) +
    stat_summary(data = subset(Res,Group == "0.98"),fun = "mean", geom = "point", size = 1, colour = "blue", na.rm = T) +
    scale_x_discrete(guide = guide_axis(angle = 45), limit = c(levels(as.factor(Res_temp$Taxon))[((cut/2)+1):cut])) +
    scale_y_continuous(breaks = seq(1, 10, 1))
  
  #### Venn Diagram for compare OTUs presence from two groups ####
  
  x <- split(Res$Taxon,Res$Group)
  
  ggVennDiagram(x, label_alpha = 0)+ 
    scale_fill_gradient(low="white",high = "#0073C2FF")#, label = "count", show_intersect = TRUE,fill = "white") #+ ggplot2::scale_fill_gradient(low="white",high = "blue")
  ggsave(path = Images_path, filename = "VennDiagram_byTaxon.svg", width = 5, height = 5)
  
}

Violinplot <- function(Res_melt)
{
  ggplot(subset(Res_melt,Taxon != "unknown"), aes(x= "alignment.length" , y= alignment.length)) +
    geom_violin() +
    geom_boxplot(width = 0.2,outlier.shape = NA) +
    geom_hline(yintercept=160, linetype = "dashed", color = "red") + theme(legend.position="none") #+ ylab("Similarity")
  
  ggplot(subset(Res_melt,Taxon != "unknown"), aes(x= "bit.score" , y= bit.score)) +
    geom_violin() +
    geom_boxplot(width = 0.2,outlier.shape = NA) +
    #geom_hline(yintercept=0.964, linetype = "dashed", color = "red") + geom_hline(yintercept=0.882, color = "red") +
    theme(legend.position="none") #+ ylab("Similarity")
  
  ggplot(subset(Res_melt,Taxon != "unknown"), aes(x= "%ID" , y= X.ID)) +
    geom_violin() +
    geom_boxplot(width = 0.2,outlier.shape = NA) +
    #geom_hline(yintercept=0.964, linetype = "dashed", color = "red") + geom_hline(yintercept=0.882, color = "red") +
    theme(legend.position="none") #+ ylab("Similarity")
  #plot(subset(Res_melt,Taxon != "unknown")$X.ID, subset(Res_melt,Taxon != "unknown")$bit.score, ylim = c(50,350), pch = 20, cex = 0.7)
  
  ggplot(subset(Res_melt,Taxon != "unknown unknown"), aes(x= "alignment.length" , y= alignment.length)) +
    geom_violin() +
    geom_boxplot(width = 0.2,outlier.shape = NA) +
    geom_jitter(alpha = 0.1) + 
    geom_hline(yintercept=160, linetype = "dashed", color = "red") + 
    theme_bw() + theme(legend.position="none") #+ ylab("Similarity")
  
  ggplot(subset(Res_melt,Taxon != "unknown unknown"), aes(x= "bit.score" , y= bit.score)) +
    geom_violin() +
    geom_boxplot(width = 0.2,outlier.shape = NA) +
    geom_jitter(alpha = 0.1) + 
    #geom_hline(yintercept=0.964, linetype = "dashed", color = "red") + geom_hline(yintercept=0.882, color = "red") +
    theme_bw() + theme(legend.position="none") #+ ylab("Similarity")
  
  ggplot(subset(Res_melt,Taxon != "unknown unknown"), aes(x= "%ID" , y= X.ID)) +
    geom_violin() +
    geom_boxplot(width = 0.2,outlier.shape = NA) +
    geom_jitter(alpha = 0.1) + 
    #geom_hline(yintercept=0.964, linetype = "dashed", color = "red") + geom_hline(yintercept=0.882, color = "red") +
    theme_bw() + theme(legend.position="none") #+ ylab("Similarity")
  
}

Fouille <- function(df)
{
  Images_path <- paste0(Images_path,"/PCA/")
  #### Analyses en Composante Principales ####
  #http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/

  library("FactoMineR")
  library("factoextra")
  #Note that, by default, the function PCA() [in FactoMineR], standardizes the data automatically during the PCA; 
  #so you don’t need do this transformation before the PCA
  
  res.pca <- PCA(df, scale.unit = TRUE, ncp = 5, graph = TRUE)
  eig.val <- get_eigenvalue(res.pca)
  eig.val
  var <- get_pca_var(res.pca)
  var$contrib
  var$cor
  scree_plot <- fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))
  ggsave(scree_plot,path=Images_path,filename="ScreePlot.svg")#, width = 30, height = 20)
  
  fviz_pca_var(res.pca, col.var = "black")
  cos2 <- fviz_cos2(res.pca, choice = "var", axes = 1:2)
  ggsave(cos2,path=Images_path,filename="Cos2.svg")#, width = 30, height = 20)
  # Color by cos2 values: quality on the factor map
  fviz_pca_var(res.pca, col.var = "cos2",
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
               repel = TRUE # Avoid text overlapping
  )
  # Contributions of variables to PC1
  contrib1 <- fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
  ggsave(contrib1,path=Images_path,filename="ContribVarToDim1.svg")#, width = 30, height = 20)
  # Contributions of variables to PC2
  contrib2 <- fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)
  ggsave(contrib2,path=Images_path,filename="ContribVarToDim2.svg")#, width = 30, height = 20)
  
  # Create a grouping variable using kmeans
  # Create 3 groups of variables (centers = 3)
  set.seed(123)
  res.km <- kmeans(var$coord, centers = 3, nstart = 25)
  grp <- as.factor(res.km$cluster)
  # Color variables by groups
  test <- fviz_pca_var(res.pca, col.var = grp, 
               palette = c("#0073C2FF", "#EFC000FF", "#868686FF"),
               legend.title = "Cluster")
  test
  #ggsave(test,path="ResultsR/test")
  
  ### Graph of individuals ### 
  ind <- get_pca_ind(res.pca)
  
  res.pca.ind <- fviz_pca_ind(res.pca, col.ind = "cos2", 
                              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                              #repel = TRUE # Avoid text overlapping (slow if many points)
                              )
  res.pca.ind
  #Res[which(Res$Num_clusters == "548"),]
  ggsave(res.pca.ind,path=Images_path,filename="PCA_Indiv.svg")#, width = 30, height = 20)
  
  fviz_pca_biplot(res.pca, 
                  palette = "jco", 
                  addEllipses = TRUE, label = "var",
                  col.var = "black", repel = TRUE,
                  legend.title = "Species") 
  
  biplot <- fviz_pca_biplot(res.pca, 
                  col.ind = as.factor(df[,"Time"]),#, palette = "jco", 
                  #col.ind.sup = "blue",
                  addEllipses = TRUE,
                  label = "var",
                  col.var = "black", repel = F,
                  legend.title = "Time")
  
  ggsave(biplot, path=Images_path, filename="PCA_Biplot.svg")#, width = 30, height = 20)
  
  # Total contribution on PC1 and PC2
  ind_cos2 <- fviz_cos2(res.pca, choice = "ind")
  contrib <- ind_cos2$data
  contrib <- contrib[order(contrib$cos2, decreasing = T),]
  ind_contrib <- fviz_contrib(res.pca, choice = "ind", axes = 1:2)
  ggsave(ind_contrib, path=Images_path, filename="Indiv_contribution.svg", width = 40, height = 10)
  
  
}

#### Hierarchical Ascending Classification #####
HA <- function(Tab, file_path, line){
  TabN <- My_normalize(Tab)
  #colSums(TabN) #each should be equal to 100
  
  clusturing <- hclust(dist(t(TabN)))
  
  {plot.new()
    svg(file_path,  width = 15, height = 15)
    plot(as.dendrogram(clusturing), horiz = T)
    #abline(v=line, col = 'red')
    dev.off()
  }
}

#### HeatMap ####
# x = x (nb_reads)
# y = Taxon 
# facet.grid = ~Bio_rep
Heatmap <- function(Tab, path, file_name, x){
  ggplot(as.data.frame(TabNMelt)) +
    geom_tile(aes(x= x , y= Taxon, fill = log(Nb_reads))) +
    #scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    facet_grid(~Sampling_replicates, scales = "free_x") +
    scale_fill_gradient(low = "#55aff4", high = "#203850")
  ggsave(path = path, filename = file_name, width = 14, height = 7, device='tiff', dpi=700)
}

#### Slider ####
Slider <- function(df){
  # Load the required libraries
  library(shiny)
  
  df <- melt(Tab)
  names(df) <- c("OTUs","Samples","Nb_reads")
  df <- df[which(df$Nb_reads != 0 ) , ]
  
  
  # Define the user interface for the Shiny app
  ui <- fluidPage(
    # Add a title to the page
    titlePanel("Interactive Table with Slider"),
    
    # Add a text input field
    textInput("max_value", "Enter the maximum of the slider you want:", value = max(df$Nb_reads)),
    
    # Define the sidebar layout
    sidebarLayout(
      # Define the sidebar
      sidebarPanel(
        # Add a slider input
        sliderInput("slider", "Select a number:", min = 0, max =  5000, value = 2000)
      ),
      
      # Define the main panel
      mainPanel(
        # Add a text output
        textOutput("percentage"),
        # Add a table output
        tableOutput("table")  
        
      )
    )
  )
  
  # Define the server function for the Shiny app
  server <- function(input, output, session) {
    
    
    # Create a reactive expression that updates the maximum value of the slider
    observeEvent(input$max_value, {
      updateSliderInput(session, "slider", max = as.numeric(input$max_value))
    })
    
    # Create a reactive expression that filters the data based on the value of the slider
    filtered_data <- reactive({
      df[df$Nb_reads < input$slider, ]
    })
    # Sort the filtered data in ascending order by the OTUs column
    sorted_data <- reactive({
      filtered_data()[order(filtered_data()$OTUs), ]
    })
    
    
    # Calculate the percentage of names displayed on the total number of names
    percentage <- reactive({
      nrow(unique(sorted_data())) / nrow(unique(df)) * 100
    })
    
    # Render the percentage as a text output
    output$percentage <- renderText({
      paste0(round(percentage(), 2), "%")
    })
    
    
    # Render the data frame as a table in the table output
    output$table <- renderTable({
      unique(sorted_data()$OTUs)
    })
    
  }
  
  # Run the Shiny app
  shinyApp(ui, server)
}

#### Phyloseq ####

Palette <- function(physeq, tax_level)
{
  set.seed(1531564)
  # species Palette 
  RColorBrewer::display.brewer.all(colorblindFriendly = TRUE)
  getPalette = colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))
  speciesList = unique(phyloseq::tax_table(physeq)[,tax_level])
  speciesPalette = getPalette(length(speciesList))
  names(speciesPalette) = speciesList
  speciesPalette["Hominidae"] <- "#444444"
  speciesPalette["Others"] <- "#444444"
  return(speciesPalette)
}



Phyloseq <- function(Res, rarefied, sampleDBool){
  #https://joey711.github.io/phyloseq/import-data.html
  #https://micca.readthedocs.-io/en/latest/phyloseq.html
  
  # OTUs object for phyloseq
  OTUs <- phyloseq::otu_table(as.data.frame(Res[,c(20:dim(Res)[2])]), taxa_are_rows = T)
  
  # TAX objct for phyloseq
  TAX = cbind(Family = Res$Family, Genus = Res$Genus, Species = Res$Species)
  rownames(TAX) <- rownames(OTUs)
  TAX <- cbind(TAX, OTUs = ifelse(is.na(TAX[,"Genus"]),TAX[,"Family"],paste(TAX[,"Genus"], TAX[,"Species"])))
  TAX <- phyloseq::tax_table(TAX)
  
  physeq <- phyloseq::phyloseq(OTUs, TAX)
  
  # Sample_data object for phyloseq
  if (sampleDBool){
    rownames(metadatas) <- metadatas$Run_Barcod
    Sample_data <- phyloseq::sample_data(metadatas)
    physeq <- phyloseq::merge_phyloseq(physeq, Sample_data)
  }
  
  if (rarefied){
    # Rarefy the samples without replacement.
    # Rarefaction is used to simulate even number of reads per sample.
    # In this example, the rarefaction depth chosen is the 90% of the minimum sample depth in the dataset
    physeq = phyloseq::rarefy_even_depth(physeq, rngseed=1, sample.size=0.9*min(phyloseq::sample_sums(physeq)), replace=F)
  }
  
  saveRDS(physeq, "ResultsR/physeq.RDS")
  return(physeq)
}



Phyloseq_Class <- function(Tab,order,filename){
  library("phyloseq")
  library("vegan")
  
  #https://micca.readthedocs.-io/en/latest/phyloseq.html
  OTUs <- otu_table(Tab, taxa_are_rows = T)
  
  #https://joey711.github.io/phyloseq/import-data.html
  
  TAX = matrix(str_split(rownames(OTUs), "_", simplify = T),nrow = nrow(OTUs), ncol = 4)
  rownames(TAX) <- rownames(OTUs)
  colnames(TAX) <- c("Class", "Family", "Genus", "Species")
  TAX <- tax_table(TAX)
  
  physeq <- phyloseq(OTUs, TAX)
  #sample_data(physeq)
  
  # Rarefy the samples without replacement.
  # Rarefaction is used to simulate even number of reads per sample.
  # In this example, the rarefaction depth chosen is the 90% of the minimum sample depth in the dataset
  # (in this case 459 reads per sample). NO
  
  physeq.rarefied = rarefy_even_depth(physeq, rngseed=1, sample.size=0.9*min(sample_sums(physeq)), replace=F)
  
  #return(physeq)
  #return(physeq.rarefied)
  
  p <- plot_bar(physeq.rarefied, fill = "Class") + scale_color_manual(na.value = "grey50")
  p$data$Family <- factor(p$data$Family, levels = rev(unique(c("unknown", as.vector(p$data$Family[order(p$data$Abundance, decreasing = T)])))))
  p$data$Sample <- factor(p$data$Sample, levels = order)
  p$data[which(p$data$Class != "Actinopterygii"),c("OTU", "Class", "Family", "Genus","Species")] <- NA
  
  return(p)
  
}



BetaDiv <- function(physeq, x)
{
  #https://joey711.github.io/phyloseq/import-data.html
  
  library("ape")
  random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
  # plot(random_tree) # pas intéressant pour moi mais il faut le faire quand même pour la suite 
  
  physeq = merge_phyloseq(random_tree, physeq)
  
  # PCoA plot using the unweighted UniFrac as distance
  wunifrac_dist = phyloseq::distance(physeq, method="unifrac", weighted=F)
  ordination = ordinate(physeq, method="PCoA", distance=wunifrac_dist)
  #plot_ordination(physeq, ordination, color = x, shape = "Size") + theme(aspect.ratio=1)
  #plot_ordination(physeq, ordination, color = x) + theme(aspect.ratio=1)
  p <- plot_ordination(physeq, ordination) + theme(aspect.ratio=1)
  
  print(p)
  
  ggsave(path = Images_path, filename = "BetaDiv.svg")
}


Tree <- function(ps, color, shape, label)
{
  #https://joey711.github.io/phyloseq/plot_tree-examples.html
  
  
  # ps = object phyloseq
  # label = "Species" 
  # color = 1 condition décrivant l'échantillon ou NULL
  # shape = 1 autre ou NULL
  
  
  library(ape)
  random_tree = rtree(ntaxa(ps), rooted = T, tip.label = taxa_names(ps))

  
  if(!missing(color) & !missing(shape) & !missing(label))
  {
    ps_tree = merge_phyloseq(ps, random_tree)
    sample_data(ps_tree)
    return(plot_tree(ps_tree, color = color, shape = shape, label.tips=label, ladderize="left", text.size = 5))
  }
  
  return(random_tree)

}

#### Accumulation Curves ####
AccPlot <- function(Tab){
  R1 <- length(unique(Tab$Taxon[which(Tab$Sample.ID == "1")]))
  R2 <- R1 + length(setdiff(unique(Tab$Taxon[which(Tab$Sample.ID == "2")]),unique(Tab$Taxon[which(Tab$Sample.ID == "1")])))
  R3 <- R2 + length(setdiff(unique(Tab$Taxon[which(Tab$Sample.ID == "3")]),unique(Tab$Taxon[which(Tab$Sample.ID == "1" | Tab$Sample.ID == "2")])))
  R4 <- R3 + length(setdiff(unique(Tab$Taxon[which(Tab$Sample.ID == "4")]),unique(Tab$Taxon[which(Tab$Sample.ID == "1" | Tab$Sample.ID == "2"| Tab$Sample.ID == "3")])))
  R5 <- R4 + length(setdiff(unique(Tab$Taxon[which(Tab$Sample.ID == "5")]),unique(Tab$Taxon[which(Tab$Sample.ID == "1" | Tab$Sample.ID == "2"| Tab$Sample.ID == "3"| Tab$Sample.ID == "4")])))
  R6 <- R5 + length(setdiff(unique(Tab$Taxon[which(Tab$Sample.ID == "6")]),unique(Tab$Taxon[which(Tab$Sample.ID == "1" | Tab$Sample.ID == "2"| Tab$Sample.ID == "3"| Tab$Sample.ID == "4"| Tab$Sample.ID == "5")])))
  R7 <- R6 + length(setdiff(unique(Tab$Taxon[which(Tab$Sample.ID == "7")]),unique(Tab$Taxon[which(Tab$Sample.ID == "1" | Tab$Sample.ID == "2"| Tab$Sample.ID == "3"| Tab$Sample.ID == "4"| Tab$Sample.ID == "5" | Tab$Sample.ID == "6")])))
  R8 <- R7 + length(setdiff(unique(Tab$Taxon[which(Tab$Sample.ID == "8")]),unique(Tab$Taxon[which(Tab$Sample.ID != "8")])))
  
  x= c(1,2,3,4,5,6,7,8)
  y = c(R1,R2,R3,R4,R5,R6,R7,R8)
  
  p <- plot( x= c(0,x),
             y = c(0,y), 
             type = "b",
             xlab = "Replicates", 
             ylab = "OTUs Richness", 
             #xaxt = "n",
             main = "Accumulation curve", )
  print(p)
  
  data <- data.frame(x=c(0.05,x),y=c(0,y))
  data <- data.frame(x=x,y=y)

  p <- ggplot(data, aes(x=x,y=y)) +
    geom_point() +
    labs(title = "Accumulation Curve", x = "Replicates", y = "OTUs Richness") +
    #xlim(0,25) + ylim(0, max(y) + 25) + 
    stat_smooth(method="lm", formula=y~log(x),fullrange = T)
  
  print(p)
  return(p)
}

AccPlot2 <- function(Tab){
  
  data <- data.frame()

  for (s in unique(Tab$Site)){
    data_temp <- data.frame()
    Tab_temp <- subset(TabNMelt,Site == s)
    R1 <- length(Tab_temp$OTUs[which(Tab_temp$Rep == "1")])
    R2 <- R1 + length(setdiff(Tab_temp$OTUs[which(Tab_temp$Rep == "2")],Tab_temp$OTUs[which(Tab_temp$Rep == "1")]))
    R3 <- R2 + length(setdiff(Tab_temp$OTUs[which(Tab_temp$Rep == "3")],unique(Tab_temp$OTUs[which(Tab_temp$Rep == "1" | Tab_temp$Rep == "2")])))
    R4 <- R3 + length(setdiff(Tab_temp$OTUs[which(Tab_temp$Rep == "4")],unique(Tab_temp$OTUs[which(Tab_temp$Rep == "1" | Tab_temp$Rep == "2"| Tab_temp$Rep == "3")])))
    R5 <- R4 + length(setdiff(Tab_temp$OTUs[which(Tab_temp$Rep == "5")],unique(Tab_temp$OTUs[which(Tab_temp$Rep == "1" | Tab_temp$Rep == "2"| Tab_temp$Rep == "3"| Tab_temp$Rep == "4")])))
    
    data_temp <- data.frame( x = c(seq(1,5)), y  = c(R1,R2,R3,R4,R5))
    data_temp$Site <- s 
    data <- rbind(data,data_temp)
  }
  rm(data_temp,Tab_temp)
  
  data %>%
    ggplot(aes(x=x,y=y, color = Site))+
    geom_point()+
    labs(title = "Accumulation Curve", x = "Replicates", y = "OTUs Richness") +
    xlim(1,30)+ ylim(0, 110) +
    stat_smooth(method="lm", formula=y~log(x), fullrange = T)
  ggsave(path= Images_path, filename = "AccumulationCurve_PlotModel.svg")
  
  
}

NbOTUsBarplot <- function(Tab){
  # Tab : An Aggregate Table
  N1 <- length(Tab$Taxon[which(Tab$Sample == "1")])
  N2 <- length(Tab$Taxon[which(Tab$Sample == "2")])
  N3 <- length(Tab$Taxon[which(Tab$Sample == "3")])
  N4 <- length(Tab$Taxon[which(Tab$Sample == "4")])
  N5 <- length(Tab$Taxon[which(Tab$Sample == "5")])
  N6 <- length(Tab$Taxon[which(Tab$Sample == "6")])
  N7 <- length(Tab$Taxon[which(Tab$Sample == "7")])
  N8 <- length(Tab$Taxon[which(Tab$Sample == "8")])
  Ntot <- length(unique(Tab$Taxon[which(Tab$Sample != "-" & Tab$Sample != "+" & Tab$Sample != "TN")]))
  
  data <- data.frame(
    name=c("Total",seq(1:8)) ,  
    value=c(Ntot,N1,N2,N3,N4,N5,N6,N7,N8)
  )
  
  p <- ggplot(data, aes(x=name, y=value)) + 
    geom_bar(stat = "identity")
  
  p
  
  return(p)
  
  
}


AccBarplot <- function(Tab){
  # Tab : An Aggregate Table
  
  library(plotly)
  
  N1 <- unique(Tab$Taxon[which(Tab$Sample == "1")])
  N2 <- setdiff(Tab$Taxon[which(Tab$Sample == "2")],unique(Tab$Taxon[which(Tab$Sample == "1")]))
  N3 <- setdiff(Tab$Taxon[which(Tab$Sample == "3")],unique(Tab$Taxon[which(Tab$Sample == "1" | Tab$Sample == "2")]))
  N4 <- setdiff(Tab$Taxon[which(Tab$Sample == "4")],unique(Tab$Taxon[which(Tab$Sample == "1" | Tab$Sample == "2"| Tab$Sample == "3")]))
  N5 <- setdiff(Tab$Taxon[which(Tab$Sample == "5")],unique(Tab$Taxon[which(Tab$Sample == "1" | Tab$Sample == "2"| Tab$Sample == "3"| Tab$Sample == "4")]))
  N6 <- setdiff(Tab$Taxon[which(Tab$Sample == "6")],unique(Tab$Taxon[which(Tab$Sample == "1" | Tab$Sample == "2"| Tab$Sample == "3"| Tab$Sample == "4"| Tab$Sample == "5")]))
  N7 <- setdiff(Tab$Taxon[which(Tab$Sample == "7")],unique(Tab$Taxon[which(Tab$Sample == "1" | Tab$Sample == "2"| Tab$Sample == "3"| Tab$Sample == "4"| Tab$Sample == "5" | Tab$Sample == "6")]))
  N8 <- setdiff(Tab$Taxon[which(Tab$Sample == "8")],unique(Tab$Taxon[which(Tab$Sample != "8")]))
  NTN <- setdiff(Tab$Taxon[which(Tab$Sample == "TN")],unique(Tab$Taxon[which(Tab$Sample == "1" | Tab$Sample == "2"| Tab$Sample == "3"| Tab$Sample == "4"| Tab$Sample == "5" | Tab$Sample == "6"| Tab$Sample == "7"| Tab$Sample == "8")]))
  Npos <- setdiff(Tab$Taxon[which(Tab$Sample == "+")],unique(Tab$Taxon[which(Tab$Sample == "1" | Tab$Sample == "2"| Tab$Sample == "3"| Tab$Sample == "4"| Tab$Sample == "5" | Tab$Sample == "6"| Tab$Sample == "7"| Tab$Sample == "8"| Tab$Sample == "TN")]))
  Nneg <- setdiff(Tab$Taxon[which(Tab$Sample == "-")],unique(Tab$Taxon[which(Tab$Sample == "1" | Tab$Sample == "2"| Tab$Sample == "3"| Tab$Sample == "4"| Tab$Sample == "5" | Tab$Sample == "6"| Tab$Sample == "7"| Tab$Sample == "8"| Tab$Sample == "TN"| Tab$Sample == "+")]))
  
  
  Cumul <- seq(1:8)
  R1 <- c(rep(0,0),rep(length(N1),8))
  R2 <- c(rep(0,1),rep(length(N2),7))
  R3 <- c(rep(0,2),rep(length(N3),6))
  R4 <- c(rep(0,3),rep(length(N4),5))
  R5 <- c(rep(0,4),rep(length(N5),4))
  R6 <- c(rep(0,5),rep(length(N6),3))
  R7 <- c(rep(0,6),rep(length(N7),2))
  R8 <- c(rep(0,7),rep(length(N8),1))
  
  data <- data.frame(Cumul, R1,R2,R3,R4,R5,R6,R7,R8)
  
  fig <- plot_ly(data, x = ~Cumul, y = ~R1, type = 'bar', name = '1') %>%
    add_trace(y = ~R2, name = '2') %>%
    add_trace(y = ~R3, name = '3') %>%
    add_trace(y = ~R4, name = '4') %>%
    add_trace(y = ~R5, name = '5') %>% #toString(N5)) %>%
    add_trace(y = ~R6, name = '6') %>% #toString(N6)) %>%
    add_trace(y = ~R7, name = '7') %>% #toString(N7)) %>%
    add_trace(y = ~R8, name = '8') %>% #toString(N8)) %>%
    layout(yaxis = list(title = 'Number of OTUs'),
           barmode = 'stack',
           legend = list(orientation = 'h'))
  
  
  ###### Autre option
  
  df <- data.frame(
    Cumul = (c(rep(1,8),rep(2,8),rep(3,8),rep(4,8),rep(5,8),rep(6,8),rep(7,8),rep(8,8))),
    Replicates = (paste0("R",rep(1:8,8))), 
    OTUs = c(rep(0,7),length(N1),
              rep(0,6),length(N2),length(N1),
              rep(0,5),length(N3),length(N2),length(N1),
              rep(0,4),length(N4),length(N3),length(N2),length(N1),
              0,0,0,length(N5),length(N4),length(N3),length(N2),length(N1),
              0,0,length(N6),length(N5),length(N4),length(N3),length(N2),length(N1),
              0,length(N7),length(N6),length(N5),length(N4),length(N3),length(N2),length(N1),
              length(N8),length(N7),length(N6),length(N5),length(N4),length(N3),length(N2),length(N1)))
  
  ggplot(data=df, aes(x=Cumul, y=OTUs, fill=Replicates)) +
    geom_bar(stat="identity") + 
    #theme(legend.position="bottom")+ 
    scale_fill_discrete(labels=c("8","7","6","5","4","3","2","1"))
  ggsave(path = Images_path , filename = "AccumulationCurves_Barplot_revGS.svg")
  
    return(fig)

}


#### Fantaxtic (lucie Version) ####
# Je change x = "Sample" par x = x_value (à passer en argument)
# Je supprime l'angle de l'axe x également
# je change le theme ligth par bw

plot_nested_bar_Lucie <- function (ps_obj, top_level, nested_level, top_merged_label = "Other", x_value,
          nested_merged_label = "Other <tax>", palette = NULL, base_clr = "#008CF0", 
          merged_clr = "grey90", include_rank = T, na_taxon_label = "<tax> (<rank>)", 
          asv_as_id = F, duplicate_taxon_label = "<tax> <id>", relative_abundances = T, 
          sample_order = NULL, ...) 
{
  library(dplyr)
  #detach("package:fantaxtic", unload = T)
  ps_tmp <- ps_obj %>% fantaxtic::name_na_taxa(include_rank = include_rank, 
                                    na_label = na_taxon_label)
  ps_tmp <- ps_tmp %>% fantaxtic::label_duplicate_taxa(tax_level = nested_level, 
                                            asv_as_id = asv_as_id, duplicate_label = duplicate_taxon_label)
  pal <- fantaxtic::taxon_colours(ps_tmp, tax_level = top_level, merged_label = top_merged_label, 
                       merged_clr = merged_clr, palette = palette, base_clr = base_clr)
  psdf <- phyloseq::psmelt(ps_tmp)
  psdf <- fantaxtic::move_label(psdf = psdf, col_name = top_level, label = top_merged_label, 
                     pos = 0)
  psdf <- fantaxtic::move_nested_labels(psdf, top_level = top_level, 
                             nested_level = nested_level, top_merged_label = top_merged_label, 
                             nested_label = gsub("<tax>", "", nested_merged_label), 
                             pos = Inf)
  if (!is.null(sample_order)) {
    if (all(sample_order %in% unique(psdf$Sample))) {
      psdf <- psdf %>% mutate(Sample = factor(Sample, 
                                              levels = sample_order))
    }
    else {
      stop("Error: not all(sample_order %in% sample_names(ps_obj)).")
    }
  }
  p <- ggnested::ggnested(psdf, aes_string(main_group = top_level, sub_group = nested_level, 
                                 x = x_value, y = "Abundance"), ..., main_palette = pal) + 
    scale_y_continuous(expand = c(0, 0)) + ggnested::theme_nested(theme_bw) + guides(fill=guide_legend(ncol=2)) #+ 
    #theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90))
  if (relative_abundances) {
    p <- p + geom_col(position = position_fill())
  }
  else {
    p <- p + geom_col()
  }
  return(p)
}
