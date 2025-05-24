Melting <- function(Res, metadatas_selected_col)
{
  Res_melt <- reshape2::melt(Res[,c(1,15:dim(Res)[2])], id = "clusters.id", variable.name = "Run_Barcod", value.name = "Nb.reads")
  Res_melt <- Res_melt[Res_melt$Nb.reads != 0,]
  Res_melt <- merge(Res_melt,Res[,c(1:14)], by = "clusters.id", all = T)
  Res_melt <- merge(Res_melt, metadatas[,metadatas_selected_col], by = "Run_Barcod", all = T)
  # Remove samples that haven't any result
  Res_melt <- Res_melt[!is.na(Res_melt$clusters.id),]
  # Creating Taxon column
  Res_melt[which(is.na(Res_melt$Family)),c("Family","Genus","Species")] <- "unknown"
  Res_melt$Taxon <- ifelse(is.na(Res_melt$Genus), Res_melt$Family,paste(Res_melt$Genus, Res_melt$Species))
  unique(Res_melt$Taxon)
  return(Res_melt)
}

#### Fantaxtic (lucie Version) ####
#Replace x = "Sample" with x = x_value (to be passed as an argument)
#Remove the x-axis angle
#Change the theme from theme_light to theme_bw

plot_nested_bar_Lucie <- function (ps_obj, top_level, nested_level, top_merged_label = "Other", x_value,
          nested_merged_label = "Other <tax>", palette = NULL, base_clr = "#008CF0", 
          merged_clr = "grey90", include_rank = T, na_taxon_label = "<tax> (<rank>)", 
          asv_as_id = F, duplicate_taxon_label = "<tax> <id>", relative_abundances = T, 
          sample_order = NULL, ...) 
{
  library(dplyr)
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
