library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

#msigdbr_species()
#msigdbr_collections()

# For mouse-specific gene sets, use:
#mm_gene_sets <- msigdbr(species = "Mus musculus")

# To see available categories for mouse:
#unique(mm_gene_sets$gs_cat)


WP_OR_GSEA <- function(instance_markers,name,path,ctgry){
  instance_markers <- instance_markers[order(-instance_markers$avg_log2FC),]
  gene_symbols <- instance_markers$Names
  entrez_ids <- mapIds(org.Hs.eg.db, keys = gene_symbols, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  wp_result <- enrichWP(gene = entrez_ids, organism = "Homo sapiens", pvalueCutoff = 0.05)
  wp_result <- setReadable(wp_result, 'org.Hs.eg.db', 'ENTREZID')
  
  if (is.null(wp_result) || nrow(wp_result@result) == 0) {
    message(paste("No significant WIKI pathways found for", name, "- skipping file."))
    return(NULL)  # Skip further processing
  }
  #gene_symbols <- instance_markers$Names
  gene_list <- instance_markers$avg_log2FC
  names(gene_list) <- rownames(instance_markers)
  plot2 <- dotplot(wp_result)
  ggsave(paste0(path,"WP_",name,".png"),plot=plot2)
  #barplot2 <- barplot(wp_result, showCategory=ctgry)
  #ggsave(paste0(path,"Bar_pot_","WP_",name,".png"),plot=barplot2)
  heatmap2 <- heatplot(wp_result, foldChange=gene_list, showCategory=ctgry)
  ggsave(paste0(path,"Heatmap_pot_","WP_",name,".jpeg"),plot=heatmap2)
  
  edox <- setReadable(wp_result, 'org.Hs.eg.db', 'ENTREZID')
  p2 <- cnetplot(edox, color.params = list(foldChange=gene_list))
  ggsave(paste0(path,"Cnet_plot_","WP_",name,".jpeg"),plot=p2)
  
  #edo_heat2 <- enrichDGN(entrez_ids)
  #edox_heat2 <- setReadable(edo_heat2, 'org.Hs.eg.db', 'ENTREZID')
  #p_heat2 <- heatplot(edox_heat2, foldChange=gene_list11, showCategory=ctgry)
  #ggsave(paste0(path,"Real_Heatmap_","WP_",name,".jpeg"),plot=p_heat2)
  
  write.xlsx(as.data.frame(wp_result),paste0(path,"WikiPathway-",name,"-Complete_DataFrame.xlsx"))
}

Reactome_OR_GSEA <- function(instance_markers,name,path,ctgry){
  library(ReactomePA)
  instance_markers <- instance_markers[order(-instance_markers$avg_log2FC),]
  gene_symbols <- instance_markers$Names
  entrez_ids <- mapIds(org.Hs.eg.db, keys = gene_symbols, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  
  reac_result <- enrichPathway(gene = entrez_ids, organism = "Homo sapiens", pvalueCutoff = 0.05)
  
  if (is.null(reac_result) || nrow(reac_result@result) == 0) {
    message(paste("No significant Reactome pathways found for", name, "- skipping file."))
    return(NULL)  # Skip further processing
  }
  reac_result <- setReadable(reac_result, 'org.Hs.eg.db', 'ENTREZID')
  #gene_symbols <- instance_markers$Names
  gene_list <- instance_markers$avg_log2FC
  names(gene_list) <- rownames(instance_markers)
  plot2 <- dotplot(reac_result)
  ggsave(paste0(path,"Reactome_",name,".png"),plot=plot2)
  print("plot 1 is saved")
  #barplot2 <- barplot(reac_result, showCategory=ctgry)
  #ggsave(paste0(path,"Bar_pot_","Reactome_",name,".png"),plot=barplot2)
  print("plot 2 bar is saved")
  heatmap2 <- heatplot(reac_result, foldChange=gene_list, showCategory=ctgry)
  ggsave(paste0(path,"Heatmap_pot_","Reactome_",name,".jpeg"),plot=heatmap2)
  
  edox <- setReadable(reac_result, 'org.Hs.eg.db', 'ENTREZID')
  p2 <- cnetplot(edox, color.params = list(foldChange=gene_list))
  ggsave(paste0(path,"Cnet_plot_","Reactome_",name,".jpeg"),plot=p2)
  
  #edo_heat2 <- enrichDGN(entrez_ids)
  #edox_heat2 <- setReadable(edo_heat2, 'org.Hs.eg.db', 'ENTREZID')
  #p_heat2 <- heatplot(edox_heat2, foldChange=gene_list11, showCategory=ctgry)
  #ggsave(paste0(path,"Real_Heatmap_","Reactome_",name,".jpeg"),plot=p_heat2)
  
  write.xlsx(as.data.frame(reac_result),paste0(path,"ReactomePathway-",name,"-Complete_DataFrame.xlsx"))
}
kegg_OR_GSEA <- function(instance_markers,name,path,ctgry){
  instance_markers <- instance_markers[order(-instance_markers$avg_log2FC),]
  gene_symbols <- instance_markers$Names
  entrez_ids <- mapIds(org.Hs.eg.db, keys = gene_symbols, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  kegg_result <- enrichKEGG(gene = entrez_ids, organism = 'human', pvalueCutoff = 0.05)
  kegg_result <- setReadable(kegg_result, 'org.Hs.eg.db', 'ENTREZID')
  
  if (is.null(kegg_result) || nrow(kegg_result@result) == 0) {
    message(paste("No significant KEGG pathways found for", name, "- skipping file."))
    return(NULL)  # Skip further processing
  }
  #gene_symbols <- instance_markers$Names
  gene_list <- instance_markers$avg_log2FC
  names(gene_list) <- rownames(instance_markers)
  plot2 <- dotplot(kegg_result)
  ggsave(paste0(path,"KEGG",name,".png"),plot=plot2)
  #barplot2 <- barplot(kegg_result, showCategory=ctgry)
  #ggsave(paste0(path,"Bar_pot_","KEGG_",name,".png"),plot=barplot2)
  heatmap2 <- heatplot(kegg_result, foldChange=gene_list, showCategory=ctgry)
  ggsave(paste0(path,"Heatmap_pot_","KEGG_",name,".jpeg"),plot=heatmap2)
  
  edox <- setReadable(kegg_result, 'org.Hs.eg.db', 'ENTREZID')
  p2 <- cnetplot(edox, color.params = list(foldChange=gene_list))
  ggsave(paste0(path,"Cnet_plot_","KEGG_",name,".jpeg"),plot=p2)
  
  #edo_heat2 <- enrichDGN(entrez_ids)
  #edox_heat2 <- setReadable(edo_heat2, 'org.Hs.eg.db', 'ENTREZID')
  #p_heat2 <- heatplot(edox_heat2, foldChange=gene_list11, showCategory=ctgry)
  #ggsave(paste0(path,"Real_Heatmap_","KEGG_",name,".jpeg"),plot=p_heat2)
  
  write.xlsx(as.data.frame(kegg_result),paste0(path,"KeggPathway-",name,"-Complete_DataFrame.xlsx"))
}

TA_cell <- read.csv("D:/Data_Sets/GSE180286/New/Immune/outputs/DEGs New/B_cells/B_Cell-TLN 1-Tumor.csv")
TA_cell <-  gsea_preprocess("D:/Data_Sets/GSE180286/New/Immune/outputs/DEGs New/B_cells/B_Cell-TLN 1-Tumor.csv")
rownames(TA_cell) <- TA_cell$Names
names(TA_cell)
WP_OR_GSEA(instance_markers = TA_cell,name = "B_cell",path = "D:/Data_Sets/GSE180286/New/Immune/outputs/DEGs New/B_cells/test/",ctgry = 20)

Reactome_OR_GSEA(instance_markers = TA_cell,name = "B_cell",path = "D:/Data_Sets/GSE180286/New/Immune/outputs/DEGs New/B_cells/test/",ctgry = 20)

wp_simp_plotter_human2(instance_markers = TA_cell,name = "B_cell",path = "D:/Data_Sets/GSE180286/New/Immune/outputs/DEGs New/B_cells/test/",ctgry = 20)
