library(seurat)





standard_processing_function <- function(GSESeeuratObject,Resolution=0.05 ){
  
  GSESeeuratObject <- NormalizeData(object =GSESeeuratObject, verbose = FALSE,normalization.method = "LogNormalize",  scale.factor = 10000)
  GSESeeuratObject <- FindVariableFeatures(object =GSESeeuratObject, nfeatures = 2000, verbose = FALSE, selection.method = 'vst')
  GSESeeuratObject  <- ScaleData(GSESeeuratObject, verbose = FALSE)
  GSESeeuratObject <- RunPCA(GSESeeuratObject, features = VariableFeatures(object = GSESeeuratObject))
  GSESeeuratObject <- FindNeighbors(GSESeeuratObject, dims = 1:30)
  GSESeeuratObject <- FindClusters(GSESeeuratObject, resolution = Resolution)
  GSESeeuratObject <- RunUMAP(GSESeeuratObject, dims = 1:30)
  DimPlot(GSESeeuratObject, reduction = "umap", label=TRUE)
  return(GSESeeuratObject)
}

standard_processing_function <- function(GSESeeuratObject ){
  
  GSESeeuratObject <- NormalizeData(object =GSESeeuratObject, verbose = FALSE,normalization.method = "LogNormalize",  scale.factor = 10000)
  GSESeeuratObject <- FindVariableFeatures(object =GSESeeuratObject, nfeatures = 3000, verbose = FALSE, selection.method = 'vst')
  GSESeeuratObject  <- ScaleData(GSESeeuratObject, verbose = FALSE)
  GSESeeuratObject <- RunPCA(GSESeeuratObject, features = VariableFeatures(object = GSESeeuratObject))
  GSESeeuratObject <- FindNeighbors(GSESeeuratObject, dims = 1:20)
  GSESeeuratObject <- FindClusters(GSESeeuratObject, resolution = 0.8)
  #GSESeeuratObject <- RunUMAP(GSESeeuratObject, dims = 1:20)
  #DimPlot(GSESeeuratObject, reduction = "umap", label=TRUE)
  harmony_integrated <- RunHarmony(GSESeeuratObject1,group.by.vars="condition",plot_convergance=TRUE)
  harmony_integrated <- RunUMAP(harmony_integrated,reduction="harmony",dims=1:20)
  harmony_integrated <- FindNeighbors(harmony_integrated,reduction="harmony",dims=1:20)
  harmony_integrated <- FindClusters(harmony_integrated,resolution=0.8)
  DimPlot(harmony_integrated,reduction = "umap",label = TRUE)
  
  harmony_integrated <- RunTSNE(harmony_integrated)
  DimPlot(harmony_integrated,reduction = "tsne",label = TRUE)
  return(GSESeeuratObject)
}

cell_calculations <- function(seurat_integrated,name,path){
  conditions <- c("Normal","Tumour","TLN")
  final <- NULL
  for(i in 1:length(conditions)){
    print(i)
    print(conditions[i])
    Idents(seurat_integrated) <- seurat_integrated@meta.data$condition
    rorpos <- subset(seurat_integrated, idents=conditions[i])
    Idents(rorpos) <- rorpos@meta.data$main_annotation
    n_cells_rorpos <- FetchData(rorpos, vars = "ident") %>%
      dplyr::count(ident) %>%
      tidyr::spread(ident, n)
    
    
    new <- as.data.frame( t(n_cells_rorpos))
    
    col<- rownames(new)
    col2 <- new[,1]
    
    new[,1] <- col
    new[,2] <- col2
    
    colnames(new) <- c("Cell_Type", "Counts")
    sumz <- sum(new$Counts)
    new$percentage <- (new$Counts/sumz)*100
    
    write.xlsx(new,paste0(path,name,"_",conditions[i],"_","counts",".xlsx"))
    #return(new)#
    pie(new$Counts, labels = new$Cell_Type, main = paste0("Pie Char of Cell Counts ",conditions[i]))
  }
  
}


number_overall_degs_finder <- function(ROR2SeeuratObject, cluster,lfc=0.5,path){
  Idents(ROR2SeeuratObject) <- ROR2SeeuratObject$seurat_clusters
  marker <- FindMarkers(ROR2SeeuratObject, ident.1 = cluster, min.pct = 0.25, logfc.threshold = lfc)
  marker <- marker[marker$p_val<0.05,]
  marker$Names <- rownames(marker)
  marker <- marker[order(-marker$avg_log2FC), ]
  write.csv(marker, paste0(path,as.character(cluster),"-",".csv"))
  return(marker)
}

labeled_overall_degs_finder <- function(ROR2SeeuratObject, cluster,lfc=0.5,path,the_factor){
  Idents(ROR2SeeuratObject) <- the_factor
  marker <- FindMarkers(ROR2SeeuratObject, ident.1 = cluster, min.pct = 0.25, logfc.threshold = lfc)
  marker <- marker[marker$p_val<0.05,]
  marker$Names <- rownames(marker)
  marker <- marker[order(-marker$avg_log2FC), ]
  write.csv(marker, paste0(path,as.character(cluster),"-",".csv"))
  return(marker)
}


five_plots <- function(seurat.obj){
  #Idents(seurat.obj) <- seurat.obj$condition
  #seurat.obj <- subset(seurat.obj,idents=c("ROR-"))
  
  s.genes <- "Ccnd1, Cdk2,Cdk4,Cdk6, Cdk7,Ccne1,Ccne2, Myc,Orc1,Mcm6 ,Pcna,Ccna2,Ccnh"
  s.genes <-  unlist(strsplit(s.genes, split = ",", fixed = TRUE))
  
  g2m.genes <- "Wee1,Cdk1,Ccnf,Nusap1,Aurka,Ccna2,Ccnb2,Mki67,Top2a,Ccnb1"
  g2m.genes <- unlist(strsplit(g2m.genes, split = ",", fixed = TRUE))
  
  # Calculate cell cycle scores
  seurat.obj <- CellCycleScoring(seurat.obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  
  Idents(seurat.obj) <- seurat.obj$seurat_clusters
  metrix <- c("nCount_RNA","nFeature_RNA","percent.mt","S.Score","G2M.Score")
  
  FeaturePlot(seurat.obj, reduction = "umap", features = metrix,pt.size = 0.4,order = TRUE, 
              min.cutoff ='q10', label = TRUE)
  #Idents(seurat_integrated) <- seurat_integrated$customclassifions
  #DimPlot(seurat_integrated,reduction = "umap",label = TRUE, repel = TRUE, split.by = "condition")
  
  #Idents(seurat.obj) <- seurat.obj$seurat_clusters
  #DimPlot(seurat.obj,reduction = "umap",label = TRUE, repel = TRUE)
}


marker_to_heatmap_function <- function(s.obj,annotations){
  Idents(s.obj) <- annotations
  markers <- FindAllMarkers(s.obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  top_10 <- markers %>%group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)
  hmap <- DoHeatmap(s.obj, features = top_10$gene)
  return(hmap)
}

cell_calculations2 <- function(seurat_integrated,name,path,conditions,annotations){
  #conditions <- c("Normal","Tumour","TLN")
  final <- NULL
  conditions_count <- unique(conditions)
  df_list <- list()
  for(i in 1:length(conditions_count)){
    print(i)
    print(conditions_count[i])
    Idents(seurat_integrated) <- conditions
    rorpos <- subset(seurat_integrated, idents=conditions[i])
    Idents(rorpos) <- annotations
    n_cells_rorpos <- FetchData(rorpos, vars = "ident") %>%
      dplyr::count(ident) %>%
      tidyr::spread(ident, n)
    
    
    new <- as.data.frame( t(n_cells_rorpos))
    
    col<- rownames(new)
    col2 <- new[,1]
    
    new[,1] <- col
    new[,2] <- col2
    new[,3] <- conditions_count[i]
    colnames(new) <- c("Cell_Type", "Counts","Conditions")
    sumz <- sum(new$Counts)
    new$percentage <- round((new$Counts/sumz)*100)
    
    write.xlsx(new,paste0(path,name,"_",conditions_count[i],"_","counts",".xlsx"))
    #return(new)#
    p1 <- ggplot(new, aes(x = "", y = Counts, fill = Cell_Type)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar("y") +
      theme_void() +
      scale_fill_brewer(palette = "Set3") +
      labs(title = paste(Name, "Pie Chart")) +
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            legend.title = element_blank()) +
      geom_text(aes(label = paste0(Cell_Type, ": ", format(percentage, nsmall = 1), "%")), 
                position = position_stack(vjust = 0.5), size = 5, color = "white")
    p1
    ggsave(paste0("Pie Char of Cell Counts ",conditions_count[i],".jpeg"),p1)
    df_list[[i]] <- new
    
  }
  combined_df <- do.call(rbind, df_list)
  write.xlsx(combined_df,paste0(path,name,"_ Combine","_","counts",".xlsx"))
  return(combined_df)
}

#seurat.obj <- read_rds("D:/Li-Porject Data/Li-Single Cell Data/outputs/Master Markers/combined_ROR.Rds")

#name_vector is the one you want to subset the object on
#the_factor is the annotation column. e.g: seurat.obj$annotations
#name is the name you want to save the file with
#conditions are the conditions based on which you want to test the data
trajectory_making_function <- function(seurat.obj,name,name_vector,the_factor,conditions=c("ROR2+","ROR2-")){
  
  for (cond in conditions){
    Idents(seurat.obj) <- seurat.obj@meta.data$condition
    seurat_object <- subset(x=seurat.obj,idents=c(cond))
    seurat_object$customclassifions <- seurat_object@meta.data$the_factor
    
    Idents(seurat_object) <- seurat_object$customclassifions
    seurat_object <- subset(seurat_object, idents =  name_vector)
    
    seurat_object <- NormalizeData(seurat_object)
    seurat_object <- FindVariableFeatures(seurat_object)
    seurat_object <- ScaleData(seurat_object)
    seurat_object <- RunPCA(seurat_object)
    seurat_object <- FindNeighbors(seurat_object, dims = 1:30)
    seurat_object <- FindClusters(seurat_object, resolution = 0.9)
    seurat_object <- RunUMAP(seurat_object, dims = 1:30, n.neighbors = 50)
    
    
    lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
    
    db_ = "D:/Li-Porject Data/Li-Single Cell Data/outputs/T-Cell-subset/T-Subset-Annotation.xlsx"
    tissue = c("Immune system")
    
    gs_list = gene_sets_prepare(db_, tissue)
    
    
    es.max = sctype_score(scRNAseqData = seurat_object[["RNA"]]$scale.data, scaled = TRUE, 
                          gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
    
    
    cL_resutls = do.call("rbind", lapply(unique(seurat_object@meta.data$seurat_clusters), function(cl){
      es.max.cl = sort(rowSums(es.max[ ,rownames(seurat_object@meta.data[seurat_object@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
      head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_object@meta.data$seurat_clusters==cl)), 10)
    }))
    
    sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
    
    sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/5] = "Unknown"
    print(sctype_scores[,1:3])
    #write.xlsx(sctype_scores,"D:/Li-Porject Data/Li-Single Cell Data/outputs/sctype_scores1.xlsx")
    
    
    
    seurat_object@meta.data$customclassif = ""
    for(j in unique(sctype_scores$cluster)){
      cl_type = sctype_scores[sctype_scores$cluster==j,]; 
      seurat_object@meta.data$customclassifions[seurat_object@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
    }
    
    a<- DimPlot(seurat_object, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassifions')
    a1 <- DimPlot(seurat_object, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassifions')
    
    
    a1 <- DimPlot(seurat_object, reduction = 'umap', group.by = 'customclassifions', label = T)
    a2 <- DimPlot(seurat_object, reduction = 'umap', group.by = 'seurat_clusters', label = T)
    a1|a2
    ggsave2(path,name,"1",cond,".png", height = 14, width = 14)
    
    cds <- as.cell_data_set(seurat_object)
    fData(cds)$gene_short_name <- rownames(fData(cds))
    
    
    
    reacreate.partition <- c(rep(1,length(cds@colData@rownames)))
    names(reacreate.partition) <- cds@colData@rownames
    reacreate.partition <- as.factor(reacreate.partition)
    
    cds@clusters$UMAP$partitions <- reacreate.partition
    
    
    list_cluster <- seurat_object@active.ident
    cds@clusters$UMAP$clusters <- list_cluster
    
    
    cds@int_colData@listData$reducedDims$UMAP <- seurat_object@reductions$umap@cell.embeddings
    
    
    
    # plot
    
    #color_palette <- brewer.pal(115, "Set1")
    
    color_palette <- viridis(115)
    cluster.before.trajectory <- plot_cells(cds,
                                            color_cells_by = 'cluster',
                                            label_groups_by_cluster = FALSE,
                                            group_label_size = 5) +
      theme(legend.position = "right")
    
    cluster.names <- plot_cells(cds,
                                color_cells_by = "customclassifions",
                                label_groups_by_cluster = FALSE,
                                group_label_size = 5) +
      scale_color_manual(values = color_palette) +
      theme(legend.position = "right")
    
    cluster.before.trajectory | cluster.names
    ggsave2(path,name,"1",cond,".png", height = 14, width = 14)
    
    ####
    cds <- learn_graph(cds, use_partition = FALSE)
    
    plot_cells(cds,
               color_cells_by = 'customclassifions',
               label_groups_by_cluster = FALSE,
               label_branch_points = FALSE,
               label_roots = FALSE,
               label_leaves = FALSE,
               group_label_size = 5)
    ggsave2(path,name,"3",cond,".png",height = 14, width = 14)
    
    ####
    
    
    cds <- order_cells(cds, reduction_method = 'UMAP')
    
    plot_cells(cds,
               color_cells_by = 'pseudotime',
               label_groups_by_cluster = FALSE,
               label_branch_points = FALSE,
               label_roots = FALSE,
               label_leaves = FALSE)
    
    ggsave2(path,name,"4",cond,".png",height = 14, width = 14)
    
    pseudotime(cds)
    cds$monocle3_pseudotime <- pseudotime(cds)
    data.pseudo <- as.data.frame(colData(cds))
    
    ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(customclassifions, monocle3_pseudotime, median), fill = customclassifions)) +
      geom_boxplot()
    
    ggsave2(path,name,"5",cond,".png",height = 14, width = 14)
  }
}


cell_cell_communcation <- function(seurat_object,the_factor,signal_type){
  data.input <- GetAssayData(seurat_object, assay = "RNA", slot = "data") 
  
  Idents(seurat_object) <- the_factor
  labels <- Idents(seurat_object)
  
  meta <- data.frame(group = labels, row.names = names(labels)) 
  
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
  CellChatDB.human <- CellChatDB.human
  showDatabaseCategory(CellChatDB.human)
  CellChatDB.use <- CellChatDB.human
  CellChatDB.use <- subsetDB(CellChatDB.human, search = signal_type)
  
  interaction <- CellChatDB.human$interaction
  complex <- CellChatDB.human$complex
  cofactor <- CellChatDB.human$cofactor
  geneInfo <- CellChatDB.human$geneInfo
  
  options(stringsAsFactors = FALSE)
  CellChatDB.use <- subsetDB(CellChatDB.human, search = signal_type)
  
  
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat)
  
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  cellchat <- projectData(cellchat, PPI.human)
  
  cellchat <- computeCommunProb(cellchat) 
  
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  
  
  
  cellchat <- aggregateNet(cellchat)
  cellchat@net$count
  cellchat@net$weight
  
  return(cellchat)
}


subsetter_function <- function(merged.obj,the_factor,immune_types,epithelial_types){
  Idents(merged.obj) <- the_factor
  merged.obj <- subset(merged.obj,idents=c(immune_types,epithelial_types))
  return(merged.obj)
}
simple_circ_function <- function(cellchat,pathways.show){
  groupSize <- as.numeric(table(cellchat@idents))
  par(mfrow = c(1, 2), xpd=TRUE)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                   weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                   weight.scale = T, label.edge= T, title.name = "Interaction weights/strength")
  
  heatmap <- netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
  heatmap
  bubble <- netVisual_bubble(cellchat)
  bubble
  print(heatmap)
  print(bubble)
  plotGeneExpression(cellchat, signaling = pathways.show)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  gg1 <- netAnalysis_signalingRole_scatter(cellchat)
  gg1
  
  ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
  ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
  ht1 + ht2
}

cell_counts_ploting <- function(seurat_integrated, name, path, conditions, annotations) {
  # Initialize empty list to store dataframes
  df_list <- list()
  
  # Get unique conditions
  conditions_count <- unique(conditions)
  
  for (i in seq_along(conditions_count)) {
    print(i)
    print(conditions_count[i])
    
    # Set identities and subset data
    Idents(seurat_integrated) <- conditions
    rorpos <- subset(seurat_integrated, idents = conditions_count[i])
    Idents(rorpos) <- annotations
    
    # Fetch and process data
    n_cells_rorpos <- FetchData(rorpos, vars = "ident") %>%
      dplyr::count(ident) %>%
      tidyr::spread(ident, n)
    
    new <- as.data.frame(t(n_cells_rorpos))
    new$Cell_Type <- rownames(new)
    new$Counts <- new[,1]
    new$Conditions <- conditions_count[i]
    
    new <- new[, c("Cell_Type", "Counts", "Conditions")]
    
    sumz <- sum(as.numeric(new$Counts))
    new$percentage <- round((as.numeric(new$Counts) / sumz) * 100, 1)
    
    write.xlsx(new, file.path(path, paste0(name, "_", conditions_count[i], "_counts.xlsx")))
    
    # Create pie chart
    p1 <- ggplot(new, aes(x = "", y = Counts, fill = Cell_Type)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar("y") +
      theme_void() +
      scale_fill_brewer(palette = "Set3") +
      labs(title = paste(name,conditions[i], "Pie Chart")) +
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            legend.title = element_blank()) +
      geom_text(aes(label = paste0(Cell_Type, ": ", format(percentage, nsmall = 1), "%")), 
                position = position_stack(vjust = 0.5), size = 5, color = "white")
    
    file_name <- file.path(path, paste0("Pie_Char_of_", name, "_Cell_Counts_", conditions_count[i], ".jpeg"))
    ggsave(file_name, plot = p1, width = 8, height = 6, dpi = 300)
    
    
    df_list[[i]] <- new
  }
  
  combined_df <- do.call(rbind, df_list)
  
  write.xlsx(combined_df, file.path(path, paste0(name, "_Combined_counts.xlsx")))
  sbar <- ggplot(combined_df, aes(x = Cell_Type, y = Counts, fill = Conditions)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(x = "Cell_types", y = "Counts", fill = "Conditions", title = "Stacked Bar Chart") +
    theme(plot.title = element_text(hjust = 0.5))
  
  file_name <- file.path(path, paste0("Stacked_Bar_Char_of_", name, "_Cell_Counts_", conditions_count[i], ".jpeg"))
  ggsave(file_name, plot = sbar, width = 8, height = 6, dpi = 300)
  
  # I have to load the function marker_to_heatmap_function
  hmap <- marker_to_heatmap_function(seurat_integrated,annotations = annotations)
  file_name <- file.path(path, paste0("Heatmap_of_", name, ".jpeg"))
  ggsave(file_name, plot = hmap, width = 12, height = 8, dpi = 300)
  return(combined_df)
}

specific_labeled_degs_finder <- function(ROR2SeeuratObject, cluster,lfc=0.5,path,the_factor,idn1,idn2,merger=F,conditions){
  if(merger==T){
    ROR2SeeuratObject$joined <- paste0(the_factor,"_",conditions)
    Idents(ROR2SeeuratObject) <- ROR2SeeuratObject$joined
    print(unique(ROR2SeeuratObject$joined))
    cluster1 <- paste0(cluster,"_",idn1)
    cluster2 <- paste0(cluster,"_",idn2)
    marker <- FindMarkers(ROR2SeeuratObject, ident.1 = cluster1, ident.2 = cluster2, min.pct = 0.25, logfc.threshold = lfc)
  }else{
    ROR2SeeuratObject <- ROR2SeeuratObject
    Idents(ROR2SeeuratObject) <- the_factor
    cluster1 <- paste0(cluster,"_",idn1)
    cluster2 <- paste0(cluster,"_",idn2)
    marker <- FindMarkers(ROR2SeeuratObject, ident.1 = cluster1,ident.2 =cluster2, min.pct = 0.25, logfc.threshold = lfc)
  }
  marker <- marker[marker$p_val<0.05,]
  marker$Names <- rownames(marker)
  marker <- marker[order(-marker$avg_log2FC), ]
  write.csv(marker, paste0(path,as.character(cluster),"-",idn1,"-",idn2,".csv"))
  return(marker)
}
specific_labeled_degs_finder_2 <- function(ROR2SeeuratObject, cluster, lfc=0.5, path, the_factor, idn1, idn2, merger=FALSE, conditions) {
  if(merger == TRUE) {
    ROR2SeeuratObject$joined <- paste0(the_factor, "_", conditions)
    Idents(ROR2SeeuratObject) <- ROR2SeeuratObject$joined
    print(unique(ROR2SeeuratObject$joined))
    cluster1 <- paste0(cluster, "_", idn1)
    cluster2 <- paste0(cluster, "_", idn2)
  } else {
    Idents(ROR2SeeuratObject) <- the_factor
    cluster1 <- paste0(cluster, "_", idn1)
    cluster2 <- paste0(cluster, "_", idn2)
  }
  
  # Check if both clusters have at least 3 cells
  cells_cluster1 <- WhichCells(ROR2SeeuratObject, idents = cluster1)
  cells_cluster2 <- WhichCells(ROR2SeeuratObject, idents = cluster2)
  
  if (length(cells_cluster1) < 3 || length(cells_cluster2) < 3) {
    warning(paste("Skipping", cluster1, "vs", cluster2, "- fewer than 3 cells in one of the groups"))
    return(NULL)
  }
  
  # Perform differential expression analysis if cell count condition is met
  marker <- FindMarkers(
    ROR2SeeuratObject, ident.1 = cluster1, ident.2 = cluster2,
    min.pct = 0.25, logfc.threshold = lfc
  )
  
  # Filter markers by p-value and sort by avg_log2FC
  marker <- marker[marker$p_val < 0.05, ]
  marker$Names <- rownames(marker)
  marker <- marker[order(-marker$avg_log2FC), ]
  
  # Save results to CSV
  write.csv(marker, paste0(path, as.character(cluster), "-", idn1, "-", idn2, ".csv"))
  
  return(marker)
}
gene_exp_heatmap <- function(s.obj,condition_inf,path,name,row.title="Conditions"){
  expression_data <- GetAssayData(s.obj, assay = "RNA", layer = "data")
  #condition_inf <- pre_B$condition
  #selected_genes <- unlist(CC_list)
  #
  selected_genes <- genes_to_plot
  expression_filtered <- expression_data[rownames(expression_data)%in% selected_genes,]
  matrix_grouped <- aggregate(t(expression_filtered), by = list(condition_inf), FUN = mean)
  rownames(matrix_grouped) <- matrix_grouped$Group.1
  matrix_grouped <- matrix_grouped[,-1]
  hmap <- Heatmap(matrix_grouped, 
                  name = "Expression", 
                  heatmap_legend_param = list(title = "Expression"),
                  column_title = "Genes",
                  row_title = row.title)
  ggsave2(paste0(path,name,".jpeg"))
}

gsea_preprocess <- function(path){
  degs <- read_csv(path)
  degs <- na.omit(degs)
  rownames(degs) <- degs$Names
  degs <- degs[order(-degs$avg_log2FC), ]
  
  return(degs)
}

DEGs_function <- function(degs_folder_path,out_main_folder_path){
  file_paths <-  list.files(path = degs_folder_path, full.names = TRUE)
  for (file in file_paths) {
    tryCatch({
      Degs <- gsea_preprocess(file)
      initial <- tail(unlist(strsplit(file, "/")), n=1)
      Name <- sub("\\.csv$", "", initial)
      main_Name <- tail(unlist(strsplit(degs_folder_path, "/")), n=1)
      
      print(paste("Running", Name, "KEGG"))
      dir.create(paste0(out_main_folder_path, main_Name), showWarnings = FALSE)
      
      kegg_simp_plotter_human2(
        instance_markers = Degs,
        name = Name,
        path = paste0(out_main_folder_path, main_Name),
        ctgry = 13
      )
      
      print(paste("Running", Name, "Wiki"))
      wp_simp_plotter_human2(
        instance_markers = Degs,
        name = Name,
        path = paste0(out_main_folder_path, main_Name),
        ctgry = 13
      )
    }, error = function(e) {
      message(paste("Error processing file:", file, "\nSkipping to next file."))
      message("Error details:", e)
    })
  }
}


DEGs_function_OR <- function(degs_folder_path, out_main_folder_path) {
  file_paths <- list.files(path = degs_folder_path, full.names = TRUE)
  
  for (file in file_paths) {
    tryCatch({
      Degs <- gsea_preprocess(file)
      initial <- tail(unlist(strsplit(file, "/")), n = 1)
      Name <- sub("\\.csv$", "", initial)
      main_Name <- tail(unlist(strsplit(degs_folder_path, "/")), n = 1)
      
      print(paste("Running", Name, "KEGG"))
      dir.create(paste0(out_main_folder_path, main_Name), showWarnings = FALSE)
      
      kegg_OR_GSEA(
        instance_markers = Degs,
        name = Name,
        path = paste0(out_main_folder_path, main_Name,"/"),
        ctgry = 13
      )
      
      print(paste("Running", Name, "Wiki"))
      wp_result <- WP_OR_GSEA(
        instance_markers = Degs,
        name = Name,
        path = paste0(out_main_folder_path, main_Name,"/"),
        ctgry = 13
      )
      Reactome_OR_GSEA(instance_markers = Degs,
                       name = Name,
                       path = paste0(out_main_folder_path, main_Name,"/"),
                       ctgry = 13)
      
      # Ensure there are results before proceeding
      if (nrow(wp_result@result) > 0) {
        # Proceed with operations on the results
      } else {
        message(paste("No significant pathways found for file:", file))
      }
      
    }, error = function(e) {
      message(paste("Error processing file:", file, "\nSkipping to next file."))
      message("Error details:", e)
    })
  }
}



DEGs_function_OR <- function(degs_folder_path, out_main_folder_path) {
  file_paths <- list.files(path = degs_folder_path, full.names = TRUE)
  
  for (file in file_paths) {
    tryCatch({
      Degs <- gsea_preprocess(file)
      initial <- tail(unlist(strsplit(file, "/")), n = 1)
      Name <- sub("\\.csv$", "", initial)
      main_Name <- tail(unlist(strsplit(degs_folder_path, "/")), n = 1)
      
      print(paste("Running", Name, "KEGG"))
      dir.create(paste0(out_main_folder_path, main_Name), showWarnings = FALSE)
      
      kegg_simp_plotter_human(
        instance_markers = Degs,
        name = Name,
        path = paste0(out_main_folder_path, main_Name),
        ctgry = 13
      )
      
      print(paste("Running", Name, "Wiki"))
      wp_result <- wp_simp_plotter_human(
        instance_markers = Degs,
        name = Name,
        path = paste0(out_main_folder_path, main_Name),
        ctgry = 13
      )
      Reactome_OR_GSEA(instance_markers = Degs,
                       name = Name,
                       path = paste0(out_main_folder_path, main_Name),
                       ctgry = 13)
      # Ensure there are results before proceeding
      if (nrow(wp_result@result) > 0) {
        # Proceed with operations on the results
      } else {
        message(paste("No significant pathways found for file:", file))
      }
      
    }, error = function(e) {
      message(paste("Error processing file:", file, "\nSkipping to next file."))
      message("Error details:", e)
    })
  }
}


degs_scatterplot <- function(data,threshold=1.5){
  data <- na.omit(data)
  data$p_val[data$p_val == 0] <- 1e-300
  data$logP <- -log10(data$p_val)
  
  data$significance <- ifelse(data$p_val < 0.05 & abs(data$avg_log2FC) > threshold, 
                              ifelse(data$avg_log2FC > 0, "Upregulated", "Downregulated"), 
                              "Not Significant")
  significant_data <- data[data$p_val < 0.05 & abs(data$avg_log2FC) > threshold, ]
  top_Namess <- significant_data[order(significant_data$p_val), ]
  # Select top 10 Names for labeling
  top_10_Namess <- top_Namess[1:min(10, nrow(top_Namess)), ]  # Ensure not to exceed available rows
  
  # Plotting
  gplot<-ggplot(data, aes(x = avg_log2FC, y = logP)) +
    geom_point(aes(color = significance), alpha = 0.7, size = 3) +  # Plot points
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
    theme_minimal() +
    labs(x = "log2 Fold Change", y = "-log10(p-value)", title = "Differential Volcano Plot") +
    geom_text_repel(data = top_10_Namess, aes(label = Names), size = 4, box.padding = 0.5) +  # Label top 10 Names
    theme(legend.position = "top") +
    geom_vline(xintercept = c(-threshold, threshold), linetype = "dashed") +  # Vertical lines for logFC threshold
    geom_hline(yintercept = -log10(0.05), linetype = "dashed")
  print(gplot)
  return(gplot)
}


gene_expression_bar_plotter <- function(s.obj, annotatations, path, gene_list,Name){
  table(s.obj$sub_annotations)
  Idents(s.obj) <- s.obj$sub_annotations
  
  genes <- as.character(unlist(gene_list))
  
  # df to plot
  df.plot <- as.data.frame(t(as.matrix(GetAssayData(s.obj, slot = "count"))[genes, ]))
  head(df.plot)
  df.plot$group <- Idents(s.obj)
  
  ## plot
  for (i in 1:length(genes)){
    ggbarplot(df.plot, x = "group", y=genes[1], add="mean_se",fill="group",  position=position_dodge(0.8),
              palette = pal_nejm()(6)) +
      xlab("Cell type") +
      ylab("Expression") +
      ggtitle(genes[i]) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
      guides(fill = F)
    ggsave(paste0(path,Name,"_",genes[i],"_barplot.pdf"), width = 4, height = 4)
  }
}

########################################################################################################
#END OF FUNCTIONS
########################################################################################################

#######################################################################################################
#START OF GSEA FUNCTIONS
######################################################################################################

wp_simp_plotter_human2 <- function(instance_markers, name, path, ctgry) {
  gene_list <- instance_markers$avg_log2FC
  names(gene_list) <- rownames(instance_markers)
  gene_symbols <- instance_markers$Names
  entrez_ids <- mapIds(org.Hs.eg.db, keys = gene_symbols, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  
  unamed_entrez_ids <- unname(entrez_ids)
  newglist <- setNames(instance_markers$avg_log2FC, unamed_entrez_ids)
  
  wp_result <- gseWP(newglist, organism = "Homo sapiens", pvalueCutoff = 0.05)
  
  if (is.null(wp_result) || nrow(wp_result@result) == 0) {
    message(paste("No significant pathways found for", name, "- skipping file."))
    return(NULL)
  }
  
  # Set readable names for the pathways
  wp_result <- setReadable(wp_result, 'org.Hs.eg.db', 'ENTREZID')
  
  # Save the complete result to CSV
  write.csv(as.data.frame(wp_result), paste0(path, "/WikiPathways-", name, "-Complete.csv"))
  
  dp <- dotplot(wp_result, showCategory = ctgry, font.size = 14)
  print(dp)
  cnp <- cnetplot(wp_result, showCategory = ctgry, foldChange = newglist) + 
    scale_color_gradient2(name = 'associated data', low = 'darkgreen', high = 'firebrick')
  print(cnp)
  hp <- heatplot(wp_result, showCategory = ctgry, foldChange = newglist)
  print(hp)
  
  ggsave2(paste0(path, "/WikiPathways-", name, "-DotPlot.jpeg"), plot = dp)
  ggsave2(paste0(path, "/WikiPathways-", name, "-CNetPlot.jpeg"), plot = cnp)
  ggsave2(paste0(path, "/WikiPathways-", name, "-HeatMap.jpeg"), plot = hp)
  
  dfr <- wp_result@result
  dfr <- dfr[dfr$qvalue < 0.25, ]
  write.csv(dfr, paste0(path, "/WikiPathways-", name, "-qvalue-filtered.csv"))
  
  if (nrow(dfr) > 0) {
    barp <- ggplot(dfr, aes(x = reorder(Description, -qvalue), y = qvalue, fill = qvalue)) +
      geom_bar(stat = "identity") +
      labs(x = "Pathways", y = "q-value", title = "Pathway Analysis") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_fill_gradient(low = "blue", high = "red")
    
    ggsave2(paste0(path, "/WikiPathways-", name, "-BarPlot.jpeg"), plot = barp)
  }
}
rec_simp_plotter_human2 <- function(instance_markers,name,path,ctgry){
  gene_list <- instance_markers$avg_log2FC
  names(gene_list) <- rownames(instance_markers)
  gene_symbols <- instance_markers$Names
  entrez_ids <- mapIds(org.Hs.eg.db, keys = gene_symbols, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  
  unamed_entrez_ids <- unname(entrez_ids)
  newglist <- setNames(instance_markers$avg_log2FC,unamed_entrez_ids)
  
  rec_result <- gsePathway(newglist, organism = "human",pvalueCutoff = 0.05)
  # Check if there are any enriched pathways
  if (is.null(rec_result) || nrow(rec_result@result) == 0) {
    message(paste("No significant pathways found for", name, "- skipping file."))
    return(NULL)  # Skip further processing
  }
  rec_result <- setReadable(rec_result, 'org.Hs.eg.db', 'ENTREZID')
  
  write.csv(as.data.frame(rec_result@result),paste0(path,"Reactome-",name,"-Complete.csv"))
  
  dp <- dotplot(rec_result, showCategory = ctgry, font.size = 14)
  print(dp)
  #bp <- barplot(rec_result, showCategory = ctgry, font.size = 14)
  #print(bp)
  cnp <- cnetplot(rec_result, showCategory=ctgry, foldChange=newglist) + 
    scale_color_gradient2(name='associated data', low='darkgreen', high='firebrick')
  print(cnp)
  hp <- heatplot(rec_result,  showCategory=ctgry, foldChange=newglist)
  print(hp)
  
  ggsave2(paste0(path,"Reactome-",name,"-DotPlot.jpeg"),plot = dp)
  #ggsave2(paste0(path,"Reactome-",name,"-BarPlot.jpeg"),plot = bp)
  ggsave2(paste0(path,"Reactome-",name,"-CNetPlot.jpeg"),plot = cnp)
  ggsave2(paste0(path,"Reactome-",name,"-HeatMap.jpeg"),plot = hp)
  write.csv(as.data.frame(rec_result),paste0(path,"Reactome-",name,"-Filtered.csv"))
  
  dfr <- rec_result@result
  dfr <- dfr[dfr$qvalue<0.25,]
  barp <- ggplot(dfr, aes(x = reorder(Description, -qvalue), y = qvalue, fill = qvalue)) +
    geom_bar(stat = "identity") +
    labs(x = "Pathways", y = "q-value", title = "Pathway Analysis") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_gradient(low = "blue", high = "red")
  
  write.csv(as.data.frame(rec_result),paste0(path,"Reactome-",name,"-qvalue-filtered.csv"))
  ggsave2(paste0(path,"Reactome-",name,"-qvlaue.jpeg"),plot = barp)
}



kegg_simp_plotter_human2 <- function(instance_markers, name, path, ctgry) {
  gene_list <- instance_markers$avg_log2FC
  names(gene_list) <- rownames(instance_markers)
  gene_symbols <- instance_markers$Names
  entrez_ids <- mapIds(org.Hs.eg.db, keys = gene_symbols, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  
  unamed_entrez_ids <- unname(entrez_ids)
  newglist <- setNames(instance_markers$avg_log2FC, unamed_entrez_ids)
  

  kegg_result <- gseKEGG(gene = newglist, organism = 'human', pvalueCutoff = 0.05)
  
  if (is.null(kegg_result) || nrow(kegg_result@result) == 0) {
    message(paste("No significant KEGG pathways found for", name, "- skipping file."))
    return(NULL)  
  }
  
  kegg_result <- setReadable(kegg_result, 'org.Hs.eg.db', 'ENTREZID')
  
  write.csv(as.data.frame(kegg_result@result), paste0(path, "/KEGG-", name, "-Complete.csv"))
  
  dp <- dotplot(kegg_result, showCategory = ctgry, font.size = 14)
  print(dp)
  cnp <- cnetplot(kegg_result, showCategory = ctgry, foldChange = newglist) + 
    scale_color_gradient2(name = 'associated data', low = 'darkgreen', high = 'firebrick')
  print(cnp)
  hp <- heatplot(kegg_result, showCategory = ctgry, foldChange = newglist)
  print(hp)
  
  ggsave2(paste0(path, "/KEGG-", name, "-DotPlot.jpeg"), plot = dp)
  ggsave2(paste0(path, "/KEGG-", name, "-CNetPlot.jpeg"), plot = cnp)
  ggsave2(paste0(path, "/KEGG-", name, "-HeatMap.jpeg"), plot = hp)
  
  # Filter results based on q-value and save filtered results
  dfr <- kegg_result@result
  dfr <- dfr[dfr$qvalue < 0.25, ]
  write.csv(dfr, paste0(path, "KeggPathway-", name, "-qvalue-filtered.csv"))
  
  # Generate and save bar plot of filtered results
  if (nrow(dfr) > 0) {
    barp <- ggplot(dfr, aes(x = reorder(Description, -qvalue), y = qvalue, fill = qvalue)) +
      geom_bar(stat = "identity") +
      labs(x = "Pathways", y = "q-value", title = "Pathway Analysis") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_fill_gradient(low = "blue", high = "red")
    
    ggsave2(paste0(path, "KeggPathway-", name, "-qvalue.jpeg"), plot = barp)
  }
}

#Over representation Analysis
#################################################################################################################

kegg_OR_GSEA <- function(instance_markers,name,path,ctgry){
  instance_markers <- instance_markers[order(-instance_markers$avg_log2FC),]
  gene_symbols <- instance_markers$Names
  entrez_ids <- mapIds(org.Hs.eg.db, keys = gene_symbols, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  kegg_result <- enrichKEGG(gene = entrez_ids, organism = 'human', pvalueCutoff = 0.05)
  kegg_result <- setReadable(kegg_result, 'org.Hs.eg.db', 'ENTREZID')
  View()
  if (is.null(kegg_result) || nrow(kegg_result@result) == 0) {
    message(paste("No significant KEGG pathways found for", name, "- skipping file."))
    return(NULL)  
  }
  #gene_symbols <- instance_markers$Names
  gene_list <- instance_markers$avg_log2FC
  names(gene_list) <- rownames(instance_markers)
  plot2 <- dotplot(kegg_result,showCategory=ctgry)
  ggsave(paste0(path,"KEGG_OR_",name,".png"),plot=plot2)
  #barplot2 <- barplot(kegg_result, showCategory=ctgry)
  #ggsave(paste0(path,"Bar_pot_","KEGG_",name,".png"),plot=barplot2)
  heatmap2 <- heatplot(kegg_result, foldChange=gene_list, showCategory=ctgry)
  ggsave(paste0(path,"Heatmap_pot_OR_","KEGG_",name,".jpeg"),plot=heatmap2)
  
  edox <- setReadable(kegg_result, 'org.Hs.eg.db', 'ENTREZID')
  p2 <- cnetplot(edox, color.params = list(foldChange=gene_list),showCategory=ctgry)
  ggsave(paste0(path,"Cnet_plot_","KEGG_OR_",name,".jpeg"),plot=p2)
  
  enrich_df <- as.data.frame(kegg_result)
  enrich_df <- enrich_df[enrich_df$p.adjust<0.05,]
  #enrich_df <- enrich_df[enrich_df$Description %in% ctgry]
  dp2 <- ggplot(enrich_df, aes(x = GeneRatio, y = reorder(Description, GeneRatio), 
                               size = Count, color = -log10(p.adjust))) +
    geom_point() +
    scale_size_continuous(range = c(1, 4)) + 
    scale_color_gradientn(colors = c("blue", "yellow", "red")) + 
    theme_minimal() +
    labs(size = "Number", color = "-log10(Pvalue)", x = "GeneRatio", y = "Terms", 
         title = paste0("KEGG_OR_",name)) +
    theme(axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 12, face = "bold"))
  ggsave(paste0(path,"New_Dot_pot_OR_","KEGG_",name,".jpeg"),plot=dp2)  
  #edo_heat2 <- enrichDGN(entrez_ids)
  #edox_heat2 <- setReadable(edo_heat2, 'org.Hs.eg.db', 'ENTREZID')
  #p_heat2 <- heatplot(edox_heat2, foldChange=gene_list11, showCategory=ctgry)
  #ggsave(paste0(path,"Real_Heatmap_","KEGG_",name,".jpeg"),plot=p_heat2)
  
  write.xlsx(as.data.frame(kegg_result),paste0(path,"KeggPathway_OR_",name,"-Complete_DataFrame.xlsx"))
}


gsea_preprocess <- function(path){
  degs <- read_csv(path)
  degs <- na.omit(degs)
  rownames(degs) <- degs$Names
  degs <- degs[order(-degs$avg_log2FC), ]
  
  return(degs)
}



Reactome_OR_GSEA <- function(instance_markers,name,path,ctgry){
  instance_markers <- instance_markers[order(-instance_markers$avg_log2FC),]
  gene_symbols <- instance_markers$Names
  entrez_ids <- mapIds(org.Hs.eg.db, keys = gene_symbols, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  
  reac_result <- enrichPathway(gene = entrez_ids, organism = "human", pvalueCutoff = 0.05)
  
  if (is.null(reac_result) || nrow(reac_result@result) == 0) {
    message(paste("No significant Reactome pathways found for", name, "- skipping file."))
    return(NULL)  # Skip further processing
  }
  reac_result <- setReadable(reac_result, 'org.Hs.eg.db', 'ENTREZID')
  #gene_symbols <- instance_markers$Names
  gene_list <- instance_markers$avg_log2FC
  names(gene_list) <- rownames(instance_markers)
  plot2 <- dotplot(reac_result)
  ggsave(paste0(path,"Reactome_OR_",name,".png"),plot=plot2)
  print("plot 1 is saved")
  #barplot2 <- barplot(reac_result, showCategory=ctgry)
  #ggsave(paste0(path,"Bar_pot_","Reactome_",name,".png"),plot=barplot2)
  print("plot 2 bar is saved")
  heatmap2 <- heatplot(reac_result, foldChange=gene_list, showCategory=ctgry)
  ggsave(paste0(path,"Heatmap_pot_","Reactome_OR_",name,".jpeg"),plot=heatmap2)
  
  edox <- setReadable(reac_result, 'org.Hs.eg.db', 'ENTREZID')
  p2 <- cnetplot(edox, color.params = list(foldChange=gene_list))
  ggsave(paste0(path,"Cnet_plot_","Reactome_OR_",name,".jpeg"),plot=p2)
  
  enrich_df <- as.data.frame(reac_result)
  enrich_df <- enrich_df[enrich_df$p.adjust<0.05,]
  #enrich_df <- enrich_df[enrich_df$Description %in% ctgry]
  dp2 <- ggplot(enrich_df, aes(x = GeneRatio, y = reorder(Description, GeneRatio), 
                               size = Count, color = -log10(p.adjust))) +
    geom_point() +
    scale_size_continuous(range = c(1, 4)) + 
    scale_color_gradientn(colors = c("blue", "yellow", "red")) + 
    theme_minimal() +
    labs(size = "Number", color = "-log10(Pvalue)", x = "GeneRatio", y = "Terms", 
         title = paste0("Reactome_OR_",name)) +
    theme(axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 12, face = "bold"))
  ggsave(paste0(path,"New_Dot_pot_OR_","Reactome_",name,".jpeg"),plot=dp2) 
  #edo_heat2 <- enrichDGN(entrez_ids)
  #edox_heat2 <- setReadable(edo_heat2, 'org.Hs.eg.db', 'ENTREZID')
  #p_heat2 <- heatplot(edox_heat2, foldChange=gene_list11, showCategory=ctgry)
  #ggsave(paste0(path,"Real_Heatmap_","Reactome_",name,".jpeg"),plot=p_heat2)
  
  write.xlsx(as.data.frame(reac_result),paste0(path,"ReactomePathway-",name,"-Complete_DataFrame.xlsx"))
}




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
  ggsave(paste0(path,"WP_OR_",name,".png"),plot=plot2)
  #barplot2 <- barplot(wp_result, showCategory=ctgry)
  #ggsave(paste0(path,"Bar_pot_","WP_",name,".png"),plot=barplot2)
  heatmap2 <- heatplot(wp_result, foldChange=gene_list, showCategory=ctgry)
  ggsave(paste0(path,"Heatmap_pot_","WP_OR_",name,".jpeg"),plot=heatmap2)
  
  edox <- setReadable(wp_result, 'org.Hs.eg.db', 'ENTREZID')
  p2 <- cnetplot(edox, color.params = list(foldChange=gene_list))
  ggsave(paste0(path,"Cnet_plot_","WP_OR_",name,".jpeg"),plot=p2)
  
  #edo_heat2 <- enrichDGN(entrez_ids)
  
  #edox_heat2 <- setReadable(edo_heat2, 'org.Hs.eg.db', 'ENTREZID')
  #p_heat2 <- heatplot(edox_heat2, foldChange=gene_list11, showCategory=ctgry)
  #ggsave(paste0(path,"Real_Heatmap_","WP_",name,".jpeg"),plot=p_heat2)
  
  enrich_df <- as.data.frame(wp_result)
  enrich_df <- enrich_df[enrich_df$p.adjust<0.05,]
  #enrich_df <- enrich_df[enrich_df$Description %in% ctgry]
  dp2 <- ggplot(enrich_df, aes(x = GeneRatio, y = reorder(Description, GeneRatio), 
                               size = Count, color = -log10(p.adjust))) +
    geom_point() +
    scale_size_continuous(range = c(1, 4)) + 
    scale_color_gradientn(colors = c("blue", "yellow", "red")) + 
    theme_minimal() +
    labs(size = "Number", color = "-log10(Pvalue)", x = "GeneRatio", y = "Terms", 
         title = paste0("Wiki_OR_",name)) +
    theme(axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 12, face = "bold"))
  ggsave(paste0(path,"New_Dot_pot_OR_","Wiki_",name,".jpeg"),plot=dp2) 
  
  write.xlsx(as.data.frame(wp_result),paste0(path,"WikiPathway_OR_",name,"-Complete_DataFrame.xlsx"))
  
}

##################################################################################################################
######################Pahse Implemenataion codes
######################################################################################################################
DEGs_function <- function(degs_folder_path,out_main_folder_path){
  file_paths <-  list.files(path = degs_folder_path, full.names = TRUE)
  for (file in file_paths) {
    tryCatch({
      Degs <- gsea_preprocess(file)
      initial <- tail(unlist(strsplit(file, "/")), n=1)
      Name <- sub("\\.csv$", "", initial)
      main_Name <- tail(unlist(strsplit(degs_folder_path, "/")), n=1)
      
      print(paste("Running", Name, "KEGG"))
      dir.create(paste0(out_main_folder_path, main_Name), showWarnings = FALSE)
      
      kegg_simp_plotter_human2(
        instance_markers = Degs,
        name = Name,
        path = paste0(out_main_folder_path, main_Name),
        ctgry = 13
      )
      
      print(paste("Running", Name, "Wiki"))
      wp_simp_plotter_human2(
        instance_markers = Degs,
        name = Name,
        path = paste0(out_main_folder_path, main_Name),
        ctgry = 13
      )
      print(paste("Running", Name, "Reactome"))
      rec_simp_plotter_human2(
        instance_markers = Degs,
        name = Name,
        path = paste0(out_main_folder_path, main_Name),
        ctgry = 13
      )
    }, error = function(e) {
      message(paste("Error processing file:", file, "\nSkipping to next file."))
      message("Error details:", e)
    })
  }
}


DEGs_function_OR <- function(degs_folder_path, out_main_folder_path) {
  file_paths <- list.files(path = degs_folder_path, full.names = TRUE)
  
  for (file in file_paths) {
    tryCatch({
      Degs <- gsea_preprocess(file)
      initial <- tail(unlist(strsplit(file, "/")), n = 1)
      Name <- sub("\\.csv$", "", initial)
      main_Name <- tail(unlist(strsplit(degs_folder_path, "/")), n = 1)
      
      print(paste("Running", Name, "KEGG"))
      dir.create(paste0(out_main_folder_path, main_Name), showWarnings = FALSE)
      kegg_OR_GSEA(
        instance_markers = Degs,
        name = Name,
        path = paste0(out_main_folder_path, main_Name,"/"),
        ctgry = 13
      )
      print(paste("Running", Name, "Wiki"))
      WP_OR_GSEA(
        instance_markers = Degs,
        name = Name,
        path = paste0(out_main_folder_path, main_Name,"/"),
        ctgry = 13
      )
      Reactome_OR_GSEA(instance_markers = Degs,
                       name = Name,
                       path = paste0(out_main_folder_path, main_Name,"/"),
                       ctgry = 13)
      
      if (nrow(wp_result@result) > 0) {
      } else {
        message(paste("No significant pathways found for file:", file))
      }
      
    }, error = function(e) {
      message(paste("Error processing file:", file, "\nSkipping to next file."))
      message("Error details:", e)
    })
  }
}