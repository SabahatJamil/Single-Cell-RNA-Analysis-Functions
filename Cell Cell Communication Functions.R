#Idents(GSE11_immune) <- GSE11_immune$condition
#tln <- subset(GSE11_immune,idents=c("TLN"))
#tumor <- subset(GSE11_immune,idents=c("Tumor"))

#cellchat_pos <- cell_cell_communcation(tln,the_factor = tln$Manual_Annotations)
#cellchat_neg <- cell_cell_communcation(tumor,the_factor = tumor$Manual_Annotations)

#cellchat_pos <- updateCellChat(cellchat_pos)
#cellchat_neg <- updateCellChat(cellchat_neg)
#object.list <- list(TLN = cellchat_pos, TUMOR = cellchat_neg)
#cellchat <- mergeCellChat(object.list, add.names = names(object.list))

cellchat_comparison1 <- function(cellchat){
  gg1 <- netVisual_heatmap(cellchat)
  gg2 <- netVisual_heatmap(cellchat, measure = "weight")
  print(gg1 + gg2)
  
  gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
  gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
  gg1 + gg2
}
cellchat2_creator<- function(s.obj, id1,id2,the_factor){
  Idents(s.obj) <- the_factor
  tln <- subset(s.obj,idents=c(id1))
  tumor <- subset(s.obj,idents=c(id2))
  
  cellchat_pos <- cell_cell_communcation(tln,the_factor = tln$sub_annotations,signal_type = "Cell-Cell Contact")
  cellchat_neg <- cell_cell_communcation(tumor,the_factor = tumor$sub_annotations,signal_type = "Cell-Cell Contact")
  
  cellchat_pos <- updateCellChat(cellchat_pos)
  cellchat_neg <- updateCellChat(cellchat_neg)
  object.list <- list(TLN = cellchat_pos, TUMOR = cellchat_neg)
  cellchat <- mergeCellChat(object.list, add.names = names(object.list))
  
  weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(object.list)) {
    netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
  }
  return(cellchat)
}



cellchat_visiulization1 <-function(s_cellchat){
  groupSize <- as.numeric(table(s_cellchat@idents))
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_circle(s_cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(s_cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  
  
  par(mfrow=c(1,1))
  netVisual_aggregate(s_cellchat, signaling = s_cellchat@netP$pathways, layout = "chord")
  netVisual_aggregate(s_cellchat, signaling = s_cellchat@netP$pathways[1], layout = "chord")
  netVisual_aggregate(s_cellchat, signaling = s_cellchat@netP$pathways[2], layout = "chord")
  netVisual_aggregate(s_cellchat, signaling = s_cellchat@netP$pathways[3], layout = "chord")
  
  
  netAnalysis_contribution(s_cellchat, signaling = s_cellchat@netP$pathways)
  
  plotGeneExpression(s_cellchat, signaling = s_cellchat@netP$pathways[1])
  plotGeneExpression(s_cellchat, signaling = s_cellchat@netP$pathways[2])
  plotGeneExpression(s_cellchat, signaling = s_cellchat@netP$pathways[3])
}
dysregulated_interactions <- function(cellchat,condition1,condition2){
  pos.dataset = condition1
  # define a char name used for storing the results of differential expression analysis
  features.name = pos.dataset
  # perform differential expression analysis
  cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
  #> Use the joint cell labels from the merged CellChat object
  # map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
  net <- netMappingDEG(cellchat, features.name = features.name)
  # extract the ligand-receptor pairs with upregulated ligands in LS
  net.up <- subsetCommunication(cellchat, net = net, datasets = condition1,ligand.logFC = 0.2, receptor.logFC = NULL)
  # extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
  net.down <- subsetCommunication(cellchat, net = net, datasets = condition2,ligand.logFC = -0.1, receptor.logFC = -0.1)
  
  gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
  gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
  
  pairLR.use.up = net.up[, "interaction_name", drop = F]
  gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 4, targets.use = c(2:4), comparison = c(1, 2),  angle.x = 90,
                          remove.isolate = T,title.name = paste0("Up-regulated signaling in ", condition1))
  #> Comparing communications on a merged object
  pairLR.use.down = net.down[, "interaction_name", drop = F]
  gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 4, targets.use = c(2:4), comparison = c(1, 2),  angle.x = 90, 
                          remove.isolate = T,title.name = paste0("Down-regulated signaling in ", condition1))
  #> Comparing communications on a merged object
  pplot <- gg1 + gg2
  print(pplot)
  
}
dysregulated_interactions(cellchat_BT2,condition1 = "TLN", condition2 = "TUMOR")

dotplot_ST_plotter <- function(cellchat,sources, targets, path){
  for (source in sources){
    p<- netVisual_bubble(cellchat, sources.use = source, targets.use = targets)
    ggsave2(paste0(path,source,".jpeg"),plot = p)
  }
}

dotplot_ST_plotter2 <- function(cellchat,sources, targets, path){
  for (source in sources){
    p<- netVisual_bubble(cellchat,comparison = c(1, 2), sources.use = source, targets.use = targets)
    ggsave2(paste0(path,source,"_combined.jpeg"),plot = p)
  }
}
####################################################################
library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
#################################################################
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)

imm_ecm <-cell_cell_communcation(GSE11_immune,the_factor = GSE11_immune$Manual_Annotations,signal_type = "ECM-Receptor")


cellchat_pos1 <- cell_cell_communcation(tln,the_factor = tln$Manual_Annotations,signal_type = "ECM-Receptor")
cellchat_neg1 <- cell_cell_communcation(tumor,the_factor = tumor$Manual_Annotations,signal_type = "ECM-Receptor")

cellchat_pos1 <- updateCellChat(cellchat_pos1)
cellchat_neg1 <- updateCellChat(cellchat_neg1)

object.list1 <- list(TLN = cellchat_pos1, Tumor = cellchat_neg1)
cellchat_ecm <- mergeCellChat(object.list1, add.names = names(object.list))

imm_ccc <-cell_cell_communcation(GSE11_immune,the_factor = GSE11_immune$Manual_Annotations,signal_type = "Cell-Cell Contact")


cellchat_pos2 <- cell_cell_communcation(tln,the_factor = tln$Manual_Annotations,signal_type = "Cell-Cell Contact")
cellchat_neg2 <- cell_cell_communcation(tumor,the_factor = tumor$Manual_Annotations,signal_type = "Cell-Cell Contact")

cellchat_pos2 <- updateCellChat(cellchat_pos2)
cellchat_neg2 <- updateCellChat(cellchat_neg2)

object.list2 <- list(TLN = cellchat_pos2, Tumor = cellchat_neg2)
cellchat_ccc <- mergeCellChat(object.list2, add.names = names(object.list))

immune_cco3 <- cell_cell_communcation(GSE11_immune,GSE11_immune$Manual_Annotations, signal_type = "Secreted Signaling")


cellchat_comparison1(cellchat_ecm)
cellchat_comparison1(cellchat_ccc)
cellchat_visiulization1(imm_ccc)
cellchat_visiulization1(immune_cco)


##$$$$$$$$$$$$$$$$$$$$$

par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = 4, targets.use = c(5:11), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = 4, targets.use = c(5:11), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))


i <- 2

single_interactions <- function(cellchat,name){
  mat <- cellchat@net$weight
  par(mfrow = c(1,1), xpd=TRUE)
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[name, ] <- mat[name, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[name],label.edge= T)
}

single_interactions(cellchat = cellchat,"B_Memory")



#here instead of b_ids we can put a single name instead and it will work fine
for (namez in b_idns){
  weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(object.list)) {
    netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= T, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]),
                     sources.use =namez)
  }
  
}


cellchat_comparison_full <- function(s.obj, id1,id2,the_factor, sources, signal_type){
  #types "Cell-Cell Contact", "ECM-Receptor", "Secreted Signaling"
  Idents(s.obj) <- the_factor
  tln <- subset(s.obj,idents=c(id1))
  tumor <- subset(s.obj,idents=c(id2))
  
  cellchat_pos <- cell_cell_communcation(tln,the_factor = tln$sub_annotations,signal_type = signal_type)
  cellchat_neg <- cell_cell_communcation(tumor,the_factor = tumor$sub_annotations,signal_type = signal_type)
  
  cellchat_pos <- updateCellChat(cellchat_pos)
  cellchat_neg <- updateCellChat(cellchat_neg)
  object.list <- list(TLN = cellchat_pos, TUMOR = cellchat_neg)
  cellchat.obj <- mergeCellChat(object.list, add.names = names(object.list))
  
  for (namez in sources){
    weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
    par(mfrow = c(1,2), xpd=TRUE)
    for (i in 1:length(object.list)) {
      cplot <- netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= T, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]),
                                sources.use =namez)
      print(cplot)
    }
    
  }
  
  ####################################
  pos.dataset = names(object.list)[1]
  # define a char name used for storing the results of differential expression analysis
  features.name = pos.dataset
  # perform differential expression analysis
  cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
  #> Use the joint cell labels from the merged CellChat object
  # map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
  net <- netMappingDEG(cellchat, features.name = features.name)
  # extract the ligand-receptor pairs with upregulated ligands in LS
  net.up <- subsetCommunication(cellchat, net = net, datasets = names(object.list)[1],ligand.logFC = 0.2, receptor.logFC = NULL)
  # extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
  net.down <- subsetCommunication(cellchat, net = net, datasets = condition2,ligand.logFC = -0.1, receptor.logFC = -0.1)
  
  gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
  gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
  
  pairLR.use.up = net.up[, "interaction_name", drop = F]
  gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 4, targets.use = c(2:4), comparison = c(1, 2),  angle.x = 90,
                          remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[1]))
  #> Comparing communications on a merged object
  pairLR.use.down = net.down[, "interaction_name", drop = F]
  gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 4, targets.use = c(2:4), comparison = c(1, 2),  angle.x = 90, 
                          remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[1]))
  #> Comparing communications on a merged object
  pplot <- gg1 + gg2
  print(pplot)
  
  #########################
  par(mfrow = c(1, 2), xpd=TRUE)
  # compare all the interactions sending from Inflam.FIB to DC cells
  for (i in 1:length(object.list)) {
    netVisual_chord_gene(object.list[[i]], sources.use = sources,  lab.cex = 0.5, title.name = paste0("Chord Signaling - ", names(object.list)[i]))
  }
  S
  netVisual_heatmap(cellchat_BT, measure = "weight")
  netVisual_heatmap(cellchat_BT)
  
  groupSize <- as.numeric(table(cellchat_BT@idents))
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_circle(cellchat_BT@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Number of interactions")
  netVisual_circle(cellchat_BT@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Interaction weights/strength")
  
  netVisual_bubble(cellchat_BT, remove.isolate = FALSE)
  netVisual_bubble(cellchat, remove.isolate = FALSE)
  
  gg1 <- rankNet(cellchat_BT, mode = "comparison", stacked = T, do.stat = TRUE)
  gg2 <- rankNet(cellchat_BT, mode = "comparison", stacked = F, do.stat = TRUE)
  gg1 + gg2
  
  netVisual_bubble(cellchat_BT,comparison = c(1, 2))
  
  gg1 <- compareInteractions(cellchat_BT, show.legend = F, group = c(1,2))
  gg2 <- compareInteractions(cellchat_BT, show.legend = F, group = c(1,2), measure = "weight")
  gg1 + gg2
  
  netVisual_bubble(cellchat_BT,comparison = c(1, 2), sources.use = b_idns)
  return(cellchat)
}






cellchat.obj_comparison_full <- function(s.obj, id1,id2,the_factor, sources, targets,signal_type){
  
  
  #types "Cell-Cell Contact", "ECM-Receptor", "Secreted Signaling"
  Idents(s.obj) <- the_factor
  tln <- subset(s.obj,idents=c(id1))
  tumor <- subset(s.obj,idents=c(id2))
  
  cellchat.obj_pos <- cell_cell_communcation(tln,the_factor = tln$sub_annotations,signal_type = signal_type)
  cellchat.obj_neg <- cell_cell_communcation(tumor,the_factor = tumor$sub_annotations,signal_type = signal_type)
  
  cellchat.obj_pos <- updatecellchat.obj(cellchat.obj_pos)
  cellchat.obj_neg <- updatecellchat.obj(cellchat.obj_neg)
  object.list <- list(TLN = cellchat.obj_pos, TUMOR = cellchat.obj_neg)
  cellchat.obj <- mergecellchat.obj(object.list, add.names = names(object.list))
  
  for (namez in sources){
    weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
    par(mfrow = c(1,2), xpd=TRUE)
    for (i in 1:length(object.list)) {
      cplot <- netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= T, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]),
                                sources.use =namez)
      print(cplot)
    }
    
  }
  
  for (namez in targets){
    weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
    par(mfrow = c(1,2), xpd=TRUE)
    for (i in 1:length(object.list)) {
      cplot <- netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= T, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]),
                                sources.use =namez)
      print(cplot)
    }
    
  }
  
  ####################################
  pos.dataset = names(object.list)[1]
  # define a char name used for storing the results of differential expression analysis
  features.name = pos.dataset
  # perform differential expression analysis
  cellchat.obj <- identifyOverExpressedGenes(cellchat.obj, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
  #> Use the joint cell labels from the merged cellchat.obj object
  # map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
  net <- netMappingDEG(cellchat.obj, features.name = features.name)
  # extract the ligand-receptor pairs with upregulated ligands in LS
  net.up <- subsetCommunication(cellchat.obj, net = net, datasets = names(object.list)[1],ligand.logFC = 0.2, receptor.logFC = NULL)
  # extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
  net.down <- subsetCommunication(cellchat.obj, net = net, datasets = condition2,ligand.logFC = -0.1, receptor.logFC = -0.1)
  
  gene.up <- extractGeneSubsetFromPair(net.up, cellchat.obj)
  gene.down <- extractGeneSubsetFromPair(net.down, cellchat.obj)
  
  pairLR.use.up = net.up[, "interaction_name", drop = F]
  gg1 <- netVisual_bubble(cellchat.obj, pairLR.use = pairLR.use.up, sources.use = 4, targets.use = c(2:4), comparison = c(1, 2),  angle.x = 90,
                          remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[1]))
  #> Comparing communications on a merged object
  pairLR.use.down = net.down[, "interaction_name", drop = F]
  gg2 <- netVisual_bubble(cellchat.obj, pairLR.use = pairLR.use.down, sources.use = 4, targets.use = c(2:4), comparison = c(1, 2),  angle.x = 90, 
                          remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[1]))
  #> Comparing communications on a merged object
  pplot <- gg1 + gg2
  print(pplot)
  
  #########################
  par(mfrow = c(1, 2), xpd=TRUE)
  # compare all the interactions sending from Inflam.FIB to DC cells
  for (i in 1:length(object.list)) {
    netVisual_chord_gene(object.list[[i]], sources.use = sources,  lab.cex = 0.5, title.name = paste0("Chord Signaling - ", names(object.list)[i]))
  }
  
  netVisual_heatmap(cellchat.obj, measure = "weight")
  netVisual_heatmap(cellchat.obj)
  
  groupSize <- as.numeric(table(cellchat.obj_BT@idents))
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_circle(cellchat.obj@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Number of interactions")
  netVisual_circle(cellchat.obj@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Interaction weights/strength")
  
  netVisual_bubble(cellchat.obj, remove.isolate = FALSE)
  netVisual_bubble(cellchat.obj, remove.isolate = FALSE)
  
  gg1 <- rankNet(cellchat.obj, mode = "comparison", stacked = T, do.stat = TRUE)
  gg2 <- rankNet(cellchat.obj, mode = "comparison", stacked = F, do.stat = TRUE)
  gg1 + gg2
  
  netVisual_bubble(cellchat.obj,comparison = c(1, 2))
  
  gg1 <- compareInteractions(cellchat.obj, show.legend = F, group = c(1,2))
  gg2 <- compareInteractions(cellchat.obj, show.legend = F, group = c(1,2), measure = "weight")
  gg1 + gg2
  
  netVisual_bubble(cellchat.obj,comparison = c(1, 2), sources.use = b_idns)
  return(cellchat.obj)
}

sig_cc <- cellchat.obj_comparison_full(comb, id1 = "Tumor", id2 = "TLN",the_factor = comb$updated_condtion,
                                       sources = b_idns, targets=t_idns, signal_type = "Secreted Signaling")

#########################################################################################################################################
###$##
###############################################################################################################
cellchat.obj_comparison_full <- function(s.obj, id1,id2,the_factor, sources, targets,signal_type){
  
  Idents(s.obj) <- the_factor
  tln <- subset(s.obj,idents=c(id1))
  tumor <- subset(s.obj,idents=c(id2))
  
  cellchat.obj_pos <- cell_cell_communcation(tln,the_factor = tln$sub_annotations,signal_type = signal_type)
  cellchat.obj_neg <- cell_cell_communcation(tumor,the_factor = tumor$sub_annotations,signal_type = signal_type)
  
  cellchat.obj_pos <- updateCellChat(cellchat.obj_pos)
  cellchat.obj_neg <- updateCellChat(cellchat.obj_neg)
  object.list <- list(TLN = cellchat.obj_pos, TUMOR = cellchat.obj_neg)
  cellchat.obj <- mergeCellChat(object.list, add.names = names(object.list))
  
  for (namez in sources){
    weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
    par(mfrow = c(1,2), xpd=TRUE)
    for (i in 1:length(object.list)) {
      cplot <- netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= T, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]),
                                sources.use =namez)
      print(cplot)
    }
    
  }
  
  for (namez in targets){
    weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
    par(mfrow = c(1,2), xpd=TRUE)
    for (i in 1:length(object.list)) {
      cplot <- netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= T, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]),
                                sources.use =namez)
      print(cplot)
    }
    
  }
  
  ####################################
  pos.dataset = names(object.list)[1]
  # define a char name used for storing the results of differential expression analysis
  features.name = pos.dataset
  # perform differential expression analysis
  cellchat.obj <- identifyOverExpressedGenes(cellchat.obj, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
  #> Use the joint cell labels from the merged cellchat.obj object
  # map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
  net <- netMappingDEG(cellchat.obj, features.name = features.name)
  # extract the ligand-receptor pairs with upregulated ligands in LS
  net.up <- subsetCommunication(cellchat.obj, net = net, datasets = names(object.list)[1],ligand.logFC = 0.2, receptor.logFC = NULL)
  # extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
  net.down <- subsetCommunication(cellchat.obj, net = net, datasets =names(object.list)[2],ligand.logFC = -0.1, receptor.logFC = -0.1)
  
  gene.up <- extractGeneSubsetFromPair(net.up, cellchat.obj)
  gene.down <- extractGeneSubsetFromPair(net.down, cellchat.obj)
  
  pairLR.use.up = net.up[, "interaction_name", drop = F]
  gg1 <- netVisual_bubble(cellchat.obj, pairLR.use = pairLR.use.up, sources.use = 4, targets.use = c(2:4), comparison = c(1, 2),  angle.x = 90,
                          remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[1]))
  #> Comparing communications on a merged object
  pairLR.use.down = net.down[, "interaction_name", drop = F]
  gg2 <- netVisual_bubble(cellchat.obj, pairLR.use = pairLR.use.down, sources.use = 4, targets.use = c(2:4), comparison = c(1, 2),  angle.x = 90, 
                          remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[1]))
  #> Comparing communications on a merged object
  pplot <- gg1 + gg2
  print(pplot)
  
  
  
  ###
  i = 1
  object.list[[1]] <- netAnalysis_computeCentrality(object.list[[1]])
  object.list[[2]] <- netAnalysis_computeCentrality(object.list[[2]])
  # combining all the identified signaling pathways from different datasets 
  pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
  ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
  ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
  draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
  
  ##########################################
  
  
  netVisual_heatmap(cellchat.obj, measure = "weight")
  netVisual_heatmap(cellchat.obj)
  
  groupSize <- as.numeric(table(cellchat.obj@idents))
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_circle(cellchat.obj@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Number of interactions")
  netVisual_circle(cellchat.obj@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Interaction weights/strength")
  
  netVisual_bubble(cellchat.obj, remove.isolate = FALSE)
  netVisual_bubble(cellchat.obj, remove.isolate = FALSE)
  
  gg1 <- rankNet(cellchat.obj, mode = "comparison", stacked = T, do.stat = TRUE)
  gg2 <- rankNet(cellchat.obj, mode = "comparison", stacked = F, do.stat = TRUE)
  gg1 + gg2
  
  netVisual_bubble(cellchat.obj,comparison = c(1, 2))
  
  gg1 <- compareInteractions(cellchat.obj, show.legend = F, group = c(1,2))
  gg2 <- compareInteractions(cellchat.obj, show.legend = F, group = c(1,2), measure = "weight")
  gg1 + gg2
  
  netVisual_bubble(cellchat.obj,comparison = c(1, 2), sources.use = b_idns)
  
  
}
