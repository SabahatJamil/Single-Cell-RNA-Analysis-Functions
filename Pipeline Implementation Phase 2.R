library(Seurat)
library(dplyr)
library(tidyverse)
library(SeuratData)
library(SeuratObject)
library(SeuratWrappers)
library(ggplot2)
library(openxlsx)

GSE1180286_immune <- read_rds("D:/Data_Sets/GSE18/New/Immune/GSE18_Immnune_N_Manual_Annotated.Rds")
GSE1180286_immune <- JoinLayers(GSE1180286_immune)
GSE1180286_B <- read_rds("D:/Data_Sets/GSE18/New/Immune/GSE1180286_B_Normalized_ScType_Annotated.Rds")
GSE1180286_B <- JoinLayers(GSE1180286_B)
GSE1180286_NK <- read_rds("D:/Data_Sets/GSE18/New/Immune/GSE1180286_NK_Normalized_ScType_Annotated.Rds")
GSE1180286_NK <- JoinLayers(GSE1180286_NK)

GSE1180286_T <- read_rds("D:/Data_Sets/GSE18/New/Immune/GSE1180286_T_Normalized_ScType_Annotated.Rds")
GSE1180286_T <- JoinLayers(GSE1180286_T)

GSE29 <- read_rds("D:/GSE29/Corrected/GSE16_Immune_N.Rds")
GSE29 <- JoinLayers(GSE29)

GSE_immune29 <- read_rds("D:/GSE29/Corrected/Immune/GSE16_Immune_Normalized_ScType_Annotated.Rds")
GSE_immune29  <- JoinLayers(GSE_immune29 )

GSE19 <- read_rds("D:/GSE19/Immune_Cells/Immune_All/GSE19_immune_N_Manual_Annotated.Rds")
GSE19  <- JoinLayers(GSE19)

GSE19Immune_Cells <-read_rds("D:/GSE19/Immune_Cells/Immune_All/GSE19_immune_N_Manual_Annotated.Rds")


GSE18Immune_Cells <- read_rds("D:/Data_Sets/GSE18/New/Immune/GSE18_Immnune_N_Manual_Annotated.Rds")
GSE18Immune_Cells <- JoinLayers(GSE18Immune_Cells)

GSE29_Immune_Cells <- read_rds("D:/GSE29/Corrected/Immune/GSE16_Immune_Normalized_ScType_Annotated.Rds")
GSE29_Immune_Cells <- JoinLayers(GSE29_Immune_Cells)
##################################################################################

comb_df <- cell_calculations4(GSE1180286_immune,name = "Immune",path = 
                                "D:/Data_Sets/GSE18/New/Immune/test/",
                              conditions = GSE1180286_immune$new_condition,annotations = GSE1180286_immune$Manual_Annotations)


#remeber this is subannotation so Manual_annotations factor can't be passed here instead the subannotation factor is passed here
cell_counts_ploting(GSE1180286_B,name = "B_cells",path = 
                      "D:/Data_Sets/GSE18/New/Immune/test/",
                    conditions = GSE1180286_B$new_condition,annotations = GSE1180286_B$sub_annotations)

cell_counts_ploting(GSE1180286_NK,name = "NK_cells",path = 
                      "D:/Data_Sets/GSE18/New/Immune/test/",
                    conditions = GSE1180286_NK$new_condition,annotations = GSE1180286_NK$sub_annotations)


cell_counts_ploting(GSE1180286_T,name = "T_cells",path = 
                      "D:/Data_Sets/GSE18/New/Immune/test/",
                    conditions = GSE1180286_T$new_condition,annotations = GSE1180286_T$sub_annotations)


#new object on main annotations
cell_counts_ploting(GSE29,name = "Immune_N",path = 
                      "D:/GSE29/Corrected/Cell_Calculations/",
                    conditions = GSE29$condition,annotations = GSE29$singlr_labels)

#new object on main annotations here too
cell_counts_ploting(GSE_immune29 ,name = "Immune_N",path = 
                      "D:/GSE29/Corrected/Immune/Cell Quantifications",
                    conditions = GSE_immune29$condition,annotations = GSE_immune29$sub_annotations)
#another new here on main annotations
cell_counts_ploting(GSE19 ,name = "Immune",path = 
                      "D:/GSE19/Immune_Cells/Immune_All/Cell Quantifications",
                    conditions = GSE19$condition,annotations = GSE19$Manual_Annotations)

#############################################################################################
#DEGs Section
############################################################################################
# for bigger combinations we can use these pretubrations meaning different combinations but caviet is ident 1 and 2
elements <- unique(GSE1180286_B$new_condition)
combinations <- combn(elements, 2)
combinations

######################

#but for controled and small we will use this code
# we just have to chnage the idn1 for most of  time if needed idn2 too, rest the function will do its work
for (i in unique(GSE1180286_immune$Manual_Annotations)){
  specific_labeled_degs_finder(GSE1180286_immune, cluster = i, path = "D:/Data_Sets/GSE18/New/Immune/outputs/DEGs New/",
                               the_factor = GSE1180286_immune$Manual_Annotations, idn1 = "TLN 2",idn2 = "Tumor",merger = T,
                               conditions =GSE1180286_immune$new_condition )
  
}


######################
#sub annotations DEGs

#change the path as it would be easier for later use
for (condition1 in c("TLN 1", "TLN 2")){
  for (i in unique(GSE1180286_B$sub_annotations)){
    specific_labeled_degs_finder(GSE1180286_B, cluster = i, path = "csub/",
                                 the_factor = GSE1180286_B$sub_annotations, idn1 = condition1,idn2 = "Tumor",merger = T,
                                 conditions =GSE1180286_B$new_condition )
  }
}

for (condition1 in c("TLN 1", "TLN 2")){
  for (i in unique(GSE1180286_NK$sub_annotations)){
    specific_labeled_degs_finder(GSE1180286_NK, cluster = i, path = "D:/Data_Sets/GSE18/New/Immune/outputs/DEGs New/sub/",
                                 the_factor = GSE1180286_NK$sub_annotations, idn1 = condition1,idn2 = "Tumor",merger = T,
                                 conditions =GSE1180286_NK$new_condition )
  }
}

for (condition1 in c("TLN 1", "TLN 2")){
  for (i in unique(GSE1180286_T$sub_annotations)){
    specific_labeled_degs_finder(GSE1180286_T, cluster = i, path = "D:/Data_Sets/GSE18/New/Immune/outputs/DEGs New/sub/",
                                 the_factor = GSE1180286_T$sub_annotations, idn1 = condition1,idn2 = "Tumor",merger = T,
                                 conditions =GSE1180286_T$new_condition )
  }
}
###################################new object
idn1 <- "TLN " 
idn2 <- "Tumour "
for (i in unique(GSE_immune29$sub_annotations)){
  specific_labeled_degs_finder(GSE_immune29, cluster = i, path = "D:/GSE29/Corrected/Immune/DEGs/",
                               the_factor = GSE_immune29$sub_annotations, idn1 = idn1,idn2 = idn2,merger = T,
                               conditions =GSE_immune29$condition )
  
}

############another new
idn1 <- "IDC-LNM"
idn2 <- "IDC-Tumour"
for (i in unique(GSE19$Manual_Annotations)){
  specific_labeled_degs_finder(GSE19, cluster = i, path = "D:/GSE19/Immune_Cells/Immune_All/DEGs/",
                               the_factor = GSE19$Manual_Annotations, idn1 = idn1,idn2 = idn2,merger = T,
                               conditions =GSE19$condition )
  
}

idn1 <- "Normal"
idn2 <- "DCIS-Tumour"
for (i in unique(GSE19$Manual_Annotations)){
  specific_labeled_degs_finder(GSE19, cluster = i, path = "D:/GSE19/Immune_Cells/Immune_All/DEGs/Normal DC Tumor/",
                               the_factor = GSE19$Manual_Annotations, idn1 = idn1,idn2 = idn2,merger = T,
                               conditions =GSE19$condition )
  
}



idn1 <- "IDC-Tumour"
idn2 <- "IDC-LNM"
for (i in unique(GSE19$Manual_Annotations)){
  specific_labeled_degs_finder(GSE19, cluster = i, path = "D:/GSE19/Immune_Cells/Immune_All/DEGs/IDC Tumor IDC LNM/",
                               the_factor = GSE19$Manual_Annotations, idn1 = idn1,idn2 = idn2,merger = T,
                               conditions =GSE19$condition )
  
}

idn1 <- "DCIS-Tumour"  
idn2 <- "IDC-LNM"
for (i in unique(GSE19$Manual_Annotations)){
  specific_labeled_degs_finder(GSE19, cluster = i, path = "D:/GSE19/Immune_Cells/Immune_All/DEGs/DCIS Tumor IDC LNM/",
                               the_factor = GSE19$Manual_Annotations, idn1 = idn1,idn2 = idn2,merger = T,
                               conditions =GSE19$condition )
  
}

###another

idn1 <- "TLN"  
idn2 <- "Tumor"
for (i in unique(GSE18Immune_Cells$Manual_Annotations)){
  specific_labeled_degs_finder(GSE18Immune_Cells, cluster = i, path = "D:/Data_Sets/GSE18/New/Immune/TLN-Tumor/",
                               the_factor = GSE18Immune_Cells$Manual_Annotations, idn1 = idn1,idn2 = idn2,merger = T,
                               conditions =GSE18Immune_Cells$condition )
  
}

idn1 <- "Tumor" 
idn2 <- "TLN"
for (i in unique(GSE18Immune_Cells$Manual_Annotations)){
  specific_labeled_degs_finder(GSE18Immune_Cells, cluster = i, path = "D:/Data_Sets/GSE18/New/Immune/DEGs/",
                               the_factor = GSE18Immune_Cells$Manual_Annotations, idn1 = idn1,idn2 = idn2,merger = T,
                               conditions =GSE18Immune_Cells$condition )
  
}


idn1 <- "TLN 1"
idn2 <- "Tumor"
for (i in unique(GSE18Immune_Cells$Manual_Annotations)){
  specific_labeled_degs_finder(GSE18Immune_Cells, cluster = i, path = "D:/Data_Sets/GSE18/New/Immune/DEGs/",
                               the_factor = GSE18Immune_Cells$Manual_Annotations, idn1 = idn1,idn2 = idn2,merger = T,
                               conditions =GSE18Immune_Cells$new_condition )
  
}

idn1 <- "TLN 2"
idn2 <- "Tumor"
for (i in unique(GSE18Immune_Cells$Manual_Annotations)){
  specific_labeled_degs_finder(GSE18Immune_Cells, cluster = i, path = "D:/Data_Sets/GSE18/New/Immune/DEGs/",
                               the_factor = GSE18Immune_Cells$Manual_Annotations, idn1 = idn1,idn2 = idn2,merger = T,
                               conditions =GSE18Immune_Cells$new_condition )
  
}

idn1 <- "Normal"
idn2 <- "Tumour "
for (i in unique(GSE29_Immune_Cells$sub_annotations)){
  specific_labeled_degs_finder(GSE29_Immune_Cells, cluster = i, path = "D:/GSE29/Corrected/Immune/New_DEGs/",
                               the_factor = GSE29_Immune_Cells$sub_annotations, idn1 = idn1,idn2 = idn2,merger = T,
                               conditions =GSE29_Immune_Cells$condition )
  
}

idn1 <- "Tumour "
idn2 <- "TLN "
for (i in unique(GSE29_Immune_Cells$sub_annotations)){
  specific_labeled_degs_finder(GSE29_Immune_Cells, cluster = i, path = "D:/GSE29/Corrected/Immune/New_DEGs/",
                               the_factor = GSE29_Immune_Cells$sub_annotations, idn1 = idn1,idn2 = idn2,merger = T,
                               conditions =GSE29_Immune_Cells$condition )
  
}

idn1 <- "TLN "
idn2 <-  "Tumour "
for (i in unique(GSE29_Immune_Cells$sub_annotations)){
  specific_labeled_degs_finder(GSE29_Immune_Cells, cluster = i, path = "D:/GSE29/Corrected/Immune/New_DEGs/",
                               the_factor = GSE29_Immune_Cells$sub_annotations, idn1 = idn1,idn2 = idn2,merger = T,
                               conditions =GSE29_Immune_Cells$condition )
  
}

#################new combined B and T DEGs
idn1 <- "TLN"
idn2 <-  "Tumor"
Idents(comb) <- comb$updated_condtion
for (i in unique(comb$sub_annotations)){
  specific_labeled_degs_finder(comb, cluster = i, path = "D:/Data_Sets/Combined Object B and T/Annotated/Combined_B_T_DEGs/",
                               the_factor = comb$sub_annotations, idn1 = idn1,idn2 = idn2,merger = T,
                               conditions =comb$updated_condtion )
  
}
#########################################################################################

#GSEA Section
#############################################################################################################
gsea_preprocess <- function(path){
  degs <- read_csv(path)
  degs <- na.omit(degs)
  rownames(degs) <- degs$Names
  degs <- degs[order(-degs$avg_log2FC), ]
  
  return(degs)
}

idn1 <- "TLN 1"
idn2 <- "Tumor"
for (i in unique(GSE18Immune_Cells$Manual_Annotations)){
  specific_labeled_degs_finder(GSE18Immune_Cells, cluster = i, path = "D:/Data_Sets/GSE18/New/Immune/DEGs/",
                               the_factor = GSE18Immune_Cells$Manual_Annotations, idn1 = idn1,idn2 = idn2,merger = T,
                               conditions =GSE18Immune_Cells$new_condition )
  
}idn1 <- "Normal"
idn2 <- "IDC-LNM"
for (i in unique(GSE19$Manual_Annotations)){
  specific_labeled_degs_finder(GSE19, cluster = i, path = "D:/GSE19/Immune_Cells/Immune_All/DEGs/Normal IDC-LNM/",
                               the_factor = GSE19$Manual_Annotations, idn1 = idn1,idn2 = idn2,merger = T,
                               conditions =GSE19$condition )
  
}kegg_simp_plotter_human(instance_markers = B_cell,name = "B-Cells",path = "D:/Data_Sets/GSE18/New/Immune/outputs/",ctgry = 13
)
wp_simp_plotter_human(instance_markers = B_cell,name = "B-Cells",path = "D:/Data_Sets/GSE18/New/Immune/outputs/",ctgry = 13
)
go_simp_plotter_human(instance_markers = B_cell,name = "B-Cells",path = "D:/Data_Sets/GSE18/New/Immune/outputs/",ctgry = 13
)

file_paths <- list.files(path = "D:/Data_Sets/GSE18/New/Immune/outputs/DEGs New/B_cells/", full.names = TRUE)
for (file in file_paths){
  Degs <- gsea_preprocess(file)
  initial <- tail(unlist(strsplit(file, "/")),n=1)
  Name <- sub("\\.csv$", "", initial)
  print(paste("Running", Name, "KEGG"))
  dir.create(paste0("D:/Data_Sets/GSE18/New/Immune/outputs/GSEA/",Name))
  kegg_simp_plotter_human(instance_markers = Degs,name = Name,path = paste0("D:/Data_Sets/GSE18/New/Immune/outputs/GSEA/",Name),ctgry = 13
  )
  print(paste("Running", Name, "Wiki"))
  wp_simp_plotter_human(instance_markers = Degs,name = Name,path = paste0("D:/Data_Sets/GSE18/New/Immune/outputs/GSEA/",Name),ctgry = 13
  )
}


file_paths <- list.files(path = "D:/Data_Sets/GSE18/New/Immune/outputs/DEGs New/T_cells/", full.names = TRUE)
path = "D:/Data_Sets/GSE18/New/Immune/outputs/DEGs New/T_cells/"
for (file in file_paths){
  Degs <- gsea_preprocess(file)
  initial <- tail(unlist(strsplit(file, "/")),n=1)
  Name <- sub("\\.csv$", "", initial)
  main_Name <- tail(unlist(strsplit(path, "/")),n=1)
  print(paste("Running", Name, "KEGG"))
  dir.create(paste0("D:/Data_Sets/GSE18/New/Immune/outputs/GSEA/",main_Name))
  kegg_simp_plotter_human(instance_markers = Degs,name = Name,path = paste0("D:/Data_Sets/GSE18/New/Immune/outputs/GSEA/",main_Name),ctgry = 13
  )
  print(paste("Running", Name, "Wiki"))
  wp_simp_plotter_human(instance_markers = Degs,name = Name,path = paste0("D:/Data_Sets/GSE18/New/Immune/outputs/GSEA/",main_Name),ctgry = 13
  )
}

# for running this function we need to load two other functions i) gsea_preprocessing and all the human enrichment functio
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
      
      kegg_results <- kegg_simp_plotter_human2(
        instance_markers = Degs,
        name = Name,
        path = paste0(out_main_folder_path, main_Name),
        ctgry = 13
      )
      
      print(paste("Running", Name, "Wiki"))
      wp_result <- WP_OR_GSEA(
        instance_markers = Degs,
        name = Name,
        path = paste0(out_main_folder_path, main_Name),
        ctgry = 13
      )
      rec_results <- Reactome_OR_GSEA (instance_markers = Degs,
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

DEGs_function(degs_folder_path = "D:/Data_Sets/GSE18/New/Immune/outputs/DEGs New/B_cells/",
              out_main_folder_path = "D:/Data_Sets/GSE18/New/Immune/outputs/GSEA/")



DEGs_function(degs_folder_path = "D:/Data_Sets/GSE18/New/Immune/outputs/DEGs New/NK_cells/",
              out_main_folder_path = "D:/Data_Sets/GSE18/New/Immune/outputs/GSEA/")

#new objects
DEGs_function(degs_folder_path = "D:/GSE29/Corrected/Immune/DEGs/",
              out_main_folder_path = "D:/GSE29/Corrected/Immune/GSEA/")

DEGs_function_OR(degs_folder_path = "D:/GSE29/Corrected/Immune/DEGs/",
                 out_main_folder_path = "D:/GSE29/Corrected/Immune/GSEA_OR/")

#Another new
DEGs_function(degs_folder_path = "D:/GSE19/Immune_Cells/Immune_All/DEGs/",
              out_main_folder_path = "D:/GSE19/Immune_Cells/Immune_All/GSEA/")

DEGs_function_OR(degs_folder_path = "D:/GSE19/Immune_Cells/Immune_All/DEGs/",
                 out_main_folder_path = "D:/GSE19/Immune_Cells/Immune_All/GSEA_OR/")
#old

DEGs_function_OR(degs_folder_path = "D:/Data_Sets/GSE18/New/Immune/outputs/DEGs New/NK_cells/",
                 out_main_folder_path = "D:/Data_Sets/GSE18/New/Immune/outputs/GSEA_OR/")

DEGs_function_OR(degs_folder_path = "D:/Data_Sets/GSE18/New/Immune/outputs/DEGs New/B_cells/",
                 out_main_folder_path = "D:/Data_Sets/GSE18/New/Immune/outputs/GSEA_OR/")

DEGs_function_OR(degs_folder_path = "D:/Data_Sets/GSE18/New/Immune/outputs/DEGs New/T_cells/",
                 out_main_folder_path = "D:/Data_Sets/GSE18/New/Immune/outputs/GSEA_OR/")

# more new

DEGs_function_OR(degs_folder_path = "D:/GSE19/Immune_Cells/Immune_All/DEGs/Normal DC Tumor/",
                 out_main_folder_path = "D:/GSE19/Immune_Cells/Immune_All/GSEA_OR/Normal DC Tumor/")

DEGs_function_OR(degs_folder_path = "D:/GSE19/Immune_Cells/Immune_All/DEGs/Normal DC Tumor/",
                 out_main_folder_path = "D:/Data_Sets/GSE18/New/Immune/DEGs/")
#######################################################################################
##FOR GSE29
####################################################################################################

DEGs_function(degs_folder_path = "D:/Data_Sets/Combined Object B and T/Annotated/Combined_B_T_DEGs/",
              out_main_folder_path = "D:/Data_Sets/Combined Object B and T/Annotated/")
DEGs_function_OR(degs_folder_path = "D:/Data_Sets/Combined Object B and T/Annotated/Combined_B_T_DEGs/",
                 out_main_folder_path = "D:/Data_Sets/Combined Object B and T/Annotated/GSEA_OR/")
