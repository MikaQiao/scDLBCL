rm(list=ls())
library(DoubletFinder)
library(Seurat)
library(SeuratData)
library(harmony)
library(readr)
library(tidyverse)
library(data.table)
library(patchwork)
library(clustree)
library(ggplotify)
library(ComplexHeatmap)
library(Scillus)
library(openxlsx)
library(clusterProfiler)
library(lisi)
library(stringr)


# library(monocle)
# library(DT)
# library(circlize)
# library(org.Hs.eg.db)
# library(future)
# library(scPred)
plan("multiprocess", workers = 20)


in_dir <- "~/project/scYZX/scDLBCL/data/rename/"
out_dir <- "~/project/scYZX/scDLBCL/analysis/output/"





################################################################################################################
# 0.1 load raw input data ---
s_list <- c("DLBCL02_ZLQ", "DLBCL04_ZLQ", 
            "DLBCL03_ZLQ", "DLBCL03_ZLZ", "DLBCL03_FF", "DLBCL06_ZLQ", "DLBCL06_ZLZ", "DLBCL06_FF", "DLBCL07_ZLQ", "DLBCL07_ZLZ")

for (i in s_list){
  expr <- sprintf("%s.data <- Read10X(data.dir = paste0(in_dir, '%s'))", i, i)
  eval(parse(text=expr))
}


pdf(paste0(out_dir, "OFig.VlnPlot.pdf"), width=10, height=5)
for (i in s_list){
  expr2 <- paste0(sprintf("CreateSeuratObject(counts = %s.data, project = '%s', min.cells = 3, min.features = 200)", i, i), " %>% PercentageFeatureSet(pattern = '^MT-', col.name = 'percent_mt') %>% PercentageFeatureSet(pattern = 'RP[SL]', col.name = 'percent_ribo') %>% PercentageFeatureSet(pattern = '^HB[^(P)]', col.name = 'percent_hb') %>% VlnPlot(features = c('nFeature_RNA', 'nCount_RNA', 'percent_mt', 'percent_ribo', 'percent_hb'), ncol = 5)")
  print(eval(parse(text=expr2)))
}
dev.off()



for (i in s_list){
  expr2 <- paste0(sprintf("%s <- CreateSeuratObject(counts = %s.data, project = '%s', min.cells = 3, min.features = 200)", i, i, i), " %>% PercentageFeatureSet(pattern = '^MT-', col.name = 'percent_mt') %>% PercentageFeatureSet(pattern = 'RP[SL]', col.name = 'percent_ribo') %>% PercentageFeatureSet(pattern = '^HB[^(P)]', col.name = 'percent_hb') %>% subset(subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 60000 & percent_mt < 20 & percent_hb < 15)")
  eval(parse(text=expr2))
}


# remove doublets: Assuming 10% doublet formation rate---
DLBCL02_ZLQ.sub <- DLBCL02_ZLQ[, my_doublet_v2(DLBCL02_ZLQ, my_rate = 0.1)@meta.data %>% filter(DF_hi.lo == "Singlet") %>% rownames(.)]
DLBCL04_ZLQ.sub <- DLBCL04_ZLQ[, my_doublet_v2(DLBCL04_ZLQ, my_rate = 0.1)@meta.data %>% filter(DF_hi.lo == "Singlet") %>% rownames(.)]
DLBCL03_ZLQ.sub <- DLBCL03_ZLQ[, my_doublet_v2(DLBCL03_ZLQ, my_rate = 0.1)@meta.data %>% filter(DF_hi.lo == "Singlet") %>% rownames(.)]
DLBCL03_ZLZ.sub <- DLBCL03_ZLZ[, my_doublet_v2(DLBCL03_ZLZ, my_rate = 0.1)@meta.data %>% filter(DF_hi.lo == "Singlet") %>% rownames(.)]
DLBCL03_FF.sub <- DLBCL03_FF[, my_doublet_v2(DLBCL03_FF, my_rate = 0.1)@meta.data %>% filter(DF_hi.lo == "Singlet") %>% rownames(.)]
DLBCL06_ZLQ.sub <- DLBCL06_ZLQ[, my_doublet_v2(DLBCL06_ZLQ, my_rate = 0.1)@meta.data %>% filter(DF_hi.lo == "Singlet") %>% rownames(.)]
DLBCL06_ZLZ.sub <- DLBCL06_ZLZ[, my_doublet_v2(DLBCL06_ZLZ, my_rate = 0.1)@meta.data %>% filter(DF_hi.lo == "Singlet") %>% rownames(.)]
DLBCL06_FF.sub <- DLBCL06_FF[, my_doublet_v2(DLBCL06_FF, my_rate = 0.1)@meta.data %>% filter(DF_hi.lo == "Singlet") %>% rownames(.)]
DLBCL07_ZLQ.sub <- DLBCL07_ZLQ[, my_doublet_v2(DLBCL07_ZLQ, my_rate = 0.1)@meta.data %>% filter(DF_hi.lo == "Singlet") %>% rownames(.)]
DLBCL07_ZLZ.sub <- DLBCL07_ZLZ[, my_doublet_v2(DLBCL07_ZLZ, my_rate = 0.1)@meta.data %>% filter(DF_hi.lo == "Singlet") %>% rownames(.)]


# set group info --
DLBCL02_ZLQ.sub[["group"]] <- "CR"
DLBCL04_ZLQ.sub[["group"]] <- "CR"
DLBCL03_ZLQ.sub[["group"]] <- "PD"
DLBCL03_ZLZ.sub[["group"]] <- "PD"
DLBCL03_FF.sub[["group"]] <- "PD"
DLBCL06_ZLQ.sub[["group"]] <- "PD"
DLBCL06_ZLZ.sub[["group"]] <- "PD"
DLBCL06_FF.sub[["group"]] <- "PD"
DLBCL07_ZLQ.sub[["group"]] <- "PD"
DLBCL07_ZLZ.sub[["group"]] <- "PD"

DLBCL02_ZLQ.sub[["time"]] <- "ZLQ"
DLBCL04_ZLQ.sub[["time"]] <- "ZLQ"
DLBCL03_ZLQ.sub[["time"]] <- "ZLQ"
DLBCL03_ZLZ.sub[["time"]] <- "ZLZ"
DLBCL03_FF.sub[["time"]] <- "FF"
DLBCL06_ZLQ.sub[["time"]] <- "ZLQ"
DLBCL06_ZLZ.sub[["time"]] <- "ZLZ"
DLBCL06_FF.sub[["time"]] <- "FF"
DLBCL07_ZLQ.sub[["time"]] <- "ZLQ"
DLBCL07_ZLZ.sub[["time"]] <- "ZLZ"


g_list <- c("CR", "PD")
t_list <- c("ZLQ", "ZLZ", "FF")
s_list <- c(c("DLBCL02_ZLQ", "DLBCL04_ZLQ", "DLBCL03_ZLQ", "DLBCL03_ZLZ", "DLBCL03_FF", 
            "DLBCL06_ZLQ", "DLBCL06_ZLZ", "DLBCL06_FF", "DLBCL07_ZLQ", "DLBCL07_ZLZ"))


g_color <- c("#61A2CA", "#FB8072")
t_color <- c("#F0907F", "#C7E8FA", "#C28DBD")

my_color <- c("#91B897", "#FCDD8D", "#E2CDDA", "#AFC9C6", "#E06C3E", "#CDD695", "#9188A4", "#FBA95E", "#D58E8E", "#8EB052", 
              "#A6DDC2", "#46A39B", "#FBAD82", "#DFB7A9", "#D966A4", "#ACD9E6", "#E2AED4", "#A4B8B3", "#B8CCE2", "#F8D7AA", 
              "#74ACD1", "#99DD8A", "#D2D5BE", "#F7B5CE", "#A6D297", "#C69D97", "#61A2CA", "#FB8072", "#C28DBD")


scDLBCL.list <- list(DLBCL02_ZLQ = DLBCL02_ZLQ.sub, DLBCL04_ZLQ = DLBCL04_ZLQ.sub, 
                     DLBCL03_ZLQ = DLBCL03_ZLQ.sub, DLBCL03_ZLZ = DLBCL03_ZLZ.sub, DLBCL03_FF = DLBCL03_FF.sub, 
                     DLBCL06_ZLQ = DLBCL06_ZLQ.sub, DLBCL06_ZLZ = DLBCL06_ZLZ.sub, DLBCL06_FF = DLBCL06_FF.sub, 
                     DLBCL07_ZLQ = DLBCL07_ZLQ.sub, DLBCL07_ZLZ = DLBCL07_ZLZ.sub)
# write_rds(scDLBCL.list, file = paste0(out_dir, "Table0.scDLBCL.list.rds"))
scDLBCL.list <- read_rds(paste0("~/project/scYZX/scDLBCL/analysis/out_10s_20220822/", "Table0.scDLBCL.list.rds"))






#  Preprocess ##############################################################################################################
## 1.0 merge data: with batch----
scDLBCL.merge <- merge(scDLBCL.list[[1]], y = c(scDLBCL.list[2:length(scDLBCL.list)]), 
                project = "scDLBCL", merge.data = T, add.cell.ids = c(names(scDLBCL.list)))

idx_black <- read_rds(paste0("~/project/scYZX/scDLBCL/analysis/out_10s_20220822/", "Table0.idx_black.rds"))

# scDLBCL <- scDLBCL.merge[, setdiff(colnames(scDLBCL.merge), idx_black)]




## 1.1.1 intergrated data ----
options(future.globals.maxSize = 10000 * 1024^2)   # 10000M







#  Harmony ################################################################################################################
## 1.3 remove batch by Harmony ----
# SCTransform is an alternative to the NormalizeData, FindVariableFeatures, ScaleData workflow.
scDLBCL <- scDLBCL.merge[, setdiff(colnames(scDLBCL.merge), idx_black)]


## Cellcycle Score
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
scDLBCL <- CellCycleScoring(scDLBCL, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)

## HSP scores
hsp_features <- list(rownames(scDLBCL)[grepl("^HSP", rownames(scDLBCL))])
scDLBCL <- AddModuleScore(object = scDLBCL, features = hsp_features, name = "HSP.Score")
scDLBCL@meta.data <- scDLBCL@meta.data %>% dplyr::rename("HSP.Score" = "HSP.Score1")



## pipeline
scDLBCL <- scDLBCL %>% 
            NormalizeData() %>% 
            FindVariableFeatures(nfeatures = 2500) %>% 
            ScaleData(vars.to.regress = c("nCount_RNA", "percent_mt", "S.Score", "G2M.Score")) %>% 
            RunPCA(verbose = F)
scDLBCL <- RunHarmony(scDLBCL, group.by.vars = "orig.ident", theta = 2, lambda = 2)

use.pcs = 1:30
scDLBCL <- scDLBCL %>% 
            RunUMAP(reduction = "harmony", dims = use.pcs) %>% 
            FindNeighbors(reduction = "harmony", dims = use.pcs) 


scDLBCL <- scDLBCL %>% FindClusters(resolution = 0.8)
scDLBCL@meta.data <- scDLBCL@meta.data %>% 
        mutate(orig.ident = factor(orig.ident, levels = s_list),
               group = factor(group, levels = g_list), 
               time = factor(time, levels = t_list))

p1 <- DimPlot(scDLBCL, reduction = "umap", group.by = "seurat_clusters", label = T, raster = T, cols = my_color)
p2 <- DimPlot(scDLBCL, reduction = "umap", group.by = "orig.ident", label = F, raster = T, cols = my_color) 
p3 <- DimPlot(scDLBCL %>% subset(time == "ZLQ"), reduction = "umap", group.by = "group", label = T, raster = T, cols = g_color) 
p4 <- DimPlot(scDLBCL %>% subset(group != "CR"), reduction = "umap", group.by = "time", label = T, raster = T, cols = t_color) 

pdf(paste0(out_dir, "Fig2b_1.pdf"), width = 9, height = 7)
p1 + p2 + p3 + p4 + plot_layout(ncol = 2)
dev.off()


# write_rds(scDLBCL, file = paste0(out_dir, "Table1.scDLBCL_1raw.rds"))
# scDLBCL1 <- read_rds(paste0("~/project/scYZX/scDLBCL/analysis/out_10s_20221030/", "Table1.scDLBCL_2harmony.rds"))



pdf(paste0(out_dir, "SFig3a.pdf"), width = 10, height = 13)
FeaturePlot(scDLBCL, features = c("PLVAP", "PECAM1", "VWF", 
                                  "CALD1", "ACTA2", "DCN",
                                  "CST3", "LYZ", 
                                  "CD3D", "CD8A", "CD4", "NKG7", 
                                  "CD79A", "MS4A1", "JCHAIN", "MKI67"),
            ncol = 4, pt.size = 0.00001, order = T, reduction = "umap", raster = T) & 
            scale_color_gradient2(midpoint = 1, low = "#313695", mid = "#FFFFBF", high = "#A50026") & 
            theme_void() & theme(panel.border = element_rect(fill = NA, colour = "black", size = 1.5), 
                                 axis.line = element_blank(),
                                 plot.title = element_text(size = 14, vjust = -4), legend.position = "top",
                                 plot.margin = unit(c(1, 1, 1, 1), "mm"))
dev.off()


####### re cluster 9
scDLBCL.sel <- subset(scDLBCL, idents = c("9")) %>% 
                    NormalizeData() %>% 
                    FindVariableFeatures(nfeatures = 2000) %>% 
                    ScaleData(vars.to.regress = c("nCount_RNA", "percent_mt", "percent_ribo")) %>% 
                    RunPCA() %>% 
                    RunHarmony(group.by.vars = "orig.ident") %>% 
                    RunUMAP(reduction = "harmony", dims = use.pcs) %>% 
                    FindNeighbors(reduction = "harmony", dims = use.pcs) %>% 
                    FindClusters(resolution = 0.8)

## re-calssification
p1 <- DimPlot(scDLBCL.sel, reduction = "umap", group.by = "seurat_clusters", label = T, cols = my_color) + theme(legend.position = "right")
p2 <- DimPlot(scDLBCL.sel, reduction = "umap", group.by = "orig.ident", label = F, cols = my_color) 
p3 <- DimPlot(scDLBCL.sel %>% subset(time == "ZLQ"), reduction = "umap", group.by = "group", label = T, cols = g_color) 
p4 <- DimPlot(scDLBCL.sel %>% subset(group == "PD"), reduction = "umap", group.by = "time", label = T, cols = t_color) 

p5 <- DimPlot(scDLBCL.sel, group.by = "predicted.cell_type1", label = T, repel = T, cols = my_color)
p6 <- DimPlot(scDLBCL.sel, group.by = "predicted.cell_type2", label = T, repel = T, cols = my_color)
p7 <- DimPlot(scDLBCL.sel, reduction = "umap", group.by = "Phase", , label = T, repel = T, cols = my_color)

pdf(paste0(out_dir, "OFig.cluster.9_harmony.pdf"), width = 20, height = 7.5)
p1 + p2 + p3 + p4 + plot_spacer() + p6 + p7 + p5 + plot_layout(ncol = 4)
dev.off()

Idents(scDLBCL.sel) <- "seurat_clusters"
cluster.dot <- DotPlot(object = scDLBCL.sel, assay = "RNA", features = rev(unique(index)[22:41]), 
                       cols = c("white", "red"), col.min = -0.1, col.max = 3) + 
                  theme_bw() + theme(axis.text.x = element_text(angle=45, hjust = 1)) + coord_flip()

pdf(paste0(out_dir, "OFig.cluster.9_dot.pdf"), width = 6, height = 6)
cluster.dot
dev.off()



# write_rds(scDLBCL.sel, file = paste0(out_dir, "Table1.scDLBCL.sel.rds"))
# scDLBCL.sel <- read_rds(paste0("~/project/scYZX/scDLBCL/analysis/out_10s_20221030/", "Table1.scDLBCL.sel.rds"))




## 1.5 annotation ----
my_levels <- c("B", "T/NK", "Myeloid", "Fibroblast", "Endothelia")
Idents(scDLBCL) <- "seurat_clusters"
scDLBCL.anno <- RenameIdents(scDLBCL, 
        "0" = "B", "2" = "B", "4" = "B", "5" = "B", "6" = "B", "7" = "B", "8" = "B", "10" = "B", 
        "13" = "B", "15" = "B", "18" = "B", "19" = "B", "20" = "B", "22" = "B","23" = "B",
        "1" = "T/NK", "3" = "T/NK", "11" = "T/NK", "16" = "T/NK", "24" = "T/NK", "25" = "T/NK", "27" = "T/NK", 
        "14" = "Fibroblast", "17" = "Endothelia", "26" = "Endothelia", 
        "12" = "Myeloid", "21" = "Myeloid")

# for cluster 9
cells.use.B <- WhichCells(scDLBCL.sel, idents = c("0", "3", "5", "6", "9"))
cells.use.T <- WhichCells(scDLBCL.sel, idents = c("1", "2", "4", "7", "8", "10"))
scDLBCL.anno <- SetIdent(scDLBCL.anno, cells = cells.use.B, value = "B")
scDLBCL.anno <- SetIdent(scDLBCL.anno, cells = cells.use.T, value = "T/NK")


levels(scDLBCL.anno) <- my_levels
scDLBCL.anno@meta.data <- scDLBCL.anno@meta.data %>% 
            mutate(cell_type = Idents(scDLBCL.anno), 
                   group = factor(group, levels = g_list),
                   seurat_clusters = cell_type)

color_re <- c("#FB8072", "#61A2CA", "#8DD3C7", "#BEBADA", "#B2DF8A", "#FDBF6F", "#FCCDE5", 
              "#A6CEE3", "#FF7F00", "#FFFF99")

pdf(paste0(out_dir, "Fig2b_2.pdf"), width = 6, height = 5)
DimPlot(scDLBCL.anno, reduction = "umap", group.by = "cell_type", label = T, raster = T, cols = color_re)
dev.off()

# write_rds(scDLBCL.anno, file = paste0(out_dir, "Table1.scDLBCL_3harmony_anno.rds"))







# T cell ######################################################################################################
## T cell sub cluster analysis----
cells.T <- c(WhichCells(scDLBCL, idents = c("1", "3", "11", "16", "24", "25", "27")), 
             WhichCells(scDLBCL.sel, idents = c("1", "2", "4", "7", "8", "10"))
            )
use.pcs = 1:40
scDLBCL.T <- scDLBCL[, cells.T] %>% 
                    NormalizeData() %>% 
                    FindVariableFeatures(nfeatures = 3000) %>% 
                    ScaleData(vars.to.regress = c("nCount_RNA", "percent_mt", "percent_ribo", "S.Score", "G2M.Score")) %>% 
                    RunPCA() %>% 
                    RunHarmony(group.by.vars = "orig.ident") %>% 
                    RunUMAP(reduction = "harmony", dims = use.pcs) %>% 
                    FindNeighbors(reduction = "harmony", dims = use.pcs) %>% 
                    FindClusters(resolution = 1.5)

## re-calssification
p1 <- DimPlot(scDLBCL.T, reduction = "umap", group.by = "seurat_clusters", label = T, cols = my_color) + theme(legend.position = "right")

pdf(paste0(out_dir, "OFig.T_harmony_raw.pdf"), width = 6, height = 5)
p1
dev.off()




# write_rds(scDLBCL.T, file = paste0(out_dir, "Table3.scDLBCL.T_1raw.rds"))
# scDLBCL.T <- read_rds(paste0(out_dir, "v3_T_3000HVG_40PC/Table3.scDLBCL.T_1raw.rds"))



## T cell annotation ----
color_T <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#E31A1C", "#FB9A99", "#FDBF6F", 
             "#FF7F00", "#CAB2D6", "#6A3D9A", "#FCDD8D", "#F7B5CE", "#61A2CA", "#FBAD82", "#B15928", "#BFBFBF")
my_levels <- c("NK", "CD8_Teff", "CD8_Tex_1", "CD8_Tex_2", "CD8_Tex_3", "CD8_Tex_4", "CD8_Tex_5", "CD8_Tex_6",
               "CD8_Prolif_1", "CD8_Prolif_2", "CD8_Prolif_3", "CD4_Naive", "CD4_Treg", "CD4_Mem", "CD4_Th1_like", "LowQual")

scDLBCL.T.anno <- RenameIdents(scDLBCL.T, 
                     "15" = "NK", "11" = "NK", "20" = "NK", 
                     "5" = "CD8_Teff", 
                     "10" = "CD8_Tex_1", "18" = "CD8_Tex_1", 
                     "3" = "CD8_Tex_2", "13" = "CD8_Tex_2", 
                     "9" = "CD8_Tex_3", 
                     "12" = "CD8_Tex_4", "7" = "CD8_Tex_4", 
                     "6" = "CD8_Tex_5", 
                     "2" = "CD8_Tex_6", "4" = "CD8_Tex_6", 
                     "17" = "CD8_Prolif_1", "16" = "CD8_Prolif_2", "19" = "CD8_Prolif_3", 
                     "8" = "CD4_Naive", "1" = "CD4_Treg", "0" = "CD4_Mem", "21" = "CD4_Th1_like", 
                     "14" = "LowQual", "22" = "LowQual"
                     )

scDLBCL.T.anno@meta.data <- scDLBCL.T.anno@meta.data %>%
        mutate(cell_type = Idents(scDLBCL.T.anno), 
               seurat_clusters = cell_type,
               orig.ident = factor(orig.ident, levels = c(s_list, "reference")))


pdf(paste0(out_dir, "Fig2e.pdf"), width = 5.5, height = 4)
DimPlot(scDLBCL.T.anno, reduction = "umap", label = F, repel = F, cols = color_T)
dev.off()



## 3.4 stat ----
p1 <- plot_stat(scDLBCL.T.anno %>% subset(time == "ZLQ"), plot_type = "prop_fill", group_by = "group", pal_setup = color_T) + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + labs(x = "")
p2 <- plot_stat(scDLBCL.T.anno %>% subset(group == "PD"), plot_type = "prop_fill", group_by = "time", pal_setup = color_T) + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + labs(x = "")
p3 <- plot_stat(scDLBCL.T.anno, plot_type = "prop_fill", group_by = "orig.ident", pal_setup = color_T) + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
          axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
pdf(paste0(out_dir, "Fig2f.pdf"), width = 9, height = 5)
p1 + p2 + p3 + plot_layout(widths = c(0.4, 0.55, 1.7), guides = "collect")
dev.off()


## 3.5 calculate dot for each cluster----
index.T <-  c(
            "CD3D", # "CD3E", "CD3G", # T cells
            "CD8A", # "CD8B", #CD8 
            "CD4", "CD40LG", # CD4
            "SELL", "CCR7", "TCF7", # Naive
            "FGFBP2", "GZMH", "GZMB", "KLRD1", "ANXA1", # Effector
            "GZMK", "CCL5", "NKG7", 
            "CXCL13", "LAG3", "HAVCR2", "TIGIT", "CTLA4", "PDCD1",# exhaustion
            "KLRB1", "GNLY", "FCGR3A", # cytotoxic
            "FOXP3", "IL2RA", "TNFRSF4", # Treg
            "IL7R", 
            "LTB", "FLT3LG", "AQP3",
            "TYMS", "MKI67", "TOP2A", "UBE2C" # Cell cycle
            )

pdf(paste0(out_dir, "SFig3d.pdf"), width = 6, height = 8)
StackedVlnPlot(obj = scDLBCL.T.anno, features = unique(index.T), cols = color_T) 
dev.off()


pdf(paste0(out_dir, "SFig3e.pdf"), width = 12, height = 14)
FeaturePlot(scDLBCL.T.anno, features = c("CD3D", "CD8A", "CD4", "SELL", "CCR7", "FOXP3", "IL2RA", "FGFBP2", "GZMH", "GZMK", "LAG3", "HAVCR2", "PDCD1", "GNLY", "FCGR3A", "MKI67"),
            ncol = 4, pt.size = 0.01, reduction = "umap", order = T, raster = T) & 
            scale_color_gradient2(midpoint = 1, low = "#313695", mid = "#FFFFBF", high = "#A50026") & 
            theme_void() & theme(panel.border = element_rect(fill = NA, colour = "black", size = 1.5), 
                                 axis.line = element_blank(),
                                 plot.title = element_text(size = 14, vjust = -4), legend.position = "top",
                                 plot.margin = unit(c(1, 1, 1, 1), "mm"))
dev.off()




# 3.5. signature scores violin: AddModuleScore ----
cd_features <- list(
    Immunosuppressive_signature_score = c("CD244", "CD160", "CTLA4", "PDCD1", "TIGIT", "LAYN", "LAG3", "HAVCR2", "CD274", "CD47", "CD96", "ENTPD1", "VSIR", "BTLA", "EBI3", "IL2RB", "IL2RA", "IL2RG"), 
    Dysfunctional_signature_score = c("LAG3", "HAVCR2", "PDCD1", "PTMS", "FAM3C", "IFNG", "AKAP5", "CD7", "PHLDA1", "ENTPD1", "SNAP47", "TNS3", "CXCL13", "RDH10", "DGKH", "KIR2DL4", "LYST", "MIR155HG", "RAB27A", "CSF1", "TNFRSF9", "CTLA4", "CD27", "CCL3", "ITGAE", "PAG1", "TNFRSF1B", "GALNT1", "GBP2", "MYO7A", "TIGIT")
    )
sig_score <- AddModuleScore(object = scDLBCL.T.anno, features = cd_features, assay = "RNA", nbin = 20, ctrl = 100, name = names(cd_features))


# 2.6. signature scores heatmap: AddModuleScore ----
data2 <- rbind(
            sig_score@meta.data %>% group_by(cell_type) %>% summarise(value = median(Immunosuppressive_signature_score1)) %>% 
                mutate(sig_score = "Immunosuppressive_signature_score"),
            sig_score@meta.data %>% group_by(cell_type) %>% summarise(value = median(Dysfunctional_signature_score2)) %>% 
                mutate(sig_score = "Dysfunctional signature score")
            ) %>% spread(cell_type, value) %>% column_to_rownames("sig_score")

pdf(paste0(out_dir, "Fig4a.pdf"), width=7, height=4)
ht <- Heatmap(t(scale(t(data2))), show_column_names = T, 
        row_split = c(1, 2), 
        heatmap_legend_param = list(title = "Scaled module score", direction = "horizontal", title_position = "lefttop"),
        cluster_rows = F, cluster_columns = F,
        width = ncol(data2)*unit(5, "mm"), height = nrow(data2)*unit(8, "mm"),
        rect_gp = gpar(col = "black"), row_names_side = "left", 
        col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)
        )
draw(ht, heatmap_legend_side="top")
dev.off()



### 
my_list <- c("PD_ZLQ", "PD_ZLZ", "PD_FF")
my_c <- c("CD8_Teff", "CD8_Tex_1", "CD8_Tex_2", "CD8_Tex_3", "CD8_Tex_4", "CD8_Tex_5", "CD8_Tex_6")

data <- sig_score@meta.data %>% filter(cell_type %in% my_c) %>% 
        filter(group %in% "PD") %>% 
        mutate(cell_type = factor(cell_type, levels = my_c), 
               my_group = factor(paste0(group, "_", time), levels = my_list))
res <- data.frame()
# my_comparisons = combn(my_list, 2, simplify = F)[c(1, 4, 6, 5)]
my_comparisons = combn(my_list, 2, simplify = F)[c(1, 3, 2)]
pdf(paste0(out_dir, "Fig4b.pdf"), width = 9, height = 3)
for (i in paste0(names(cd_features), 1:2)[1:2]){
    stat.s.data1 <- data %>% group_by(cell_type) %>% rstatix::wilcox_test(eval(parse(text = paste0(i, " ~ my_group"))), comparisons = my_comparisons) %>% 
          rstatix::add_xy_position(x = "my_group", dodge = 1)

    res <- rbind(res, stat.s.data1)
    p <- ggplot() +
          geom_violin(data = data, aes_string(x = "my_group", y = i, fill = "my_group"),
                       position = position_dodge(1), scale = "width") + 
          geom_boxplot(data = data, aes_string(x = "my_group", y = i, fill = "my_group"), 
                       position = position_dodge(1), width = 0.3, outlier.shape = NA) + 
          stat_summary(data = data, aes_string(x = "my_group", y = i, fill = "my_group"), fun = "median", 
                       geom = "point", shape = 16, size = 2.2, color = "white", position = position_dodge(1)) + 
          facet_wrap(~ cell_type, ncol = 7) + 
          # scale_fill_manual(values = c("#91B897", "#F0907F", "#C7E8FA", "#C28DBD")) + 
          scale_fill_manual(values = c("#F0907F", "#C7E8FA", "#C28DBD")) + 
          ggpubr::stat_pvalue_manual(stat.s.data1, label = "p.adj.signif", 
                        tip.length = 0.01, bracket.nudge.y = 0, size = 2.8) +
          labs(x = "", y = "Module score", title = i, fill = "") + ggthemes::theme_few(base_size = 10) + 
          theme(plot.title = element_text(hjust = 0.5, size = 10), 
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                legend.position = "right")
    print(p)
}
dev.off()






## anno merge ----
scDLBCL.T.anno2 <- RenameIdents(scDLBCL.T, 
                     "15" = "NK", "11" = "NK", "20" = "NK", 
                     "5" = "CD8_Teff", 
                     "10" = "CD8_Tex", "18" = "CD8_Tex", 
                     "3" = "CD8_Tex", "13" = "CD8_Tex", 
                     "9" = "CD8_Tex", 
                     "12" = "CD8_Tex", "7" = "CD8_Tex", 
                     "6" = "CD8_Tex", 
                     "2" = "CD8_Tex", "4" = "CD8_Tex", 
                     "17" = "CD8_Prolif", "16" = "CD8_Prolif", "19" = "CD8_Prolif", 
                     "8" = "CD4_Naive", "1" = "CD4_Treg", "0" = "CD4_Mem", "21" = "CD4_Th1_like", 
                     "14" = "LowQual", "22" = "LowQual"
                     )

scDLBCL.T.anno2@meta.data <- scDLBCL.T.anno2@meta.data %>%
        mutate(cell_type = Idents(scDLBCL.T.anno2), 
               seurat_clusters = cell_type,
               orig.ident = factor(orig.ident, levels = s_list))

## CD8_Teff
my.index.T <- c("IFNGR1", "TNFSF14", "KLRG1", "CD28", "CXCR3", "CCL3L1", "KLRB1", "GZMK", "CD40", 
                "KLF13", "HMGN1", "KLF12", "KLF3", "TFDP2", "ZFP36L2", "HOPX", "HSF1", "CEBPD")
my.dot.T <- DotPlot(object = scDLBCL.T.anno2 %>% subset(group == "PD") %>% subset(idents = c("CD8_Teff")), 
    assay = "RNA", features = my.index.T, cols = c("white", "red"), col.min = -0.2, group.by = "time")
heatdat.T <- my.dot.T$data %>% dplyr::select(avg.exp, features.plot, id) %>% spread(id, avg.exp) %>% 
           column_to_rownames("features.plot") %>% arrange(desc(row_number()))
pdf(paste0(out_dir, "Fig4f_1.pdf"), width=3, height=5)
Heatmap(t(scale(t(heatdat.T[my.index.T, ]))), name = "Z-score", show_column_names = T,
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),
        cluster_rows = F, cluster_columns = F, 
        width = ncol(heatdat.T)*unit(5, "mm"), height = nrow(heatdat.T)*unit(5, "mm"),
        rect_gp = gpar(col = NA), row_names_side = "right", 
        row_split = rep(c("Effector molecular", "Transcription factors"), c(9, 9)), 
        row_title_gp = gpar(fill = c("#9ECAE1", "#F1B6DA"), font = 1, border = F, fontsize = 9), 
        col = circlize::colorRamp2(c(-1.5, 0, 1.5), c("#4575b4","white","#d73027"))
        )
dev.off()


## CD8_Tex
my.index.T <- c("CCL3", "PHLDA1", "IFNG", "LYST", "TIGIT", "TNFRSF1B", "CTLA4", "HAVCR2", "PDCD1", "LAG3", "LAYN",
                "TSC22D1", "PRDM1", "IRF7", "IRF9", "TOX", "RBPJ", "TRPS1", "TOX2", "RFX5", "BATF")
my.dot.T <- DotPlot(object = scDLBCL.T.anno2 %>% subset(group == "PD") %>% subset(idents = c("CD8_Tex")), 
    assay = "RNA", features = my.index.T, cols = c("white", "red"), col.min = -0.2, group.by = "time")
heatdat.T <- my.dot.T$data %>% dplyr::select(avg.exp, features.plot, id) %>% spread(id, avg.exp) %>% 
           column_to_rownames("features.plot") %>% arrange(desc(row_number()))
pdf(paste0(out_dir, "Fig4f_2.pdf"), width=3, height=5)
Heatmap(t(scale(t(heatdat.T[my.index.T, ]))), name = "Z-score", show_column_names = T,
        row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),
        cluster_rows = F, cluster_columns = F, 
        width = ncol(heatdat.T)*unit(5, "mm"), height = nrow(heatdat.T)*unit(5, "mm"),
        rect_gp = gpar(col = NA), row_names_side = "right", 
        row_split = rep(c("Inhibitor molecular", "Transcription factors"), c(11, 10)), 
        row_title_gp = gpar(fill = c("#91D1C2", "#F1B6DA"), font = 1, border = F, fontsize = 9), 
        col = circlize::colorRamp2(c(-1.5, 0, 1.5), c("#4575b4","white","#d73027"))
        )
dev.off()


 
 



# 3.7. monocle2 ----
library(monocle)
library(ggridges)
library(RColorBrewer)


   

# ## CD8 ----
scDLBCL.T.anno.CD8 <- scDLBCL.T.anno %>% subset(cell_type %in% c("CD8_Teff", "CD8_Tex_1", "CD8_Tex_2", "CD8_Tex_3", 
                                                "CD8_Tex_4", "CD8_Tex_5", "CD8_Tex_6", "CD8_Prolif_1", "CD8_Prolif_2", "CD8_Prolif_3")) %>%
    subset(group == "PD")
cds.T.anno <- my_monocle2(scDLBCL.T.anno.CD8)



pdf(paste0(out_dir, "OFig.trajectory.T_State.pdf"), width = 11, height = 3)
plot_cell_trajectory(cds.T.anno, color_by = "State", show_branch_points = F,
    show_tree = T, cell_size = 0.5) + facet_wrap(~State, nrow = 1)
dev.off()



pData(cds.T.anno)$my_group <- factor(paste0(pData(cds.T.anno)$group, "_", pData(cds.T.anno)$time),
                              levels = c("PD_ZLQ", "PD_ZLZ", "PD_FF"))
pData(cds.T.anno)$Pseudotime <- max(pData(cds.T.anno)$Pseudotime) - pData(cds.T.anno)$Pseudotime
pdf(paste0(out_dir, "Fig4d.pdf"), width = 4, height = 3)
p1 <- plot_cell_trajectory(cds.T.anno, color_by = "Pseudotime", show_branch_points = F, 
        show_tree = T, cell_size = 0.5) + theme(legend.position = "right") + scale_color_viridis_c()
p1
dev.off()


pdf(paste0(out_dir, "Fig4e_1.pdf"), width = 12.5, height = 3.5)
p1 <- plot_cell_trajectory(cds.T.anno, color_by = "seurat_clusters", show_branch_points = F, 
        show_tree = T, cell_size = 0.5) + theme(legend.position = "right") + 
        theme(legend.position = "none") + 
        scale_colour_manual(values = color_T[1:11]) 
p2 <- plot_cell_trajectory(cds.T.anno, 
        color_by = "Ident", show_branch_points = F, 
        show_tree = T, cell_size = 0.5) + theme(legend.position = "right") + facet_wrap(~my_group, nrow = 1) + 
        scale_colour_manual(values = color_T[1:11])
p1 + p2 + plot_layout(widths = c(0.5, 1.55), guides = "collect")
dev.off()



## Branch state statistics
# CD8
df_cds.T.anno <- cds.T.anno@phenoData@data %>% 
        mutate(Branch = ifelse(State %in% c("1", "2", "9"), "State 3",
                            ifelse(State %in% c("3", "4", "5"), "State 2", "State 1")))

mydata_g1 <- df_cds.T.anno %>% group_by(seurat_clusters, Branch) %>% summarise(N=n()) %>%  
              group_by(seurat_clusters) %>% mutate(Freq = (N/sum(N)) * 100, Cumsum = cumsum(Freq))
mydata_g2 <- df_cds.T.anno %>% 
              group_by(my_group, Branch) %>% summarise(N=n()) %>%  
              group_by(my_group) %>% mutate(Freq = (N/sum(N)) * 100, Cumsum = cumsum(Freq))

p1 <- ggplot(mydata_g1, aes(x = Branch, y = Freq/100, fill = seurat_clusters)) + 
      geom_col(position = position_dodge2(width = 0.8, preserve = "single")) + 
      scale_y_continuous(limits = c(0,1), labels = scales::percent) + theme_classic() + labs(x = "", y = "Percentage") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      scale_fill_manual(values = color_T[1:11])
p2 <- ggplot(mydata_g2, aes(x = Branch, y = Freq/100, fill = my_group)) + 
      geom_col(position = position_dodge2(width = 0.8, preserve = "single")) + 
      scale_y_continuous(limits = c(0,1), labels = scales::percent) + theme_classic() + labs(x = "", y = "Percentage") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      scale_fill_manual(values = c("#F0907F", "#C7E8FA", "#C28DBD"))

pdf(paste0(out_dir, "Fig4e_2.pdf"), width = 11, height = 4)
p1 + p2 + plot_layout(widths = c(2.7, 1), guides = "collect")
dev.off()








# Myeloid cell ######################################################################################################
## Myeloid cell sub cluster analysis----
scDLBCL.M <- subset(scDLBCL, idents = c("12", "21")) %>% 
                    NormalizeData() %>% 
                    FindVariableFeatures(nfeatures = 2000) %>% 
                    ScaleData(vars.to.regress = c("nCount_RNA", "percent_mt")) %>% 
                    RunPCA() %>% 
                    RunHarmony(group.by.vars = "orig.ident") %>% 
                    RunUMAP(reduction = "harmony", dims = use.pcs) %>% 
                    FindNeighbors(reduction = "harmony", dims = use.pcs) %>% 
                    FindClusters(resolution = 1)

# scDLBCL.M <- read_rds(paste0("~/project/scYZX/scDLBCL/analysis/out_10s_20221030/v1_M_2000HVG_30PCs/", "Table4.scDLBCL.M_1raw.rds"))



## 4.4 Myeloid cell annotation ----
color_M <- c("#8DD3C7", "#CCEBC5", "#BEBADA", "#FB8072", "#61A2CA", "#FBAD82", "#BFBFBF")
scDLBCL.M.anno <- RenameIdents(scDLBCL.M, 
                    "11" = "M_c01_cDC1_CLEC9A", "9" = "M_c02_cDC2_CLEC10A", 
                    "14" = "M_c03_DC_LAMP3", "5" = "M_c04_pDC_CLEC4C", 
                    "1" = "M_c05_MP_IL1B", 
                    "0" = "M_c06_MP_C1QB", "4" = "M_c06_MP_C1QB", "6" = "M_c06_MP_C1QB", "12" = "M_c06_MP_C1QB",
                    "3" = "M_c06_MP_C1QB", "7" = "M_c06_MP_C1QB",
                    "2" = "M_c07_LowQual", "8" = "M_c07_LowQual", "10" = "M_c07_LowQual", "13" = "M_c07_LowQual")

scDLBCL.M.anno@meta.data <- scDLBCL.M.anno@meta.data %>%
        mutate(cell_type = as.character(Idents(scDLBCL.M.anno)), 
               seurat_clusters = cell_type,
               orig.ident = factor(orig.ident, levels = s_list))


pdf(paste0(out_dir, "Fig2c.pdf"), width = 6.5, height = 4)
DimPlot(scDLBCL.M.anno, reduction = "umap", label = F, repel = T, cols = color_M)
dev.off()



## 4.5 calculate dot for each cluster----
index.M <-  c(
            "CST3", "LYZ", # Pan-Myeloid
            "CLEC9A", "XCR1", # cDC1
            "CLEC10A", "CD1C", # cDC2            
            "LAMP3", "IDO1", "CLEC4A", "CCL19", # LAMP3 cDC
            "CLEC4C", "IRF8", "TCF4", "IRF7", # pDC            
            "IL1B", "CD86", # M1
            "S100A8", "S100A9", "FCN1", 
            "C1QB", "C1QA", "APOE", "CD68", "CD163", "CD14", # M2
            "GPNMB"
)

cluster.dot.M <- DotPlot(object = scDLBCL.M.anno, assay = "RNA", features = rev(index.M), 
                       cols = c("white", "red"), col.min = -0.2) + 
                    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + coord_flip()
pdf(paste0(out_dir, "SFig3b.pdf"), width = 5, height = 7)
cluster.dot.M
dev.off()


pdf(paste0(out_dir, "SFig3c.pdf"), width = 12, height = 10.5)
FeaturePlot(scDLBCL.M.anno, features = c("CLEC9A", "XCR1", "CLEC10A", "CD1C", "LAMP3", "IDO1", "CLEC4C", "IRF8", "IL1B", "C1QB", "APOE", "GPNMB"),
            ncol = 4, pt.size = 0.01, order = T, reduction = "umap", raster = F) & 
            scale_color_gradient2(midpoint = 0.8, low = "#313695", mid = "#FFFFBF", high = "#A50026") & 
            theme_void() & theme(panel.border = element_rect(fill = NA, colour = "black", size = 1.5), 
                                 axis.line = element_blank(),
                                 plot.title = element_text(size = 14, vjust = -4), legend.position = "top",
                                 plot.margin = unit(c(1, 1, 1, 1), "mm"))
dev.off()







# cellchat ######################################################################################################
library(CellChat)
library(Seurat)
library(patchwork)
library(tidyverse)
library(readr)
library(doMC)
registerDoMC(cores=30)
options(stringsAsFactors = F)


g_color <- c("#91B897", "#F0907F", "#C7E8FA", "#C28DBD")
my_interaction <- read.xlsx(paste0("~/project/scYZX/scDLBCL/analysis/out_10s_20221030/cellchat/", "cellchat_localDB_CellChatDB.human.xlsx")) %>% 
    mutate(ID = interaction_name) %>% column_to_rownames("ID")


CellChatDB.local <- CellChatDB.human
CellChatDB.local$interaction <- my_interaction


# 1. CellChat object & merge together----
my_anno <- rbind(
    Idents(scDLBCL.M.anno) %>% data.frame() %>% setNames("cell_type"),
    Idents(scDLBCL.T.anno2) %>% data.frame() %>% setNames("cell_type")
) %>% rownames_to_column("ID") %>%
    dplyr::filter(!cell_type %in% c("M_c07_LowQual", "LowQual"))


scDLBCL.sub <- subset(scDLBCL, cells = my_anno$ID)
scDLBCL.sub@meta.data <- scDLBCL.sub@meta.data %>% rownames_to_column("ID") %>% left_join(., my_anno, by = "ID") %>% column_to_rownames("ID")
Idents(object = scDLBCL.sub) <- 'cell_type'
levels(scDLBCL.sub)
scDLBCL.sub$cell_type <- factor(scDLBCL.sub$cell_type, levels = levels(scDLBCL.sub))

scDLBCL.sub %>% Idents() %>% table()


scDLBCL.CR_ZLQ <- scDLBCL.sub %>% subset(group == "CR" & time == "ZLQ")
cc.CR_ZLQ <- my_cellchat(object = scDLBCL.CR_ZLQ, group.by = "cell_type", CellChatDB = CellChatDB.local)

scDLBCL.PD_ZLQ <- scDLBCL.sub %>% subset(group == "PD" & time == "ZLQ")
cc.PD_ZLQ <- my_cellchat(object = scDLBCL.PD_ZLQ, group.by = "cell_type", CellChatDB = CellChatDB.local)

scDLBCL.PD_ZLZ <- scDLBCL.sub %>% subset(group == "PD" & time == "ZLZ")
cc.PD_ZLZ <- my_cellchat(object = scDLBCL.PD_ZLZ, group.by = "cell_type", CellChatDB = CellChatDB.local)

scDLBCL.PD_FF <- scDLBCL.sub %>% subset(group == "PD" & time == "FF")
cc.PD_FF <- my_cellchat(object = scDLBCL.PD_FF, group.by = "cell_type", CellChatDB = CellChatDB.local)

object.list <- list(CR_ZLQ = cc.CR_ZLQ, PD_ZLQ = cc.PD_ZLQ, PD_ZLZ = cc.PD_ZLZ, PD_FF = cc.PD_FF)
cellchat.merge <- mergeCellChat(object.list, add.names = names(object.list))


# cellchat.merge <- read_rds(paste0("~/project/scYZX/scDLBCL/analysis/out_10s_20221030/cellchat/", "cellchat_scDLBCL.rds"))


# 2. Compare the total number of interactions and interaction strength ----
pdf(paste0(out_dir, "OFig.compareInteractions.pdf"), height = 3)
gg1 <- compareInteractions(cellchat.merge, show.legend = F, group = c(1, 2, 3, 4)) + scale_fill_manual(values = g_color)
gg2 <- compareInteractions(cellchat.merge, show.legend = F, group = c(1, 2, 3, 4), measure = "weight") + scale_fill_manual(values = g_color)
gg1 + gg2
dev.off()


# 2. interactions and interaction strength among different cell populations ----
pdf(paste0(out_dir, "OFig.netVisual_circle.pdf"), width = 14, height = 12)
weight.max <- getMaxWeight(object.list, attribute = c("idents", "count"))
par(mfrow = c(2, 2), xpd = T, mai = c(0, 0.2, 0.5, 0.2))
for (i in 1:length(object.list)) {
  groupSize <- as.numeric(table(cellchat.merge@idents[[i]]))
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], 
    vertex.weight = groupSize, edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()



cellchat <- cellchat.merge

cellchat <- computeCommunProb(cellchat,type = "truncatedMean", trim = 0.1)



# 5.Identify the upgulated and down-regulated signaling ligand-receptor pairs ----
##  selected pairs
my_pairs <- c("APOE_SORL1", "ARF1_FAS", "ARF1_HMMR", "GSTP1_HMMR", # "HSP90B1_LDLR", 
              "LGALS9_CD44", "LGALS9_CD45", "LGALS9_HAVCR2", "CD86_CTLA4", "NECTIN2_TIGIT", "SIGLEC1_SPN", "CD48_CD244A", 
              "CXCL9_CXCR3", "CXCL10_CXCR3", "IL10_IL10RA_IL10RB", 
              "MIF_CD74_CXCR4", "HBEGF_CD44", "MMP9_CD44")
pairLR.use <- data.frame(interaction_name = my_pairs)

pdf(paste0(out_dir, "Fig5b.pdf"), width = 4, height = 5)
p2 <- netVisual_bubble(cellchat, sources.use = c(6), targets.use = c(8:10), 
  pairLR.use = pairLR.use, max.quantile = 0.85, 
  comparison = c(2, 3, 4), angle.x = 90) + 
  scale_y_discrete(labels = rev(my_pairs)) 
p2
dev.off()


pdf(paste0(out_dir, "OFig.netVisual_bubble.pdf"), width = 12, height = 20)
p1 <- netVisual_bubble(cellchat, sources.use = c(5), targets.use = c(7:14), comparison = c(1, 2, 3, 4), angle.x = 90)
p2 <- netVisual_bubble(cellchat, sources.use = c(6), targets.use = c(7:14), comparison = c(1, 2, 3, 4), angle.x = 90)
p1 + p2 + plot_layout(guides = "collect")
dev.off()







