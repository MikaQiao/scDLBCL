#######################################################################################
library(DoubletFinder)
my_doublet_v2 <- function(Pool, my_rate){
    ## Pre-process Seurat object (sctransform) -----------------------------------------------------------------------------------
    use.pcs = 1:30
    Pool <- Pool %>% SCTransform() %>% RunPCA() %>% RunUMAP(dims = use.pcs)

    # Find best PK
    sweep.res.list <- paramSweep_v3(Pool, PCs = use.pcs, sct = TRUE)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))

    # Find best nExp
    annotations <- Pool@meta.data$seurat_clusters
    homotypic.prop <- modelHomotypic(annotations)
    nExp_poi <- round(my_rate * nrow(Pool@meta.data))  ## Assuming 10% doublet formation rate - tailor for your dataset
    nExp_poi.adj <- round(nExp_poi*(1 - homotypic.prop))

    # Find Doublet
    Pool <- doubletFinder_v3(Pool, PCs = use.pcs, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
    table(Pool@meta.data[grepl("DF.classifications", colnames(Pool@meta.data))])
    
    my.reuse.pANN <- colnames(Pool@meta.data)[grepl("pANN_", colnames(Pool@meta.data))]
    Pool <- doubletFinder_v3(Pool, PCs = use.pcs, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = my.reuse.pANN, sct = TRUE)
    my.df.class <- colnames(Pool@meta.data)[grepl("DF.classifications", colnames(Pool@meta.data))][1]
    print(table(Pool@meta.data[my.df.class]))

    # # visualization
    Pool@meta.data <- Pool@meta.data %>% mutate(DF_hi.lo = .[[my.df.class]])
    print(table(Pool@meta.data$DF_hi.lo))
    # DimPlot(Pool, group.by = "DF_hi.lo")
    return(Pool)
}




#######################################################################################
library(Seurat)
library(patchwork)
library(ggplot2)

## remove the x-aNxis text and tick 
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot <- function(obj, feature, cols, pt.size = 0, 
                           plot.margin = unit(c(-1, 0, -1, 0), "cm"),
                           ...) {
  p <- VlnPlot(obj, features = feature, cols = cols, pt.size = pt.size, ... )  + 
              labs(x = "", y = feature, title = NULL) + 
              theme(legend.position = "none", 
                    # panel.border = element_rect(fill = NA, colour = "black", size = 1.5),
                    # axis.line = element_blank(), 
                    axis.title.x = element_blank(), 
                    axis.text.x = element_blank(), 
                    axis.ticks.x = element_blank(),
                    axis.title.y = element_text(size = rel(1), vjust = 0.5, angle = 0), 
                    axis.text.y = element_blank(), 
                    # axis.ticks.y = element_blank(),
                    plot.margin = plot.margin) + 
              aes(color = obj$seurat_clusters) + 
              scale_color_manual(values = cols)
  return(p)
}

## extract the max value of the y axis
extract_max <- function(p){
  ymax <- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function: v1
StackedVlnPlot <- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-1, 0, -1, 0), "cm"),
                          ...) {
  
  plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj, feature = x, plot.margin = unit(c(0, 0, 0, 0), "mm"), ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]] <- plot_list[[length(plot_list)]] +
    theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1), axis.ticks.x = element_line())
  
  # # # change the y-axis tick to only max value 
  # ymaxs <- purrr::map_dbl(plot_list, extract_max)
  # plot_list <- purrr::map2(plot_list, ymaxs, function(x,y) x + 
  #                           scale_y_continuous(breaks = c(y)) + 
  #                           expand_limits(y = y))
  
  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}





#######################################################################################
my_monocle2 <- function(object){
    ### read input
    data <- as(object@assays$RNA@counts, "sparseMatrix")
    fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
    pd <- new("AnnotatedDataFrame", data = cbind(object@meta.data, Ident = object@active.ident))
    fd <- new("AnnotatedDataFrame", data = fData)

    # if(all(data == floor(data))) {
    #     expressionFamily <- negbinomial.size() # choose
    # } else if(any(data < 0)){
    #     expressionFamily <- uninormal()
    # } else {
    #     expressionFamily <- tobit()
    # }
    expressionFamily <- negbinomial.size()

    my_cds <- newCellDataSet(data, phenoData = pd, featureData = fd, lowerDetectionLimit = 0.5, 
                          expressionFamily = expressionFamily) 
    my_cds <- estimateSizeFactors(my_cds)
    my_cds <- estimateDispersions(my_cds)

    ### choose definition genes
    ordering_genes <- object[["RNA"]]@var.features      # object@assays$RNA@var.features
    my_cds <- setOrderingFilter(my_cds, ordering_genes)
    # plot_ordering_genes(my_cds)
    dim(exprs(my_cds))

    ### reduce & sort & visual
    my_cds <- reduceDimension(my_cds, reduction_method = "DDRTree", max_components = 4)
    my_cds <- orderCells(my_cds)

    return(my_cds)
}






#######################################################################################
my_cellchat <- function(object, group.by, CellChatDB){
  # 1. Create a CellChat object
  cellchat <- createCellChat(object = object, group.by = group.by)

  ## Set the ligand-receptor interaction database: Secreted Signaling, ECM-receptor, Cell-Cell contact
  # CellChatDB <- CellChatDB.human
  # unique(CellChatDB$interaction$annotation)
  cellchat@DB <- CellChatDB 

  ## Preprocessing the expression data for cell-cell communication analysis
  cellchat <- subsetData(cellchat) # necessary even if using the whole database: Secreted Signaling, Cell-Cell Contact, ECM-Receptor

  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)  # project gene expression to PPI network (optional)


  # 2. Inference of cell-cell communication network ----
  ## Compute the communication probability and infer cellular communication network
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  ## Infer the cell-cell communication at a signaling pathway level
  cellchat <- computeCommunProbPathway(cellchat)
  ## Calculate the aggregated cell-cell communication network
  cellchat <- aggregateNet(cellchat)

  # Compute the network centrality scores
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

  return(cellchat)
}

