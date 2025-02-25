
load_h5 = function(samples, metadata, join_layers = TRUE, 
                   qc_mitopattern = qc_mitopattern,
                   qc_ribopattern = qc_ribopattern){
  
  names(samples) = gsub(".*\\/","",gsub("\\.h5","",samples))
  
  sample_list = lapply(samples, function(sample){
    CreateSeuratObject(Read10X_h5(sample),
                       project =  gsub(".*\\/","",gsub("\\.h5","",sample)))
  })
  
  return(combine_objects(sample_list, metadata, join_layers) %>%
           add_qc(qc_mitopattern = qc_mitopattern,
                  qc_ribopattern = qc_ribopattern))
  
}

combine_objects = function(sample_list, metadata, join_layers = TRUE){
  
  if(length(sample_list)>1){
    obj = merge(sample_list[[1]], sample_list[-1],
                    add.cell.ids = names(sample_list)
    )
  }else{
    obj = sample_list[[1]]
  }
  
  if(join_layers){
    obj = JoinLayers(obj)
  }
  
  obj$orig.ident = 
    factor(obj$orig.ident, levels = names(sample_list))
  
  cellmetadata = metadata[match(obj$orig.ident, metadata$Sample),]
  rownames(cellmetadata) = Cells(obj)
  obj = AddMetaData(obj, cellmetadata)
  
  nCells = as.data.frame(table(obj$orig.ident))
  
  obj$nCells = nCells$Freq[match(obj$orig.ident, nCells$Var1)]
  
  return(obj)
}

read_metadata = function(metadatapath, required = c("Sample")){
  
  metadata = read.csv(metadatapath)
  
  if(!is.null(required)){
    if(!all(required %in% colnames(metadata))){
      stop(paste0(paste0(required[! required %in% colnames(metadata)], 
                         collapse = ","), 
                  " is a (are) required parameter(s) but cannot be found in ",
                  metadatapath))
    }
  }
  metadata
}

add_qc = function(alldata, qc_mitopattern, qc_ribopattern){
  alldata$percent.mito = 
    PercentageFeatureSet(alldata, pattern = qc_mitopattern)
  alldata$percent.ribo =
    PercentageFeatureSet(alldata, pattern = qc_ribopattern)
  alldata
}


qc_vln = function(alldata, qc_groupby){
  VlnPlot(alldata, 
          features = c("nFeature_RNA","nCount_RNA",
                       "percent.mito", "percent.ribo"),
          group.by = qc_groupby,
          pt.size = 0, ncol = 2) & violintheme()
}

filter_obj = function(obj, qc_mitoMax, qc_riboMin, qc_nFeatureMin,
                      qc_nCountMax,qc_minCells){
  obj = subset(obj, subset = percent.mito < qc_mitoMax &
                     percent.ribo > qc_riboMin &
                     nFeature_RNA > qc_nFeatureMin &
                     nCount_RNA < qc_nCountMax,
                   features = rownames(obj)[
                     Matrix::rowSums(obj@assays$RNA@layers$counts>0) >
                       qc_minCells])
}


checkTopExpressedGenes = function(obj){
  C <- FetchData(obj, assay = "RNA", layer = "counts",
                 Features(obj))
  C <- Matrix::t(C/Matrix::rowSums(C)) * 100

  most_expressed <- order(apply(C, 1, median), decreasing = T)[20:1]

  expr_df = as.data.frame(t(as.matrix(C[most_expressed, ]))) %>%
    pivot_longer(1:20, names_to = "Gene")
  expr_df$Gene = factor(expr_df$Gene, levels = rownames(C)[most_expressed])
  
  ggplot(expr_df, aes(x = value, y = Gene)) +
    geom_boxplot(outlier.size = 0.1, staplewidth = 0.5, fill = "lightcoral") +
    theme_classic() +
    labs(x = "% total count per cell")
  
}



