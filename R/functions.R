
load_h5 = function(samples, metadata, join_layers = TRUE, 
                   qc_perc_patterns){
  
  names(samples) = gsub(".*\\/","",gsub("\\.h5","",samples))
  
  sample_list = lapply(samples, function(sample){
    CreateSeuratObject(Read10X_h5(sample),
                       project =  gsub(".*\\/","",gsub("\\.h5","",sample)))
  })
  
  return(combine_objects(sample_list, metadata, join_layers) %>%
           add_qc(qc_perc_patterns))
  
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

add_qc = function(obj, qc_perc_patterns){
  for(n in names(qc_perc_patterns)){
    obj = AddMetaData(obj,
                      PercentageFeatureSet(obj, pattern = qc_perc_patterns[[n]]),
                      n)
  }
  obj
}


qc_vln = function(obj, qc_groupby, qc_plotvars, thresholds = list(), pt.size = 0,...){
  vp = VlnPlot(obj, 
          features = qc_plotvars,
          group.by = qc_groupby, layer = "counts",
          pt.size = pt.size, ncol = 2, ...) & violintheme()
  if(length(thresholds)>0){
    variables = sapply(vp, function(x){colnames(x$data)[1]})
    
    for(thresh in names(thresholds)){
      if(thresh %in% variables){
        vp[[match(thresh, variables)]] = 
          vp[[match(thresh, variables)]] & 
          geom_hline(yintercept = thresholds[[thresh]], linetype = 2)
      }else{
        warning(paste0(thresh, " is not a variable in the violin plots"))
      }
    }
  }
  
  addColors(vp, scale_type = "fill")
}

filter_obj = function(obj, qc_filt_rules,qc_minCells){
  
  em = FetchData(obj, Features(obj))>0
  
  ruledf = as.data.frame(lapply(qc_filt_rules, function(rule){
    do.call(rule[[2]], list(FetchData(obj, rule[[1]]), rule[[3]]))
  }))
  
  obj = subset(obj, cells = Cells(obj)[apply(ruledf, 1, min)>0],
               features = rownames(obj)[
                 Matrix::colSums(em) >
                   qc_minCells])
}

filter_obj_old = function(obj, qc_mitoMax, qc_riboMin, qc_nFeatureMin,
                      qc_nCountMax,qc_minCells){
  
  em = FetchData(obj, Features(obj))>0
  
  obj = subset(obj, subset = percent.mito < qc_mitoMax &
                 percent.ribo > qc_riboMin &
                 nFeature_RNA > qc_nFeatureMin &
                 nCount_RNA < qc_nCountMax,
               features = rownames(obj)[
                 Matrix::colSums(em) >
                   qc_minCells])
}



remove_genes = function(obj, patterns, genes){
  
  patterns = collate_gene_patterns(patterns, genes)
  
  if(length(patterns)>0){
    allgenes = unique(unlist(lapply(obj@assays, Features)))
    
    toremove = unique(unlist(lapply(patterns, function(pattern){
      Features(obj)[grepl(pattern,Features(obj))]
    })))
    goodfeatures = allgenes[!allgenes %in% toremove]
    
    obj = subset(obj, 
                 features = goodfeatures)
  }
  obj
}


collate_gene_patterns = function(patterns, genes){
  genes = sapply(genes, function(g){
    paste0("^",g,"$")
  })
  
  c(patterns, genes)
}

checkTopExpressedGenes = function(obj, maxCellsPerSample = NULL){
  if(!is.null(maxCellsPerSample)){
    obj = subset(obj, downsample = maxCellsPerSample)
  }
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

read_param_file = function(file){
  if(file.exists(file)){
    file_contents = read.csv(file, header = FALSE)[,1]
  }else{
    file_contents = c()
  }
  file_contents
}

cellCycleScores = function(obj, slist, g2mlist){
  if(length(slist)==0){
    slist = Seurat::cc.genes.updated.2019$s.genes
  }
  if(length(g2mlist) == 0){
    g2mlist = Seurat::cc.genes.updated.2019$g2m.genes
  }
  
  if(is.null(obj@assays$RNA@layers$counts)){
    obj_joined = JoinLayers(obj)
    
    obj_joined = runCC(obj_joined, slist, g2mlist)
    
    obj = AddMetaData(obj, obj_joined@meta.data)
  }else{
    obj = runCC(obj, slist, g2mlist)
  }
  obj
}

runCC = function(obj, slist, g2mlist){
  obj = NormalizeData(obj)
  
  obj = CellCycleScoring(object = obj,
                         g2m.features = g2mlist,
                         s.features = slist)
}

plot_CC = function(obj, group.by){
  addColors(VlnPlot(obj, features = c("S.Score", "G2M.Score"), group.by = group.by,
          pt.size = 0) & violintheme(), scale_type = "fill")
}



findDoublets = function(obj){
  
  if(! "expectedCells" %in% colnames(obj@meta.data)){
    stop('For doublet identification to be run, metadata.csv must contain the column "expectedCells", indicating for each sample how many cells were expected based on the number of cells added to the well.')
  }
  
  if(! "gemWell" %in% colnames(obj@meta.data)){
    obj$gemWell = obj$orig.ident
  }
  
  df.list = future_lapply(SplitObject(obj, split.by = "gemWell"), function(x){
    x = NormalizeData(x)
    x = FindVariableFeatures(x)
    x = ScaleData(x, vars.to.regress = c("nFeature_RNA", "percent.mito"))
    x = RunPCA(x, npcs = 20)
    x = RunUMAP(x, dims = 1:10)
    
    
    totalExpected = 
      sum(unique(x@meta.data[,c("orig.ident","expectedCells")])$expectedCells)
    expProc = c(0.004,0.008,0.016,0.024,0.032,0.04,0.048,0.056,0.064,0.072,0.08)[
      which.min(abs(totalExpected - c(500,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000)))]
    
    sweep.res <- paramSweep(x)
    sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    
    mypK = as.numeric(as.character(bcmvn$pK[
      which.max(bcmvn$BCmetric)]))
    x <- doubletFinder(x, pN = 0.25, pK = mypK, nExp = expProc*ncol(x), PCs = 1:10)
    DF.name = colnames(x@meta.data)[grepl("DF.classification", colnames(x@meta.data))]
    df_meta = x@meta.data[, DF.name, drop = FALSE]
    colnames(df_meta) = "DF"
    
    df_meta
    
  })
  
  names(df.list) = NULL
  dfmetadata = do.call("rbind",df.list)
  obj = AddMetaData(obj, dfmetadata)
  
  rm(df.list)
    
  obj
}

plotDF = function(obj, group.by, DFrun){
  if(DFrun){
    
    obj = NormalizeData(obj)
    obj = FindVariableFeatures(obj)
    obj = ScaleData(obj, vars.to.regress = c("nFeature_RNA", "percent.mito"))
    obj = RunPCA(obj, npcs = 20)
    obj = RunUMAP(obj, dims = 1:10)
    
    list(
      dimplot = DimPlot(obj, group.by = "DF"),
      vlnplot = VlnPlot(obj, "nFeature_RNA", pt.size = 0, split.by = "DF",
                        group.by = group.by)
    )
  }else{
    list(
      dimplot = NULL,
      vlnplot = NULL
    )
  }
}


removeDoublets = function(obj){
  subset(obj, subset = DF != "Doublet")
}


integrate_obj = function(obj, join_layers, split_group){
  if(! split_group %in% colnames(obj@meta.data)){
    stop("split_group must be a column name in obj@meta.data")
  }
  if(length(unique(obj@meta.data[[split_group]])) == 1){
    stop(paste0("obj has only one value of split_group parameter ", split_group, ". Either change split_group to a column name with multiple values, or skip integration by setting run_int = FALSE"))
  }
  
  if(join_layers){
    obj = SplitObject(obj, split_group)
    
    obj = lapply(obj, function(x){
      x <- NormalizeData(x, verbose = FALSE)
      x <- FindVariableFeatures(x, selection.method = "vst",
                                nfeatures = 2000, verbose = FALSE)
    })
    
    data.anchors <- FindIntegrationAnchors(object.list = obj, 
                                           dims = 1:30, reduction = "cca")
    rm(obj)
    gc()
    
    obj <- IntegrateData(anchorset = data.anchors, dims = 1:30, new.assay.name = "CCA")
    
    rm(data.anchors)
    gc()
    
  }else{
    
    obj = JoinLayers(obj)

    obj[[obj@active.assay]] <- split(obj[[obj@active.assay]], f = obj@meta.data[[split_group]])
    
    obj <- NormalizeData(obj)
    obj <- FindVariableFeatures(obj)
    obj <- ScaleData(obj)
    obj <- RunPCA(obj)
    
    obj <- IntegrateLayers(
      object = obj, method = CCAIntegration,
      orig.reduction = "pca", new.reduction = "integrated.cca",
      verbose = FALSE
    )
  }
  obj
}


moveReduction = function(obj, oldname, newname, newkey){
  oldred = obj@reductions[[oldname]]
  obj@reductions[[newname]] = CreateDimReducObject(embeddings = oldred@cell.embeddings,
                       loadings = oldred@feature.loadings,
                       projected = oldred@feature.loadings.projected,
                       assay = oldred@assay.used,
                       key = newkey
                       )
  obj@reductions[[oldname]] = NULL
  obj
}


runClustering = function(obj,reduction, clusteringres, dims){
  obj <- FindNeighbors(obj, reduction = reduction,
                       dims = dims)
  for(res in clusteringres){
    obj <- FindClusters(obj, resolution = res,
                        cluster.name = paste0("clusters_", reduction,"_",res))
  }
  obj
}


runClusterDE = function(obj, group, assay, ...){
  obj = SetIdent(obj, value = group)
  
  maxn = min(table(obj@active.ident))
  
  message(paste0("Max n for cluster DE is ", maxn))
  
  obj = subset(obj, downsample = maxn)
  
  obj@active.assay = assay
  obj = JoinLayers(obj)
  
  list(de_res = FindAllMarkers(obj, ...),
       maxn = maxn)

  
}

