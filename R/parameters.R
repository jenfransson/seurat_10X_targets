getparams = function(){
  list(
  
    # Parameters for loading data
    tar_target(
      load_joinlayers,
      TRUE
    ),
    
    # Parameters for QC and filtering
    
    tar_target(
      name = qc_groupby,
      command = "Sample"
    ),
    tar_target(
      qc_perc_patterns,
      list(
        percent.mito = "^MT-",
        percent.ribo = "^RP[SL]",
        percent.hb = "^HB[^(P|E|S)]"
      )
    ),
    tar_target(
      qc_plotvars,
      
      c("nFeature_RNA",
        "nCount_RNA",
        names(qc_perc_patterns))
    ),
    
    
    tar_target(
      # List of rules to use when filtering cells.
      # Each rule is composed of a list, with the variable of interest (should be accessible using FetchData),
      # the operator to apply, and the threshold. An optional fourth value will be used as an explanation in the report.
      qc_filt_rules,
      list(
        list("nFeature_RNA",">",200,"Min number of unique features per cell"),
        list("nCount_RNA","<",1000000, "Max number of reads (unique UMIs) per cell"),
        list("percent.mito","<",20,
             "Max % mitochondrial reads per cell"),
        list("percent.ribo",">",5, "Min % ribosomal reads per cell")
      )
    ),
    
    
    
    tar_target(
      name = qc_minCells,
      command = 10
    ),
    tar_target(
      name = qc_maxCellsTopExpr,
      command = NULL
    ),
    tar_target(
      name = qc_genesToRemove,
      command = c("MALAT1")
    ),
    tar_target(
      name = qc_genePatternsToRemove,
      command = c("^RP[SL]", "^MT-")
    ),
    tar_target(
      name = qc_runDF,
      command = TRUE
    ),
    
    # Dimensional reduction parameters
    tar_target(
      # Should normalization be performed with SCTransform?
      dimred_sct,
      FALSE
    ),
    tar_target(
      dimred_varstoregress,
      c("nFeature_RNA")
    ),
    tar_target(
      dimred_nHVG,
      2000
      ## Seurat default: 2000
    ),
    tar_target(
      dimred_umap_nn,
      30
      ## Seurat default: 30
    ),
    tar_target(
      dimred_umap_ncomp,
      2
      ## Seurat default: 2
    ),
    tar_target(
      dimred_umap_mindist,
      0.3
      ## Seurat default: 0.3
    ),
    tar_target(
      dimred_umap_dims,
      1:20
    ),
    
    # Integration parameters
    
    tar_target(
      run_int,
      TRUE
    ),
    tar_target(
      int_split,
      "Sample"
    ),
    tar_target(
      int_umap_nn,
      30
      ## Seurat default: 30
    ),
    tar_target(
      int_umap_ncomp,
      2
      ## Seurat default: 2
    ),
    tar_target(
      int_umap_mindist,
      0.3
      ## Seurat default: 0.3
    ),
    tar_target(
      int_umap_dims,
      1:20
    ),
    
    # Clustering parameters
    
    tar_target(
      # Which reduction should be used for clustering?
      clus_reduction,
      "pca" 
    ),
    tar_target(
      # Which resolution(s) should be used for clustering?
      clus_resolutions,
      c(0.25,0.5,1,1.5)
    ),
    tar_target(
      # Which dimensions of the reduction should be used for clustering?
      clus_dims,
      1:30
    ),
    tar_target(
      # Which resolution should be used for cluster DE?
      clus_de_res,
      0.25
    ),
    tar_target(
      clus_de_logfc.threshold,
      0.1
      # Seurat default 0.1
    ),
    tar_target(
      clus_de_test.use,
      "wilcox"
      # Seurat default "wilcox"
    ),
    tar_target(
      clus_de_min.pct,
      0.01
      # Seurat default 0.01
    ),
    tar_target(
      clus_de_min.diff.pct,
      -Inf
      # Seurat default -Inf
    ),
    tar_target(
      clus_de_only.pos,
      TRUE
      # Seurat default FALSE
    )
  )
}