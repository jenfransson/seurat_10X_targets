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
      name = qc_mitopattern,
      command = "^MT-"
    ),
    tar_target(
      name = qc_ribopattern,
      command = "^RP[SL]"
    ),
    tar_target(
      name = qc_mitoMax,
      command = 20
    ),
    tar_target(
      name = qc_riboMin,
      command = 5
    ),
    tar_target(
      name = qc_nFeatureMin,
      command = 200
    ),
    tar_target(
      name = qc_nCountMax,
      command = 1000000
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
    )
  )
}