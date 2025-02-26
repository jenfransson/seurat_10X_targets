getparams = function(){
  list(
  
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
      name = qc_genesToRemove,
      command = c("MALAT1")
    ),
    tar_target(
      name = qc_genePatternsToRemove,
      command = c("^RP[SL]", "^MT-")
    ),
    tar_target(
      name = qc_runDF,
      command = FALSE
    )
  )
}