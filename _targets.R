# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed.






# Set target options:
tar_option_set(
  packages = c("tibble", "Seurat", "ggplot2", "dplyr", "tidyr","future.apply", "DoubletFinder"),  # Packages that your targets need for their tasks.
  controller = crew::crew_controller_local(workers = 2, seconds_idle = 60),
  trust_timestamps = TRUE
)

options(future.globals.maxSize= 10^10)

# Run the R scripts in the R/ folder with your custom functions:
tar_source()
# tar_source("other_functions.R") # Source other scripts as needed.


# Replace the target list below with your own:
c(getparams(),
  list(
  tar_target(
    name = data_list,
    command = list.files("data/h5", pattern = "h5", full.names = TRUE),
    format = "file"
  ),
  tar_target(
    name = metadata_file,
    command = "data/metadata.csv",
    format = "file"
  ),
  tar_target(
    name = metadata,
    command = read_metadata(metadata_file, required = c("Sample"))
  ),
  tar_target(
    name = obj_orig,
    command = load_h5(data_list, metadata, join_layers = load_joinlayers,
                      qc_mitopattern = qc_mitopattern, 
                      qc_ribopattern = qc_ribopattern)
  ),
  tar_target(
    name = obj_orig_meta,
    command = obj_orig@meta.data
  ),
  
  tar_target(
    name = qc_unfiltered_vln,
    command = qc_vln(obj_orig, qc_groupby = qc_groupby, 
                     thresholds = list(percent.mito = qc_mitoMax,
                                       percent.ribo = qc_riboMin,
                                       nFeature_RNA = qc_nFeatureMin,
                                       nCount_RNA = qc_nCountMax))
  ),
  tar_target(
    name = obj_filt,
    command = filter_obj(obj_orig, qc_mitoMax, qc_riboMin, qc_nFeatureMin,
                         qc_nCountMax,qc_minCells)
  ),
  tar_target(
    name = qc_filtered_vln,
    command = qc_vln(obj_filt, qc_groupby = qc_groupby)
  ),
  tar_target(
    top_expr_boxplot,
    checkTopExpressedGenes(obj_filt, qc_maxCellsTopExpr)
  ),
  tar_target(
    obj_filt_meta,
    obj_filt@meta.data
  ),
  tar_target(
    obj_filt2,
    remove_genes(obj_filt,patterns = qc_genePatternsToRemove,genes = qc_genesToRemove)
  ),
  tar_target(
    obj_filt2_meta,
    obj_filt2@meta.data
  ),
  tar_target(
    CC_s_genes,
    read_param_file("parameters/cc_S_genes.csv"),
    format = "file"
  ),
  tar_target(
    CC_g2m_genes,
    read_param_file("parameters/cc_G2M_genes.csv"),
    format = "file"
  ),
  tar_target(
    obj_cc,
    cellCycleScores(obj_filt2, slist = CC_s_genes, g2mlist = CC_g2m_genes)
  ),
  tar_target(
    cc_violin,
    plot_CC(obj_cc, qc_groupby)
  ),
  tar_target(
    obj_df,
    if(qc_runDF){
      findDoublets(obj_cc)
    }else{
      obj_cc
    }
  ),
  tar_target(
    doublet_plots,
    if(qc_runDF){
      plotDF(obj_df, group.by = qc_groupby, DFrun = qc_runDF)
    }else{
      obj_cc
    }
  ),
  tar_target(
    obj_filt_final,
    if(qc_runDF){
      removeDoublets(obj_df)
    }else{
      obj_cc
    }
  ), 
  tar_target(
    obj_normalized,
    if(dimred_sct){
      SCTransform(obj_filt_final, vars.to.regress = dimred_varstoregress, 
                  variable.features.n = dimred_nHVG)
    }else{
      NormalizeData(obj_filt_final) %>%
        FindVariableFeatures(nfeatures = dimred_nHVG) %>%
        ScaleData(vars.to.regress = dimred_varstoregress)
    }
  ),
  tar_target(
    obj_pca,
    RunPCA(obj_normalized)
  ),
  tar_target(
    obj_pca_umap,
    RunUMAP(obj_pca, n.neighbors = dimred_umap_nn,
            n.components= dimred_umap_ncomp,
            min.dist = dimred_umap_mindist,
            dims = dimred_umap_dims)
  )
  
))
