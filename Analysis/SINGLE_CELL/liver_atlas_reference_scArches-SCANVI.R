library(Seurat)
library(dplyr)
library(anndata)

GSE115469_RNA <- readRDS('/home/ubuntu/single_cell/Liver_sc_atlas/GSE115469_Normalized_with_annotation.RDS')
GSE149614_RNA <- readRDS('/home/ubuntu/single_cell/Liver_sc_atlas/GSE149614_Normalized_with_annotation.RDS')

#In order to get shared genes, we also import the query datasets:
for (i in c('GSE242889', 'GSE151530', 'GSE217235', 'GSE234241', 'GSE140228', 'GSE136103',
            'GSE202642', 'GSE125449', 'GSE162616', 'GSE182159', 'GSE169446', 'GSE236382', 'GSE124395')) {
  
  assign(paste0(i, '_RNA'), readRDS(paste0('/home/ubuntu/single_cell/Liver_sc_atlas/', i, '_Normalized.RDS')))
  print(i)
  
}



all_files = ls()[grep('_RNA', ls())]

shared_genes <- NULL

for (i in all_files[1:length(all_files)]) {
  tmp <- get(i)
  DefaultAssay(tmp) = 'RNA'
  if(length(shared_genes) == 0){
    shared_genes = rownames(tmp)
  }else{
    shared_genes <- intersect(shared_genes, rownames(tmp))
  }
}

write.csv(shared_genes, '/home/ubuntu/single_cell/Liver_sc_atlas/results/reference_combined_15_datasets_RNA_shared_genes_20231009.csv')


#Save RNA count data with only overlapped genes across training sets:
for (i in all_files) {
  tmp <- get(i)
  DefaultAssay(tmp) = 'RNA'
  assign(i, subset(tmp, features = shared_genes))
}

#combine all reference RNA count matrices:
RNA <- merge(GSE115469_RNA, GSE149614_RNA)

save.image("/home/ubuntu/single_cell/Liver_sc_atlas/results/combined_15_datasets_shared_genes_reference_RNA_count_20231009.RData")
#


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Continue with SCANVI embedding:
#

SCANVI_embed_ref <- function(save_reference_latent_file = '/home/ubuntu/single_cell/more_pro/impute_protein_pipeline/train_reference_datasets/by_cell_type/datasets/PBMC_GBM_BM_ref_RNA_latent_representations.csv',
                             save_visulization_plot_file_by_celltype = '/home/ubuntu/single_cell/more_pro/impute_protein_pipeline/train_reference_datasets/by_cell_type/datasets/visualize_PBMC_GBM_BM_SCANVI_by_celltype.pdf',
                             save_visulization_plot_file_by_patient = '/home/ubuntu/single_cell/more_pro/impute_protein_pipeline/train_reference_datasets/by_cell_type/datasets/visualize_PBMC_GBM_BM_SCANVI_by_study.pdf',
                             save_visulization_plot_file_by_study = '/home/ubuntu/single_cell/more_pro/impute_protein_pipeline/train_reference_datasets/by_cell_type/datasets/visualize_PBMC_GBM_BM_SCANVI_by_study.pdf',
                             ref_model_path = 'demo_scanvi/PBMC_GBM_BM/',
                             save_top_2000_HVGs_file = '/home/ubuntu/single_cell/more_pro/impute_protein_pipeline/train_reference_datasets/by_cell_type/datasets/top_1000_HVGs_PBMC_GBM_BM.csv'
){
  RNA <- NormalizeData(RNA, normalization.method = "LogNormalize", scale.factor = 10000)
  RNA <- FindVariableFeatures(RNA, selection.method = "vst", assay = "RNA", nfeatures = 2000)
  top2000 <- head(VariableFeatures(RNA, assay = "RNA"), 2000)
  RNA <- RNA[top2000]
  print(RNA)
  
  library(reticulate)
  use_python('/home/ubuntu/chenlab_deeplearning/chenlab_deeplearning_V2/anaconda3/envs/SCANVI/bin/python', required = T)
  sc <- import('scanpy', convert = FALSE)
  scvi <- import('scvi', convert = FALSE)
  #setwd('/home/ubuntu/single_cell/scarches/')
  setwd('/home/ubuntu/single_cell/scarches-0.4-2.0/')
  sca <- import('scarches', convert = FALSE)
  torch <- import('torch', convert = FALSE)
  remove_sparsity <- import('scarches.dataset.trvae.data_handling', convert = FALSE)
  plt <- import('matplotlib.pyplot', convert = FALSE)
  np <- import('numpy', convert = FALSE)
  os <- import('os', convert = FALSE)
  scvi$settings$progress_bar_style = 'tqdm'
  adata <- sc$AnnData(
    X   = t(as.matrix(GetAssayData(RNA,slot='counts',assay = "RNA"))), #scVI requires raw counts
    obs = RNA[[]],
    var = GetAssay(RNA)[[]]
  )
  print(adata)
  
  #SCANVI model setting:
  sc$settings$set_figure_params(dpi=as.integer(200), frameon=FALSE)
  sc$set_figure_params(dpi=as.integer(200))
  sc$set_figure_params(figsize=c(4, 4))
  torch$set_printoptions(precision=as.integer(3), sci_mode=FALSE, edgeitems=as.integer(7))
  
  condition_key = 'orig.ident'
  cell_type_key = 'cell_type'
  
  #Create SCANVI model and train it on fully labelled reference dataset:
  sca$dataset$setup_anndata(adata, batch_key=condition_key, labels_key=cell_type_key)
  
  vae = sca$models$SCVI(
    adata,
    n_layers=as.integer(2),
    encode_covariates=TRUE,
    deeply_inject_covariates=FALSE,
    use_layer_norm="both",
    use_batch_norm="none")
  
  vae$train()
  
  scanvae = sca$models$SCANVI$from_scvi_model(vae, unlabeled_category = "Unknown")
  print(paste0("Labelled Indices: ", length(scanvae$`_labeled_indices`)))
  print(paste0("Unlabelled Indices: ", length(scanvae$`_unlabeled_indices`)))
  
  scanvae$train(max_epochs=as.integer(20))
  
  #Getting the latent representation and visualization
  reference_latent = scanvae$get_latent_representation()
  
  reference_latent = as.matrix(reference_latent)
  rownames(reference_latent) = colnames(RNA)
  
  RNA[['scvi']] <- CreateDimReducObject(embeddings = reference_latent, key = "scvi_", assay = DefaultAssay(RNA))
  # Find clusters, then run UMAP, and visualize
  RNA <- FindNeighbors(RNA, dims = 1:10, reduction = 'scvi')
  RNA <- FindClusters(RNA, resolution =1)
  
  RNA <- RunUMAP(RNA, dims = 1:10, reduction = 'scvi', n.components = 2)
  
  pdf(paste0(save_visulization_plot_file_by_celltype), width = 12, height = 10)
  p = DimPlot(RNA, reduction = "umap",label=T,group.by='cell_type')
  print(p)
  dev.off()
  
  pdf(paste0(save_visulization_plot_file_by_patient), width = 12, height = 10)
  p = DimPlot(RNA, reduction = "umap",label=T,group.by='orig.ident')
  print(p)
  dev.off()
  
  RNA$orig.ident[grep('HCC', RNA$orig.ident)] = 'GSE149614'
  RNA$orig.ident[grep('TLH', RNA$orig.ident)] = 'GSE115469'
  
  pdf(paste0(save_visulization_plot_file_by_study), width = 12, height = 10)
  p = DimPlot(RNA, reduction = "umap",label=T,group.by='orig.ident')
  print(p)
  dev.off()
  
  setwd('/home/ubuntu/single_cell/Liver_sc_atlas/results/')
  
  scanvae$save(ref_model_path, overwrite=TRUE)
  
  write.csv(reference_latent, save_reference_latent_file)
  
  write.csv(as.character(py_to_r(adata$var$index)), save_top_2000_HVGs_file)
  
}

#Run function to embed ref data sets:
SCANVI_embed_ref(save_reference_latent_file = '/home/ubuntu/single_cell/Liver_sc_atlas/results/reference_combined_2_datasets_RNA_scArches-SCANVI_latent_representations_20231010.csv',
                 save_visulization_plot_file_by_celltype = '/home/ubuntu/single_cell/Liver_sc_atlas/results/visualize_reference_combined_2_datasets_RNA_scArches-SCANVI_latent_by_celltype_20231010.pdf',
                 save_visulization_plot_file_by_patient = '/home/ubuntu/single_cell/Liver_sc_atlas/results/visualize_reference_combined_2_datasets_RNA_scArches-SCANVI_latent_by_patient_20231010.pdf',
                 save_visulization_plot_file_by_study = '/home/ubuntu/single_cell/Liver_sc_atlas/results/visualize_reference_combined_2_datasets_RNA_scArches-SCANVI_latent_by_study_20231010.pdf',
                 ref_model_path = 'reference_combined_2_datasets_RNA_scArches-SCANVI_latent_20231010/',
                 save_top_2000_HVGs_file = '/home/ubuntu/single_cell/Liver_sc_atlas/results/top_2000_HVGs_reference_combined_2_datasets_RNA_scArches-SCANVI_latent_20231010.csv')

