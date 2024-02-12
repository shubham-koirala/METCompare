library(Seurat)
library(dplyr)
library(anndata)

library(reticulate)
use_python('/home/ubuntu/chenlab_deeplearning/chenlab_deeplearning_V2/anaconda3/envs/SCANVI/bin/python', required = T)
sc <- import('scanpy', convert = FALSE)
scvi <- import('scvi', convert = FALSE)
setwd('/home/ubuntu/single_cell/scarches-0.4-2.0/')
sca <- import('scarches', convert = FALSE)
torch <- import('torch', convert = FALSE)
remove_sparsity <- import('scarches.dataset.trvae.data_handling', convert = FALSE)
plt <- import('matplotlib.pyplot', convert = FALSE)
np <- import('numpy', convert = FALSE)
os <- import('os', convert = FALSE)

load("/home/ubuntu/single_cell/Liver_sc_atlas/results/combined_15_datasets_shared_genes_reference_RNA_count_20231009.RData")

RNA_ref <- RNA
RNA_ref[['study']] = ''
RNA_ref$study[grep('HCC', RNA_ref$orig.ident)] = 'GSE149614'
RNA_ref$study[grep('TLH', RNA_ref$orig.ident)] = 'GSE115469'

#Merge query datasets:
RNA <- NULL

for (i in c('GSE242889', 'GSE151530', 'GSE217235', 'GSE234241', 'GSE140228', 'GSE136103',
            'GSE202642', 'GSE125449', 'GSE162616', 'GSE182159', 'GSE169446', 'GSE236382', 'GSE124395')) {
  
   tmp <- get(paste0(i, '_RNA'))
   tmp[['study']] = i
   rownames(tmp@meta.data) <- colnames(tmp@assays$RNA@counts)
    
   if(length(RNA) == 0){
      RNA <- tmp
    }else{
      RNA <- merge(RNA, tmp)
    }
    print(i)
  
}

#--------------------------------------------------------------------------------------------------
#changeable parameters:
#HCC, combined tumor and control group:
save_prediction_query_set_cell_type_file = '/home/ubuntu/single_cell/Liver_sc_atlas/results/liver_combined_13_query_scArches-SCANVI_latent_pred_celltype_reference_2_datasets_20231010.csv'
save_query_set_cell_latent_representations_file = '/home/ubuntu/single_cell/Liver_sc_atlas/results/liver_combined_13_query_RNA_latent_representations_scArches-SCANVI_reference_2_datasets_20231010.csv'
reference_latent_representations_file = '/home/ubuntu/single_cell/Liver_sc_atlas/results/reference_combined_2_datasets_RNA_scArches-SCANVI_latent_representations_20231010.csv'
reference_top_2000_HVGs_file = '/home/ubuntu/single_cell/Liver_sc_atlas/results/top_2000_HVGs_reference_combined_2_datasets_RNA_scArches-SCANVI_latent_20231010.csv'
ref_model_path = '/home/ubuntu/single_cell/Liver_sc_atlas/results/reference_combined_2_datasets_RNA_scArches-SCANVI_latent_20231010/'
save_query_model_path = '/home/ubuntu/single_cell/Liver_sc_atlas/results/reference_combined_2_datasets_query_combined_13_RNA_scArches-SCANVI_latent_20231010/'
save_visualize_ref_query_sets_by_cell_type_file = '/home/ubuntu/single_cell/Liver_sc_atlas/results/visualize_reference_combined_2_datasets_query_combined_13_RNA_scArches-SCANVI_by_celltype_20231010.pdf'
save_visualize_ref_query_sets_by_patient_file = '/home/ubuntu/single_cell/Liver_sc_atlas/results/visualize_reference_combined_2_datasets_query_combined_13_RNA_scArches-SCANVI_by_patient_20231010.pdf'
save_visualize_ref_query_sets_by_study_file = '/home/ubuntu/single_cell/Liver_sc_atlas/results/visualize_reference_combined_2_datasets_query_combined_13_RNA_scArches-SCANVI_by_study_20231010.pdf'
save_visualize_query_sets_by_cell_type_file = '/home/ubuntu/single_cell/Liver_sc_atlas/results/visualize_query_combined_13_RNA_scArches-SCANVI_by_celltype_20231010.pdf'
save_visualize_query_sets_by_study_file = '/home/ubuntu/single_cell/Liver_sc_atlas/results/visualize_query_combined_13_RNA_scArches-SCANVI_by_study_20231010.pdf'
save_visualize_query_sets_by_patient_file = '/home/ubuntu/single_cell/Liver_sc_atlas/results/visualize_query_combined_13_RNA_scArches-SCANVI_by_patient_20231010.pdf'

#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
#Below are the common code for all external validation sets:
#But the code lines for QC by mitochodrial DNAs need to be manually revised according to each data set.
save_top_2000_HVGs = read.csv(reference_top_2000_HVGs_file, stringsAsFactors = F, row.names = 1)$x

RNA[['cell_type']] = 'Unknown'

RNA <- NormalizeData(RNA, normalization.method = "LogNormalize", scale.factor = 10000)

RNA <- RNA[save_top_2000_HVGs] #If there are HVGs that are not existed in the query data set's genes, use zeros to impute.
print(RNA)


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

setwd('/home/ubuntu/single_cell/Liver_sc_atlas/results/')
ref_model_path = ref_model_path

#Perform surgery on reference model and train on query dataset without cell type labels
model = sca$models$SCANVI$load_query_data(
  adata,
  ref_model_path,
  freeze_dropout = TRUE
)

model$`_unlabeled_indices` = np$arange(adata$n_obs)
model$`_labeled_indices` = list()

print(paste0("Labelled Indices: ", length(model$`_labeled_indices`)))
print(paste0("Unlabelled Indices: ", length(model$`_unlabeled_indices`)))

model$train(
  max_epochs=as.integer(100),
  plan_kwargs=dict(weight_decay=0.0),
  check_val_every_n_epoch=as.integer(10)
)

query_latent = sc$AnnData(model$get_latent_representation())
query_latent$obs$predictions = model$predict()
write.csv(as.character(py_to_r(query_latent$obs$predictions)), save_prediction_query_set_cell_type_file)

query_latent = model$get_latent_representation()
query_latent = as.matrix(query_latent)
rownames(query_latent) = colnames(RNA)
write.csv(query_latent, save_query_set_cell_latent_representations_file)

reference_latent <- read.csv(reference_latent_representations_file, stringsAsFactors = F, row.names = 1)
ref_query_latent <- rbind(reference_latent, query_latent)

RNA$cell_type = as.character(py_to_r(model$predict()))
RNA_ref_query <- merge(RNA_ref, RNA)

rownames(ref_query_latent) = colnames(RNA_ref_query)

RNA_ref_query[['scvi']] <- CreateDimReducObject(embeddings = as.matrix(ref_query_latent), key = "scvi_", assay = DefaultAssay(RNA_ref_query))
# Find clusters, then run UMAP, and visualize
RNA_ref_query <- FindNeighbors(RNA_ref_query, dims = 1:10, reduction = 'scvi')
RNA_ref_query <- FindClusters(RNA_ref_query, resolution =1)

RNA_ref_query <- RunUMAP(RNA_ref_query, dims = 1:10, reduction = 'scvi', n.components = 2)

pdf(save_visualize_ref_query_sets_by_cell_type_file, width = 12, height = 10)
p = DimPlot(RNA_ref_query, reduction = "umap",label=T, repel = T, group.by='cell_type', raster=FALSE)
print(p)
dev.off()

pdf(save_visualize_ref_query_sets_by_patient_file, width = 12, height = 10)
p = DimPlot(RNA_ref_query, reduction = "umap",label=F, repel = T, group.by='orig.ident', raster=FALSE)
print(p)
dev.off()


pdf(save_visualize_ref_query_sets_by_study_file, width = 12, height = 10)
p = DimPlot(RNA_ref_query, reduction = "umap",label=F, repel = T, group.by='study', raster=FALSE)
print(p)
dev.off()

model$save(save_query_model_path, overwrite=TRUE)

#visualize query RNA alone:
RNA[['scvi']] <- CreateDimReducObject(embeddings = as.matrix(query_latent), key = "scvi_", assay = DefaultAssay(RNA))
# Find clusters, then run UMAP, and visualize
RNA <- FindNeighbors(RNA, dims = 1:10, reduction = 'scvi')
RNA <- FindClusters(RNA, resolution =1)

RNA <- RunUMAP(RNA, dims = 1:10, reduction = 'scvi', n.components = 2)

pdf(save_visualize_query_sets_by_cell_type_file, width = 12, height = 10)
p = DimPlot(RNA, reduction = "umap",label=F, repel = T, group.by='cell_type', raster=FALSE)
print(p)
dev.off()

pdf(save_visualize_query_sets_by_patient_file, width = 12, height = 10)
p = DimPlot(RNA, reduction = "umap",label=F, repel = T, group.by='orig.ident', raster=FALSE)
print(p)
dev.off()

pdf(save_visualize_query_sets_by_study_file, width = 12, height = 10)
p = DimPlot(RNA, reduction = "umap",label=F, repel = T, group.by='study', raster=FALSE)
print(p)
dev.off()

#---------------------------------------------------------------------------------------







