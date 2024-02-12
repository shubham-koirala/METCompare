load("/egr/research-aidd/chenruo4/infercnv/IC_input_epithelial_malignant_cells_20231221.RData")
library(Matrix)

#subset11_ann = data.frame(V2=subset11_ann$V2, stringsAsFactors = F, row.names = rownames(subset11_ann))
cell_countmatrix_query <- cell_countmatrix[, rownames(cell_metadata)[cell_metadata$V2 != 'CD8T']]
cell_countmatrix_ref <- cell_countmatrix[, rownames(cell_metadata)[cell_metadata$V2 == 'CD8T']]

cell_metadata_query <- data.frame(V2 = cell_metadata[colnames(cell_countmatrix_query), ], stringsAsFactors = F, row.names = colnames(cell_countmatrix_query))
cell_metadata_ref <- data.frame(V2 = cell_metadata[colnames(cell_countmatrix_ref), ], stringsAsFactors = F, row.names = colnames(cell_countmatrix_ref))

for(i in 1:213){

ice_obj <- infercnv::CreateInfercnvObject(#raw_counts_matrix=subset22_cm,
  raw_counts_matrix=cbind(cell_countmatrix_query[, ((i-1)*1000+1): (i*1000+1000)], cell_countmatrix_ref),
  gene_order_file=gene_ann,
  annotations_file = data.frame(V2= c(cell_metadata_query[((i-1)*1000+1): (i*1000), 'V2'], cell_metadata_ref$V2), stringsAsFactors = F, row.names = c(rownames(cell_metadata_query)[((i-1)*1000+1): (i*1000)], rownames(cell_metadata_ref))),
  ref_group_names=c("CD8T"))

ice_obj<-infercnv::run(ice_obj,
                       cutoff=0.1,
                       out_dir= paste0("output_epithelial_malignant_cells_20231221_batch", i),
                       cluster_by_group=T,
                       denoise=T,
                       HMM=T,
                       save_rds = T,
                       save_final_rds = T
)

}