library('octad')

files = Sys.glob(file.path('C:\\Users\\Dmitry Leshchiner\\Desktop\\final_analysis\\signaturescut\\', 'res_*.csv', fsep="\\")) 

for (i in files){
  data = read.csv(i)
  data$log2FoldChange <- data$log2FoldChange*-1
  nom = strsplit(i, 'res_')[[1]][2]
  
  sRGES = runsRGES(data, max_gene_size=50, permutations=10000)
  
  sRGESf = sRGES[sRGES$sRGES < -0.2,]
  sRGESf = sRGESf[sRGESf$n > 3,]
  sRGESf = sRGESf[!grepl("BRD-", sRGESf$pert_iname),]
  
  if (length(sRGESf$pert_iname) > 0){
    srges_nom = sprintf('C:\\Users\\Dmitry Leshchiner\\Desktop\\final_analysis\\sRGES_cut\\sRGES_%s', nom)
    srgesf_nom = sprintf('C:\\Users\\Dmitry Leshchiner\\Desktop\\final_analysis\\sRGES_cut\\sRGESf_%s', nom)
    
    write.csv(sRGES, srges_nom)
    write.csv(sRGESf, srgesf_nom)
    
    #skip_to_next <- FALSE
    #output = 'C:\\Users\\Dmitry Leshchiner\\Desktop\\Sigs_24dec\\enrich\\'
    #tryCatch(octadDrugEnrichment(sRGES=sRGESf, target = c('chembl_targets','mesh','ChemCluster'), 
    #                             outputFolder = output, enrichFolder = strsplit(nom, '.csv')[[1]]), error = function(e) { skip_to_next <<- TRUE})
    #if(skip_to_next) { next }     
  #} else {print(nom)
  }
}

# lincs_signatures <- get_ExperimentHub_data("EH7271")
# lincs_sig_info <- get_ExperimentHub_data("EH7270")  