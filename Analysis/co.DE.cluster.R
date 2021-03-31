library(org.Hs.eg.db)
require(dplyr)
require(bedr)
require(foreach)

gene_with_protein_product <- read.delim("~/Project/BreastCancerMetaPotenial/client-side/Data/HGNC/gene_with_protein_product.txt", stringsAsFactors=FALSE)
mapping.df                <- gene_with_protein_product[,c('ensembl_gene_id','symbol')]
mapping.df                <- mapping.df[mapping.df$ensembl_gene_id != '',]
rownames(mapping.df)      <- mapping.df$ensembl_gene_id

load('client-side/output/DE.breast.cancer.R.output/DE.breast.cancer.RData')
load('client-side/output/DE.colorectal.cancer.R.output/DE.colorectal.cancer.RData')
load('client-side/output/DE.NET.pancreatic.cancer.R.output/DE.NET.pancreatic.cancer.RData')
load('client-side/output/DE.NET.si.cancer.R.output/DE.NET.si.cancer.RData')
load('client-side/output/DE.prostate.cancer.R.output/DE.prostate.cancer.RData')
DE.rs.list        <- list(BRCA.Basal.DE.rs, BRCA.Her2.DE.rs, BRCA.LumB.DE.rs, PRAD.DE.rs,COAD.DE.rs,NET.PAAD.DE.rs,NET.SI.DE.rs)
names(DE.rs.list) <- c('BRCA.Basal', 'BRCA.Her2', 'BRCA.LumB','PRAD', 'COAD', 'PNET', 'SINET')

load('client-side/output/hg19.gene.info.R.output/hg19.gene.info.RData')
distance <- 60000
get.co.up.DE.cluster <- function(DE.rs) {
    up.gene.symbol    <- mapping.df[DE.rs$tumor.intrinsic.DE.gene.rs$up.gene,'symbol']
    up.gene.symbol    <- up.gene.symbol[is.na(up.gene.symbol) == FALSE]
    idx               <- match(up.gene.symbol,hg19.gene.info$genename)
    up.gene.coor.df   <- hg19.gene.info[idx,]     
    up.gene.coor.df   <- up.gene.coor.df[complete.cases(up.gene.coor.df),]

    bed.region        <- paste("chr",up.gene.coor.df$chrom,':',up.gene.coor.df$start,"-",up.gene.coor.df$end,sep='')
    names(bed.region) <- up.gene.coor.df$genename

    bed.sorted.region <- bedr.sort.region(bed.region)
    merged.region     <- bedr.merge.region(x=bed.sorted.region,distance=distance,check.chr = TRUE,check.sort = TRUE)

    
    co.up.clusters    <- setdiff(merged.region,bed.sorted.region)
    if(length(co.up.clusters) == 0){
        return(NULL)
    }
    region.size       <- bedr:::size.region(co.up.clusters)
    df                <- data.frame(co.up.clusters=co.up.clusters,region.size=region.size)
    df$co.up.clusters <- as.character(df$co.up.clusters)

    df$gene.name <- foreach(x= df$co.up.clusters %>% as.character(),.combine='c') %do% {
        f <- bed.sorted.region %in.region% x
        s <- bed.sorted.region[f]
        paste(names(bed.region)[match(s,bed.region)],collapse =":")
    }

    df$gene.number <- foreach(x= df$co.up.clusters %>% as.character(),.combine='c') %do% {
        f <- bed.sorted.region %in.region% x
        sum(f)
    }

    df <- df[order(df$gene.number,decreasing = TRUE),]
    df
}

get.co.dn.DE.cluster <- function(DE.rs) {
  dn.gene.symbol    <- mapping.df[DE.rs$tumor.intrinsic.DE.gene.rs$dn.gene,'symbol']
  dn.gene.symbol    <- dn.gene.symbol[is.na(dn.gene.symbol) == FALSE]
  idx               <- match(dn.gene.symbol,hg19.gene.info$genename)
  dn.gene.coor.df   <- hg19.gene.info[idx,]     
  dn.gene.coor.df   <- dn.gene.coor.df[complete.cases(dn.gene.coor.df),]
  
  bed.region        <- paste("chr",dn.gene.coor.df$chrom,':',dn.gene.coor.df$start,"-",dn.gene.coor.df$end,sep='')
  names(bed.region) <- dn.gene.coor.df$genename
  
  bed.sorted.region <- bedr.sort.region(bed.region)
  merged.region     <- bedr.merge.region(x=bed.sorted.region,distance=distance,check.chr = TRUE,check.sort = TRUE)
  
  co.dn.clusters    <- setdiff(merged.region,bed.sorted.region)
  if(length(co.dn.clusters) == 0){
      return(NULL)
  }
  region.size       <- bedr:::size.region(co.dn.clusters)
  df                <- data.frame(co.dn.clusters=co.dn.clusters,region.size=region.size)
  df$co.dn.clusters <- as.character(df$co.dn.clusters)
  
  df$gene.name <- foreach(x= df$co.dn.clusters %>% as.character(),.combine='c') %do% {
    f <- bed.sorted.region %in.region% x
    s <- bed.sorted.region[f]
    paste(names(bed.region)[match(s,bed.region)],collapse =":")
  }
  
  df$gene.number <- foreach(x= df$co.dn.clusters %>% as.character(),.combine='c') %do% {
    f <- bed.sorted.region %in.region% x
    sum(f)
  }
  
  df <- df[order(df$gene.number,decreasing = TRUE),]
  df
}

co.up.DE.cluster.list <- lapply(DE.rs.list,get.co.up.DE.cluster)
co.dn.DE.cluster.list <- lapply(DE.rs.list,get.co.dn.DE.cluster)
names(co.up.DE.cluster.list) <- names(DE.rs.list)
names(co.dn.DE.cluster.list) <- names(DE.rs.list)


pooled.up.DE.cluster.df <- foreach(i = 1: length(co.up.DE.cluster.list),.combine='rbind') %do% {
    df             <-   co.up.DE.cluster.list[[i]]
    if( is.null(nrow(df)) == FALSE){
        df$cancer.type <- names(co.up.DE.cluster.list)[i]
        df
    }else{
        NULL  
    }
}

pooled.dn.DE.cluster.df <- foreach(i = 1: length(co.dn.DE.cluster.list),.combine='rbind') %do% {
  df             <-   co.dn.DE.cluster.list[[i]]
  if(is.null(nrow(df)) == FALSE){
    df$cancer.type <- names(co.dn.DE.cluster.list)[i]
    df
  }else{
    NULL  
  }
}

BRCA.LumB.fake.cluster.rs <- foreach(i = 1:1000, .combine='rbind') %do% {
    fake.DE.rs <- list()
    fake.DE.rs$tumor.intrinsic.DE.gene.rs$up.gene <- sample(mapping.df$ensembl_gene_id,length(BRCA.LumB.DE.rs$tumor.intrinsic.DE.gene.rs$up.gene))
    get.co.up.DE.cluster(fake.DE.rs)
}

BRCA.Basal.fake.cluster.rs <- foreach(i = 1:1000, .combine='rbind') %do% {
    fake.DE.rs <- list()
    fake.DE.rs$tumor.intrinsic.DE.gene.rs$up.gene <- sample(mapping.df$ensembl_gene_id,length(BRCA.Basal.DE.rs$tumor.intrinsic.DE.gene.rs$up.gene))
    get.co.up.DE.cluster(fake.DE.rs)
}


save(file = 'client-side/output/co.DE.cluster.R.output/co.DE.cluster.RData',list=c('pooled.up.DE.cluster.df','pooled.dn.DE.cluster.df','BRCA.LumB.fake.cluster.rs','BRCA.Basal.fake.cluster.rs'))


start <- 13001942 #start of gene GCDH
end   <- 15560762 # end of  gene WIZ



flag              <- hg19.gene.info$chrom == '19' & hg19.gene.info$start >= start & hg19.gene.info$end <= end
chr19.p13.12.gene <- hg19.gene.info[flag,'genename'] %>% as.character()
flag              <- grepl(chr19.p13.12.gene,pattern = 'MIR') | grepl(chr19.p13.12.gene,pattern = 'LINC') | grepl(chr19.p13.12.gene,pattern = 'LOC') |grepl(chr19.p13.12.gene,pattern = 'SNORA')
chr19.p13.12.gene <- chr19.p13.12.gene[!flag]

basal.up.gene <- mapping.df[BRCA.Basal.DE.rs$tumor.intrinsic.DE.gene.rs$up.gene,'symbol'] %>% unique
basal.dn.gene <- mapping.df[BRCA.Basal.DE.rs$tumor.intrinsic.DE.gene.rs$dn.gene,'symbol'] %>% unique

c.gene           <- intersect(basal.up.gene,chr19.p13.12.gene)
c.gene.df       <- hg19.gene.info[hg19.gene.info$genename %in% c.gene,]
c.gene.df       <- c.gene.df[order(c.gene.df$start),]
c.gene.df$chrom <- paste('chr',c.gene.df$chrom,sep='')
write.table(x=c.gene.df[,c('chrom','start','end','genename')],file='client-side/output/co.DE.cluster.R.output//basal.up.chr19.gene.bed',quote=FALSE,row.names=FALSE)



get.id <- function(x){
    l <- strsplit(x = x,split='\\.')  %>% unlist 
    l[1]
}
MET500.seg.data                         <- read.csv("~/Project/Cancer2CellLine/client-side/Data/CNV/cnv_v4.csv", stringsAsFactors=FALSE)
MET500.seg.data                         <- MET500.seg.data[,c('Pipeline_ID','Chr','Start','End','Log2_Coverage_Ratio')]
colnames(MET500.seg.data)               <- c('ID','chrom','loc.start','loc.end','seg.mean')
MET500.seg.data$ID                      <- sapply(MET500.seg.data$ID,get.id)

load('~/Project/Cancer2CellLine/server-side/RData/MET500.RData')
load('client-side/output/Select.pure.sample.breast.cancer.R.output/Select.pure.sample.breast.cancer.RData')

pure.MET.breast.cancer.Basal.patient <- MET500.sample.meta$MET500.id[match(pure.MET.breast.cancer.Basal.sample,MET500.sample.meta$Run)]
flag <-(MET500.seg.data$ID %in% pure.MET.breast.cancer.Basal.patient) & (MET500.seg.data$chrom == 19) &
       (MET500.seg.data$loc.end >= start) & (MET500.seg.data$loc.start <= end)
bed.data <- MET500.seg.data[flag,]  
bed.data$genename <- bed.data$ID
bed.data <- bed.data[order(bed.data$genename),]
write.table(x=bed.data[,c('chrom','loc.start','loc.end','genename')],file='client-side/output/co.DE.cluster.R.output//MET500.basal.cnv.call.bed',quote=FALSE,row.names=FALSE)


