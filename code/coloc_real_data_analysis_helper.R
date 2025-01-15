get_top_mts_coloc <- function(colocf,group){
  coloc_res <- readRDS(colocf)
  #coloc_res <- coloc_res[coloc_res$PP4>0.8,]
  coloc_res$id <- paste0(coloc_res$id,"|",group)
  rownames(coloc_res) <- NULL
  return(coloc_res)
}

get_top_genes_coloc <- function(colocf){
  coloc_res <- readRDS(colocf)
  rownames(coloc_res) <- NULL
  return(coloc_res)
}

get_top_mts_ctwas <- function(ctwasf,regions){
  finemap_res <- readRDS(ctwasf)
  finemap_res <- finemap_res$finemap_res
  finemap_res <- finemap_res[finemap_res$region_id %in% regions,]
  finemap_res <- finemap_res[(finemap_res$susie_pip>0.8) & (finemap_res$type!="SNP") & !is.na((finemap_res$cs)),]
  return(finemap_res)
}

get_top_eQTL_genes_ctwas <- function(gwas_id,tissue,regions,mapping_table){
  ctwasf <- paste0("/project/xinhe/shengqian/single_tissue_screen/processed_weights/expression_weights/",gwas_id,"/",gwas_id,"_",tissue,".finemap_regions_res.RDS")
  susie_alpha_res <- readRDS(ctwasf)
  susie_alpha_res <- susie_alpha_res$susie_alpha_res
  susie_alpha_res <- susie_alpha_res[susie_alpha_res$region_id %in% regions, ]
  susie_alpha_res <- anno_susie_alpha_res(susie_alpha_res,
                                          mapping_table = mapping_predictdb,
                                          map_by = "molecular_id",
                                          drop_unmapped = TRUE)
  
  combined_res <- combine_gene_pips(susie_alpha_res,
                                    group_by = "gene_name",
                                    by = "type",
                                    method = "combine_cs",
                                    filter_cs = TRUE)
  
  combined_res <- combined_res[combined_res$combined_pip>0.8, c("gene_name","combined_pip")]
  
  return(combined_res)
}

get_top_sQTL_genes_ctwas <- function(gwas_id,tissue,regions,mapping_table){
  ctwasf <- paste0("/project/xinhe/shengqian/single_tissue_screen/processed_weights/splicing_weights/",gwas_id,"/",gwas_id,"_",tissue,".finemap_regions_res.RDS")
  susie_alpha_res <- readRDS(ctwasf)
  susie_alpha_res <- susie_alpha_res$susie_alpha_res
  susie_alpha_res <- susie_alpha_res[susie_alpha_res$region_id %in% regions, ]
  susie_alpha_res <- anno_susie_alpha_res(susie_alpha_res,
                                          mapping_table = mapping_predictdb,
                                          map_by = "molecular_id",
                                          drop_unmapped = TRUE)
  
  combined_res <- combine_gene_pips(susie_alpha_res,
                                    group_by = "gene_name",
                                    by = "type",
                                    method = "combine_cs",
                                    filter_cs = TRUE)
  
  combined_res <- combined_res[combined_res$combined_pip>0.8, c("gene_name","combined_pip")]
  
  return(combined_res)
}

