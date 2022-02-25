# Projet d'annotation fonctionnel avec R Shiny 

# Auteurs :

# MATHE Meije : meije.mathe@univ-rouen.fr
# PETY Solene : solene.pety@etu.univ-rouen.fr
# LETERRIER Bryce : bryce.leterrier@univ-rouen.fr
# OLLIVIER Louis : louis.ollivier@etu.univ-rouen.fr

# M2.1 BIMS - Univ. Rouen Normandie 

# 2021 - 2022

get_gene_list <- function(data){
  # select log2 fold change 
  original_gene_list <- data$log2FC
  # set gene names
  names(original_gene_list) <- data$ID
  # omit any NA values 
  gene_list<-na.omit(original_gene_list)
  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)
  return(gene_list)
}

get_kegg_gene_list <- function(data, organism){
  original_gene_list <- data$log2FC
  names(original_gene_list) <- data$ID
  ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)
  # remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
  dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
  # Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
  data2 = data[data$ID %in% dedup_ids$ENSEMBL,]
  # Create a new column in df2 with the corresponding ENTREZ IDs
  data2$Y = dedup_ids$ENTREZID
  # Create a vector of the gene unuiverse
  kegg_gene_list <- data2$log2FC
  # Name vector with ENTREZ ids
  names(kegg_gene_list) <- data2$Y
  # omit any NA values 
  kegg_gene_list<-na.omit(kegg_gene_list)
  # sort the list in decreasing order (required for clusterProfiler)
  kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
  l <- list(data2, kegg_gene_list)
  return(l)
}

get_ego <- function(data, organism){
  # fonction enrichGO()
  # input :  vector of genes
  # output : enrichment in GO categories after FDR control.
  ego <- enrichGO(gene = data$ID,
                  OrgDb  = organism,
                  keyType = 'ENSEMBL',
                  ont = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05)
  return(ego)
}

get_gsego <- function(gene_list, organism){
  gsego <- gseGO(geneList=gene_list, 
                  ont ="ALL", 
                  keyType = "ENSEMBL", 
                  nPerm = 10000, 
                  minGSSize = 3, 
                  maxGSSize = 800, 
                  pvalueCutoff = 0.05, 
                  verbose = TRUE, 
                  OrgDb = organism, 
                  pAdjustMethod = "none")
  return(gsego)
}

get_ekk <- function(data){
  ekk <- enrichKEGG(gene = data$Y,
                    organism = 'mmu',
                    pvalueCutoff = 0.05)
  return(ekk)
}

get_gsekk <- function(kegg_gene_list){
  gsekk <- gseKEGG(geneList = kegg_gene_list,
                 organism = 'mmu',
                 pvalueCutoff = 0.05,
                 verbose = FALSE)
  return(gsekk)
}





