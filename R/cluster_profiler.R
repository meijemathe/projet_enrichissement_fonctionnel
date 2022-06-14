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

get_ego <- function(data, organism, process, pvalue){
  # fonction enrichGO()
  # input :  vector of genes
  # output : enrichment in GO categories after FDR control.
  ego <- enrichGO(gene = data$ID,
                  OrgDb  = organism,
                  keyType = 'ENSEMBL',
                  ont = process,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = pvalue)
  # qvalueCutoff  = 0.05)
  return(ego)
}

get_gsego <- function(gene_list, organism, process, pvalue){
  gsego <- gseGO(geneList = gene_list, 
                 ont = process, 
                 keyType = "ENSEMBL", 
                 nPerm = 10000, 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = pvalue, 
                 verbose = TRUE, 
                 OrgDb = organism, 
                 pAdjustMethod = "BH")
  return(gsego)
}

db_to_organism <- function(db) {
  # if this doesn't work 
  # list of OrgDb packages is available here :
  # https://www.bioconductor.org/packages/release/BiocViews.html#___OrgDb
  # kegg friendly genome 3 letters codes are available here :
  # https://www.genome.jp/dbget-bin/www_bfind?genome
  if (grepl('Mm.', db)) {
    # mouse 
    org = 'mmu'
  } else if (grepl('Hs.', db)) {
    # human 
    org = 'hsa'
  } else if (grepl('Rn.', db)) {
    # rat 
    org = 'rno'
  } else if (grepl('Dm.', db)) {
    # fly 
    org = 'dme'
  } else if (grepl('At.', db)){
    # arabidopsis
    org = 'ath'
  } else if (grepl('Sc.', db)) {
    # yeast
    org = 'sce'
  } else if (grepl('Dr.', db)) {
    # zebrafish
    org = 'dre'
  } else if (grepl('Ce.', db)) {
    # worm
    org = 'cel'
  } else if (grepl('Bt.', db)) {
    # bovine
    org = 'bta'
  } else if (grepl('Ss.', db)) {
    # pig
    org = 'ssc'
  } else if (grepl('Gg.', db)) {
    # chicken
    org = 'gga'
  } else if (grepl('Mmu.', db)) {
    # rhesus
    org = 'mcc'
  } else if (grepl('Cf.', db)) {
    # canine
    org = 'cfa'
  } else if (grepl('EcK12.', db)) {
    # e coli strain K12
    org = 'eco'
  } else if (grepl('Xl.', db)) {
    # xenopus
    org = 'xla'
  } else if (grepl('Ag.', db)) {
    # anopheles
    org = 'aga'
  } else if (grepl('Pt.', db)) {
    # chimpanzee
    org = 'ptr'
  } else if (grepl('EcSakai.', db)) {
    # e coli strain sakai
    org = 'ecs'
  } else if (grepl('Mxanthus.', db)) {
    # myxococcus xanthus DK 1622
    org = 'mxa'
  }
  
  return(org)
}

get_ekk <- function(data, pvalue, organism){
  ekk <- enrichKEGG(gene = data$Y,
                    organism = db_to_organism(organism),
                    pvalueCutoff = pvalue,
                    pAdjustMethod = "BH")
  return(ekk)
}

get_gsekk <- function(kegg_gene_list, pvalue, organism){
  gsekk <- gseKEGG(geneList = kegg_gene_list,
                   organism = db_to_organism(organism),
                   pvalueCutoff = pvalue,
                   pAdjustMethod = "BH",
                   verbose = TRUE)
  return(gsekk)
}





