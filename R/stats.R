# Projet d'annotation fonctionnel avec R Shiny 

# Auteurs :

# MATHE Meije : meije.mathe@univ-rouen.fr
# PETY Solene : solene.pety@etu.univ-rouen.fr
# LETERRIER Bryce : bryce.leterrier@univ-rouen.fr
# OLLIVIER Louis : louis.ollivier@etu.univ-rouen.fr

# M2.1 BIMS - Univ. Rouen Normandie 

# 2021 - 2022

# DOMAINS OVER REPRESENTATION ANALYSIS

org_to_ensembldb <- function(organism){
  # to get all the available datasets run : 
  # > library(biomaRt)
  # > ensembl <- useEnsembl(biomart = "genes")
  # > datasets <- listDatasets(ensembl)
  # > head(datasets)
  # plants :
  # > listEnsemblGenomes()
  # > ensembl_plants <- useEnsemblGenomes(biomart = "plants_mart")
  # > datasets_plants <- listDatasets(ensembl_plants)
  # > head(datasets_plants)
  
  if(organism == 'org.Mm.eg.db'){
    dataset = "mmusculus_gene_ensembl"
  } else if(organism == 'org.Hs.eg.db'){
    dataset = "hsapiens_gene_ensembl"
  } else if(organism == 'org.Rn.eg.db'){
    dataset = "rnorvegicus_gene_ensembl"
  } else if(organism == 'org.Dm.eg.db'){
    dataset = "dmelanogaster_gene_ensembl"
  } else if(organism == 'org.At.tair.db'){
    dataset = "athaliana_eg_gene"
  } else if(organism == 'org.Sc.sgd.db'){
    dataset = "scerevisiae_gene_ensembl"
  } else if(organism == 'org.Dr.eg.db'){
    dataset = "drerio_gene_ensembl"
  } else if(organism == 'org.Ce.eg.db'){
    dataset = "celegans_gene_ensembl"
  } else if(organism == 'org.Bt.eg.db'){
    dataset = "btaurus_gene_ensembl"
  } else if(organism == 'org.Ss.eg.db'){
    dataset = "sscrofa_gene_ensembl"
  } else if(organism == 'org.Gg.eg.db'){
    dataset = "ggallus_gene_ensembl"
  } else if(organism == 'org.Mmu.eg.db'){
    dataset = "mmulatta_gene_ensembl"
  } else if(organism == 'org.Cf.eg.db'){
    dataset = "clfamiliaris_gene_ensembl"
  } else if(organism == 'org.EcK12.eg.db'){
    dataset = "escherichia_coli_str_k_12_substr_mg1655_gca_000005845"
  } else if(organism == 'org.Xl.eg.db'){
    dataset = "xtropicalis_gene_ensembl"
  } else if(organism == 'org.Ag.eg.db'){
    dataset = "anopheles_gambiae"
  } else if(organism == 'org.Pt.eg.db'){
    dataset = "ptroglodytes_gene_ensembl"
  } else if(organism == 'org.EcSakai.eg.db'){
    dataset = "escherichia_coli_o157_h7_str_sakai_gca_000008865"
  } else if(organism == 'org.Mxanthus.db'){
    dataset = "myxococcus_xanthus_dk_1622_gca_000012685 "
  }
  return(dataset)
}

get_table_ORA_domains <- function(data, organism, pvalue_lim){
  ## =======================================================================
  ## DATA PREPARATION 
  ## =======================================================================
  # select log2 fold change 
  original_gene_list <- data$log2FC
  # set gene names
  names(original_gene_list) <- data$ID
  # omit any NA values 
  gene_list<-na.omit(original_gene_list)
  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)
  ## =======================================================================
  ## GET REFSEQ IDS
  ## =======================================================================
  ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "REFSEQ", OrgDb=organism)
  # remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
  dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
  
  # Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
  data3 = data[data$ID %in% dedup_ids$ENSEMBL,]
  
  # Create a new column in df2 with the corresponding ENTREZ IDs
  data3$Y = dedup_ids$REFSEQ
  
  # Create a vector of the gene universe
  domain_gene_list <- data3$log2FC
  
  # Name vector with ENTREZ ids
  names(domain_gene_list) <- data3$Y
  
  # omit any NA values 
  domain_gene_list<-na.omit(domain_gene_list)
  
  # sort the list in decreasing order (required for clusterProfiler)
  domain_gene_list = sort(domain_gene_list, decreasing = TRUE)
  # DOMAINS ANNOTATION FOR INTEREST LIST
  dataset = org_to_ensembldb(organism)
  ensembl = useMart(biomart = "ensembl", dataset = dataset)
  
  refseqids = names(domain_gene_list)
  domains = getBM(attributes = c("entrezgene_id", "refseq_mrna", "interpro", "interpro_description"), 
                  filters = "refseq_mrna",
                  values = refseqids, 
                  mart = ensembl)
  # REFERENCE GENE LIST
  gene_ref <- keys(get(organism), "ENTREZID")
  
  # DOMAINS ANNOTATION FOR REFERENCE LIST
  domain_ref_id <- keys(org.Mm.eg.db, "REFSEQ")
  domain_ref = getBM(attributes = c("entrezgene_id", "refseq_mrna", "interpro", "interpro_description"), 
                     filters = "refseq_mrna",
                     values = domain_ref_id, 
                     mart = ensembl)
  gene_ref <- keys(get(organism), "ENTREZID")
  # DATATABLE : Interpro ID, m, X, BgRatio, GeneRatio, pvalue, adjusted pvalue
  table <- data.frame(interproID = unique(domain_ref$interpro), domain = unique(domain_ref$interpro_description))
  k = length(gene_list) #total nb of genes in the interest list
  n = length(gene_ref) #total nb of genes in the reference list
  table$m <- table(domain_ref$interpro)[table$interpro] # nb of annotated genes in the reference list
  table$X <- table(domains$interpro)[table$interpro] # nb of annotated genes in the interest list
  table <- na.omit(table) # remove genes that are in the reference list but not in the interest one
  table$BgRatio <- signif(100*table$m/(table$m+n), 3) # compute background ratio for each domain
  table$GeneRatio <- signif(100*table$X/k, 3) # compute gene ratio for each domain
  table$pvalue <- signif(phyper(table$X-1, table$m, n, k, lower.tail = FALSE), digits = 6) # compute p-value for each domain
  table$padjust = signif(p.adjust(table$pvalue, method = "hochberg"), digits = 6) # compute adjusted p-value for each domain
  # FINAL RESULTS
  res <- table[c("interproID","pvalue","padjust","BgRatio","GeneRatio","X","domain")]
  res.signif <- res[res$padjust <= pvalue_lim ,]
  res.signif <- res.signif[order(res.signif$padjust),]
  
  return(res.signif)
}

domains_ORA_barplot <- function(res.signif){
  top10 <- res.signif[1:10,]
  top10$X <- as.numeric(top10$X)
  p<-ggplot(data=top10, aes(x = X, y = domain, fill = padjust)) +
    xlab("Count") +
    geom_bar(stat="identity") +
    scale_fill_gradient2(high='red', mid='blue', space='Lab') + 
    theme_minimal()
  return(p)
}


domains_ORA_dotplot <- function(res.signif){
  top10 <- res.signif[1:10,]
  top10$X <- as.numeric(top10$X)
  size1 <- c()
  for(i in 1:length(top10$X)){
    if(top10$X[i] > 400){
      size1 <- c(size1, 2)
    }else if(top10$X[i] > 200){
      size1 <- c(size1, 1)
    }else{
      size1 <- c(size1, 0.7)
    }
  }
  
  p<-ggplot(top10, aes(x=GeneRatio, y=domain)) + 
    geom_dotplot(binaxis='y', stackdir='center') +
    geom_point(aes(size=size1, color = padjust)) +
    scale_colour_gradient2(mid = "blue", high = "red") + 
    
    theme_minimal()
  return(p)
}

get_table_GSEA_domains <- function(data, organism, pvalue_lim){
  ## =======================================================================
  ## DATA PREPARATION 
  ## =======================================================================
  # select log2 fold change 
  original_gene_list <- data$log2FC
  # set gene names
  names(original_gene_list) <- data$ID
  # omit any NA values 
  gene_list<-na.omit(original_gene_list)
  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)
  ## =======================================================================
  ## GET REFSEQ IDS
  ## =======================================================================
  ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "REFSEQ", OrgDb=organism)
  # remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
  dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
  
  # Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
  data3 = data[data$ID %in% dedup_ids$ENSEMBL,]
  
  # Create a new column in df2 with the corresponding ENTREZ IDs
  data3$Y = dedup_ids$REFSEQ
  
  # Create a vector of the gene universe
  domain_gene_list <- data3$log2FC
  
  # Name vector with ENTREZ ids
  names(domain_gene_list) <- data3$Y
  
  # omit any NA values 
  domain_gene_list<-na.omit(domain_gene_list)
  
  # sort the list in decreasing order (required for clusterProfiler)
  domain_gene_list = sort(domain_gene_list, decreasing = TRUE)
  # DOMAINS ANNOTATION FOR INTEREST LIST
  dataset = org_to_ensembldb(organism)
  ensembl = useMart(biomart = "ensembl", dataset = dataset)
  
  refseqids = names(domain_gene_list)
  domains = getBM(attributes = c("entrezgene_id", "refseq_mrna", "interpro", "interpro_description"), 
                  filters = "refseq_mrna",
                  values = refseqids, 
                  mart = ensembl)
  
  # DOMAINS ANNOTATION FOR REFERENCE LIST
  domain_ref_id <- keys(get(organism), "REFSEQ")
  domain_ref = getBM(attributes = c("entrezgene_id", "refseq_mrna", "interpro", "interpro_description"), 
                     filters = "refseq_mrna",
                     values = domain_ref_id, 
                     mart = ensembl)
  
  data_gsea <- left_join(domain_ref, data3[,c("log2FC", "Y")], by = c("refseq_mrna" = "Y"))
  
  domain2gene <- data_gsea[, c("interpro","entrezgene_id")]
  domain2gene <- domain2gene[!is.na(domain2gene$entrezgene_id),]
  domain2gene$entrezgene_id <- as.character(domain2gene$entrezgene_id)
  
  domain2name <- data_gsea[, c("interpro","interpro_description")]
  domain2name$entrezgene_id <- as.character(domain2name$interpro)
  geneList2 = data_gsea$log2FC[!is.na(data_gsea$log2FC)]
  names(geneList2) = data_gsea$entrezgene_id[!is.na(data_gsea$log2FC)]
 
  geneList2 <- geneList2[order(unlist(geneList2), decreasing = TRUE)]
  
  y = GSEA(geneList2, 
           pvalueCutoff = 0.1,
           TERM2GENE = domain2gene, 
           TERM2NAME = domain2name)
  table <- as.data.frame(summary(y))
  table <- table[, -c(10,11,12)]

  res.signif <- table[table$p.adjust <= pvalue_lim ,]
  res.signif <- res.signif[order(res.signif$p.adjust),]
  return(list(res.signif, y))
}








