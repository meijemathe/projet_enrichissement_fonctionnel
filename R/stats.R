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
  # > head(datasets)
  if(organism == 'org.Mm.eg.db'){
    dataset = "mmusculus_gene_ensembl"
  }
  # btaurus_gene_ensembl
  # celegans_gene_ensembl
  # clfamiliaris_gene_ensembl
  # dmelanogaster_gene_ensembl
  # drerio_gene_ensembl
  # ggallus_gene_ensembl 
  # hsapiens_gene_ensembl 
  # mmulatta_gene_ensembl
  # ptroglodytes_gene_ensembl
  # scerevisiae_gene_ensembl
  # sscrofa_gene_ensembl
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
  gene_ref <- keys(org.Mm.eg.db, "ENTREZID")
  
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




