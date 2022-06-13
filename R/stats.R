# Projet d'annotation fonctionnel avec R Shiny 

# Auteurs :

# MATHE Meije : meije.mathe@univ-rouen.fr
# PETY Solene : solene.pety@etu.univ-rouen.fr
# LETERRIER Bryce : bryce.leterrier@univ-rouen.fr
# OLLIVIER Louis : louis.ollivier@etu.univ-rouen.fr

# M2.1 BIMS - Univ. Rouen Normandie 

# 2021 - 2022

# DOMAINS OVER REPRESENTATION ANALYSIS

get_domains_gene_list <- function(data, organism){
    original_gene_list <- data$log2FC
    names(original_gene_list) <- data$ID
    # GET REFSEQ IDS
    ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "REFSEQ", OrgDb=organism)
    # remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
    dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
    # Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
    data3 = data[data$ID %in% dedup_ids$ENSEMBL,]
    # Create a new column in df2 with the corresponding ENTREZ IDs
    data3$Y = dedup_ids$REFSEQ
    # Create a vector of the gene unuiverse
    domain_gene_list <- data3$log2FC
    # Name vector with ENTREZ ids
    names(domain_gene_list) <- data3$Y
    # omit any NA values 
    domain_gene_list<-na.omit(domain_gene_list)
    # sort the list in decreasing order (required for clusterProfiler)
    domain_gene_list = sort(domain_gene_list, decreasing = TRUE)
    return(domain_gene_list)
}

get_domains_interest <- function(domain_gene_list, organism){
    # TODO : gene_esembl pas en dur : ajouter une colonne au tableau d'organismes et
    ensembl = useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
    refseqids = names(domain_gene_list)
    domains = getBM(attributes = c("entrezgene_id", "refseq_mrna", "interpro", "interpro_description"), 
                  filters = "refseq_mrna",
                  values = refseqids, 
                  mart = ensembl)
    return(domains)
}

get_genes_refrence <- function(organism){
  # REFERENCE GENE LIST
  gene_ref <- keys(organism, "ENTREZID")
  return(gene_ref)
}

get_domains_reference <- function(organism){
    # DOMAINS ANNOTATION FOR REFERENCE LIST
    # TODO : gene_esembl pas en dur : ajouter une colonne au tableau d'organismes et
    ensembl = useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
    domain_ref_id <- keys(organism, "REFSEQ")
    domain_ref = getBM(attributes = c("entrezgene_id", "refseq_mrna", "interpro", "interpro_description"), 
                   filters = "refseq_mrna",
                   values = domain_ref_id, 
                   mart = ensembl)
    return(domain_ref)
}

get_domains_ORA_datatable <- function(domain_ref, gene_interest, gene_ref){
    # DATATABLE : Interpro ID, m, X, BgRatio, GeneRatio, pvalue, adjusted pvalue
    table <- data.frame(interproID = unique(domain_ref$interpro), domain = unique(domain_ref$interpro_description))
    k = length(gene_interest) #total nb of genes in the interest list
    n = length(gene_ref) #total nb of genes in the reference list
    table$m <- table(domain_ref$interpro)[table$interpro] # nb of annotated genes in the reference list
    table$X <- table(domains$interpro)[table$interpro] # nb of annotated genes in the interest list
    table <- na.omit(table) # remove genes that are in the reference list but not in the interest one
    table$BgRatio <- signif(100*table$m/(table$m+n), 3) # compute background ratio for each domain
    table$GeneRatio <- signif(100*table$X/k, 3) # compute gene ratio for each domain
    table$pvalue <- signif(phyper(table$X-1, table$m, n, k, lower.tail = FALSE), digits = 6) # compute p-value for each domain
    table$p.adjust = signif(p.adjust(table$pvalue, method = "hochberg"), digits = 6) # compute adjusted p-value for each domain
    # FINAL RESULTS
    domain_ORA <- table[c("interproID","pvalue","p.adjust","BgRatio","GeneRatio","X","domain")]
    return(domain_ORA)
    
}

domains_barplot <- function(domain_ORA, pvalue){
    res.signif <- domain_ORA[domain_ORA$p.adjust <= pvalue ,]
    res.signif <- res.signif[order(res.signif$p.adjust),]
    top10 <- res.signif[1:10,]
    #plot only first 10 signif genes on plot
    barplot <- ggplot(data = top10, aes(x=X, y=domain, fill = p.adjust)) +
        geom_bar(stat = "identity") +
        scale_fill_gradient2(high = 'red', mid='blue', space='Lab') + 
        theme_minimal()
    return(barplot)
}






