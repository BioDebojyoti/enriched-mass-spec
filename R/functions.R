packages <- c(
  "this.path",
  "comprehenr",
  "dplyr",
  "readr",
  "stringr",
  "tidyverse",
  "tidyr",
  "EnhancedVolcano",
  "ggplot2",
  "plotly",
  "org.Hs.eg.db",
  "clusterProfiler",
  "enrichR",
  "enrichplot",
  "shiny",
  "shinyBS",
  "bslib",
  "DT"
)

for(pkg in packages){
  library(pkg, character.only = TRUE)
}

ont_list <- c("MF", "CC", "BP")
names(ont_list) <- c("Molecular Function", "Cellular Component",  "Biological Process")

rename_column <- function(s) {
  if (grepl("^PG.",s)){
    s <- gsub("PG.","",s)
    return(s)
  } else {
    p1 <- stringr::str_split(s, " ng ")[[1]][2]
    p2 <- stringr::str_split(p1, "DIA ")[[1]][2]
    s_new <- stringr::str_split(p2, "_")[[1]][1]
    return(s_new)
  }
}

clean_string <- function(x){
  first_char <- substr(x, 1, 1)
  rest_str <- substr(x, 2, nchar(x))
  if(nchar(x) == 1 || grepl("[a-zA-Z0-9]", first_char)){
    return(x)
  } else {
    return(clean_string(rest_str))
  }
}


first_entry <- function(s){
  s = clean_string(s)
  return(gsub("_HUMAN","",str_split(s,";")[[1]][1]))
}

make_numeric_if_possible <- function(df) {
  df[] <- lapply(df, function(col) {
    col_cleaned <- gsub(",", ".", col)
    if (all(suppressWarnings(!is.na(as.numeric(col_cleaned))))) {
      return(as.numeric(col_cleaned))
    } else {
      return(col)
    }
  })
  return(df)
}

read_mass_spec_tsv <- function(filename, 
                               log2FC_cutoff = 1.5, 
                               pval_cutoff = 0.05,
                               logfc_col2use = "log2FC",
                               pval_col2use = "Pvalue",
                               protein_name_col = "ProteinNames",
                               gene_name_col = "Genes",
                               uniprotid_col = "UniProtId",
                               protein_description_col = "ProteinDescriptions"
                               ){
  
  raw_data <- read_tsv(file = filename,
                       col_names = TRUE,
                       show_col_types = FALSE) 
  
  raw_data <- make_numeric_if_possible(raw_data) %>%
    dplyr::select(
      all_of(c(logfc_col2use, pval_col2use, protein_name_col, gene_name_col,uniprotid_col, protein_description_col))
    )
  
  column_dictionary <- data.frame(
    label = c("log2FC", "Pvalue", "ProteinName", "GeneName", "UniProtId", "ProteinDescriptions"),
    name = c(logfc_col2use, pval_col2use, protein_name_col, gene_name_col,uniprotid_col, protein_description_col)
  )
  
  # print(names(raw_data))
  
  for(icd in seq(1,nrow(column_dictionary))){
  column2use <- to_vec(for(i in seq_along(names(raw_data))) if(names(raw_data)[i] == column_dictionary[["name"]][icd]) i)
  names(raw_data)[column2use] <- column_dictionary[["label"]][icd]
  }
   
  # print(names(raw_data))
  
  raw_data <- raw_data %>%
    as.data.frame() %>%
    dplyr::mutate(
      ProteinName = sapply(ProteinName, first_entry),
      GeneName = sapply(GeneName, first_entry),
      UniProtId = sapply(UniProtId, first_entry)
    ) %>%
    dplyr::mutate(Threshold = ifelse(
      abs(log2FC) > log2FC_cutoff & Pvalue < pval_cutoff,
      "Significant", "Not significant"))
  
  for(c in c("ProteinName", "GeneName", "UniProtId")){
    attr(raw_data[[c]], "names") <- NULL
  }
  
  return(raw_data)
}

script_path <- this.path::this.path()
reviewed_uniprotids <- read.delim(
  paste0(dirname(script_path),"/uniprotkb_AND_reviewed_true_AND_model_o_2024_12_04.tsv"))
reviewed_uniprotids <- (reviewed_uniprotids %>% dplyr::select(Entry))$Entry

duplicated_gene_rows <- function(data2check, 
                                 x = "log2FC"){
  
  duplicated_genes <- (data2check %>% 
                         dplyr::count(GeneName) %>% 
                         dplyr::filter(n!=1))$GeneName
  return(data2check %>% 
           dplyr::filter(GeneName %in% duplicated_genes) %>%
           dplyr::mutate(
             Not_reviewed =  ifelse(!(GeneName %in% duplicated_genes), 
                                    FALSE, 
                                    ifelse((!(UniProtId %in% reviewed_uniprotids)),
                                           TRUE, FALSE))
           ) %>%
           dplyr::arrange(GeneName)
  )
}

filter_duplicate_genes <- function(data2clean, 
                                   x = "log2FC",
                                   y = "Pvalue",
                                   log2FC_cutoff=1.5, 
                                   pval_cutoff=0.05){
  
  duplicated_genes <- (data2clean %>% 
                         dplyr::count(GeneName) %>% 
                         dplyr::filter(n!=1))$GeneName
  
  clean_data <- data2clean %>%
    dplyr::mutate(
      Threshold = ifelse(
        abs(log2FC) > log2FC_cutoff & Pvalue < pval_cutoff, "Significant", "Not significant")
    ) %>%
    dplyr::filter(GeneName != "NaN") %>%
    dplyr::group_by(GeneName) %>%
    dplyr::arrange(GeneName, desc(Threshold)) %>%
    dplyr::filter(!(grepl("fragment", ProteinDescriptions, ignore.case = TRUE) & (GeneName %in% duplicated_genes)  & (n() > 2) & (Threshold =="Not significant"))) %>%
    dplyr::filter(!(!(UniProtId %in% reviewed_uniprotids)  & (GeneName %in% duplicated_genes) & (n() > 2)  & (Threshold =="Not significant"))) %>%      
    dplyr::filter(row_number() == 1) %>%
    ungroup()
  
  return (clean_data)
}

volcano_plot <- function(data2plot, 
                         x = "log2FC", 
                         y = "Pvalue",
                         log2FC_cutoff=1.5, 
                         pval_cutoff=0.05, 
                         label_genes = NULL,
                         title = ""){
  
  column2use <- to_vec(
    for(i in seq_along(names(data2plot))) 
      if(names(data2plot)[i] == x) i)
  
  names(data2plot)[column2use] <-"log2FC"
  
  data2plot <- filter_duplicate_genes(data2clean = data2plot,
                                      x = "log2FC",
                                      log2FC_cutoff = log2FC_cutoff) %>%
    as.data.frame()
  
  metrics <- data2plot %>% dplyr::select(log2FC) %>% summary()
  # print(metrics)
  x1 = floor(as.numeric(gsub(" ", "", str_split(metrics[1], pattern = ":")[[1]][2])))
  x2 = ceiling(as.numeric(gsub(" ", "", str_split(metrics[6], pattern = ":")[[1]][2])))
  greater <- max(abs(x1), x2)
  xlim = c(-greater, greater)
  
  ylim = c(0,ceiling(-log10(min(data2plot[y], na.rm = TRUE))))
  
  row.names(data2plot) <- data2plot$GeneName
  
  enhanced_volcano_plot <- EnhancedVolcano(data2plot,
                                           lab = data2plot$GeneName,
                                           x = "log2FC",
                                           y = y,
                                           pCutoff = pval_cutoff,
                                           FCcutoff = log2FC_cutoff,
                                           selectLab = label_genes,
                                           pointSize = 3.0,
                                           labSize = 6.0,
                                           boxedLabels = TRUE,
                                           labCol = 'black',
                                           legendLabSize = 16,
                                           legendIconSize = 5.0,
                                           drawConnectors = TRUE,
                                           widthConnectors = 1.0,
                                           title = title,
                                           subtitle = "",
                                           xlim = xlim
  )
  
  enhanced_volcano_plot <- enhanced_volcano_plot +
    xlab("log2FC") + 
    ylab(paste0("-log10[",y,"]")) +
    scale_y_continuous(limits = ylim, expand = c(0, 0))
  
  return(enhanced_volcano_plot)
}

# GO over-representation analysis
enrichGO_pathway <- function(data2use, 
                             ont = "MF",
                             x = "log2FC", 
                             y_pval = "Pvalue",
                             pval_cutoff = 0.05,
                             ego_pval_threshold = 0.05,
                             log2FC_cutoff = 1.5,
                             gene_set = "all"){
  
  
  column2use <- to_vec(
    for(i in seq_along(names(data2use))) 
      if(names(data2use)[i] == x) i)
  
  
  names(data2use)[column2use] <-"log2FC"
  
  data2use <- filter_duplicate_genes(
    data2clean = data2use,
    x = "log2FC",
    y = y_pval,
    log2FC_cutoff = log2FC_cutoff,
    pval_cutoff = pval_cutoff
  )
  
  geneList <- data2use$log2FC
  names(geneList) <- data2use$GeneName
  
  hs <- org.Hs.eg.db
  my.symbols <- names(geneList)
  entrezid <- AnnotationDbi::select(hs, 
                                    keys = my.symbols,
                                    columns = c("ENTREZID", "SYMBOL"),
                                    keytype = "SYMBOL", 
                                    multiVals = "first")
  
  entrezid <- entrezid[!duplicated(entrezid$SYMBOL),]
  
  names(geneList) <- entrezid$ENTREZID
  
  geneList <- sort(geneList, decreasing = TRUE)
  
  if (gene_set == "all") {
    gene <- data2use %>%
      dplyr::filter((abs(log2FC) > log2FC_cutoff) &
                      ( Pvalue < pval_cutoff)) %>%
      dplyr::pull(GeneName)
  } else if (gene_set == "up") {
    gene <- data2use %>%
      dplyr::filter((log2FC > log2FC_cutoff) &
                      (Pvalue < pval_cutoff)) %>%
      dplyr::pull(GeneName)
  } else if (gene_set == "down") {
    gene <- data2use %>%
      dplyr::filter((-log2FC > log2FC_cutoff) &
                      (Pvalue < pval_cutoff)) %>%
      dplyr::pull(GeneName)
  } else {
    stop(paste0("gene set \"", gene_set, "\" is not allowed!"))
  }
  
  gene <- comprehenr::to_vec(for (g in gene)
    entrezid %>% dplyr::filter(SYMBOL == g) %>% dplyr::pull(ENTREZID))
  
  print(paste0("Number of genes to test: ", as.character(length(gene))))
  print(paste0("Number of background genes: ", as.character(length(geneList))))
  
  ego <- clusterProfiler::enrichGO(gene = gene,
                                   universe = names(geneList),
                                   OrgDb = 'org.Hs.eg.db',
                                   ont = ont, 
                                   pvalueCutoff = ego_pval_threshold, 
                                   readable = TRUE)
  
  names(geneList) <- entrezid$SYMBOL[match(names(geneList), entrezid$ENTREZID)]
  
  
  if(nrow(ego@result %>% dplyr::filter(p.adjust < ego_pval_threshold)) != 0){
    ego_heatmap <- heatplot(ego, foldChange=geneList, showCategory = 10)  +
      coord_flip()
    
    
    ego_gene_list <- ego_heatmap$data$Gene
    num_genes <- length(unique(ego_gene_list))
    
    if(num_genes > 50){
      ego_heatmap <- ego_heatmap + 
        theme(axis.text.y = element_text(
          size = 5, hjust = 1))
    }
    
    ego_barplot <- barplot(ego, order = TRUE)
    ego_dotplot <- dotplot(ego)
    ego_cnetplot <- cnetplot(ego)
    ego2 <- enrichplot::pairwise_termsim(ego)
    ego2_nrow <- nrow(
      ego2@result %>% 
        dplyr::filter(p.adjust < ego_pval_threshold)
    )
    if(!is.null(ego2_nrow) && ego2_nrow > 2){
      ego2_treeplot <- enrichplot::treeplot(
        ego2,
        cluster.params = list(n = min(5,ego2_nrow))
      )  
    } else {
      ego2_treeplot <- NULL
    } 
    
    
  } else {
    ego_heatmap <- NULL
    ego_barplot <- NULL
    ego_dotplot <- NULL
    ego_cnetplot <- NULL
    ego2 <- NULL
    ego2_treeplot <- NULL  
  }
  
  
  return(list(
    ego@result, 
    ego_heatmap,
    ego_barplot,
    ego_dotplot,
    ego_cnetplot,
    ego2,
    ego2_treeplot,
    ego,
    geneList
  ))
}

# Modified sampling function
sample_with_priority <- function(x, expected_genelist, sample_size) {
  # Find genes in x that match the expected_genelist
  prioritized_genes <- intersect(x, expected_genelist)
  # Exclude prioritized genes from random sampling
  remaining_genes <- setdiff(x, prioritized_genes)
  # Total number of genes to sample
  n_to_sample <- min(sample_size, length(x))
  # Combine prioritized genes with random sampling
  sampled_genes <- c(
    prioritized_genes,
    sample(remaining_genes, size = max(0, n_to_sample - length(prioritized_genes)))
  )
  return(sampled_genes) 
}

# reduced genes displayed on heatplot
clean_heatplot <- function(
    e_obj,
    sample_size = 10,
    genelist = NULL,
    show_category = 10,
    expected_genelist=NULL,
    seed = 123,
    pval_enrich = 0.05
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Perform random sampling with prioritization
  random_core_genes <- sapply(
    lapply(
      str_split(as.data.frame(e_obj)[, "geneID"], "/"), 
      function(x) sample_with_priority(x, expected_genelist, sample_size)
    ), 
    function(y) base::paste0(y, collapse = "/")
  )
  
  # random_core_genes <- sapply(lapply(str_split(as.data.frame(e_obj)[, "geneID"], "/"), 
  # function(x)
  # x[c(1:min(sample_size, length(x)))]), function(y) base::paste0(y, collapse = "/"))
  # sample(x, size = min(sample_size, length(x)))), function(y) base::paste0(y, collapse = "/"))
  
  enriched_indices <- which(e_obj@result$p.adjust < pval_enrich)
  
  e_obj@result$geneID[enriched_indices] <- random_core_genes
  clean_heat_plot <- enrichplot::heatplot(e_obj, foldChange = genelist, showCategory = show_category)
  
  clean_heat_plot <- clean_heat_plot + ggplot2::coord_flip() 
  return(clean_heat_plot)
}

# KEGG Over-Representation Analysis
enrichKEGG_pathway <- function(data2use,
                               x = "log2FC",
                               y_pval = "Pvalue",
                               pval_cutoff = 0.05,
                               ekegg_pval_threshold = 0.05,
                               log2FC_cutoff = 1.5,
                               gene_set = "all") {
  
  
  column2use <- to_vec(
    for(i in seq_along(names(data2use))) 
      if(names(data2use)[i] == x) i)
  
  
  names(data2use)[column2use] <-"log2FC"
  
  data2use <- filter_duplicate_genes(
    data2clean = data2use,
    x = "log2FC",
    y = y_pval,
    log2FC_cutoff = log2FC_cutoff,
    pval_cutoff = pval_cutoff
  )
  
  geneList <- data2use$log2FC
  names(geneList) <- data2use$GeneName
  
  hs <- org.Hs.eg.db
  my.symbols <- names(geneList)
  entrezid <- AnnotationDbi::select(hs, 
                                    keys = my.symbols,
                                    columns = c("ENTREZID", "SYMBOL"),
                                    keytype = "SYMBOL", 
                                    multiVals = "first")
  
  entrezid <- entrezid[!duplicated(entrezid$SYMBOL),]
  
  names(geneList) <- entrezid$ENTREZID
  
  geneList <- sort(geneList, decreasing = TRUE)
  
  
  if (gene_set == "all") {
    gene <- data2use %>%
      dplyr::filter((abs(log2FC) > log2FC_cutoff) &
                      (Pvalue < pval_cutoff)) %>%
      dplyr::pull(GeneName)
  } else if (gene_set == "up") {
    gene <- data2use %>%
      dplyr::filter((log2FC > log2FC_cutoff) &
                      (Pvalue < pval_cutoff)) %>%
      dplyr::pull(GeneName)
  } else if (gene_set == "down") {
    gene <- data2use %>%
      dplyr::filter((-log2FC > log2FC_cutoff) &
                      (Pvalue < pval_cutoff)) %>%
      dplyr::pull(GeneName)
  } else {
    stop(paste0("gene set \"", gene_set, "\" is not allowed!"))
  }
  
  gene <- comprehenr::to_vec(for (g in gene)
    entrezid %>% dplyr::filter(SYMBOL == g) %>% dplyr::pull(ENTREZID))
  
  print(paste0("Number of genes to test: ", as.character(length(gene))))
  print(paste0("Number of background genes: ", as.character(length(geneList))))
  
  ekegg <- clusterProfiler::enrichKEGG(gene = gene,
                                       universe = names(geneList),
                                       organism = 'hsa',
                                       pvalueCutoff = ekegg_pval_threshold)
  
  transformed_geneid <-  sapply(
    lapply(
      str_split(as.data.frame(ekegg)[, "geneID"], "/"), 
      function(x) entrezid$SYMBOL[match(x, entrezid$ENTREZID)]
    ),
    function(y) base::paste0(y, collapse = "/")
  )
  
  enriched_indices <- which(ekegg@result$p.adjust < ekegg_pval_threshold)
  ekegg@result$geneID[enriched_indices] <- transformed_geneid
  
  names(geneList) <- entrezid$SYMBOL[match(names(geneList), entrezid$ENTREZID)]
  
  if(nrow(ekegg@result %>% dplyr::filter(p.adjust < ekegg_pval_threshold)) != 0){
    ekegg_heatmap <- heatplot(ekegg, foldChange=geneList, showCategory = 10)  +
      coord_flip() +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
    
    
    ekegg_gene_list <- ekegg_heatmap$data$Gene
    num_genes <- length(unique(ekegg_gene_list))
    
    if(num_genes > 50){
      ekegg_heatmap <- ekegg_heatmap + 
        theme(axis.text.y = element_text(
          size = 5, hjust = 1))
    }
    
    ekegg_barplot <- barplot(ekegg, order = TRUE)
    ekegg_dotplot <- dotplot(ekegg)
    ekegg_cnetplot <- cnetplot(ekegg)
    
    ekegg2 <- enrichplot::pairwise_termsim(ekegg)
    ekegg2_nrow <- nrow(
      ekegg2@result %>% 
        dplyr::filter(p.adjust < ekegg_pval_threshold)
    )
    if(!is.null(ekegg2_nrow) && ekegg2_nrow > 2){
      ekegg2_treeplot <- enrichplot::treeplot(
        ekegg2,
        cluster.params = list(n = min(5,ekegg2_nrow))
      )  
    } else {
      ekegg2_treeplot <- NULL
    }    
    
  } else {
    ekegg_heatmap <- NULL
    ekegg_barplot <- NULL 
    ekegg_dotplot <- NULL 
    ekegg_cnetplot <- NULL
    ekegg2 <- NULL
    ekegg2_treeplot <- NULL
  }
  
  
  return(list(
    ekegg@result, 
    ekegg_heatmap,
    ekegg_barplot,
    ekegg_dotplot,
    ekegg_cnetplot,
    ekegg2,
    ekegg2_treeplot,
    ekegg,
    geneList
  ))
}

# Reactome Over-Representation Analysis
enrichReactome_pathway <- function(data2use,
                                   x = "log2FC",
                                   y_pval = "Pvalue",
                                   pval_cutoff = 0.05,
                                   ereactome_pval_threshold = 0.05,
                                   log2FC_cutoff = 1.5,
                                   gene_set = "all"){
  
  column2use <- to_vec(
    for(i in seq_along(names(data2use))) 
      if(names(data2use)[i] == x) i)
  
  
  names(data2use)[column2use] <-"log2FC"
  
  data2use <- filter_duplicate_genes(
    data2clean = data2use,
    x = "log2FC",
    y = y_pval,
    log2FC_cutoff = log2FC_cutoff,
    pval_cutoff = pval_cutoff
  )  
  
  geneList <- data2use$log2FC
  names(geneList) <- data2use$GeneName
  
  hs <- org.Hs.eg.db
  my.symbols <- names(geneList)
  entrezid <- AnnotationDbi::select(hs, 
                                    keys = my.symbols,
                                    columns = c("ENTREZID", "SYMBOL"),
                                    keytype = "SYMBOL", 
                                    multiVals = "first")
  
  entrezid <- entrezid[!duplicated(entrezid$SYMBOL),]
  
  names(geneList) <- entrezid$ENTREZID
  geneList <- sort(geneList, decreasing = TRUE)
  
  if (gene_set == "all") {
    gene <- data2use %>%
      dplyr::filter((abs(log2FC) > log2FC_cutoff) &
                      (Pvalue < pval_cutoff)) %>%
      dplyr::pull(GeneName)
  } else if (gene_set == "up") {
    gene <- data2use %>%
      dplyr::filter((log2FC > log2FC_cutoff) &
                      (Pvalue < pval_cutoff)) %>%
      dplyr::pull(GeneName)
  } else if (gene_set == "down") {
    gene <- data2use %>%
      dplyr::filter((-log2FC > log2FC_cutoff) &
                      (Pvalue < pval_cutoff)) %>%
      dplyr::pull(GeneName)
  } else {
    stop(paste0("gene set \"", gene_set, "\" is not allowed!"))
  }
  
  gene <- comprehenr::to_vec(for (g in gene)
    entrezid %>% dplyr::filter(SYMBOL == g) %>% dplyr::pull(ENTREZID))
  
  print(paste0("Number of genes to test: ", as.character(length(gene))))
  print(paste0("Number of background genes: ", as.character(length(geneList))))
  
  ereactome <- ReactomePA::enrichPathway(gene = gene,
                                         universe = names(geneList),
                                         organism = 'human',
                                         pvalueCutoff = ereactome_pval_threshold,
                                         readable = TRUE)
  names(geneList) <- entrezid$SYMBOL[match(names(geneList), entrezid$ENTREZID)]
  
  if(nrow(ereactome@result %>% dplyr::filter(p.adjust < ereactome_pval_threshold)) != 0){
    ereactome_heatmap <- heatplot(ereactome, foldChange=geneList, showCategory = 10)  +
      coord_flip()
    
    
    ereactome_gene_list <- ereactome_heatmap$data$Gene
    num_genes <- length(unique(ereactome_gene_list))
    
    if(num_genes > 50){
      ereactome_heatmap <- ereactome_heatmap + 
        theme(axis.text.y = element_text(
          size = 5, hjust = 1))
    }
    
    ereactome_barplot <- barplot(ereactome, order = TRUE)
    ereactome_dotplot <- dotplot(ereactome)
    ereactome_cnetplot <- cnetplot(ereactome)
    ereactome2 <- enrichplot::pairwise_termsim(ereactome)
    ereactome2_nrow <- nrow(
      ereactome2@result %>% 
        dplyr::filter(p.adjust < ereactome_pval_threshold)
    )
    if(!is.null(ereactome2_nrow) && ereactome2_nrow > 2){
      ereactome2_treeplot <- enrichplot::treeplot(
        ereactome2,
        cluster.params = list(n = min(5,ereactome2_nrow))
      )  
    } else {
      ereactome2_treeplot <- NULL
    }
    
  } else {
    ereactome_heatmap <- NULL
    ereactome_barplot <- NULL
    ereactome_dotplot <- NULL 
    ereactome_cnetplot <- NULL 
    ereactome2 <- NULL
    ereactome2_treeplot <- NULL
    
  }
  
  return(list(
    ereactome@result, 
    ereactome_heatmap,
    ereactome_barplot,
    ereactome_dotplot,
    ereactome_cnetplot,
    ereactome2,
    ereactome2_treeplot,
    geneList
  ))
}

# PANTHER / KEGG Over-Representation Analysis
enrichR_pathway <- function(data2use, 
                            database = "Panther_2016",
                            x = "log2FC",
                            y_pval = "Pvalue",
                            pval_cutoff = 0.05,
                            enrichr_pval_threshold = 0.05,
                            log2FC_cutoff = 1.5,
                            gene_set = "all"){
  
  column2use <- to_vec(
    for(i in seq_along(names(data2use))) 
      if(names(data2use)[i] == x) i)
  
  
  names(data2use)[column2use] <-"log2FC"
  
  data2use <- filter_duplicate_genes(
    data2clean = data2use,
    x = "log2FC",
    y = y_pval,
    log2FC_cutoff = log2FC_cutoff,
    pval_cutoff = pval_cutoff
  )  
  
  geneList <- data2use$log2FC
  names(geneList) <- data2use$GeneName
  
  hs <- org.Hs.eg.db
  my.symbols <- names(geneList)
  entrezid <- AnnotationDbi::select(hs, 
                                    keys = my.symbols,
                                    columns = c("ENTREZID", "SYMBOL"),
                                    keytype = "SYMBOL", 
                                    multiVals = "first")
  
  entrezid <- entrezid[!duplicated(entrezid$SYMBOL),]
  
  names(geneList) <- entrezid$ENTREZID
  
  geneList <- sort(geneList, decreasing = TRUE)
  
  if (gene_set == "all") {
    gene <- data2use %>%
      dplyr::filter((abs(log2FC) > log2FC_cutoff) &
                      (Pvalue < pval_cutoff)) %>%
      dplyr::pull(GeneName)
  } else if (gene_set == "up") {
    gene <- data2use %>%
      dplyr::filter((log2FC > log2FC_cutoff) &
                      (Pvalue < pval_cutoff)) %>%
      dplyr::pull(GeneName)
  } else if (gene_set == "down") {
    gene <- data2use %>%
      dplyr::filter((-log2FC > log2FC_cutoff) &
                      (Pvalue < pval_cutoff)) %>%
      dplyr::pull(GeneName)
  } else {
    stop(paste0("gene set \"", gene_set, "\" is not allowed!"))
  }
  
  gene <- comprehenr::to_vec(for (g in gene)
    entrezid %>% dplyr::filter(SYMBOL == g) %>% dplyr::pull(ENTREZID))
  
  print(paste0("Number of genes to test: ", as.character(length(gene))))
  print(paste0("Number of background genes: ", as.character(length(geneList))))
  
  enrich_obj <- enrichR::enrichr(gene = gene, databases = database)
  
  names(geneList) <- entrezid$SYMBOL[match(names(geneList), entrezid$ENTREZID)]
  
  if(nrow(enrich_obj[[database]] %>% dplyr::filter(Adjusted.P.value < enrichr_pval_threshold)) != 0){
    
    enrich_obj_barplot <- enrich_obj[[database]] %>%
      dplyr::filter(Adjusted.P.value < pval_threshold) %>%
      ggplot(aes(x = reorder(Term, -Adjusted.P.value), y = -log10(Adjusted.P.value))) +
      geom_bar(aes(fill = Adjusted.P.value, linewidth = Combined.Score), stat = "identity") +
      coord_flip() +
      theme(
        axis.text.y = element_text(size = 10, hjust = 1),
        plot.margin = margin(10, 10, 10, 20)
      ) +
      # theme_minimal() +
      labs(x = "Enriched Terms", y = "-log10(Adjusted P-value)") +
      scale_fill_gradient(low = "red", high = "blue", name = "Adjusted P-value")
    
    enrich_obj_dotplot <- enrich_obj[[database]] %>%
      dplyr::filter(Adjusted.P.value < pval_threshold) %>%
      ggplot(aes(x = reorder(Term, -Adjusted.P.value), y = -log10(Adjusted.P.value))) +
      geom_point(aes(fill = Adjusted.P.value, size = Combined.Score), shape = 21) +
      coord_flip() +
      theme(
        axis.text.y = element_text(size = 10, hjust = 1),
        plot.margin = margin(10, 10, 10, 20)
      ) +
      # theme_minimal() +
      labs(x = "Enriched Terms", y = "-log10(Adjusted P-value)") +
      scale_fill_gradient(low = "red", high = "blue", name = "Adjusted P-value")
    
    
    
    
    enrich_obj_gene_list <- enrich_obj_barplot$data$Gene
    num_genes <- length(unique(enrich_obj_gene_list))
    
    if(num_genes > 50){
      enrich_obj_barplot <- enrich_obj_barplot + 
        theme(axis.text.y = element_text(
          size = 5, hjust = 1)) 
      
      enrich_obj_dotplot <- enrich_obj_dotplot + 
        theme(axis.text.y = element_text(
          size = 5, hjust = 1))
    }
    
  } else {
    enrich_obj_barplot <- NULL
    enrich_obj_dotplot <- NULL
  }
  
  
  return(list(
    enrich_obj[[database]], 
    enrich_obj_barplot,
    enrich_obj_dotplot,
    enrich_obj,
    geneList
  )
  )
}


pca <- function(df){
  
  names(df) <- to_vec(
    for(cl in names(df)) 
      gsub("\\(","",
           gsub("\\)","",
                gsub("\\[","",
                     gsub("\\]","",
                          gsub("_"," ",
                               gsub("PC1","",
                                    gsub("\\([0-9][0-9].[0-9][0-9]%\\)","",cl)
                               )
                          )
                     )
                )
           )
      )
  )
  
  names_df <- names(df)
  
  names(df) <- to_vec(
    for(i in seq_along(names_df)) 
      if(i%%2 == 0) sub(" ","_",names_df[i]) else sub("  ","_",names_df[i])
  )
  
  
  df_long <- df %>%
    pivot_longer(
      cols = everything(),
      names_to = c("pc", "group"),
      names_sep = "_",
      values_to = "value"
    )
  
  
  # Print the transformed dataframe
  # Separate pc1 and pc2 into different columns and combine them
  df_pc1 <- df_long %>%
    dplyr::filter(pc == "x") %>%
    dplyr::select(-pc) %>%
    dplyr::rename(x = value)
  
  df_pc2 <- df_long %>%
    dplyr::filter(pc == "y") %>%
    dplyr::select(-pc) %>%
    dplyr::rename(y = value)
  
  df_combined <- cbind(df_pc1, df_pc2 %>% dplyr::select(-group))
  
  df_combined$group <- factor(df_combined$group)
  
  df_combined <- df_combined %>%
    dplyr::arrange(group)
  
  names(df_combined) <- c("group", "PC1", "PC2")
  
  # Create a plot with ggplot2
  return(
    list(
      df_combined,
      ggplotly(
        ggplot(df_combined, aes(x = PC1, y = PC2, color = group)) +
          geom_point(size = 3) +
          theme_minimal()
      )
    )
  )
}


combined_pca <- function(pca_file_list){
  
  pca_df_list <- pca_file_list[[1]]
  
  for(i in seq_along(pca_df_list)){
    
    curr_df = readr::read_tsv(pca_df_list[i], show_col_types = FALSE)
    curr_df_pca_results <- pca(curr_df)
    
    curr_df_pca = curr_df_pca_results[[1]]
    curr_df_pca$group <- factor(curr_df_pca$group)
    
    print(curr_df_pca)
    
    groups <- levels(curr_df_pca$group)
    
    experiment <- paste0(groups[1]," vs ", groups[2])
    
    curr_df_pca["experiment"] <- experiment
    
    if(i==1){
      total_df <- curr_df_pca
    } else {
      total_df <- rbind(total_df, curr_df_pca)
    }
    
  }
  
  total_pca <- total_df %>%
    ggplot2::ggplot(aes(PC1,PC2,color=group)) + 
    ggplot2::geom_point(size = 3) +
    ggplot2::facet_wrap(~experiment) 
  
  total_pca <- ggplotly(total_pca)
  
  merged_pca <- total_df %>%
    ggplot2::ggplot(aes(PC1, PC2, color=group, shape = experiment)) + 
    ggplot2::geom_point(size = 3)
  
  merged_pca <- ggplotly(merged_pca)
  
  return(list(total_pca, merged_pca))
}



