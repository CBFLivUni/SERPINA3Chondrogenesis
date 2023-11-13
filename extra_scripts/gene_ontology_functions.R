
ReadInGMT <- function(db)
{
  # Get filename with entrezids
  gmt_files <- paste(gmt_dir,list.files(gmt_dir)[grep(db,list.files(gmt_dir),fixed=TRUE)],sep="/")
  gmt_file <- gmt_files[grep("entrez",gmt_files)]
  
  # Read in gmt file
  pathways <- read.gmt(gmt_file)
  
  # Return
  return(pathways)
}



RunGSEA <- function(genes, database_names=".", db_list, id_col="ENTZ", stat_col="log2FoldChange", group_col="TRAIT", group_name, minimum_gs=5, id="", dir=".", alpha_thr=0.05)
{
  # Format to genelist
  gene_list <- genes[[stat_col]]; names(gene_list) <- genes[[id_col]]
  gene_list <- gene_list[is.finite(gene_list)]
  
  # Remove NAs, non-numeric & duplicated values
  gene_list <- gene_list[!is.na(gene_list) & is.numeric(gene_list)]
  gene_list <- gene_list[!is.na(names(gene_list))]
  gene_list <- gene_list[!duplicated(gene_list)]
  
  # Print number of genes
  print(paste0(">>>>> RUNNING GSEA vs geneset of size :- ",length(gene_list)))
  
  # Loop over databases and run GSEA
  gseas <- lapply(unique(database_names), function(db_name)
  {
    print(paste("   >> Running GSEA against :-",db_name))
    
    # Get db-specific terms
    sub.db_list <- db_list[grepl(db_name, names(db_list), fixed=TRUE)]
    
    # Run fGSEA
    gsea <- fgsea(pathways = sub.db_list[[db_name]], 
                  stats    = gene_list,
                  # nperm    = 1000, # Warning in fgsea(pathways = db, stats = gene_list, nperm = 1000, minSize = 5,  : // You are trying to run fgseaSimple. It is recommended to use fgseaMultilevel. To run fgseaMultilevel, you need to remove the nperm argument in the fgsea function call.
                  minSize  = minimum_gs,
                  maxSize  = 500)
    
    # Get sig GSEA
    sig_gsea <- gsea[gsea$padj < alpha_thr,]
    
    # Plot leading edges for sig pathways
    print(".....Saving Leading Edge plots.....")
    plot_dir <- paste0(dir, "plots/"); dir.create(plot_dir, recursive=FALSE, showWarnings=FALSE)
    for (term in unique(sig_gsea$pathway)) {
      
      # Get term name, adding new lines in for long titles (for title)
      term_name <- strsplit(term, split="_")[[1]]
      for (i in 1:length(term_name)) { if (i %% 10 == 0) { term_name[[i]] <- paste0(term_name[[i]],"\n") } }
      term_name <- paste0(term_name, collapse=" ")
      
      # Plot
      leadedge_plt <- plotEnrichment(sub.db_list[[db_name]][[term]], gene_list) + 
        labs(title=mgsub(term_name,c("_","REACTOME"),c(" ",""))) + 
        xlab("Rank") + ylab("Enrichment Score") + 
        theme(axis.title=element_text(size=20),
              axis.text=element_text(size=16),
              title=element_text(size=14)) 
      
      # Save
      ggsave(paste0(plot_dir,id,".",term,".leading_edge_plt.png"),  leadedge_plt, units="px", width=4000, height=2500)
 
    }
    
    # Make df
    gsea <- data.frame(gsea)
    
    if (is.null(gsea) == TRUE | nrow(gsea) == 0) 
    {
      gsea <- NULL
      
    } else {
      
      # Add ontology
      gsea["Ontology_Source"] <- db_name
      gsea[group_col] <- group_name
      
    }

    
    # Return
    return(gsea)
  })

  # Remove nulls
  gseas <- do.call("rbind",gseas[!is.null(gseas)])
  
  # Save
  save_file <- paste0(dir, group_name,".stat_is_",stat_col, ".gsea_results.tsv")
  write.table(data.frame(gseas[,!grepl("leadingEdge",colnames(gseas))]), save_file, sep="\t", row.names=FALSE, quote=FALSE) # Need to remove leading edge as it prevents saving
  
  # Return
  return(gseas)
}


