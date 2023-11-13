

# Format a PVCA (derived from )
ProcessPVCA <- function(res, method_key, source_toggle=TRUE) {
  
  # Get as df
  df <- data.frame(as.data.frame(res$label), t(as.data.frame(res$dat)), Method=method_key)
  colnames(df) <- c("Effect", "Variance", "Method")
  
  # Add effect type
  df["Type"] <- "Single"
  df[grepl("[:]",df$Effect),]["Type"] <- "Int"
  df["Int_Type"] <- ""
  
  # Add interaction type
  df[grepl("LESION",df$Effect),]["Int_Type"] <- "Lesion"
  df[grepl("KSHV_TYPE",df$Effect),]["Int_Type"] <- "KSHV Type"
  df[grepl("GENDER",df$Effect),]["Int_Type"] <- "Gender"
  if (source_toggle == TRUE) { df[grepl("SOURCE",df$Effect),]["Int_Type"] <- "Source" }
  df[grepl("BATCH",df$Effect),]["Int_Type"] <- "Batch"
  df[!grepl("[:]",df$Effect),]["Int_Type"] <- ""
  
  # Return
  return(df)
}


# Function to loop over factor of interest and run pairwise DE between treatment/control
RunGroupwiseDE <- function(txi_obj, meta, annot, group_col, contrast, design_formula, dir=".") {
  
  # Packages
  require(DESeq2)
  require(apeglm)
  require(ashr)
  
  # Loop over days
  ds2_res <- lapply(unique(meta[[group_col]]), function(group) { 
    
    # Subset metadata
    print(paste0(">> On ",group))
    sub_meta <- meta[meta[[group_col]] == group,]
    
    # Adjust txi colnames
    sub.txi_obj <- txi_obj
    sub.txi_obj$abundance <- sub.txi_obj$abundance[rowMeans(sub.txi_obj$abundance) > 0,rownames(sub_meta)]
    sub.txi_obj$counts<- sub.txi_obj$counts[rownames(sub.txi_obj$abundance),rownames(sub_meta)]
    sub.txi_obj$length <- sub.txi_obj$length[rownames(sub.txi_obj$abundance),rownames(sub_meta)]
    
    # Create dds object from Txi
    dds <- DESeqDataSetFromTximport(sub.txi_obj, colData = sub_meta[colnames(sub.txi_obj$counts),], as.formula(design_formula))
    dds <- DESeq(dds, minReplicatesForReplace=100)
    # boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
    
    # LFC shrinkake
    res  <- lfcShrink(dds, contrast = contrast, type="ashr")
    res_df <- data.frame(res); res_df <- data.frame(merge(annot, res_df, by.x="ENSG", by.y="row.names"), GROUP=group)
    
    # Save
    save_name <- paste0(dir, gsub("^[.]","",group),"-",group_col,".contrast-",paste(contrast, collapse="_"),".deseq2_results.tsv")
    print(paste0("..... Saving to : ",save_name,"....."))
    
    # Save
    write.table(res_df, save_name, quote=FALSE, row.names=FALSE,sep="\t")
    
    # Return
    return(res_df)
  })
  
  # Return
  names(ds2_res) <- unique(meta[[group_col]])
  return(ds2_res)
}




NiceVolcanoPlots <- function(res, comparison, dir, padj_thr=0.01, topn_genes=10, label_sig_thr=1e-20, n.max_ovlp=25) {
  
  # Packages
  require(ggplot2)
  require(ggrepel)
  
  # Loop over results df
  lapply(unique(res[[comparison]]), function(contrast) {
    
    # Subset by contrast
    res_df <- res[res[[comparison]] == contrast,]
    
    # Remove NA
    res_df <- res_df[!is.na(res_df$padj),]
    
    # Get sig
    sig.res_df <- res_df[res_df$padj < padj_thr & !is.na(res_df$padj),]
    
    # Add sig key
    res_df["Sig"] <- "Non-Sig"
    res_df[which(res_df$ENSG %in% sig.res_df$ENSG),]["Sig"] <- paste0("FDR<",padj_thr)
    
    # Add lfc key
    res_df["Direction"] <- "Up"
    res_df[res_df$log2FoldChange < 0,]["Direction"] <- "Down"
    
    # Add colour key
    res_df["ColorKey"] <- paste(res_df$Direction, res_df$Sig, sep=" ")
    res_df[grepl("Non-Sig",res_df$ColorKey),]["ColorKey"] <- "Non-Sig"
    
    # Order plotting
    res_df$ColorKey <- factor(res_df$ColorKey , levels=c("Non-Sig", paste0("Down FDR<",padj_thr), paste0("Up FDR<",padj_thr)))
    
    # Order by LFCs (for labelling)
    res_df <- res_df[order(res_df$log2FoldChange, decreasing=TRUE),]
    
    # Define genes to label
    sig.res_df <- res_df[res_df$padj < padj_thr & !is.na(res_df$padj),]
    #label_df <- rbind(head(sig.res_df, topn_genes), tail(sig.res_df, topn_genes), sig.res_df[sig.res_df$padj < label_sig_thr,])
    sig.res_df$pi <- sig.res_df$log2FoldChange * -log10(sig.res_df$padj)
    sig.res_df <- sig.res_df[order(sig.res_df$pi, decreasing=TRUE),]
    label_df <- rbind(head(sig.res_df, topn_genes), tail(sig.res_df, topn_genes),sig.res_df[sig.res_df$padj < label_sig_thr,])
    label_df <- label_df[!duplicated(label_df),]
    
    
    
    
    # Add size key
    res_df["SHAPE"] <- "Circle"
    res_df[which(res_df$ENSG %in% label_df$ENSG),]["SHAPE"] <- "Triangle"
    
    # Plot
    plt <- ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), fill=ColorKey, color=ColorKey)) + 
      geom_vline(xintercept=0, linetype="dashed", color="darkgrey") +
      geom_point(alpha=0.75, shape=21, size=3) + 
      geom_text_repel(data=label_df, mapping=aes(x=log2FoldChange, y=-log10(padj), label=SYMBOL), size=7,
                      show.legend = FALSE, min.segment.length = 0, max.overlaps=Inf,seed = 12345) +
      scale_fill_manual(values=c("Non-Sig"="lightgrey","Up FDR<0.01"="#BB5566", "Down FDR<0.01"="#aed2e6")) +
      scale_color_manual(values=c("Non-Sig"="lightgrey","Up FDR<0.01"="black", "Down FDR<0.01"="black")) +
      geom_hline(yintercept=-log10(padj_thr), linetype="dashed", color="darkgrey") +
      ylab("-log10 FDR") + xlab("Log2 FC") + ggtitle(mgsub::mgsub(contrast, c("-Day_","_"), c(": "," vs "))) +
      xlim(-10,10) +
      theme_bw(base_size=18) + theme(legend.title=element_blank(), legend.position="top", legend.direction="horizontal", 
                                     axis.title.x = element_text(size=26, face="bold"), axis.title.y = element_text(size=26, face="bold"),
                                     plot.title = element_text(size=40), legend.text = element_text(size=20),
                                     axis.text = element_text(size=20)) +
      guides(color = guide_legend(override.aes = list(size=10)), text=FALSE)

    # Save
    setwd(dir)
    ggsave(paste0(contrast, ".volcano_plot.padj_thr",padj_thr,".png"), plt, units="in", width=8, height=9)
  })
}

#get the cosine similarity to all the datasets
getSim <- function(dataset,foldchangeTable,metadata){
  sim <- sapply(seq_along(foldchangeTable)[-ncol(foldchangeTable)],cos.sim,foldchangeTable,dataset)
  sim <- data.frame(colnames(foldchangeTable)[-ncol(foldchangeTable)],sim)
  sim <- merge(sim,metadata,by.x=1,by.y="combined")
  sim$zscore <- scale(sim$sim)
  colnames(sim)[1] <- "ID"
  sim
  
}

cos.sim <- function(i,X,query) {
  
  A <- query[,grep("log|GeneSymbol|HGNC.symbol|gene_name",colnames(query))]
  B  <-  X[,i,drop=F]
  merged <- merge(A,B,by.x=1,by.y="row.names")
  
  A <- merged[,2]
  B <- merged[,3]
  
  return( sum(A*B,na.rm = TRUE)/sqrt(sum(A^2,na.rm = TRUE)*sum(B^2,na.rm = TRUE)) )
}

runSkeletalVis <- function(query,skeletalVisTable,metadata){
  
  query <- na.omit(query)
  query <- query[,c("SYMBOL","log2FoldChange","padj")]
  colnames(query)[1] <- "gene_name"
  similarityScores <-  getSim(query,skeletalVisTable,metadata)
  similarityScores <- similarityScores[ order(similarityScores$zscore,decreasing =TRUE),]
  colnames(similarityScores)[2] <- "cosineSimilarity"
  similarityScores <- similarityScores %>% select(ID, cosineSimilarity, zscore,everything()) %>% as.data.frame()
  
  return(similarityScores)
  
}
