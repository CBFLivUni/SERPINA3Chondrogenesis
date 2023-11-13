
RunPCA <- function(m)
{
  # Run PCA
  res <- prcomp(t(m), center=TRUE, scale=TRUE) # Scale and center values
}



Plot2DPCA <- function(meta, var, pcs=c(1,2), color_by=NA, shape_by=NA, default_color="black")
{
  # Set cell lines and axis labels
  xlabel <- paste("PC",pcs[[1]]," ~ ",round(var$variance.percent[pcs[[1]]],2),"%",sep="")
  ylabel <- paste("PC",pcs[[2]]," ~ ",round(var$variance.percent[pcs[[2]]],2),"%",sep="")
  
  if (is.na(shape_by))
  {
  
    # Plot PCA plot, coloured by color-by
    pca_plt <- ggplot(meta,
                      aes(x=.data[[paste0("PC",pcs[[1]])]],
                          y=.data[[paste0("PC",pcs[[2]])]], 
                          fill=.data[[color_by]],
                          color=.data[[color_by]])) + 
      geom_point(size=3, alpha=1, color=default_color, shape=21) + xlab(xlabel) + ylab(ylabel) + 
      geom_vline(xintercept=0,linetype="dashed",color="grey") + 
      geom_hline(yintercept=0,linetype="dashed",color="grey") +
      theme_bw(base_size=16) + theme(legend.position = "none")
  
  } else {
    
    # Plot PCA plot, coloured by color-by
    pca_plt <- ggplot(meta,
                      aes(x=.data[[paste0("PC",pcs[[1]])]],
                          y=.data[[paste0("PC",pcs[[2]])]], 
                          shape=.data[[shape_by]],
                          group=.data[[shape_by]],
                          fill=.data[[color_by]],
                          color=.data[[color_by]])) + 
      geom_point(size=3, alpha=1, color=default_color) + xlab(xlabel) + ylab(ylabel) + 
      theme_bw() + 
      geom_vline(xintercept=0,linetype="dashed",color="grey") + 
      geom_hline(yintercept=0,linetype="dashed",color="grey") +
      theme_bw(base_size=16) + theme(legend.position = "none")
    
  }
    
  # Return
  return(pca_plt)
}

