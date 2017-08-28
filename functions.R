# Function for creating a common legend for 2 ggplot2 figures.
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid.newpage()
  grid.draw(combined)
  
}

# Function to extract %GC and geneID from gff file

extract_gc_from_gff <- function(inputgff, outputFolder){
  
  for(i in inputgff){
    d1 <- read.table(i, as.is = TRUE, sep = "\t")
    d2 <- t(sapply(strsplit(d1[[9]],";"), `[`, c(1:4)))
    d <- cbind(d1[,1], d2[,c(1,3)])
    colnames(d) <- c("contig", "contig_geneID", "GC")
    d <- data.frame(d)
    
    # Remove unwanted characters
    d$GC <- gsub(d$GC, pattern = "gc_cont=", replacement = "")
    d$contig_geneID <- gsub(d$contig_geneID, pattern = "ID=", replacement = "")
    d$contig_geneID <- gsub(d$contig_geneID, pattern = ".", replacement = "",
                            fixed = TRUE)
    d$GC <- as.numeric(d$GC); d <- d[!is.na(d$GC),]
    d <- data.frame(d, Genome = rep(gsub(".*/","",i), nrow(d)))
    
    if(i!=inputgff[1]) d.tot <- rbind(d.tot, d) else d.tot <- d
    # Export data
    write.table(d.tot, paste(outputFolder,"/","seqid_GC_",gsub(".*/","",i), ".tsv", sep = ""), row.names = FALSE,
                quote = FALSE)
  }
  
  # Plot
  # p <- easyGgplot2::ggplot2.histogram(data = d.tot, xName = 'GC',
  #                   groupName = 'Genome', alpha = 0.5,
  #                   legendPosition = "top", binwidth = 0.01)
  # print(p)
}

# Function to link gc content with gene functions

gc2function <- function(seq_id_gc, gene_id_seq_id, functions, gc_thresh = 0.75,
                        output = FALSE, label){
  # Import data
  seqid_GC <- read.table(seq_id_gc,
                         header = TRUE, stringsAsFactors = FALSE, sep = " ")
  gene_oid_2_seq_id <- read.table(gene_id_seq_id,
                                  stringsAsFactors = TRUE, col.names = c("gene_oid", "seq_id"),
                                  colClasses = "character", sep = " ")
  function_file <- read.table(functions,
                              stringsAsFactors = FALSE,  sep = "\t", header = TRUE,
                              quote=NULL, comment='')
  function_file$gene_oid <- as.character(function_file$gene_oid)
  
  # Merge files
  merged_df <- dplyr::inner_join(seqid_GC, gene_oid_2_seq_id, by = c("contig_geneID" = "seq_id"))
  merged_df <- dplyr::inner_join(merged_df, function_file, by = "gene_oid")
  tot_genes <- nrow(merged_df)
  merged_df <- data.frame(merged_df, genome_id = rep(label, nrow(merged_df)))
  
  # Filter out based on threshold
  merged_df <- merged_df[merged_df$GC>gc_thresh,]
  
  # Rank from large to small
  merged_df <- merged_df[order(merged_df$GC), ]
  
  # Print brief summary
  cat(date(), ' --- There are', nrow(merged_df), "genes with >", gc_thresh, "%\n" ,sep = " ")
  cat(date(), ' --- This is', round(100*nrow(merged_df)/tot_genes,2), "% of all genes\n",sep = " ")
  cat(date(), ' --- The 10 genes with the highest GC% are:\n', sep = " ")
  output <- tail(data.frame(function_id = merged_df[, ncol(merged_df)-2], function_name = merged_df[, ncol(merged_df)-1], 
                            GC = 100*merged_df$GC), n = 10)
  print(output)
  # Write output table
  if(output == TRUE){
    write.table(merged_df, "GC_function.tsv",
                row.names = FALSE, quote = FALSE)
  }
  return(merged_df)
}
