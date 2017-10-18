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

# Merge annotation files outputted by IMG ER annotation pipeline into 
# one merged dataframe that can be used for further analysis
# Input is list of IMG annotation dirs, e.g. :
# c("./IMG_annotation/IMG_2737471681", "./IMG_annotation/IMG_2737471682")
merge_annotations <- function(paths_gen, genoid_seqid = FALSE){
  cat(" --- I will merge the annotation files from the following genomes:\n" ,sep = " ")
  print(data.frame(Genomes = paste(paths_gen)))
  sum_tmp <- 0
  for(genomes in paths_gen){
    paths_ann <- list.files(list.dirs(paste(genomes, "/IMG_Data", sep = ""))[2],
                            full.names = TRUE)
    
    # Only retain paths to functional annotation files
    paths_ann <- paths_ann[grep(pattern = "ko|pfam|tigrfam|cog|gene_oid", x = paths_ann)]
    
    # Now extract relevant information from all these files
    input_trimmed_ko <- paths_ann[grep("*ko*", paths_ann)] %>% 
      read.table(header = TRUE, sep ="\t", quote=NULL, comment='',
                 stringsAsFactors = FALSE) %>% 
      dplyr::select(one_of(c("gene_oid", "ko_id",
                             "ko_name")))
    input_trimmed_cog <- paths_ann[grep("*cog*", paths_ann)] %>% 
      read.table(header = TRUE, sep ="\t", quote=NULL, comment='',
                 stringsAsFactors = FALSE) %>% 
      dplyr::select(one_of(c("gene_oid", "cog_id",
                             "cog_name")))
    input_trimmed_pfam <- paths_ann[grep("*pfam*", paths_ann)] %>% 
      read.table(header = TRUE, sep ="\t", quote=NULL, comment='',
                 stringsAsFactors = FALSE) %>% 
      dplyr::select(one_of(c("gene_oid", "pfam_id",
                             "pfam_name")))
    input_trimmed_tigrfam <- paths_ann[grep("*tigrfam*", paths_ann)] %>% 
      read.table(header = TRUE, sep ="\t", quote=NULL, comment='',
                 stringsAsFactors = FALSE) %>% 
      dplyr::select(one_of(c("gene_oid", "tigrfam_id",
                             "tigrfam_name")))
    # Merge these files into one dataframe
    input_trimmed_ko$gene_oid <- as.character(input_trimmed_ko$gene_oid)
    input_trimmed_cog$gene_oid <- as.character(input_trimmed_cog$gene_oid)
    input_trimmed_pfam$gene_oid <- as.character(input_trimmed_pfam$gene_oid)
    input_trimmed_tigrfam$gene_oid <- as.character(input_trimmed_tigrfam$gene_oid)
    
    merged <- dplyr::full_join(input_trimmed_ko, input_trimmed_cog, by = "gene_oid")
    merged <- dplyr::full_join(merged, input_trimmed_pfam, by = "gene_oid")
    merged <- dplyr::full_join(merged, input_trimmed_tigrfam, by = "gene_oid")
    
    # This will allow the inclusion of all seqids insteads of just the annotated
    # genes
    if(genoid_seqid == TRUE){
      input_trimmed_allgenid <- paths_ann[grep("*gene_oid*", paths_ann)] %>% 
        read.table(header = FALSE, sep =" ", quote=NULL, comment='',
                   stringsAsFactors = FALSE) %>% 
        dplyr::select(one_of(c("V1")))
      colnames(input_trimmed_allgenid) <- "gene_oid"
      input_trimmed_allgenid$gene_oid <- as.character(input_trimmed_allgenid$gene_oid)
      merged <- dplyr::full_join(merged, input_trimmed_allgenid, by = "gene_oid")
      print(nrow(input_trimmed_allgenid))
    }
    
    # add genome ID
    merged$Genome_id <- gsub(genomes, pattern="^.*IMG_", replacement = "")
    sum_tmp <- length(unique(merged$gene_oid)) + sum_tmp
    print(sum_tmp)
    
    if(genomes != paths_gen[1]) {
      merged_final <- rbind(merged_final, merged)
    } else merged_final <- merged
  }
  cat(date(), ' --- Sucessfully merged files\n')
  return(merged_final)
}


# Modified igraph functions

plot_network_custom <- function (g, physeq = NULL, type = "samples", color = NULL, shape = NULL, 
                                 point_size = 4, alpha = 1, label = "value", hjust = 1.35, 
                                 line_weight = 0.5, line_color = color, line_alpha = 0.4, 
                                 layout.method = layout.fruchterman.reingold, title = NULL, label_size = 3) 
{
  if (vcount(g) < 2) {
    stop("The graph you provided, `g`, has too few vertices. \\n         Check your graph, or the output of `make_network` and try again.")
  }
  if (type %in% c("taxa", "species", "OTUs", "otus", "otu")) {
    type <- "taxa"
  }
  edgeDF <- data.frame(get.edgelist(g))
  edgeDF$id <- 1:length(edgeDF[, 1])
  vertDF <- layout.method(g)
  colnames(vertDF) <- c("x", "y")
  vertDF <- data.frame(value = get.vertex.attribute(g, "name"), 
                       vertDF)
  if (!is.null(physeq)) {
    extraData <- NULL
    if (type == "samples" & !is.null(sample_data(physeq, 
                                                 FALSE))) {
      extraData = data.frame(sample_data(physeq))[as.character(vertDF$value), 
                                                  , drop = FALSE]
    }
    else if (type == "taxa" & !is.null(tax_table(physeq, 
                                                 FALSE))) {
      extraData = data.frame(tax_table(physeq))[as.character(vertDF$value), 
                                                , drop = FALSE]
    }
    if (!is.null(extraData)) {
      vertDF <- data.frame(vertDF, extraData)
    }
  }
  graphDF <- merge(reshape2::melt(edgeDF, id = "id"), vertDF, 
                   by = "value")
  p <- ggplot(vertDF, aes(x, y))
  p <- p + theme_bw() + theme(panel.grid.major = element_blank(), 
                              panel.grid.minor = element_blank(), axis.text.x = element_blank(), 
                              axis.text.y = element_blank(), axis.title.x = element_blank(), 
                              axis.title.y = element_blank(), axis.ticks = element_blank(), 
                              panel.border = element_blank())
  p <- p + geom_point(aes_string(color = color, shape = shape), 
                      size = point_size, na.rm = TRUE)
  if (!is.null(label)) {
    p <- p + geom_text(aes_string(label = label), size = label_size, 
                       hjust = hjust, na.rm = TRUE)
  }
  p <- p + geom_line(aes_string(group = "id", color = line_color), 
                     graphDF, size = line_weight, alpha = line_alpha, na.rm = TRUE)
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}