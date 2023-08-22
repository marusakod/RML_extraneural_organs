library(tidyverse)
library(WGCNA)
library(openxlsx)
library(strex)
library(DESeq2)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(parallel)
library(clusterProfiler)
library(scico)
library(igraph)
library(readr)
library(randomcoloR)
library(ggnewscale)
library(testit)
library(biomaRt)
library(plotly)

get_feature_counts <- function(counts_dir, organ,meta, file_prefix, overwrite = FALSE){

  # make a file with feature counts for each organ if it doesn't already exist
  f <- file.path(counts_dir, organ, paste0(file_prefix, 'merged_FeatureCounts.rds'))

  if(file.exists(f) & overwrite == FALSE){
    readRDS(f)
  }else{
    # read all txt files in feature counts dir for organ
    dir <- file.path(counts_dir, organ)
    # get counts for all samples in metadata
    counts <- list.files(dir, pattern = paste(meta$SampleID, collapse = "|"), full.names = TRUE) %>%
      lapply(., function(x){
      c <- readr::read_table(x)
      colnames(c)[2] <-gsub(paste0('data/FeatureCounts/', organ, '/(.+).txt'), '\\1', x)
      c %>% column_to_rownames('Identifier') %>% as.matrix()
    })

    merged <- do.call('cbind', counts) %>%
    # adjust column order in counts to match row order in metadata
      .[, meta$SampleID]
    saveRDS(merged, f)
    print(paste0('Feature counts created for ', organ))
    return(merged)
  }

}


# find gene symbols and biotypes for all ensembl ids

get_feature_data <- function(ensembl_vec, symbols, mart_dataset, species, out_dir, overwrite = FALSE){

  f <- file.path(out_dir, paste0(species, '_feature_data.rds'))

  if(file.exists(f) & overwrite == FALSE){
    readRDS(f)
  }else{

  df <- data.frame(ensembl_gene_id = ensembl_vec)
  #
  mart <- useMart('ensembl', dataset = mart_dataset)
  genes_and_biotypes <- getBM(filters = 'ensembl_gene_id',
                              attributes = c('ensembl_gene_id', 'gene_biotype', symbols),
                              values = ensembl_vec, mart = mart)

  res <- merge(df, genes_and_biotypes, by = 'ensembl_gene_id', all.x = TRUE) %>%
    dplyr::rename('gene_symbol' = symbols) %>%
    group_by(ensembl_gene_id, gene_biotype) %>%
    dplyr::slice_head() %>%
    ungroup() %>%
    as.data.frame() %>%
    # add broad biotype categories
    mutate(gene_biotype_broad = case_when(str_detect(gene_symbol, '^mt-') ~ 'Mitochondrial RNA',
                                          str_detect(gene_biotype , '^IG_') ~ 'IG gene',
                                          str_detect(gene_biotype, '^TR_') ~ 'TR gene',
                                          str_detect(gene_biotype, 'pseudogene') ~ 'Pseudogene',
                                          str_detect(gene_biotype, 'misc_RNA|miRNA|snRNA|snoRNA|rRNA|scRNA|scaRNA') ~ 'Non-coding RNA',
                                          str_detect(gene_biotype, 'lncRNA') ~ 'Long non-coding RNA',
                                          str_detect(gene_biotype, 'protein_coding') ~ 'Protein-coding',
                                          TRUE ~ 'other'))

  rownames(res) <- res$ensembl_gene_id
  res <- res[ensembl_vec, ]
  saveRDS(res, f)
  return(res)
  }
}


# get counts in metadata

extract_group_counts <- function(organ, meta){

  all_c <- all_FC_counts_list[[organ]]
  sub <- all_c[, meta$SampleID]

  sub
}




# add qc metrics to metadata

add_qc <- function(meta, counts, organ, feature_meta, file_prefix, out_dir, overwrite = FALSE){

  f <- file.path(out_dir, paste0(file_prefix, organ, '_metadata_w_qc.rds'))

  if(file.exists(f) & overwrite == FALSE){

    readRDS(f)
  }else{
  assert("Rownamws of metadata and column names of counts don't match",
         identical(rownames(meta), colnames(counts)))

  # get library size and the number of genes detected
  meta$lib_size <- Matrix::colSums(counts)
  meta$feature_counts <- Matrix::colSums(counts > 0)

  # get pecentages of different biotypes
  biotype_counts <- sapply(unique(feature_meta$gene_biotype_broad),
                           FUN = get_pct_per_biotype,
                feature_data = feature_meta,
                counts = counts)

  final <- cbind(meta, biotype_counts)
  saveRDS(final, f)
  return(final)
  }


}



get_pct_per_biotype <- function(biotype, feature_data, counts){
  biotype_ensembls <- feature_data %>% filter(gene_biotype_broad == biotype) %>%
    pull(ensembl_gene_id)

  biotype_counts <- counts[biotype_ensembls, ]
  biotype_percent <- Matrix::colSums(biotype_counts)
  biotype_percent
}


# make percent barplot for biotypes
make_pct_biotype_bar <- function(meta_qc, feature_data){
  bios <- unique(feature_data$gene_biotype_broad)
  df <- meta_qc %>% dplyr::select(all_of(c('SampleID', bios))) %>%
    reshape2::melt()

  ggplot(df, aes(fill=variable, y=value, x=SampleID)) +
    geom_bar(position="fill", stat="identity") +
    theme_light() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(x = NULL, y ='% total', fill = 'Biotype')

}


get_library_complexity <- function(counts, meta, cols_to_add){ # from #RNAseqQC package

 lib_comp <- lapply(1: ncol(counts), FUN=function(x) {
    cts <- sort(counts[, x], decreasing=T)
    tibble(
      sample_id = SummarizedExperiment::colnames(counts)[x],
      fraction_of_counts = cumsum(cts)/sum(cts),
      fraction_of_genes = (1:length(cts))/length(cts)
    )
  })
}



###### PCA
do_pca <- function(counts){
  vst <- vst(counts, blind = TRUE)
  rv <- rowVars(vst)
  gs <- order(rv, decreasing=TRUE)[seq_len(500)]
  pca <- prcomp(t(vst[gs, ]))
  pca

}

get_pca_meta_df <- function(meta,pca){
  # meta must have a SampleID column matching column names in pca$x
  pc_df <- as.data.frame(pca$x[, 1:3]) %>%
    rownames_to_column('SampleID') %>%
    merge(., meta, by = 'SampleID')

  vars <- round((pca$sdev^2)/sum(pca$sdev^2)*100, digits = 1)
  pc_df$PC1_var <- vars[1]
  pc_df$PC2_var <- vars[2]
  pc_df$PC3_var <- vars[3]
  pc_df
}


# get top 20 PC1 and PC2 loadings

get_top_loadings <- function(pca, feature_data){
  # first component
  pc1 <- names(sort(abs(pca$rotation[, 1]), decreasing = TRUE)[1:50])
  pc2 <- names(sort(abs(pca$rotation[, 2]), decreasing = TRUE)[1:50])

  pc1u <- setdiff(pc1, pc2)
  pc2u <- setdiff(pc2, pc1)
  both <- intersect(pc1, pc2)

  top_genes <- c(pc1u,pc2u,both)

  top <- pca$rotation[top_genes, 1:2] %>%
    as.data.frame() %>%
    mutate(pc = c(rep('PC1', length(pc1u)), rep('PC2', length(pc2u)), rep('PC1_PC2', length(both)))) %>%
    rownames_to_column('ensembl_gene_id') %>%
    merge(., feature_data, by = 'ensembl_gene_id', all.x = TRUE, all.y = FALSE) %>%
    arrange(pc)

  return(top)

}

# it might work if I use df.u for plotting the parent ggplot instead of oroginal PCs

get_top_loadings2 <- function(pcobj, feature_data){
  nobs.factor <- sqrt(nrow(pcobj$x) - 1)
  d <- pcobj$sdev
  u <- sweep(pcobj$x, 2, 1 / (d * nobs.factor), FUN = '*')
  v <- pcobj$rotation

  choices = 1:2
  scale = 1
  obs.scale = 1 -scale
  var.scale = scale
  circle = FALSE
  circle.prob = 0.69
  pc.biplot = TRUE

  # Scores
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[,choices], 2, d[choices]^obs.scale, FUN='*'))

  # Directions
  v <- sweep(v, 2, d^var.scale, FUN='*')
  df.v <- as.data.frame(v[, choices])

  names(df.u) <- c('xvar', 'yvar')
  names(df.v) <- names(df.u)

  if(pc.biplot) {
    df.u <- df.u * nobs.factor
  }

  # Scale the radius of the correlation circle so that it corresponds to
  # a data ellipse for the standardized PC scores
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)

  # Scale directions
  v.scale <- rowSums(v^2)
  df.v <- r * df.v / sqrt(max(v.scale))
  colnames(df.v) <- c('PC1', 'PC2')
  df.v <- as.matrix(df.v)

  pc1 <- sort(abs(df.v[, 1]), decreasing = TRUE)[1:20]
  pc2 <- sort(abs(df.v[, 2]), decreasing = TRUE)[1:20]
  top_genes <- unique(names(c(pc1, pc2)))

  top <- df.v[top_genes, 1:2] %>%
    as.data.frame() %>%
    rownames_to_column('ensembl_gene_id') %>%
    merge(., feature_data, by = 'ensembl_gene_id', all.x = TRUE, all.y = FALSE)

  return(top)
}

# function for making a pca plot for mouse organs
make_pca_plot <- function(pca_df,
                          sample_label = FALSE){

  pca_df_RML <- pca_df %>% filter(Treatment == 'RML6')
  pca_df_NBH <- pca_df %>% filter(Treatment == 'NBH')
  p <- ggplot() +
    geom_point(data =  pca_df_RML,
               aes(x = PC1,
                   y = PC2,
                   shape = batch,
                   color = wpi),
               size = 2)  +
    scale_color_brewer(palette = 'Oranges')+
    ggnewscale::new_scale_color() +
    geom_point(data = pca_df_NBH,
               aes(x = PC1,
                   y = PC2,
                   shape = batch,
                   color = wpi),
               size = 2) +
    scale_color_brewer(palette = 'Greys') +
  # facet_wrap(~Region) +

  xlab(paste("PC1: ", unique(pca_df$PC1_var), "% variance", sep="")) +
  ylab(paste("PC2: ", unique(pca_df$PC2_var), "% variance", sep="")) +
  theme_light() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1/1,
        legend.position = 'none',
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        strip.background = element_rect(color = '#99A3A4', fill = 'white'),
        panel.background = element_rect(color = '#99A3A4'),
        strip.text = element_text(size = 9, color = 'black')
        )
  p

}



make_pca_plot2 <- function(pca_df,
                           proportional_axes = FALSE,
                          sample_label = FALSE){



  p <- ggplot(data = pca_df, aes(x =PC1, y = PC2,fill = Treatment_stage)) +
    geom_point(shape = 21) +
    xlab(paste("PC1: ", unique(pca_df$PC1_var), "% variance", sep="")) +
    ylab(paste("PC2: ", unique(pca_df$PC2_var), "% variance", sep="")) +
    theme_light() +
    facet_wrap(~ Region) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 9),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 9),
          strip.background = element_rect(color = '#99A3A4', fill = 'white'),
          panel.background = element_rect(color = '#99A3A4'),
          strip.text = element_text(size = 9, color = 'black')) +
    # scale_color_manual(values = c('#B03A2E', '#2874A6'), breaks = c('RML6', 'NBH')) +
    scale_fill_manual(values = c('#FADBD8', '#E74C3C', '#943126',
                                 '#D6EAF8', '#3498DB', '#1F618D'))
    # geom_vline(xintercept = 0, linetype = 'dashed', linewidth = 0.5, color = 'lightgrey')
    # geom_hline(yintercept = 0, linetype = 'dashed', linewidth = 0.5, color = 'lightgrey')
    # geom_segment(data =loadings_df, mapping = aes(x = 0, xend = PC1, y = 0, yend = PC2), color ='black',
    #              inherit.aes = FALSE)
    #


  if(proportional_axes == TRUE){

    p <- p + coord_fixed()
  }


  if(sample_label == TRUE){
    p <- p + ggrepel::geom_text_repel(aes(label = SampleID), size = 3, max.overlaps = Inf)
  }

  p
}


make_stage_pca <- function(pca_df, proportional_axes = FALSE,
                           sample_label = FALSE, legend.pos = 'none'){

  p <- ggplot(data = pca_df, aes(x =PC1, y = PC2, fill = Treatment)) +
    geom_point(size = 1, shape = 21, color = '#7F8C8D') +
    labs(
    x = paste("PC1: ", unique(pca_df$PC1_var), "% variance", sep=""),
    y = paste("PC2: ", unique(pca_df$PC2_var), "% variance", sep=""),
    fill = NULL) +
    theme_light() +
    facet_wrap(~stage) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = legend.pos,
          legend.direction = 'horizontal',
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 9),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 9),
          strip.background = element_rect(color = '#99A3A4', fill = 'white'),
          panel.background = element_rect(color = '#99A3A4'),
          strip.text = element_text(size = 9, color = 'black')) +
    # scale_color_manual(values = c('#B03A2E', '#2874A6'), breaks = c('RML6', 'NBH')) +
    scale_fill_manual(values = c('#D7DBDD', '#E48017'),
                      breaks = c('NBH', 'RML6'))
  # geom_vline(xintercept = 0, linetype = 'dashed', linewidth = 0.5, color = 'lightgrey')
  # geom_hline(yintercept = 0, linetype = 'dashed', linewidth = 0.5, color = 'lightgrey')
  # geom_segment(data =loadings_df, mapping = aes(x = 0, xend = PC1, y = 0, yend = PC2), color ='black',
  #              inherit.aes = FALSE)
  #


  if(proportional_axes == TRUE){

    p <- p + coord_fixed()
  }


  if(sample_label == TRUE){
    p <- p + ggrepel::geom_text_repel(aes(label = wpi), size = 3, max.overlaps = Inf)
  }

  p



}


make_3d_pca_plot <- function(pca_df){
plot_ly(pca_df, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Treatment_stage,
        colors = c('#FADBD8', '#E74C3C', '#943126',
                   '#D6EAF8', '#3498DB', '#1F618D'),
        sizes = 1) %>%
  add_markers() %>%
  layout(
    scene = list(
      xaxis = list(title = "PC1"),
      yaxis = list(title = "PC2"),
      zaxis = list(title = "PC3")
    ),
    showlegend = TRUE
  )
}


detect_bottom_outliers <- function(df, vec) {
  # Calculate the quartiles and IQR
  Q1 <- quantile(df[, vec], 0.25)
  Q3 <- quantile(df[, vec], 0.75)
  IQR <- Q3 - Q1

  # Calculate lower bound for outliers
  lower_bound <- Q1 - 1.5 * IQR

  outliers <- df %>% filter(get(vec) < lower_bound) %>%
    pull(SampleID)

  return(outliers)
}



make_pca_legends <- function(pca_df, pal, treatment, legend.text = TRUE, legend = 'color'){
  p <- ggplot(pca_df, aes(PC1, PC2, color = wpi, shape = batch)) +
                geom_point(size = 2)

  if(legend == 'shape'){
    p <- p + guides(color = 'none') +
    theme_light() +
    labs(shape = 'Batch') +
    theme(legend.text = element_text(size = 8),
            legend.title = element_text(size = 9))
    return(get_legend(p))
  }else{

  if(length(pal) == 1){
  p <- p + scale_color_brewer(palette = pal)
  }else{
  p <- p + scale_color_manual(values = pal)
  }
  p <- p + labs(color = paste0('Timepoint\n', treatment)) +
    theme_light() +
    guides(shape = 'none')

  if(legend.text == TRUE){
   p <- p + theme(legend.text = element_text(size = 8),
          legend.title = element_text(size = 9))
  }else{
    p <- p + theme(legend.text = element_blank(),
                   legend.title = element_text(size = 9))
  }

  return(get_legend(p))

  }
}


###### DESEQ2

# general function for running DESeq2

run_deseq2 <- function(counts, meta, formula, ctrl_level, main_var){
# check if rownames in metadata match colnames in counts

    if(!identical(rownames(meta), colnames(counts))){
      # reorder rows in metadata
      meta <- meta[colnames(counts), ]
    }

  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = meta,
                                design = formula)
  # keep only genes that have at lest 10 counts
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]

  dds[[main_var]] <- relevel(dds[[main_var]], ref = ctrl_level)

  dds <- DESeq(dds)
  dds
}


# # find pca outliers
# find_vst_outliers <- function(pca_dt, mad_cutoff = 3, group_var){
#
#   cols <- c('SampleID', group_var, 'PC1', 'PC2')
#   pca_dt <- pca_dt %>%
#     dplyr::select(all_of(cols)) %>%
#     rename(Condition = group_var) %>%
#     as.data.table()
#
#   # find outliers
#   outliers_dt = pca_dt %>%
#     .[, med_1   := median(PC1), by = Condition ] %>%
#     .[, med_2   := median(PC2), by = Condition] %>%
#     .[, mad_1   := mad(PC1 - med_1) ] %>%
#     .[, mad_2   := mad(PC2 - med_2) ] %>%
#     .[, outlier := (abs(PC1 - med_1) > mad_cutoff * mad_1) |
#         (abs(PC2 - med_2) > mad_cutoff * mad_2) ]
#
#   return(as.data.frame(outliers_dt) %>% dplyr::select(-PC1, -PC2))
#
#   }





# extract nice results table from DESeqDataSet object

get_nice_df_from_dds <- function(dds, feature_meta,
                                 out_dir, fname_prefix, meta_vars_to_add, overwrite = FALSE){

  file <- file.path(out_dir, paste0(fname_prefix, '_deseq2_results.nice.rds'))

  if(file.exists(file) & overwrite == FALSE){
   readRDS(file)
  }else{
  res  <- as.data.frame(results(dds))
  # res$geneSymbol <- mapIds(AnnDb_object, keys = rownames(res), keytype = 'ENSEMBL', column = 'SYMBOL')
  res <- res %>% rownames_to_column("ensembl_gene_id")

  # get biotype of gene
  # mart <- useMart('ensembl', dataset = mart_dataset)
  # genes_and_biotypes <- getBM(filters = 'ensembl_gene_id',
  #                             attributes = c('ensembl_gene_id', 'gene_biotype'),
  #                             values = res$ensembl_gene_id, mart = mart)

  res <- merge(res, feature_meta, by = 'ensembl_gene_id', all.x = TRUE, all.y = FALSE)

  # label down and upregulated genes

  res <- res %>% mutate(gene_type = case_when(log2FoldChange >= 0.5 & padj <= 0.05 ~ "Upregulated",
                                              log2FoldChange <= -0.5 & padj <= 0.05 ~ "Downregulated",
                                              TRUE ~ "Nonsignificant"))

  if(length(meta_vars_to_add) != 0){

    if(length(meta_vars_to_add) == 1){
      res[, meta_vars_to_add] <- unique(as.data.frame(colData(dds))[, meta_vars_to_add])

    }else{

    meta_to_add <- as.data.frame(colData(dds))[, meta_vars_to_add, drop = FALSE] %>%
      unique() %>% as_tibble() # selected metadata variables can only have one value

    res <- merge(res, meta_to_add, all.x = TRUE)

    }
  }
  # save files
  saveRDS(res, file)
  return(res)
  }
}



# get normalized counts for one gene

gene_norm_counts <- function(geneId, dds){
  data <- plotCounts(dds, gene = geneId, intgroup = c("sex", "Condition"), returnData = TRUE)
  data$geneID <- geneId
  data <- data %>% rownames_to_column("SampleID")
  data
}



# find common genes in multiple deseq2 results

common_genes <- function(deseq2_results){
  Reduce(intersect, sapply(deseq2_results, FUN = function(x){x$ensembl_gene_id}, simplify = FALSE))
}


# universal ORA function
doORA <- function(genes, uni, AnnDbi_obj, out_dir, file_prefix, overwrite = FALSE){

f <- file.path(out_dir, paste0(file_prefix, '_ORA_results.rds'))

if(file.exists(f) & overwrite == FALSE){
  readRDS(f)
}else{
    res <- enrichGO(gene = genes,
                    OrgDb = AnnDbi_obj,
                    keyType = "ENSEMBL",
                    ont = "BP",
                    minGSSize = 10,
                    maxGSSize = 300,
                    universe = uni)

    res <- clusterProfiler::simplify(res)

    res <- res@result %>% filter(p.adjust < 0.05)

    # res <- res %>% filter(p.adjust <= 0.05)
    # res$timestage <- stage
    #
    if(nrow(res) == 0){
      return("No significant ORA results")
    }else{
    res <- res %>%
      arrange(p.adjust)  %>%
      mutate(pcat = case_when(p.adjust < 10^-8 ~ "<10-8",
                              p.adjust >10^-8 & p.adjust <10^-7  ~ "10^-8-10^-7",
                              p.adjust >10^-7 & p.adjust <10^-6  ~ "10^-7-10^-6",
                              p.adjust >10^-6 & p.adjust <10^-5  ~ "10^-6-10^-5",
                              p.adjust >10^-5 & p.adjust <10^-4 ~ "10^-5-10^-4",
                              p.adjust >10^-4 & p.adjust <10^-3 ~ "10^-4-10^-3",
                              p.adjust >10^-3 & p.adjust <0.01 ~ "10^-3-0.01",
                              p.adjust > 0.01 & p.adjust <0.05 ~ "0.01-0.05",
                              TRUE ~ "0.05-1"))
      saveRDS(res, f)
      return(res)

    }
}
}



########## WGCNA

get_input_for_wgcna <- function(counts, sample_info, file.prefix, out_dir, overwrite = FALSE){

  file <- file.path(out_dir, paste0(file.prefix, '_wgcna_input_data.rds'))
  if(file.exists(file) & overwrite == FALSE){
    readRDS(file)
  }else{

  # check if rownames in sample info match column names in counts
  assert("Rownames in metadata don't match column names in count data",
         identical(rownames(sample_info), colnames(counts)))

  # identify outlier genes and samples
  gsg <-  goodSamplesGenes(t(counts))
  # remove genes that are detectd as outliers
  wgcna_counts <- counts[gsg$goodGenes == TRUE,]

  # Normalization
  dds_for_wgcna  <- DESeqDataSetFromMatrix(countData = wgcna_counts,
                                           colData = sample_info,
                                           design = ~ 1) # not spcifying mode
  # kepp genes that at least in 3 samples have a count of 10 or more
  dds_for_wgcna_filtered  <- dds_for_wgcna[rowSums(counts(dds_for_wgcna) >= 10) >= 3, ]

  # perform variance stabilising transformation
  dds_norm <- vst(dds_for_wgcna_filtered)

  # get normalized counts
  norm.counts <- assay(dds_norm) %>% t()

  norm.counts

  saveRDS(norm.counts, file)
  return(norm.counts)

  }

}


pick_soft_threshold <- function(wgcna_input, power_vec, file.prefix, out_dir, overwrite = FALSE){
  f <- file.path(out_dir, paste0(file.prefix, '_sft.rds'))
  if(file.exists(f) & overwrite == FALSE){
    readRDS(f)
  }else{

    sft <- pickSoftThreshold(wgcna_input,
                       powerVector = power_vec,
                       networkType = "signed hybrid")
    saveRDS(sft, f)
    return(sft)
  }

}



make_plot_for_power_selection <- function(sft){

  sft.data <-  sft$fitIndices


  ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
    geom_point() +
    geom_text(nudge_y = 0.1) +
    geom_hline(yintercept = 0.8, color = 'red') +
    labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
    theme_classic()

}


adj_tom_cluster <- function(norm.counts, thr, out_dir, file.prefix, overwrite = FALSE){

  f <- file.path(out_dir, paste0(file.prefix, '_wgcna_modules_labels.rds'))

  if(file.exists(f) & overwrite == FALSE){
    readRDS(f)
  }else{
  soft_thr <- thr

  adj <- adjacency(datExpr = norm.counts, type = "signed hybrid", power = soft_thr)

  # transform the adjacency into Topological Overlap Matrix, and calculate the corresponding dissimilarity
  TOM <- TOMsimilarity(adj, TOMType = "signed")
  dissTOM <- 1-TOM

  # cluster the genes
  geneTree <- hclust(as.dist(dissTOM), method = "average")

  # cut the clustering tree
  minModuleSize = 30

  # module identification using dynamic tree cut
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = minModuleSize)

  dynamicColors = labels2colors(dynamicMods)
  saveRDS(dynamicColors, f)
  return(dynamicColors)

  }


}



merge_similar_modules <- function(norm.counts, colors, n,
                                  file.prefix, out_dir, overwrite = FALSE){

  f <- file.path(out_dir, paste0(file.prefix, '_final_modules_and_MEs.rds'))

  if(file.exists(f) & overwrite == FALSE){
    readRDS(f)
  }else{

  ## calculate module eigengenes
  MElist <- moduleEigengenes(norm.counts, color = colors)
  MEs <- MElist$eigengenes

  if(length(unique(colors)) < n){
  names(colors) <- colnames(norm.counts)
    final <- list("moduleMEs" = MEs,
                "moduleColors" = colors)

  }else{

    # calculate dissimilarity of module eigengenes
    MEDiss = 1-cor(MEs)
    # cluster module eigengenes
    METree = hclust(as.dist(MEDiss), method = "average")


    MEDissThres = 0.3 # height cut of 0.25, corresponds to a correlation of 0.75

    # Call an automatic merging function
    merge = mergeCloseModules(norm.counts, colors, cutHeight = MEDissThres, verbose = 3)
    # The merged module colors
    mergedColors = merge$colors
    # Eigengenes of the new merged modules:
    mergedMEs = merge$newMEs
   names(mergedColors) <- colnames(norm.counts)
    # return new module eigengenes and module colors

    final <- list("moduleMEs" = mergedMEs,
                "moduleColors" = mergedColors)

  }

  saveRDS(final, f)
  return(final)


  }

}


# info_for_organ <- function(organ, Batch){
#
#   organ_info <- all_samples_info %>% filter(Region == organ) %>% dplyr::select(SampleID, Treatment, wpi)
#
#   if(organ == "spleen"){
#     return(organ_info)
#   }else{
#
#   feature_counts_name <- paste0("FeatureCount", organ)
#
#   colnames <- paste(str_to_title(organ), ".Cohort", c(".1", ".2"), sep = "")
#
#
#   batches <- vector("numeric", length = ncol(all_FC_counts_list[[feature_counts_name]]))
#   names(batches) <- str_before_last(colnames(all_FC_counts_list[[feature_counts_name]]), pattern = "_")
#
#   for(i in 1:length(batches)){
#     if(names(batches)[i] %in% blood_muscle_cohorts[, colnames[1]]){
#       batches[i] <- 1
#     }else if(names(batches)[i] %in% blood_muscle_cohorts[, colnames[2]]){
#       batches[i] <- 2
#     }else{
#       print("The sample is not in the main cohort or it is not assigned to a batch!")
#     }
#
#   }
#
#   names(batches) <- paste0(names(batches), "_", organ)
#
#
#   organ_info  <- organ_info[match(names(batches), organ_info$SampleID), ]
#   organ_info$batch <- batches
#
#   organ_batchX_info <- organ_info %>% filter(batch == Batch)
#
#   return(organ_batchX_info)
#
#   }
#
# }


#######################
#
# extract_FC <- function(organ, batch_info){
#
#   fc_name <- paste0("FeatureCount", organ)
#
#   batch_counts <- all_FC_counts_list[[fc_name]][, batch_info$SampleID]
#   batch_counts
#
# }

#########################


#
# run_DESeq2 <- function(info, counts_per_tp){
#
#   timepoints <- unique(info$wpi)
#
#   deseq2_coldata <- sapply(timepoints,
#                            FUN = function(timepoint){
#
#                              cd <-info %>% filter(wpi == timepoint) %>% dplyr::select(SampleID, Treatment) %>%
# column_to_rownames("SampleID")
#                              cd
#                            }, simplify = FALSE)
#
#
#   dds_objects <- mapply(counts = counts_per_tp,
#                         meta = deseq2_coldata,
#                         FUN = function(counts, meta){
#
#
#                           dds <- DESeqDataSetFromMatrix(countData = counts,
#                                                         colData = meta,
#                                                         design = ~ Treatment)
#
#                           keep <- rowSums(counts(dds)) >= 10  # keep genes with at least 10 counts in total
#                           dds <- dds[keep,]
#
#                           dds$Treatment <- relevel(dds$Treatment, ref = "NBH")
#                           dds
#
#                         }, SIMPLIFY = FALSE)
#
#
#   deseq2_results <- mapply(x = dds_objects,
#                            tp = timepoints,
#                            FUN = function(x, tp){
#                              dds <- DESeq(x)
#                              res <- results(dds)
#                              res2 <- as.data.frame(res) %>% mutate(timepoint = tp) %>% rownames_to_column("geneID")
#                              res2$gene_name <- mapIds(org.Mm.eg.db, keys= res2$geneID ,keytype="ENSEMBL", column="SYMBOL")
#                              res2
#                            },
#                            SIMPLIFY = FALSE)
#
# }


################################


# rld_muscle <- sapply(muscle_counts_per_tp,
#                      FUN = rlog,
#                      simplify = FALSE,
#                      USE.NAMES = TRUE,
#                      blind = TRUE)
#
# pca_500 <- function(rld){
#   rv <- rowVars(rld)
#   select <- order(rv, decreasing=TRUE)[seq_len(500)]
#   pca <- prcomp(t(rld[select,]))
#   pca
# }
#
#
# pca_muscle <- sapply(rld_muscle, FUN = pca_500, simplify = FALSE)
#
# make_PCA_plot <- function(timepoint, pca_res){
#
# pca.data <- muscle_batch1_info %>% filter(wpi == timepoint)
# pca.data$PC1 <- pca_res$x[, 1]
# pca.data$PC2 <- pca_res$x[, 2]
# pca.data$PC3 <- pca_res$x[, 3]
# pca.data$PC4 <- pca_res$x[, 4]
#
# # calculate variance explained per principal component
# pca.var <- pca_res$sdev^2
# pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
#
# # make a plot
#
# pca.plot <- ggplot(data=pca.data, aes(x=PC1, y=PC2, color = Treatment)) +
#   geom_point(size = 2) +
#   xlab(paste("PC1: ", pca.var.per[1], "% variance", sep="")) +
#   ylab(paste("PC2:", pca.var.per[2], "% variance", sep="")) +
#   theme_light() +
#   scale_color_manual(values = c("#81BF84", "#E14637"), breaks = c("NBH", "RML6")) +
#   theme(aspect.ratio = 1/1) +
#   geom_text(aes(label=SampleID),vjust=2, size = 3)
#
# return(pca.plot)
#
# }
#
# all_muscle_PCA_plots <- mapply(timepoint = timepoints,
#                                pca_res = pca_muscle,
#                                FUN = make_PCA_plot,
#                                SIMPLIFY = FALSE)
#



# make_timecourse_plot_for_module <- function(module_column){
  #
  #   ggplot(muscle_batch1_info_and_MEs, aes(x = wpi, y = get(module_column), color = Treatment)) +
  #     geom_point(aes(group = Treatment), position = position_dodge(width = 0.75)) +
  #     facet_wrap(~wpi, scales = "free_x", ncol = 8) +
  #     theme_light() +
  #     labs(x = "Time point", y =  "Module eigengene", title = paste(gsub("ME", "", module_column), "module", sep = " "), color
# = NULL )+
  #     theme(panel.grid.major = element_blank(),
  #           strip.text = element_text(size = 14),
  #           panel.grid.minor = element_blank(),
  #           panel.spacing = unit(0, "cm"),
  #           plot.title = element_text(size = 16),
  #           legend.text = element_text(size = 14),
  #           axis.title = element_text(size = 14),
  #           axis.text.x = NULL,
  #           axis.text.y = element_text(size = 12))
  # }



#########################################



##########################################

merge_MEs_and_sample_info <- function(MEs, info){

  info_and_MEs <- cbind(MEs, info)

  info_and_MEs$wpi <- factor(info_and_MEs$wpi, levels = c("4", "8", "12", "14", "16", "18", "20", "term"))
  info_and_MEs <- info_and_MEs %>% mutate(disease_stage_2 = case_when(wpi %in% c("4", "8") ~ "early",
                                                                      wpi %in% c("12", "14", "16") ~ "presymptomatic",
                                                                      TRUE ~ "symptomatic"))

  info_and_MEs

}


########################################


Treatment_ME_corr_signif <- function(module, info_and_MEs){ # module name with and ME prefix e.i. MEorange

  ME <-info_and_MEs[, c(module, "disease_stage_2", "Treatment")] %>% group_by(disease_stage_2, Treatment)

  # # test if ME values are normally distributed in both groups
  # ME_split <- group_split(ME)
  # shapiro_pvals <- sapply(ME_split, FUN = function(x){
  #   name <- paste(unique(x$disease_stage_2), unique(x$Treatment), sep = "_")
  #   res <- shapiro.test(x[,1] %>% flatten_dbl())
  #   p <- res$p.value
  #   names(p) <- name
  #   return(p)
  # },
  # simplify = TRUE)
  #
  # shapiro_pvals_df <- data.frame(disease_stage_2 = str_before_first(names(shapiro_pvals), pattern = "_"),
  #                                Treatment = str_after_first(names(shapiro_pvals), pattern = "_"),
  #                                shapiro.pval = shapiro_pvals)
  #
  #
  # all_dists <- c()
  #
  # for(i in unique(shapiro_pvals_df$disease_stage_2)){
  #  shaps <-  shapiro_pvals_df %>% filter(disease_stage_2 == i) %>% dplyr::select(shapiro.pval) %>% flatten_dbl()
  #  shaps_normal <- shaps > 0.05
  #  if(sum(shaps_normal) == 2){
  #    dist <- "normal"
  #  }else{
  #    dist <- "notNormal"
  #  }
  #   all_dists[i] <- dist
  # }
  #
  # # calculate significnce betwen groups using a t-test if a distribution is normal or a Mann Whitney U if the distribution is
# not normal

  signif <- sapply(X = unique(ME$disease_stage_2), FUN = function(x){
    stage_ME <- ME %>% filter(disease_stage_2 == x)
    stage_ME <- as.data.frame(stage_ME[, c(module, "Treatment")])
    stage_ME$Treatment <- factor(stage_ME$Treatment, levels = c("RML6", "NBH"))
    #
    #   if(y == "normal"){
    #     # test homogeneity of variances
    #
    #     res.ftest <- var.test(get(module) ~ Treatment, data = stage_ME)
    #     ftest.pval <- res.ftest$p.value
    #
    #     if(ftest.pval > 0.05){ # use a classic t-test which assumes equality of the 2 variances
    #
    #     res.ttest <- t.test(get(module) ~ Treatment, data = stage_ME, var.equal = TRUE)
    #     pval <- res.ttest$p.value
    #     stat <- res.ttest$statistic
    #     stat_type <- "t"
    #     }else{
    #     res.ttest <- t.test(get(module) ~ Treatment, data = stage_ME, var.equal = FALSE)
    #     pval <- res.ttest$p.value
    #     stat <- res.ttest$statistic
    #     stat_type <- "t"
    #     }
    #
    #
    #   }else{
    res.wilcox <- wilcox.test(get(module) ~ Treatment, data = stage_ME, exact = TRUE)
    pval <- res.wilcox$p.value
    stat <- res.wilcox$statistic
    stat_type <- "W"


    final_df <- data.frame(stat = stat, pval = pval, stat_type = stat_type, disease_stage_2 = x, module = module)

  },
  simplify  = TRUE)

  signif <- as.data.frame(t(signif)) %>% remove_rownames() %>% mutate(signif_label = case_when(pval > 0.05 ~ "ns",
                                                                                               pval <= 0.05 & pval > 0.01 ~ "*",
                                                                                               pval <= 0.01 & pval > 0.001 ~
"**",
                                                                                               pval <= 0.001 & pval > 0.0001 ~
"***",
                                                                                               TRUE ~ "****"))

  signif$disease_stage_2 <- as.character(signif$disease_stage_2)
  signif$module <- as.character(signif$module)
  signif$module <- gsub("ME", "", signif$module)
  signif$pval <- as.double(signif$pval)
  signif$stat <- as.double(signif$stat)
  signif$stat_type <- as.character(signif$stat_type)
  signif$heatmap_label <- gsub("ns", "", signif$signif_label)

  signif



}


#############################################


make_wgcna_deseq2_info_df <- function(deseq2_results, norm.counts, moduleColors){

  gene_module_assignment <- data.frame(geneID = colnames(norm.counts),
                                       module = moduleColors)

  get_module_gene_pvalues <- function(time){

    de <- deseq2_results %>% filter(timepoint == time) %>% dplyr::select(geneID, pvalue, timepoint, log2FoldChange, gene_name)

    de_m  <- merge(gene_module_assignment, de, by = "geneID", all.x = TRUE) %>%
      drop_na(pvalue)  %>% # remove genes which don't have corresponding p-values
      mutate(log_pvalue = -log10(pvalue))

    de_m

  }

  genes_modules_pvals <- sapply(X = unique(deseq2_results$timepoint),
                                FUN = get_module_gene_pvalues,
                                simplify = FALSE,
                                USE.NAMES = TRUE)

  return(bind_rows(genes_modules_pvals))


}


###############################################


get_moduleMembership <- function(MEs_df, norm.counts, moduleColors){

  gene_module_assignment <- data.frame(geneID = colnames(norm.counts),
                                       module = moduleColors)

  module.membership.measure <- cor(MEs_df, norm.counts, use = 'p')
  module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nrow(norm.counts))

  # order the genes based on module membership for each module

  order_by_module_membership <- function(module_name){

    module_genes <- gene_module_assignment %>% filter(module == gsub("ME",  "", module_name)) %>% dplyr::select(geneID) %>%
flatten_chr()

    all_genes <- module.membership.measure.pvals[module_name, ]
    all_genes <- sort(all_genes)
    all_gene_symbols <- mapIds(org.Mm.eg.db, keys = names(all_genes), keytype = "ENSEMBL", "SYMBOL")

    all_genes_ordered_df <- data.frame(geneID = names(all_genes),
                                       pval = all_genes,
                                       geneSymbol = all_gene_symbols,
                                       module = module_name) %>%
      mutate(is.module.gene = case_when(geneID %in% module_genes ~ TRUE,
                                        TRUE ~ FALSE))

    all_genes_ordered_df

  }

  modules_all_hubs <- sapply(X = rownames(module.membership.measure.pvals),
                             FUN = order_by_module_membership,
                             simplify = FALSE,
                             USE.NAMES = TRUE)

  # merge MM values for all modules

  modules_all_hubs_merged <- bind_rows(modules_all_hubs) %>% filter(is.module.gene == TRUE) %>% dplyr::rename(MM.pval = pval) %>%
dplyr::select(geneID, MM.pval)
  modules_all_hubs_merged
}


#########################################


merge_deseq2_pvals_and_MM <- function(pvals, MMs){

  merge(pvals, MMs, by = "geneID", all.x = TRUE) %>%
    mutate(log_MMpval = -log10(MM.pval))

}


##########################################


combined_MM_pval_rank <- function(df, hub_number){ # hub_number = final number of hub genes

  # check whether genes in the module are predominantly up or down-regulated

  direction_count <- table(df$log2FoldChange > 0)
  if(direction_count[2] > direction_count[1]){
    direction <- "Upregulated"
  }else{
    direction <- "Downregulated"
  }

  direction_options <- c("Upregulated", "Downregulated")

  # get number of genes in module

  module_gene_n <- length(unique(df$geneID))


  # restructure the dataframe so that p-values and L2fc values for each timepoints will be columns
  restructure <- function(time){
    one_tp <- df %>% filter(timepoint == time)

    one_tp <-  one_tp %>%
      replace_na(list(gene_name = "",
                      pvalue = 1,
                      timepoint = unique(one_tp$timepoint),
                      log_pvalue = 0)) %>%
      dplyr::select(-timepoint)


    colnames(one_tp)[c(3,4,6)] <- paste(colnames(one_tp)[c(3,4,6)], time, sep = "_")

    one_tp

  }

  df_restructured_per_tp <- sapply(unique(df$timepoint),
                                   FUN = restructure,
                                   simplify = FALSE,
                                   USE.NAMES = TRUE)

  df_restructured_per_tp_gene_n <- sapply(df_restructured_per_tp, FUN = function(x){length(unique(x$geneID))}, simplify = TRUE)
  names(df_restructured_per_tp_gene_n) <- 1:length(unique(df$timepoint))

  # reorder df_restructured_per_tp so that the df which has all the genes in the modules goes into merging first

  order <- as.numeric(names(sort(df_restructured_per_tp_gene_n, decreasing = TRUE)))

  df_restructured_per_tp <- df_restructured_per_tp[order]


  df_restructured <- Reduce(left_join, df_restructured_per_tp)

  # caluclate sum of all pvalues for all timepoints in individual time stages

  early_timepoints <- c("4", "8")
  presymptomatic_timepoints <- c("12", "14", "16")
  symptomatic_timepoints <- c("18", "20", "term")

  all_timepoints_found <- unique(df$timepoint)

  early_timpoints_found <- early_timepoints[early_timepoints %in% all_timepoints_found]
  presymptomatic_timepoints_found <- presymptomatic_timepoints[presymptomatic_timepoints %in% all_timepoints_found]
  symptomatic_timepoints_found <- symptomatic_timepoints[symptomatic_timepoints %in% all_timepoints_found]

  # sum log_values if the gene has a fold change that matched direction, if not subtract log_pval

  get_timestage_score <- function(gene, tps){

    print(gene)
   pattern <- paste(paste0("_", tps), collapse = "|")
   cols <- colnames(df_restructured)[grepl(pattern, colnames(df_restructured))]
   log2fc_cols <- cols[grepl("log2FoldChange", cols)]
   logPval_cols <- cols[grepl("log_pvalue", cols)]

   df_restructured_fcs <- as.double(df_restructured %>% filter(geneID  == gene) %>% dplyr::select(all_of(log2fc_cols)))
   fcs_type <- sapply(df_restructured_fcs, FUN = function(x){

    if(is.na(x)){
      x <- 0
    }

     if(x == 0){
       return(direction_options[which(direction_options != direction)])
     }else if(x < 0){
       return("Downregulated")
     }else{
       return("Upregulated")
     }
   })

   fcs_direction_match <- fcs_type == direction
   df_restructured_log_pvals <- as.double(df_restructured %>% filter(geneID == gene) %>% dplyr::select(all_of(logPval_cols)))

   if(any(is.na(df_restructured_log_pvals)) == 1){
     df_restructured_log_pvals[which(is.na(df_restructured_log_pvals))] <- 0
   }

   for(i in 1:length(df_restructured_log_pvals)){
     if(fcs_direction_match[i]  == 0){
       df_restructured_log_pvals[i] <- -df_restructured_log_pvals[i]
     }
   }

   final_sum <- sum(df_restructured_log_pvals)
   final_sum

  }


  # df_restructured <- df_restructured %>% filter(geneID %in% final_genes)
  df_restructured$log_pvalue_sum_early <- sapply(unique(df_restructured$geneID),
                                                 FUN = get_timestage_score,
                                                 tps = early_timpoints_found)

  df_restructured$log_pvalue_sum_presymptomatic <- sapply(unique(df_restructured$geneID),
                                                                 FUN = get_timestage_score,
                                                                 tps = presymptomatic_timepoints_found)

  df_restructured$log_pvalue_sum_symptomatic <- sapply(unique(df_restructured$geneID),
                                                       FUN = get_timestage_score,
                                                       tps = symptomatic_timepoints_found)



# get the sum of the 1st 2 time stages and sum of last 2 time stages and multiply them

  # sum1 <- df_restructured$log_pvalue_sum_presymptomatic + df_restructured$log_pvalue_sum_early
  # sum2 <- df_restructured$log_pvalue_sum_presymptomatic + df_restructured$log_pvalue_sum_symptomatic
  # final_significance_score <- sum1*sum2
  # df_restructured$log_final_significance_score <- final_significance_score
  #

  # select all columns that need to  be ranked

  cols_to_rank <- df_restructured %>% dplyr::select(all_of(c(colnames(df_restructured)[grepl("^log_.*",
colnames(df_restructured))])))
  ranked_cols <- as.data.frame(apply(cols_to_rank, MARGIN = 2, FUN = rank, simplify = TRUE))
  colnames(ranked_cols) <- paste("rank", colnames(ranked_cols), sep = "_")

  # sum the ranks for gene significance at individual time stages

  gene_significance_rank_sum <- ranked_cols$rank_log_pvalue_sum_early + ranked_cols$rank_log_pvalue_sum_presymptomatic +
ranked_cols$rank_log_pvalue_sum_symptomatic
  gene_significance_rank_sum_rank <- rank(gene_significance_rank_sum)

  full_df <- cbind(df_restructured, ranked_cols)
  full_df$gene_significance_rank_sum <- gene_significance_rank_sum
  full_df$gene_significance_rank_sum_rank <- gene_significance_rank_sum_rank

  # resolve ties

  # get all tied ranks

  tied_ranks <- gene_significance_rank_sum_rank[duplicated(gene_significance_rank_sum_rank)]

  resolve_ties <- function(tie){
    ties_df <- full_df %>% filter(gene_significance_rank_sum_rank == tie)
    ties_df$rank_log_pvalue_sum_late <- ties_df$rank_log_pvalue_sum_presymptomatic + ties_df$rank_log_pvalue_sum_symptomatic
    if(nrow(ties_df) == 1){

      return(ties_df)
    }else{

    # get original ranks
    ties_n <- nrow(ties_df)

    if(ties_n %% 2 == 0){ # if ties_n is an odd number
      # get 2 closes integers to a tie
      int1 <- floor(tie)
      int2 <- ceiling(tie)
      # get the number of remaining integers in each direction
      remaining_n <- (ties_n - 2)/2
      from <- int1 - remaining_n
      to <- int2 + remaining_n
      original_ranks <- c(from:to)
    }else{
      int1 <- tie
      remaining_n <- (ties_n -1)/2
      from <- int1 - remaining_n
      to <- int1 + remaining_n
      original_ranks <- c(from:to)
    }

    # assign original ranks back to genes based on the sum of log(pvalues) from 2 latest timepoints

    # ties_df$rank_log_pvalue_sum_late <- ties_df$rank_log_pvalue_sum_presymptomatic + ties_df$rank_log_pvalue_sum_symptomatic
    ties_df <- ties_df %>% arrange(desc(rank_log_pvalue_sum_late))
    ties_df$gene_significance_rank_sum_rank <- rev(original_ranks)
    return(ties_df)
    }
  }

 full_df_final <-  bind_rows(sapply(unique(full_df$gene_significance_rank_sum_rank),
         FUN = resolve_ties,
         simplify = FALSE))

 full_df_final <- full_df_final %>% arrange(desc(gene_significance_rank_sum_rank))

 # total_combined_score <- mean(gene_significance_rank_sum_rank, ranked_cols$rank_log_MMpval)
 # total_combined_score <- ranked_cols$rank_log_pvalue_sum_early + ranked_cols$rank_log_pvalue_sum_presymptomatic +
ranked_cols$rank_log_pvalue_sum_symptomatic + ranked_cols$rank_log_MMpval

 # total_combined_score <- ranked_cols$rank_log_final_significance_score + ranked_cols$rank_log_MMpval

  # append individual and combined ranks back to original dataframe
#
#   full_df <- cbind(df_restructured, ranked_cols)
#  full_df$gene_significance_rank_sum <- gene_significance_rank_sum
# full_df$gene_significance_rank_sum_rank <- gene_significance_rank_sum_rank
#   #full_df <- full_df %>% rowwise() %>% mutate(final_combined_score = mean(gene_significance_rank_sum_rank, rank_log_MMpval))
#   #full_df$final_combined_score <- total_combined_score
#   full_df <- full_df %>% arrange(desc(gene_significance_rank_sum_rank))

  hubs <- full_df_final$geneID[1:hub_number]
  hubSymbols <- full_df_final$gene_name[1:hub_number]
  full_df_final <- full_df_final %>% mutate(is.hub = case_when(geneID %in% hubs ~ TRUE,
                                                   TRUE ~ FALSE))

  n_empty <-  nrow(full_df_final) - hub_number
  full_df_final$hubLabel <- c(hubSymbols, rep("", n_empty))

  full_df_final



}


##############################################################


get_igraph_for_wgcna_network <- function(norm.counts, thr){

  soft_thr <- thr

  adj <- adjacency(datExpr = norm.counts, type = "signed hybrid", power = soft_thr)

  # transform the adjacency into Topological Overlap Matrix, and calculate the corresponding dissimilarity
  TOM <- TOMsimilarity(adj, TOMType = "signed")
  dissTOM <- 1-TOM

  diss_TOM_igraph <- graph_from_adjacency_matrix(dissTOM, mode = "undirected", weighted = TRUE)
  diss_TOM_igraph <- diss_TOM_igraph %>% set_vertex_attr(name = "geneID", value = colnames(norm.counts))
  diss_TOM_igraph
}



##########################################################


get_node_and_edgelist_for_module <- function(igraph, vtx_attrs, filenameinfo){
  module_name <- unique(vtx_attrs$module)
  module_size <- nrow(vtx_attrs)

  # make a new column with random ranks for hub genes

  hubs_random_ranks <- sample(c(1:20), size = 20, replace = FALSE)
  ranks_nh <- c(21:module_size)
  nonhubs_random_ranks <- sample(ranks_nh, size = length(ranks_nh), replace = FALSE)

  vtx_attrs$random_ranks <- c(hubs_random_ranks, nonhubs_random_ranks)


  all_in_module <- vtx_attrs$geneID
  module_graph <- subgraph(igraph, V(igraph)[geneID %in% all_in_module])
  module_graph_mst <- mst(module_graph)

  # reoder rows in vtx_atrrs dataframe so that they will match the order of vertices in igraph object

  geneID_igraph <- get.vertex.attribute(module_graph_mst, name = "geneID")

  vtx_attrs <- vtx_attrs[match(geneID_igraph, vtx_attrs$geneID), ]

  nodelist <- vtx_attrs
  rownames(nodelist) <- V(module_graph_mst)
  nodelist <- nodelist %>% rownames_to_column("Id") %>% dplyr::rename(Label = geneID)

  nodelist_file_name <- paste(filenameinfo, module_name, "dissTOMmst_nodelist.csv", sep = "_")

  edgelist <- as.data.frame(as_edgelist(module_graph_mst)) %>% dplyr::rename(Source = V1, Target = V2)

  edgelist_file_name <- paste(filenameinfo, module_name, "dissTOMmst_edgelist.csv", sep = "_")

  write_excel_csv(nodelist, nodelist_file_name)
  write_excel_csv(edgelist, edgelist_file_name)

}


