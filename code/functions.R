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
library(metapod)

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

extract_group_counts <- function(organ, meta, counts){

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

merge_MEs_and_meta <- function(MEs, meta){

  MEs_trans <- MEs %>% rownames_to_column('SampleID') %>%
    merge(., meta, by = 'SampleID')

  return(MEs_trans)
}


Treatment_ME_corr_signif <- function(module, meta){ # module name with and ME prefix e.i. MEorange

  ME <- meta[, c(module, "stage", "Treatment")] %>% group_by(stage, Treatment)

  signif <- sapply(X = unique(ME$stage), FUN = function(x){
    stage_ME <- ME %>% filter(stage == x)
    stage_ME <- as.data.frame(stage_ME[, c(module, "Treatment")])
    stage_ME$Treatment <- factor(stage_ME$Treatment, levels = c("RML6", "NBH"))
    res.wilcox <- wilcox.test(get(module) ~ Treatment, data = stage_ME, exact = TRUE)
    pval <- res.wilcox$p.value
    stat <- res.wilcox$statistic
    stat_type <- "W"


    final_df <- data.frame(stat = stat, pval = pval, stat_type = stat_type, stage = x, module = module)

  },
  simplify  = TRUE)

  signif <- as.data.frame(t(signif)) %>% remove_rownames() %>% mutate(signif_label = case_when(pval > 0.05 ~ "ns",
                                                                                               pval <= 0.05 & pval > 0.01 ~ "*",
                                                                                               pval <= 0.01 & pval > 0.001 ~
                                                                                                 "**",
                                                                                               pval <= 0.001 & pval > 0.0001 ~
                                                                                                 "***",
                                                                                               TRUE ~ "****"))

  signif$stage <- unlist(signif$stage)
  signif$module <- as.character(signif$module)
  signif$module <- gsub("ME", "", signif$module)
  signif$pval <- as.double(signif$pval)
  signif$stat <- as.double(signif$stat)
  signif$stat_type <- as.character(signif$stat_type)
  signif$heatmap_label <- gsub("ns", "", signif$signif_label)

  return(signif)

}


all_modules_ME_signif <- function(meta_w_MEs){ #

  # extract all module names (with ME prefix)
  all_modules <- colnames(meta_w_MEs) %>%
    .[grepl('^ME.*', .)]

  out <- lapply(all_modules, FUN = Treatment_ME_corr_signif, meta = meta_w_MEs) %>%
    bind_rows()
  return(out)

}



# make a heatmap that shows ME significance
make_ME_signif_heatmap <- function(ME_signif_df, legend = 'none'){

   p <-   pvals_heatmap <- ggplot(ME_signif_df %>% filter(module != "grey"),
                          aes(x = stage, y = module, fill = -log10(pval))) + # leave out the grey module with unassigned genes

    geom_tile() +
    theme_minimal() +
    labs(x = NULL, y = NULL, fill = bquote(-Log[10]*(p-value))) +
    theme(panel.grid.major = element_blank(),
          #panel.spacing = unit(0, "mm"),
          #strip.text.y.left = element_text(angle = 0),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = 9),
          axis.text.y = element_text(color = "black", size = 8),
          axis.ticks.y = element_blank()) +
    scico::scale_fill_scico(palette = "devon", direction = -1) +
    geom_text(aes(label = heatmap_label, color = heatmap_label), size = 4, vjust = 0.7) +
    scale_color_manual(values = c("white", "black", "black", "white", "white"), breaks = c("**", "" ,"*", "***", "****"), guide = "none")

  if(legend == 'none'){
    p <- p + theme(legend.position = 'none')
  }else{
    p <- p + theme(legend.position = "right",
                   legend.direction = "vertical",
                   legend.title = element_text(size = 9),
                   axis.text.x = element_text(size = 9),
                   legend.key.height = unit(0.5, "cm"),
                   legend.key.width = unit(0.3, "cm"))
  }

   return(p)


}


# calculate module preservation scores
get_module_pres <- function(b1_input, b2_input, b1_modules, out_dir, file.prefix, overwrite = FALSE){

  f <- file.path(out_dir, paste0(file.prefix, '_module_preservation_stats.rds'))

  if(file.exists(f) & overwrite == FALSE){
    readRDS(f)
  }else{

    setLabels <- c("Batch1", "Batch2")
    multiExpr  = list(Batch1 = list(data = b1_input), Batch2 = list(data = b2_input))
    multiColor = list(Batch1 = b1_modules$moduleColors, Batch2 = b1_modules$moduleColors)

    enableWGCNAThreads(nThreads = 30)

    mp <- modulePreservation(multiExpr,
                             multiColor,
                             referenceNetworks = c(1,2),
                             nPermutations = 1000,
                             randomSeed = 1,
                             parallelCalculation = TRUE)
    saveRDS(mp, f)
    return(mp)

  }

}



get_moduleMembership <- function(module_df, MEs_ls, norm.counts){


  module.membership.measure <- cor(MEs_ls$moduleMEs, norm.counts, use = 'p')
  module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nrow(norm.counts))

  all_module_hubs = vector("list", length = nrow(module.membership.measure.pvals)) %>%
    setNames(rownames(module.membership.measure.pvals))

  # order the genes based on module membership for each module
  for(i in rownames(module.membership.measure.pvals)){
    all_genes = module.membership.measure.pvals[i, ]
    all_genes = sort(all_genes)

    df = data.frame('ensembl_gene_id' = names(all_genes),
                    'MM.pval' = all_genes,
                    'module' = gsub('ME', '', i))
    rownames(df) <- NULL

    all_module_hubs[[i]] = df

  }

  all_module_hubs = all_module_hubs %>%
    bind_rows()

  # merge with module dataframe
  out = merge(module_df, all_module_hubs, by = c('ensembl_gene_id', 'module'), all.x = TRUE, all.y = FALSE) %>%
    arrange(MM.pval)

  return(out)


}



# plot for showing correlation of gene significance score and module membership
make_MM_deseq2_corr_plot <- function(df){

  ggplot(df,
         aes(x = -log10(deseq2_sum.pval),
             y = -log10(MM.pval))) +
    geom_point(size = 0.3) +
    theme_light() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 8),
          plot.title = element_text(size = 8, face = "bold.italic")
    ) +

    geom_smooth(method = "lm",
                color = "#220303",
                inherit.aes = TRUE,
                se = TRUE,
                linewidth = 0.2,
                fill = "#CCD1D1" ) +

    labs(x = 'Gene significance score',
         y = bquote(-Log[10]*.("MM")),
         title = paste(unique(df$module), "module", sep = " ")) +

    ggpubr::stat_cor(p.accuracy = 0.001, size = 3)

}


# plot for showing correlation between b1 and b2 gene significance scores
make_b1_b1_MM_corr_plot <- function(df, x, y, hubx, huby, module_name){

  ggplot() +
  geom_point(data =df,
             shape = 21,
             aes(x = -log10(get(x)),
                 y = -log10(get(y)),
                 color = get(huby)),
             size = 1.5,
             fill = "white")  +
  scale_color_manual(values = c("#FDFEFE", "#A569BD"), breaks = c(FALSE, TRUE)) +

  ggnewscale::new_scale_color() +

  geom_point(data = df,
             aes(x = -log10(get(x)),
                 y = -log10(get(y)),
                 color = get(hubx)),
             size = 0.3) +
  scale_color_manual(values = c("lightgrey", "#220303"), breaks = c(FALSE, TRUE)) +


  theme_light() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 8),
        # legend.position = "bottom",
        legend.position = "none",
        plot.title = element_text(size = 8, face = "bold.italic"),
        legend.direction = "horizontal",
        legend.box = "vertical",
        axis.title = element_text(size = 8)) +
  labs(x = bquote(-Log[10]*.("MM (main cohort)")),
       y =  bquote(-Log[10]*.("MM (validation cohort)")),
       title = module_name
  ) +
  ggpubr::stat_cor(p.accuracy = 0.001, inherit.aes = FALSE,
                   data = df,
                   aes(x = -log10(get(x)),
                       y = -log10(get(y))),
                   size = 3) +

  geom_smooth(method = "lm",
              color = "#220303",
              inherit.aes = FALSE,
              data = df,
              aes(x = -log10(get(x)),
                  y = -log10(get(y))),
              se = TRUE,
              linewidth = 0.2,
              fill = "#CCD1D1" )
}



# boxplots with MEs

make_disease_stage_boxplots_for_modules <- function(Module, MEs_w_meta, ME_signif){


  module_column <- paste0("ME", Module)
  # for each disease stage group determine what is the maximum ME value

  max_y <- MEs_w_meta %>%
    dplyr::select(all_of(module_column), stage) %>%
    group_by(stage) %>%
    summarise(max = max(get(module_column)))

  final_y <- max_y$max + 0.025


  lines <- tibble(stage = c("early", "presymptomatic", "symptomatic"),
                  x = c(1, 1, 1),
                  xend = c(2, 2, 2),
                  y = final_y,
                  yend = y)


  stars <- ME_signif %>% filter(module == Module) %>% dplyr::select(stage, signif_label)
  stars <- stars[match(max_y$stage, stars$stage), ]
  stars$x <- c(1.5, 1.5, 1.5)
  stars$y <- final_y + 0.025
  stars$stage <- as.character(stars$stage)


  ggplot(MEs_w_meta, aes(x = Treatment, y = get(module_column))) +
    # geom_rect(aes(fill = panel_color), xmin = -Inf,xmax = Inf,
    #           ymin = -Inf,ymax = Inf,alpha = 0.1) +
    scale_fill_manual(values = c("white", "#E48016"), breaks = c("middle", "outer"), guide = "none") +
    new_scale_fill() +
    stat_boxplot(geom = "errorbar", linetype = 1, width = 0.4, color = "#616A6B")+
    geom_boxplot(aes(fill = Treatment), notch = FALSE, outlier.size = 0, outlier.colour = "white", color = "#616A6B", width = 0.4) +
    scale_fill_manual(values = c("#D7DBDD", "#E48016"), breaks = c("NBH", "RML6")) +

    geom_jitter(aes(group = Treatment), size = 0.5, color = "black", width = 0.1, height = 0) +
    theme_light() +
    facet_wrap(~stage, scales = "fixed", strip.position = "top", ncol = 1) +
    labs(x = NULL, y = "Module eigengene", title = paste(gsub("ME", "", module_column), "module", sep = " "), fill = NULL) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.length.y = unit(0, "mm"),
          legend.position = "none",

          #legend.direction = "horizontal",
          legend.text = element_text(size = 9),
          strip.text = element_text(colour = "black", size = 8),

          strip.background = element_rect(fill = "white", color = "#B3B3B3"),
          panel.spacing = unit(0, "mm"),
          panel.background = element_rect(fill = "white", color = "#B3B3B3"),
          plot.title = element_text(size = 8, face = "bold.italic"),
          axis.title = element_text(size = 8),
          axis.text.x = element_text(size = 8)) +

    geom_segment(data = lines, aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE, color = "#220303") +

    geom_text(data = stars, aes(x = x, y = y, label = signif_label, angle = -90), inherit.aes = FALSE) +

    coord_flip()


}



#
# make_b1_b1_rank_corr_plot <- function(df, x, y, hubx, huby, module_name){
#
#   ggplot() +
#     geom_point(data =df,
#                shape = 21,
#                aes(x = get(x),
#                    y = get(y),
#                    color = get(huby)),
#                size = 1.5,
#                fill = "white")  +
#     scale_color_manual(values = c("#FDFEFE", "#A569BD"), breaks = c(FALSE, TRUE)) +
#
#     ggnewscale::new_scale_color() +
#
#     geom_point(data = df,
#                aes(x = get(x),
#                    y = get(y),
#                    color = get(hubx)),
#                size = 0.3) +
#     scale_color_manual(values = c("lightgrey", "#220303"), breaks = c(FALSE, TRUE)) +
#
#
#     theme_light() +
#     theme(panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           axis.text = element_text(size = 8),
#           # legend.position = "bottom",
#           legend.position = "none",
#           plot.title = element_text(size = 10, face = "bold.italic"),
#           legend.direction = "horizontal",
#           legend.box = "vertical",
#           axis.title = element_text(size = 9)) +
#     labs(x = "Module gene significance (main cohort)",
#          y = "Module gene significance (validation cohort)",
#          title = module_name
#     ) +
#     ggpubr::stat_cor(p.accuracy = 0.001, inherit.aes = FALSE,
#                      data = df,
#                      aes(x = get(x),
#                          y = get(y)),
#                      size = 3) +
#
#     geom_smooth(method = "lm",
#                 color = "#220303",
#                 inherit.aes = FALSE,
#                 data = df,
#                 aes(x = get(x),
#                     y = get(y)),
#                 se = TRUE,
#                 linewidth = 0.2,
#                 fill = "#CCD1D1" )
# }
#
#


# # get ranking of module genes based on per-stage deseq2 results
# get_combined_pval_rank <- function(deseq2_res){
#
#   # check whether genes in the module are predominantly up or down-regulated
#   dir <- table(deseq2_res$log2FoldChange > 0)
#   if(dir[2] > dir[1]){
#     dir_o <- "Upregulated"
#   }else{
#     dir_o <- "Downregulated"
#   }
#
#   # add fold change direction to pvalues
#   deseq2_res = deseq2_res %>%
#     mutate(log_pval = -log10(deseq2.pval)) %>%
#     mutate(log_pval_dir = case_when(log2FoldChange >= 0 ~ log_pval,
#                                        TRUE ~ log_pval*(-1)))
#
#   # get ranks in all timepoints
#   if(dir_o == 'Upregulated'){
#   deseq2_res_w_ranks = deseq2_res %>%
#     group_by(stage) %>%
#     mutate(rank = rank(log_pval_dir))
#   }else{
#     deseq2_res_w_ranks = deseq2_res %>%
#       group_by(stage) %>%
#       mutate(rank = rank(-log_pval_dir))
#   }
#   # combined rank
#   combn_rank_all_tps = deseq2_res_w_ranks %>%
#     group_by(ensembl_gene_id) %>%
#     summarise(rank_sum_all = sum(rank))
#
#   combn_rank_later_tps = deseq2_res_w_ranks %>%
#     filter(stage != 'early') %>%
#     group_by(ensembl_gene_id) %>%
#     summarise(rank_sum_last = sum(rank))
#
#   combn_rank_merged = merge(combn_rank_all_tps, combn_rank_later_tps, by = 'ensembl_gene_id')
#
#   # resolve ties and assign final ranks
#   combn_rank_merged$final_rank = rank(combn_rank_merged$rank_sum_all)
#
#   # get all tied ranks
#
#   resolve_ties <- function(tie){
#     ties_df <- combn_rank_merged %>%
#       filter(final_rank == tie)
#
#     if(nrow(ties_df) == 1){
#       return(ties_df)
#     }else{
#
#       # get original ranks
#       ties_n <- nrow(ties_df)
#
#       if(ties_n %% 2 == 0){ # if ties_n is an odd number
#         # get 2 closest integers to a tie
#         int1 <- floor(tie)
#         int2 <- ceiling(tie)
#         # get the number of remaining integers in each direction
#         remaining_n <- (ties_n - 2)/2
#         from <- int1 - remaining_n
#         to <- int2 + remaining_n
#         original_ranks <- c(from:to)
#       }else{
#         int1 <- tie
#         remaining_n <- (ties_n -1)/2
#         from <- int1 - remaining_n
#         to <- int1 + remaining_n
#         original_ranks <- c(from:to)
#       }
#
#       # assign original ranks back to genes based on the sum of log(pvalues) from 2 latest timepoints
#
#       # ties_df$rank_log_pvalue_sum_late <- ties_df$rank_log_pvalue_sum_presymptomatic + ties_df$rank_log_pvalue_sum_symptomatic
#       ties_df <- ties_df %>% arrange(desc(rank_sum_last))
#       ties_df$final_rank <- rev(original_ranks)
#       return(ties_df)
#     }
#   }
#
#   full_df_final <-  bind_rows(sapply(unique(combn_rank_merged$final_rank),
#                                      FUN = resolve_ties,
#                                      simplify = FALSE))
#
#   full_df_final <- full_df_final %>% arrange(desc(final_rank))
#
#   return(full_df_final)
#
#
# }
#



# get node and edgelists for plotting modules in gephi


get_igraph_for_wgcna_network <- function(norm.counts, thr){

  soft_thr <- thr

  adj <- adjacency(datExpr = norm.counts, type = "signed hybrid", power = soft_thr)

  # transform the adjacency into Topological Overlap Matrix, and calculate the corresponding dissimilarity
  TOM <- TOMsimilarity(adj, TOMType = "signed")
  dissTOM <- 1-TOM

  diss_TOM_igraph <- graph_from_adjacency_matrix(dissTOM, mode = "undirected", weighted = TRUE)
  diss_TOM_igraph <- diss_TOM_igraph %>% set_vertex_attr(name = "ensembl_gene_id", value = colnames(norm.counts))
  diss_TOM_igraph
}



##########################################################


get_node_and_edgelist_for_module <- function(ig, genes_df, module_name){

  genes_df = genes_df %>% filter(module == module_name) %>%
    mutate(geneID = case_when(!is.na(gene_symbol) ~ gene_symbol,
                              TRUE ~ ensembl_gene_id))

  module_size <- nrow(genes_df)


  all_in_module <- genes_df$ensembl_gene_id
  module_graph <- subgraph(ig, V(ig)[colnames(norm.counts) %in% all_in_module])
  module_graph_mst <- mst(module_graph)

  # reoder rows in vtx_atrrs dataframe so that they will match the order of vertices in igraph object

  geneID_igraph <- get.vertex.attribute(module_graph_mst, name = "ensembl_gene_id")

  genes_df <- genes_df[match(geneID_igraph, genes_df$ensembl_gene_id), ]

  nodelist <- genes_df
  rownames(nodelist) <- V(module_graph_mst)
  nodelist <- nodelist %>% rownames_to_column("Id") %>%
    arrange(MM.pval) %>%
    mutate(
      is.hub = row_number() <= 20
    ) %>%
    mutate(
      Label = case_when(is.hub == TRUE ~ geneID,
                            TRUE ~ "")
    ) %>%
    mutate(hub_rank = case_when(is.hub == TRUE ~ 2,
                                TRUE ~ 1)) %>%
    mutate(log.MM = -log10(MM.pval))


  # add random hub ranks (this is important for plotting in Gephi)
   set.seed(123435)
   hubs_random_ranks <- sample(c((module_size - 19):module_size), size = 20, replace = FALSE)
   ranks_nh <- c(1:c(module_size - 20))
   set.seed(123435)
   nonhubs_random_ranks <- sample(ranks_nh, size = length(ranks_nh), replace = FALSE)


   nodelist$random_ranks <- c(hubs_random_ranks, nonhubs_random_ranks)



  edgelist <- as.data.frame(as_edgelist(module_graph_mst)) %>% dplyr::rename(Source = V1, Target = V2)

  return(list("nodelist" = nodelist,
              "edgelist" = edgelist))

}


get_files_for_gephi <- function(norm.counts, thr, genes_df, module_name, nodelist_f, edgelist_f, overwrite = FALSE){


  if(file.exists(nodelist_f) & file.exists(edgelist_f) & overwrite == FALSE){
    edges = read.csv(edgelist_f)
    nodes = read.csv(nodelist_f)

    return(list("nodelist" = nodes, "edgelist" = edges))
  }else{

    # get graph
    gr = get_igraph_for_wgcna_network(norm.counts, thr = thr)

    # get node and edgelist
    ne = get_node_and_edgelist_for_module(ig = gr, genes_df = genes_df, module_name = module_name)

    write_excel_csv(ne[["nodelist"]], nodelist_f)
    write_excel_csv(ne[["edgelist"]], edgelist_f)

    return(ne)


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




#############################################
#
#
# make_wgcna_deseq2_info_df <- function(deseq2_results, norm.counts, moduleColors){
#
#   gene_module_assignment <- data.frame(geneID = colnames(norm.counts),
#                                        module = moduleColors)
#
#   get_module_gene_pvalues <- function(time){
#
#     de <- deseq2_results %>% filter(timepoint == time) %>% dplyr::select(geneID, pvalue, timepoint, log2FoldChange, gene_name)
#
#     de_m  <- merge(gene_module_assignment, de, by = "geneID", all.x = TRUE) %>%
#       drop_na(pvalue)  %>% # remove genes which don't have corresponding p-values
#       mutate(log_pvalue = -log10(pvalue))
#
#     de_m
#
#   }
#
#   genes_modules_pvals <- sapply(X = unique(deseq2_results$timepoint),
#                                 FUN = get_module_gene_pvalues,
#                                 simplify = FALSE,
#                                 USE.NAMES = TRUE)
#
#   return(bind_rows(genes_modules_pvals))
#
#
# }


###############################################

#
# get_moduleMembership <- function(MEs_df, norm.counts, moduleColors){
#
#   gene_module_assignment <- data.frame(geneID = colnames(norm.counts),
#                                        module = moduleColors)
#
#   module.membership.measure <- cor(MEs_df, norm.counts, use = 'p')
#   module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nrow(norm.counts))
#
#   # order the genes based on module membership for each module
#
#   order_by_module_membership <- function(module_name){
#
#     module_genes <- gene_module_assignment %>% filter(module == gsub("ME",  "", module_name)) %>% dplyr::select(geneID) %>%
# flatten_chr()
#
#     all_genes <- module.membership.measure.pvals[module_name, ]
#     all_genes <- sort(all_genes)
#     all_gene_symbols <- mapIds(org.Mm.eg.db, keys = names(all_genes), keytype = "ENSEMBL", "SYMBOL")
#
#     all_genes_ordered_df <- data.frame(geneID = names(all_genes),
#                                        pval = all_genes,
#                                        geneSymbol = all_gene_symbols,
#                                        module = module_name) %>%
#       mutate(is.module.gene = case_when(geneID %in% module_genes ~ TRUE,
#                                         TRUE ~ FALSE))
#
#     all_genes_ordered_df
#
#   }
#
#   modules_all_hubs <- sapply(X = rownames(module.membership.measure.pvals),
#                              FUN = order_by_module_membership,
#                              simplify = FALSE,
#                              USE.NAMES = TRUE)
#
#   # merge MM values for all modules
#
#   modules_all_hubs_merged <- bind_rows(modules_all_hubs) %>% filter(is.module.gene == TRUE) %>% dplyr::rename(MM.pval = pval) %>%
# dplyr::select(geneID, MM.pval)
#   modules_all_hubs_merged
# }
#
#
# #########################################
#
#
# merge_deseq2_pvals_and_MM <- function(pvals, MMs){
#
#   merge(pvals, MMs, by = "geneID", all.x = TRUE) %>%
#     mutate(log_MMpval = -log10(MM.pval))
#
# }


##########################################

#
# combined_MM_pval_rank <- function(df, hub_number){ # hub_number = final number of hub genes
#
#   # check whether genes in the module are predominantly up or down-regulated
#
#   direction_count <- table(df$log2FoldChange > 0)
#   if(direction_count[2] > direction_count[1]){
#     direction <- "Upregulated"
#   }else{
#     direction <- "Downregulated"
#   }
#
#   direction_options <- c("Upregulated", "Downregulated")
#
#   # get number of genes in module
#
#   module_gene_n <- length(unique(df$geneID))
#
#
#   # restructure the dataframe so that p-values and L2fc values for each timepoints will be columns
#   restructure <- function(time){
#     one_tp <- df %>% filter(timepoint == time)
#
#     one_tp <-  one_tp %>%
#       replace_na(list(gene_name = "",
#                       pvalue = 1,
#                       timepoint = unique(one_tp$timepoint),
#                       log_pvalue = 0)) %>%
#       dplyr::select(-timepoint)
#
#
#     colnames(one_tp)[c(3,4,6)] <- paste(colnames(one_tp)[c(3,4,6)], time, sep = "_")
#
#     one_tp
#
#   }
#
#   df_restructured_per_tp <- sapply(unique(df$timepoint),
#                                    FUN = restructure,
#                                    simplify = FALSE,
#                                    USE.NAMES = TRUE)
#
#   df_restructured_per_tp_gene_n <- sapply(df_restructured_per_tp, FUN = function(x){length(unique(x$geneID))}, simplify = TRUE)
#   names(df_restructured_per_tp_gene_n) <- 1:length(unique(df$timepoint))
#
#   # reorder df_restructured_per_tp so that the df which has all the genes in the modules goes into merging first
#
#   order <- as.numeric(names(sort(df_restructured_per_tp_gene_n, decreasing = TRUE)))
#
#   df_restructured_per_tp <- df_restructured_per_tp[order]
#
#
#   df_restructured <- Reduce(left_join, df_restructured_per_tp)
#
#   # caluclate sum of all pvalues for all timepoints in individual time stages
#
#   early_timepoints <- c("4", "8")
#   presymptomatic_timepoints <- c("12", "14", "16")
#   symptomatic_timepoints <- c("18", "20", "term")
#
#   all_timepoints_found <- unique(df$timepoint)
#
#   early_timpoints_found <- early_timepoints[early_timepoints %in% all_timepoints_found]
#   presymptomatic_timepoints_found <- presymptomatic_timepoints[presymptomatic_timepoints %in% all_timepoints_found]
#   symptomatic_timepoints_found <- symptomatic_timepoints[symptomatic_timepoints %in% all_timepoints_found]
#
#   # sum log_values if the gene has a fold change that matched direction, if not subtract log_pval
#
#   get_timestage_score <- function(gene, tps){
#
#     print(gene)
#    pattern <- paste(paste0("_", tps), collapse = "|")
#    cols <- colnames(df_restructured)[grepl(pattern, colnames(df_restructured))]
#    log2fc_cols <- cols[grepl("log2FoldChange", cols)]
#    logPval_cols <- cols[grepl("log_pvalue", cols)]
#
#    df_restructured_fcs <- as.double(df_restructured %>% filter(geneID  == gene) %>% dplyr::select(all_of(log2fc_cols)))
#    fcs_type <- sapply(df_restructured_fcs, FUN = function(x){
#
#     if(is.na(x)){
#       x <- 0
#     }
#
#      if(x == 0){
#        return(direction_options[which(direction_options != direction)])
#      }else if(x < 0){
#        return("Downregulated")
#      }else{
#        return("Upregulated")
#      }
#    })
#
#    fcs_direction_match <- fcs_type == direction
#    df_restructured_log_pvals <- as.double(df_restructured %>% filter(geneID == gene) %>% dplyr::select(all_of(logPval_cols)))
#
#    if(any(is.na(df_restructured_log_pvals)) == 1){
#      df_restructured_log_pvals[which(is.na(df_restructured_log_pvals))] <- 0
#    }
#
#    for(i in 1:length(df_restructured_log_pvals)){
#      if(fcs_direction_match[i]  == 0){
#        df_restructured_log_pvals[i] <- -df_restructured_log_pvals[i]
#      }
#    }
#
#    final_sum <- sum(df_restructured_log_pvals)
#    final_sum
#
#   }
#
#
#   # df_restructured <- df_restructured %>% filter(geneID %in% final_genes)
#   df_restructured$log_pvalue_sum_early <- sapply(unique(df_restructured$geneID),
#                                                  FUN = get_timestage_score,
#                                                  tps = early_timpoints_found)
#
#   df_restructured$log_pvalue_sum_presymptomatic <- sapply(unique(df_restructured$geneID),
#                                                                  FUN = get_timestage_score,
#                                                                  tps = presymptomatic_timepoints_found)
#
#   df_restructured$log_pvalue_sum_symptomatic <- sapply(unique(df_restructured$geneID),
#                                                        FUN = get_timestage_score,
#                                                        tps = symptomatic_timepoints_found)
#
#
#
# # get the sum of the 1st 2 time stages and sum of last 2 time stages and multiply them
#
#   # sum1 <- df_restructured$log_pvalue_sum_presymptomatic + df_restructured$log_pvalue_sum_early
#   # sum2 <- df_restructured$log_pvalue_sum_presymptomatic + df_restructured$log_pvalue_sum_symptomatic
#   # final_significance_score <- sum1*sum2
#   # df_restructured$log_final_significance_score <- final_significance_score
#   #
#
#   # select all columns that need to  be ranked
#
#   cols_to_rank <- df_restructured %>% dplyr::select(all_of(c(colnames(df_restructured)[grepl("^log_.*",
# colnames(df_restructured))])))
#   ranked_cols <- as.data.frame(apply(cols_to_rank, MARGIN = 2, FUN = rank, simplify = TRUE))
#   colnames(ranked_cols) <- paste("rank", colnames(ranked_cols), sep = "_")
#
#   # sum the ranks for gene significance at individual time stages
#
#   gene_significance_rank_sum <- ranked_cols$rank_log_pvalue_sum_early + ranked_cols$rank_log_pvalue_sum_presymptomatic +
# ranked_cols$rank_log_pvalue_sum_symptomatic
#   gene_significance_rank_sum_rank <- rank(gene_significance_rank_sum)
#
#   full_df <- cbind(df_restructured, ranked_cols)
#   full_df$gene_significance_rank_sum <- gene_significance_rank_sum
#   full_df$gene_significance_rank_sum_rank <- gene_significance_rank_sum_rank
#
#   # resolve ties
#
#   # get all tied ranks
#
#   tied_ranks <- gene_significance_rank_sum_rank[duplicated(gene_significance_rank_sum_rank)]
#
#   resolve_ties <- function(tie){
#     ties_df <- full_df %>% filter(gene_significance_rank_sum_rank == tie)
#     ties_df$rank_log_pvalue_sum_late <- ties_df$rank_log_pvalue_sum_presymptomatic + ties_df$rank_log_pvalue_sum_symptomatic
#     if(nrow(ties_df) == 1){
#
#       return(ties_df)
#     }else{
#
#     # get original ranks
#     ties_n <- nrow(ties_df)
#
#     if(ties_n %% 2 == 0){ # if ties_n is an odd number
#       # get 2 closes integers to a tie
#       int1 <- floor(tie)
#       int2 <- ceiling(tie)
#       # get the number of remaining integers in each direction
#       remaining_n <- (ties_n - 2)/2
#       from <- int1 - remaining_n
#       to <- int2 + remaining_n
#       original_ranks <- c(from:to)
#     }else{
#       int1 <- tie
#       remaining_n <- (ties_n -1)/2
#       from <- int1 - remaining_n
#       to <- int1 + remaining_n
#       original_ranks <- c(from:to)
#     }
#
#     # assign original ranks back to genes based on the sum of log(pvalues) from 2 latest timepoints
#
#     # ties_df$rank_log_pvalue_sum_late <- ties_df$rank_log_pvalue_sum_presymptomatic + ties_df$rank_log_pvalue_sum_symptomatic
#     ties_df <- ties_df %>% arrange(desc(rank_log_pvalue_sum_late))
#     ties_df$gene_significance_rank_sum_rank <- rev(original_ranks)
#     return(ties_df)
#     }
#   }
#
#  full_df_final <-  bind_rows(sapply(unique(full_df$gene_significance_rank_sum_rank),
#          FUN = resolve_ties,
#          simplify = FALSE))
#
#  full_df_final <- full_df_final %>% arrange(desc(gene_significance_rank_sum_rank))
#
#  # total_combined_score <- mean(gene_significance_rank_sum_rank, ranked_cols$rank_log_MMpval)
#  # total_combined_score <- ranked_cols$rank_log_pvalue_sum_early + ranked_cols$rank_log_pvalue_sum_presymptomatic +
# ranked_cols$rank_log_pvalue_sum_symptomatic + ranked_cols$rank_log_MMpval
#
#  # total_combined_score <- ranked_cols$rank_log_final_significance_score + ranked_cols$rank_log_MMpval
#
#   # append individual and combined ranks back to original dataframe
# #
# #   full_df <- cbind(df_restructured, ranked_cols)
# #  full_df$gene_significance_rank_sum <- gene_significance_rank_sum
# # full_df$gene_significance_rank_sum_rank <- gene_significance_rank_sum_rank
# #   #full_df <- full_df %>% rowwise() %>% mutate(final_combined_score = mean(gene_significance_rank_sum_rank, rank_log_MMpval))
# #   #full_df$final_combined_score <- total_combined_score
# #   full_df <- full_df %>% arrange(desc(gene_significance_rank_sum_rank))
#
#   hubs <- full_df_final$geneID[1:hub_number]
#   hubSymbols <- full_df_final$gene_name[1:hub_number]
#   full_df_final <- full_df_final %>% mutate(is.hub = case_when(geneID %in% hubs ~ TRUE,
#                                                    TRUE ~ FALSE))
#
#   n_empty <-  nrow(full_df_final) - hub_number
#   full_df_final$hubLabel <- c(hubSymbols, rep("", n_empty))
#
#   full_df_final
#
#
#
# }


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


