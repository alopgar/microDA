#' DESeq2 expanded results table:
#'
#' \code{microDA} function from the \code{deseq2} group (expanding \code{DESeq2} functionality).\cr
#' Builds an expanded DESeq2 results table filtering by q-value and fold-change significance thresholds, and adding average
#' group abundances. Requires \code{DESeq2} package.
#' @param dsq Deseq object.
#' @param contrast String vector with 3 elements: Name of the factor to separate data, group-1 and group-2.
#' @param siglevel q-value threshold.
#' @param lfcCutoff Fold-change threshold.
#' @param phyloseq Phyloseq object used to build the RA table and the Deseq object.
#' @param taxorder Name of taxonomic level to reorder rows from output table.
#' @param abund.grps Groups used to calculate average abundances. Defaults to \code{contrast[1]}.
#' @param RA Logical. Whether to show abundances in terms of relative abundance (TRUE) or deseq-normalized counts (FALSE). Defaults TRUE.
#' @keywords DESeq2 phyloseq
#' @export

sigtabfun <- function(dsq, contrast, siglevel, lfcCutoff, phyloseq, taxorder, abund.grps = contrast[1], RA = F){
  #dsq=dsq.m; contrast=c("Breed", "IB", "DU"); siglevel=alpha; lfcCutoff=fc.cut; phyloseq=physeq_clr; taxorder="Genus"
  #rm(dsq,contrast,siglevel,lfcCutoff,phyloseq,taxorder)
  # OTU COUNTS/RA + PREVALENCE:
  if(RA){
    otu.ra <- otu_table(transform(phyloseq, "compositional"))@.Data*100
    otu.ra.id <- merge(sample_data(phyloseq)[,abund.grps, drop = F], t(otu.ra), by = "row.names")
    otu.ra <- dplyr::select(otu.ra.id, -Row.names) %>% group_by_(abund.grps) %>% summarise_all(mean) %>%
      column_to_rownames(abund.grps) %>% t() %>% round(3)
  }else{
    #otu.ra <- otu_table(phyloseq)@.Data
    otu.ra <- counts(dsq, normalized=TRUE)
    otu.ra.id <- merge(sample_data(phyloseq)[,abund.grps, drop = F], t(otu.ra), by = "row.names")
    otu.ra <- dplyr::select(otu.ra.id, -Row.names) %>% group_by_(abund.grps) %>% summarise_all(mean) %>%
      column_to_rownames(abund.grps) %>% t() %>% round(3)
  }
  non_zero_counts <- apply(otu_table(phyloseq), 1, function(c)sum(c!=0))
  # DESEQ RESULTS:
  if(length(contrast) == 1 && !is.list(contrast)){dsq.results <- results(dsq, name = contrast, cooksCutoff = Inf)
  }else{dsq.results <- results(dsq, contrast = contrast, cooksCutoff = Inf)}
  sigtab <- dsq.results
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phyloseq)[rownames(sigtab), ], "matrix"))
  # Order by one taxonomy column:
  x = tapply(sigtab$log2FoldChange, sigtab[,taxorder], function(x) max(x))
  x = sort(x, TRUE)
  sigtab[,taxorder] = factor(as.character(sigtab[,taxorder]), levels=names(x))
  # Filtered table (qval + fc):
  sigout <- merge(otu.ra, sigtab[which(sigtab$padj <= siglevel), DAoutcols], by = "row.names")
  sigout <- sigout[abs(sigout$log2FoldChange) >= lfcCutoff,]
  sigout <- merge(non_zero_counts, sigout, by.x = "row.names", by.y = "Row.names") %>%
    rename(Prevalence = x) %>% mutate(Prevalence = round(Prevalence/nsamples(phyloseq), 2)) %>%
    column_to_rownames("Row.names")
  sigout[] <- lapply(sigout, function(x) if(is.factor(x)) as.character(x) else x)
  sigout <- sigout[order(sigout[,taxorder]),]
  # Output:
  end <- list("Full" = sigtab, "Signif" = sigout)
  return(end)
}
