#' Differential abundance plot from ALDEx2:
#'
#' \code{microDA} function from the \code{funs_aldex2} group (expanding \code{ALDEx2} functionality).\cr
#' Builds a differential abundance volcano plot-like graph with ALDEx2 objects.\cr
#' \code{ALDEx2} package required.
#' @param aldex.test Data frame from \code{\link[ALDEx2]{aldex.ttest}} function or similar.
#' @param aldex.effect Output from \code{\link[ALDEx2]{aldex.effect}} function.
#' @param taxonomy Add taxonomy if wanted. Defaults to NULL.
#' @param taxrank Taxonomy rank column corresponding to feature names from aldex.effect.
#' Defaults to "row.names" when a taxonomy data frame is given. If \code{taxonomy = NULL}, it gets NULL value aswell.
#' @param test Used statistical test to make aldex.test data frame. Defaults to "welch".
#' @param cutoff Significance p-value threshold. Defaults to 0.1.
#' @param ra_thr Relative abundance threshold to consider differences as big enough. Defaults to 1.
#' @param layer Adds predefined standard title, subtitle and captions to plot. Defaults to FALSE.
#' @param psize Size of plot points. Defaults to 3.
#' @keywords ALDEx2 volcano abundance
#' @export

aldex.ggplot <- function(aldex.test, aldex.effect, taxonomy = NULL, taxrank = NULL, test = "welch", cutoff = 0.1,
                         ra_thr = 1, layer = F, psize = 3){
  #aldex.test = t.test; aldex.effect = x.effect; taxonomy = aldetax
  #taxrank = "Phylum"; test = "welch"; ra_thr = 1; cutoff = 0.05
  #rm(aldex.test, aldex.effect, taxonomy, taxrank, test, ra_thr, cutoff)
  #rm(effect.tax, aldex.final, exp.P, ggcolor)
  if(!is.null(taxonomy)){
    effect.tax <- merge(aldex.effect, taxonomy, by = "row.names", all.x = T, sort = F) %>%
      column_to_rownames("Row.names") %>% mutate_if(is.factor, as.character)
  } else {
    effect.tax <- aldex.effect
  }

  aldex.final <- data.frame(aldex.test, effect.tax)
  if(test == "welch"){exp.P <- aldex.final$we.eBH}
  else if(test == "wilcox"){exp.P <- aldex.final$wi.eBH}
  else if(test == "glm"){exp.P <- aldex.final$glm.eBH}
  else if(test == "kruskal"){exp.P <- aldex.final$kw.eBH}

  if(is.null(taxrank)){
    aldex.final$color <- "p-val & log2 Diff"
    aldex.final[exp.P <= cutoff & abs(aldex.final$diff.btw) < ra_thr, "color"] <- "p-val"
    aldex.final[exp.P > cutoff & abs(aldex.final$diff.btw) >= ra_thr, "color"] <- "log2 Diff"
    aldex.final[exp.P > cutoff & abs(aldex.final$diff.btw) < ra_thr, "color"] <- "NS"
    aldex.final$color <- factor(aldex.final$color, levels = c("NS", "log2 Diff", "p-val", "p-val & log2 Diff"))
    ggcolor <- c("darkred", "darkgreen", "darkgrey", "black")
    names(ggcolor) <- c("p-val & log2 Diff", "p-val", "log2 Diff", "NS" )
  } else {
    ggcolor <- gg_color_hue(length(unique(aldex.final[,taxrank])))
    names(ggcolor) <- unique(aldex.final[,taxrank])
    aldex.final$color <- aldex.final[,taxrank]
    aldex.final[exp.P <= cutoff & abs(aldex.final$diff.btw) < ra_thr, "color"] <- "p-val"
    aldex.final[exp.P > cutoff & abs(aldex.final$diff.btw) >= ra_thr, "color"] <- "log2 Diff"
    aldex.final[exp.P > cutoff & abs(aldex.final$diff.btw) < ra_thr, "color"] <- "NS"
    aldex.final$color <- factor(aldex.final$color, levels = c("NS", "log2 Diff", "p-val",
                                                              intersect(unique(aldex.final$color), names(ggcolor))))
    ggcolor["p-val"] <- "darkgreen"
    ggcolor["log2 Diff"] <- "darkgrey"
    ggcolor["NS"] <- "black"
  }

  contr.lvls <- str_subset(colnames(aldex.final), "rab.win.") %>% str_replace_all("rab.win.", "")

  aldex.plot <- ggplot(aldex.final, aes(x = diff.btw, y = -log10(exp.P))) +
    geom_hline(yintercept = -log10(cutoff), linetype = "dashed") +
    geom_vline(xintercept = ra_thr, linetype = "dashed") +
    geom_vline(xintercept = -ra_thr, linetype = "dashed") +
    geom_point(aes(color = color), size = psize, alpha = 0.6) + #size = abs(effect),
    scale_color_manual(values=ggcolor) +
    scale_size(guide = FALSE) +
    theme_bw() +
    guides(colour = guide_legend(override.aes = list(size=4))) +
    xlab(expression( ~~ Log[2] ~~ "Fold Change" )) +
    ylab(expression( ~~ -Log[10] ~~ "padj" )) +
    theme(text = element_text(size = 16), axis.title = element_text(size = 14),
          legend.title = element_blank(), legend.position="top")

  layer.labs <- list(labs(title = "ALDEx2 Differential Abundance",
                          subtitle = paste(contr.lvls[1], "vs", contr.lvls[2]),
                          caption = paste("FDR cutoff =", cutoff,
                                          "\nlog2(Fold Change) cutoff = +/-", ra_thr,
                                          "\nTest =", test,
                                          "\ny < 0:", colnames(aldex.final)[6], ">", colnames(aldex.final)[7])))

  if(isTRUE(layer)){aldex.plot <- aldex.plot + layer.labs}

  outlist <- list("fulltab" = dplyr::select(aldex.final, -color),
                  "sigtab" = aldex.final[aldex.final$we.eBH <= cutoff & abs(aldex.final$diff.btw) >= ra_thr,
                                         -which(colnames(aldex.final) == "color")])

  if(!is.null(taxonomy) & !is.null(taxrank)){
    outlist$plot <- aldex.plot + labs(subtitle = paste("Taxonomic level -", taxrank))
  } else {
    outlist$plot <- aldex.plot
  }
  return(outlist)
}
