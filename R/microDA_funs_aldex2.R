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

#' Effect size calculation for GLM model:
#'
#' \code{microDA} function from the \code{funs_aldex2} group (expanding \code{ALDEx2} functionality).\cr
#' Calculates effet size when using GLM with two factors in a similar way to \code{\link[ALDEx2]{aldex.effect}}.\cr
#' \code{ALDEx2} package required. THIS IS A TEST!!!
#' @param clr Data frame transformed to CLR.
#' @param metadata Data frame with phenotypic variables including factors.
#' @param factor 2-length vector with factor names.
#' @param levels ??
#' @keywords ALDEx2 glm effect
#' @export

aldex.effect.glm <- function(clr, metadata, factor, levels){
  factordat <- as.factor(metadata[,factor])
  names(factordat) <- rownames(metadata)
  factor_both <- subset(factordat, subset = factordat %in% levels) %>% droplevels()
  factor_l1 <- subset(factordat, subset = factordat %in% levels[1]) %>% droplevels()
  factor_l2 <- subset(factordat, subset = factordat %in% levels[2]) %>% droplevels()

  glm.mci <- getMonteCarloInstances(clr)

  cl2p <- NULL; concatl1 <- NULL; concatl2 <- NULL
  for ( m in glm.mci[names(factor_both)] ) cl2p <- cbind( cl2p, m )
  for ( m in glm.mci[names(factor_l1)] ) concatl1 <- cbind( concatl1, m )
  for ( m in glm.mci[names(factor_l2)] ) concatl2 <- cbind( concatl2, m )

  # Difference within sample groups:
  concat <- list(concatl1, concatl2)
  l2d_win <- list()
  for(con in seq(concat)){
    if ( ncol(concat[[con]]) < 10000 ){sampling.size <- ncol(concat[[con]])} else {sampling.size <- 10000}
    sampl1 <- t(apply(concat[[con]], 1, function(x){sample(x, sampling.size)}))
    sampl2 <- t(apply(concat[[con]], 1, function(x){sample(x, sampling.size)}))
    l2d_win[[con]] <- abs( sampl1 - sampl2 )
  }
  rm(sampl1, sampl2)

  ncol.wanted <- min( sapply( l2d_win, ncol ) )

  l2d_win  <- lapply( l2d_win, function(arg) { arg[,1:ncol.wanted] } )

  # Difference between sample groups:
  sample.size <- min(ncol(concatl1), ncol(concatl2))
  if ( sample.size < 10000 ){sampling.size <- sample.size} else {sampling.size <- 10000}

  smpl1 <- t(apply(concatl1, 1, function(x){sample(x, sampling.size)}))
  smpl2 <- t(apply(concatl2, 1, function(x){sample(x, sampling.size)}))

  rab_all <- apply( cl2p, 1, median )
  rab_win_1 <- apply( concatl1, 1, median )
  rab_win_2 <- apply( concatl2, 1, median )
  l2d_btw <- smpl2 - smpl1

  rm(smpl1, smpl2)

  # Median effect size: diff.btw / max(diff.win)
  win.max <- matrix( 0 , nrow=numFeatures(clr) , ncol=ncol.wanted )
  l2d_effect <- matrix( 0 , nrow=numFeatures(clr) , ncol=ncol(l2d_btw) )
  rownames(l2d_effect) <- getFeatureNames(clr)

  for ( i in 1:numFeatures(clr) ) {
    win.max[i,] <- apply( ( rbind( l2d_win[[1]][i,] , l2d_win[[2]][i,] ) ) , 2 , max )
    l2d_effect[i,] <- l2d_btw[i,] / win.max[i,]
  }

  rownames(win.max)   <- getFeatureNames(clr)
  attr(l2d_win, "max") <- win.max
  rm(win.max)

  dif_win  <- apply( attr(l2d_win,"max"), 1, median )
  dif_btw <- apply( l2d_btw, 1, median )
  effect  <- apply( l2d_effect, 1, median )
  glm.effect <- as.data.frame(cbind(rab_all, rab_win_1, rab_win_2, dif_btw, dif_win, effect))
  colnames(glm.effect) <- c("rab.all", paste0("rab.win.", levels[1]), paste0("rab.win.", levels[2]), "diff.btw",
                            "diff.win", "effect")

  return(list("effect" = glm.effect,
              "conds" = list("both" = factor_both, "level1" = factor_l1, "level2" = factor_l2)))
}


#' Don't remember...
#'
#' \code{microDA} function from the \code{funs_aldex2} group (expanding \code{ALDEx2} functionality).\cr
#' Does something with RA.\cr
#' \code{ALDEx2} package required. THIS IS A TEST!!!
#' @param effect.glm ??
#' @param raw.data ??
#' @param sampletax ??
#' @param closured ??
#' @keywords ALDEx2 abundance
#' @export

aldex.ra.out <- function(effect.glm, raw.data, sampletax, closured = FALSE){
  sampledata <- raw.data
  if(isFALSE(closured)){sampledata <- apply(raw.data, 2, function(x){x/sum(x)*100})}
  level_1 <- gsub("rab_win", "", colnames(effect.glm$effect)[2])
  level_2 <- gsub("rab_win", "", colnames(effect.glm$effect)[3])
  samplediff <- effect.glm$effect[sampletax,]
  samplediff$meanRA <- rowMeans(sampledata[sampletax,aldemeta.dual$IDs])
  samplediff[,paste0("meanRA", level_1)] <- rowMeans(sampledata[sampletax,aldemeta.low$IDs])
  samplediff[,paste0("meanRA", level_2)] <- rowMeans(sampledata[sampletax,aldemeta.high$IDs])
  return(samplediff)
}
