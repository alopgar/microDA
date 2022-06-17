#' Add color palette based in hue:
#'
#' \code{microDA} function from the \code{funs_ggplot} group (expanding \code{ggplot2} functionality).\cr
#' Builds a color palette based in hues from \code{\link{hcl}} function.
#' @param n Number of colors to include.
#' @details This function uses \code{\link{hcl}} specifying luminosity 65 and chroma 100 but changing the hue, which is generated
#' sequentially from 15 to 375 according to the number of colors desired (n).
#' @keywords color hcl hue
#' @export
#' @examples
#' a <- gg_color_hue(10)
#' # Visualize colors:
#' colortools::pizza(a)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


#' Count distribution plot:
#'
#' \code{microDA} function that classifies reads according to their counts values and generates a count distribution plot for each sample.
#' @param data Matrix with counts of sequencing reads.
#' @keywords distribution counts ggplot
#' @export

count_distrib <- function(data){
  counts_dist <- NULL
  # Fundamental counts values: 0, 1, 2, 3-10, 10-100, 100-1000, >1000
  for(i in seq(ncol(data))){
    zero <- sum(data[,i] == 0)
    one <- sum(data[,i] == 1)
    two <- sum(data[,i] == 2)
    r10 <- sum(data[,i] >= 3 & data[,i] < 10)
    r100 <- sum(data[,i] >= 10 & data[,i] < 100)
    r1000 <- sum(data[,i] >= 100 & data[,i] < 1000)
    more <- sum(data[,i] >= 1000)
    row <- cbind(zero, one, two, r10, r100, r1000, more)
    counts_dist <- rbind(counts_dist, row)
    rownames(counts_dist)[i] <- colnames(data)[i]
    rm(zero, one, two, r10, r100, r1000, more, row)
  }
  counts_dist <- as.data.frame(counts_dist) %>% rownames_to_column() %>% dplyr::arrange(zero, one, two)
  colnames(counts_dist)[1] <- "samples"
  # Data frame containing order index per sample (based on number of zeros):
  counts_dist_sort <- as.data.frame(cbind("order" = as.numeric(row.names(counts_dist)), "samples" = counts_dist$samples), stringsAsFactors = F)
  # Long data format and merge with sample order index:
  counts_dist_ggdat_pre <- reshape2::melt(counts_dist, id.vars = "samples")
  counts_dist_ggdat <- merge(counts_dist_ggdat_pre, counts_dist_sort, by = "samples")
  counts_dist_ggdat$order <- as.numeric(counts_dist_ggdat$order)
  # GGPLOT:
  if (nrow(data) < 500) {syseq <- 200}
  else syseq <- 500
<<<<<<< HEAD
  ggcounts <- ggplot(counts_dist_ggdat) +
    aes(x = reorder(samples, -order), y = value, fill = variable) +
=======
  ggcounts <- ggplot(counts_dist_ggdat) + 
    aes(x = reorder(samples, -order), y = value, fill = variable) + 
>>>>>>> 23f50d155fa2b9fd1fb2ccfe1412524ecec1364f
    geom_bar(stat = "identity", width = 1, position = position_stack(reverse = T)) +
    xlab("Samples") + ylab("Number of features") +
    theme(axis.text.x = element_blank(), axis.title = element_text(size=14, face="bold")) +
    scale_y_continuous(breaks = seq(0, nrow(data), syseq), expand = expansion(mult = c(0, 0.03))) +
<<<<<<< HEAD
    scale_fill_discrete(name = "Counts per feature",
=======
    scale_fill_discrete(name = "Counts per feature", 
>>>>>>> 23f50d155fa2b9fd1fb2ccfe1412524ecec1364f
                        labels = c("Zero", "Singletons", "Doubletons", "3 to 10", "10 to 100", "100 to 1000", "1000 +"))
  return(ggcounts)
}


#' Rarefaction curves:
#'
#' \code{microDA} function from the \code{funs_ggplot} group (expanding \code{ggplot2} functionality).\cr
#' Builds a color palette based in hues from \code{\link{hcl}} function.
#' @param physeq phyloseq class object, from which abundance data are extracted.
#' @param step Step size for sample size in rarefaction curves.
#' @param color Character string. The name of the variable to map to colors in the plot. This can be a sample variable
#' (among the set returned by \code{sample_variables(physeq)}) or taxonomic rank (among the set returned by
#' \code{rank_names(physeq)}). Finally, The color scheme is chosen automatically by \code{\link[ggplot2]{ggplot}}, but it can
#' be modified afterward with an additional layer using \code{\link[ggplot2]{scale_colour_manual}}. Defaults NULL.
#' @param label Character string. The name of the variable to map to text labels on the plot. Similar to color option but for
#' plotting text. Defaults NULL.
#' @param plot Logical, whether the graphic should be plotted or not. Defaults TRUE.
#' @param parallel Whether rarefaction should be parallelized or not (using parallel framework).
#' @param se Logical, whether standard errors should be computed or not. Defaults TRUE.
#' @details This function is adapted from \code{vegan} \code{\link[vegan]{rarecurve}} function. It adds a ggplot for
#' representing rarefaction curves.
#' @keywords rarefaction vegan ggplot
#' @export

ggrare <- function(physeq, step = 10, color = NULL, label = NULL, plot = T, parallel = F, se = T) {
  x <- as(otu_table(physeq), "matrix")
  if (taxa_are_rows(physeq)) { x <- t(x) }

  ## This script is adapted from vegan `rarecurve` function
  tot <- rowSums(x)
  S <- rowSums(x > 0)
  nr <- nrow(x)

  rarefun <- function(i) {
    cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    y <- rarefy(x[i, ,drop = FALSE], n, se = se)
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
    } else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
    }
  }
  if (parallel) {
    out <- mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
  } else {
    out <- lapply(seq_len(nr), rarefun)
  }
  df <- do.call(rbind, out)

  ## Get sample data
  if (!is.null(sample_data(physeq, FALSE))) {
    sdf <- as(sample_data(physeq), "data.frame")
    sdf$Sample <- rownames(sdf)
    data <- merge(df, sdf, by = "Sample")
    labels <- data.frame(x = tot, y = S, Sample = rownames(x))
    labels <- merge(labels, sdf, by = "Sample")
  }

  ## Add, any custom-supplied plot-mapped variables
  if( length(color) > 1 ){
    data$color <- color
    names(data)[names(data)=="color"] <- deparse(substitute(color))
    color <- deparse(substitute(color))
  }
  if( length(label) > 1 ){
    labels$label <- label
    names(labels)[names(labels)=="label"] <- deparse(substitute(label))
    label <- deparse(substitute(label))
  }

  p <- ggplot(data = data, aes_string(x = "Size", y = ".S", group = "Sample", color = color))
  p <- p + labs(x = "Sample Size", y = "Species Richness")
  if (!is.null(label)) {
    p <- p + ggrepel::geom_text_repel(data = labels, aes_string(x = "x", y = "y", label = label,
                                                                color = color), size = 4, hjust = 0)
  }
  p <- p + geom_line()
  if (se) { ## add standard error if available
    p <- p + geom_ribbon(aes_string(ymin = ".S - .se", ymax = ".S + .se", color = NULL, fill = color), alpha = 0.2)
  }
  if (plot) {
    plot(p)
  }
  invisible(p)
}


#' Sunburst from data:
#'
#' \code{microDA} function from the \code{funs_ggplot} group (expanding \code{ggplot2} functionality).\cr
#' Builds a \code{sunburstR} plot from a microbiome data frame.
#' @param sunreads Feature reads data frame.
#' @param sunfeat Feature taxonomy (or other) classification.
#' @param nranks Number of ranks in taxonomic (or other) classification.
#' @param aggregate Whether glomming reads by sunburst classification or not (in case any class is repeated). Defaults TRUE.
#' @param plottype Values can be 1 (\code{\link[sunburstR]{sunburst}} function) or 2 (\code{\link[sunburstR]{sund2b}} function).
#' @keywords sunburstR sunburst
#' @export

make_sunburst <- function(sunreads, sunfeat, nranks, aggregate = T, plottype = 1){
  sunclass <- sunfeat
  sunclass <- as.data.frame(apply(sunclass, 2, function(x){gsub("-", "_", x)}))
  sunclass <- unite(sunclass, col=Ranking, seq(nranks), sep = "-")

  sunreads <- sunreads
  sunrd <- apply(sunreads, 1, mean)
  sunra <- apply(sunreads, 2, function(x){x/sum(x)*100})
  sunra <- apply(sunra, 1, mean)

  sunb <- data.frame(sunclass, sunra)
  if(isTRUE(aggregate)){
    sunb <- aggregate(sunra ~ Ranking, sunb, sum)
  }

  if(plottype == 1){sb <- sunburst(sunb, legend = F, count = TRUE, width = "100%", height = 500)}
  else if(plottype == 2){sb <- sund2b(sunb, width = "100%", height = 600)}
  else{stop()}

  return(sb)
}


#' Principal Component Analysis plot with centroids:
#'
#' \code{microDA} function from the \code{funs_ggplot} group (expanding \code{ggplot2} functionality).\cr
#' Builds a PCA plot adding centroids according to a user-defined factor.
#' @param pca PCA object from \code{\link{prcomp}}.
#' @param pdata Phenotypic data frame for factor correspondence.
#' @param fit Whether to add linear fitting (using \code{ggordiplots} \code{\link[ggordiplots]{gg_envfit}} function) of
#' variable to PCA plot or not. Defaults to NULL (no fitting).
#' @param group Factor name for centroid grouping.
#' @param axes Length=2 vector with PCA axes to represent. Defaults to c(1,2).
<<<<<<< HEAD
#' @param centr Define if group centroids might be represented. Defaults to FALSE.
=======
#' @param centr Define if group centroids might be represented. Defaults FALSE.
>>>>>>> 23f50d155fa2b9fd1fb2ccfe1412524ecec1364f
#' @keywords PCA envfit ggplot
#' @export

ggpca <- function(pca, pdata, fit = NULL, group, axes = c(1,2), centr = FALSE){
  uscores <- data.frame(pca$x)
  uscores1 <- merge(pdata, uscores, by = "row.names") %>% column_to_rownames("Row.names")

  comp1 <- uscores1[,ncol(pdata) + axes[1], drop = F]
  comp2 <- uscores1[,ncol(pdata) + axes[2], drop = F]

  selected_pcs <- cbind(comp1, comp2)
<<<<<<< HEAD

=======
  
>>>>>>> 23f50d155fa2b9fd1fb2ccfe1412524ecec1364f
  if(isTRUE(centr)){
    centroids <- aggregate(selected_pcs, by = list(uscores1[,group]), mean)
    names(centroids) <- c(group, "CenC1", "CenC2")
    uscores1 <- merge(rownames_to_column(uscores1), centroids, by = group) %>% column_to_rownames("rowname")
  }
  exp.var <- summary(pca)$importance["Proportion of Variance", axes]

  if(!is.null(fit)){
    set.seed(17)
    p <- gg_envfit(pca, dplyr::select(pdata, fit), groups = pdata[,group], choices = c(1,2), perm = 999,
                   alpha = 0.7, pt.size = 1.5, arrow.col = "blue", plot = F)
    if(isTRUE(centr)){
      p_plot <- p$plot + geom_point(data = centroids, aes(x = CenC1, y = CenC2), size = 5) + labs(color = group)
    } else{
      p_plot <- p$plot
    }
  } else{
<<<<<<< HEAD
    p_plot <- ggplot(uscores1, aes_string(x =colnames(comp1), y = colnames(comp2), col = group)) +
      geom_point(size = 3) +
      xlab(paste0(colnames(comp1), " (", round(exp.var[1] * 100, 2), "%)")) +
      ylab(paste0(colnames(comp2), " (", round(exp.var[2] * 100, 2), "%)"))
    if(isTRUE(centr)){
      p_plot <- p_plot +
        geom_point(data = centroids, aes_string(x = "CenC1", y = "CenC2"), size = 5) +
        geom_segment(aes_string(x = colnames(comp1), y = colnames(comp2), xend = "CenC1", yend = "CenC2"))
    }
=======
      p_plot <- ggplot(uscores1, aes_string(x =colnames(comp1), y = colnames(comp2), col = group)) +
        geom_point(size = 3) +
        xlab(paste0(colnames(comp1), " (", round(exp.var[1] * 100, 2), "%)")) +
        ylab(paste0(colnames(comp2), " (", round(exp.var[2] * 100, 2), "%)"))
      if(isTRUE(centr)){
        p_plot <- p_plot +
          geom_point(data = centroids, aes_string(x = "CenC1", y = "CenC2"), size = 5) +
          geom_segment(aes_string(x = colnames(comp1), y = colnames(comp2), xend = "CenC1", yend = "CenC2"))
      }
>>>>>>> 23f50d155fa2b9fd1fb2ccfe1412524ecec1364f
  }
  return(p_plot)
}


#' Scatterplot of groups means and ellipses:
#'
#' \code{microDA} function from the \code{funs_ggplot} group (expanding \code{baseplot} functionality).\cr
#' Builds an x-y scatterplot of groups means and 95% confidence ellipses classifying the data into groups.\cr
#' Function made by Michael Greenacre. \code{ellipse} package required.
#' @param x The x-variable.
#' @param y The y-variable.
#' @param group The grouping variable.
#' @param wt Set of weights on the cases (operates when ellipse=1).
#' @param varnames Vector of two labels for the axes (default x and y).
#' @param groupnames Vector of labels for the groups (default is 1, 2, etc...).
#' @param groupcols Vector a colours for the groups.
#' @param xlim Possible new limits for the plot's x axis.
#' @param ylim Possible new limits for the plot's y axis.
#' @param lwd Line width for the ellipse (default is 1).
#' @param lty Line type for the ellipse (default is 1).
#' @param add If TRUE, ellipses/intervals are added to existing plot (default=FALSE).
#' @param ellipse Set ellipse<0 for regular data-covering ellipses; ellipse=0 (default) for normal-theory confidence ellipses;
#' ellipse=1 for bootstrap confidence ellipses; ellipse=2 for confidence error bars along ellipse axes (not implemented yet!);
#' ellipse=3 for normal-theory confidence error bars lined up with axes; ellipse=4 for bootstrap confidence error bars along axes.
#' @param shade TRUE for ellipse shading (default=FALSE).
#' @param frac If value = 'proportional', part defining the width of the bars at the edges of confidence intervals
#' (for ellipse = 3 and 4).
#' @param cex Character expansion factor for group names.
#' @param alpha Ellipses confidence level (default=0.95).
#' @param shownames Shows group centroids as a label (default=TRUE).
#' @param showcentr Shows group centroids as a dot (default=FALSE).
#' @keywords PCA ellipse plot
#' @export

CIplot_biv <- function(x, y, group, wt=rep(1/length(x),length(x)), varnames=c("x","y"),
                       groupnames=sort(unique(group)), groupcols=rainbow(length(unique(group))),
                       shownames=TRUE, showcentr=FALSE, xlim=c(NA,NA), ylim=c(NA,NA), lty=1, lwd=1,
                       add=FALSE, alpha=0.95, ellipse=0, shade=FALSE, frac=0.01, cex=1.2){
  # first find min and max of all ellipses
  require(ellipse)
  groups <- sort(unique(group))
  xmin <- mean(x)
  xmax <- mean(x)
  ymin <- mean(y)
  ymax <- mean(y)
  for(k in groups) {
    points <- cbind(x,y)[group==k,]
    npts <- nrow(points)
    if(!is.matrix(points)) npts <- 1
    if(npts>=2) {
      covpoints <- var(points)
      meanpoints <- as.numeric(apply(points, 2, mean))
      rconf  <- sqrt(2 * (npts-1) * qf(alpha, 2, npts-2)/(npts*(npts-2)))
      conf.elip <- ellipse::ellipse(covpoints/npts, centre=meanpoints, level=alpha)
      if(ellipse<0) conf.elip <- ellipse(covpoints, center=meanpoints, level=alpha)
      xmin <- min(xmin, conf.elip[,1])
      xmax <- max(xmax, conf.elip[,1])
      ymin <- min(ymin, conf.elip[,2])
      ymax <- max(ymax, conf.elip[,2])
    }
  }
  par(mar=c(4.2,4.2,1,1), cex.axis=0.8)
  if(is.na(xlim[1])) xlim=c(xmin,xmax)
  if(is.na(ylim[1])) ylim=c(ymin,ymax)
  if(!add) plot(x, y, asp=1, type="n", xlab=varnames[1], ylab=varnames[2], main="", xlim=xlim, ylim=ylim, font.lab=2, cex.lab=1.5)
  groups.col <- groupcols
  index <- 0

  # -----------------------------------------------
  # ellipse negative: regular data covering regions

  if(ellipse<0) {
    for(k in groups) {
      index <- index+1
      points <- cbind(x,y)[group==k,]
      npts <- nrow(points)
      if(!is.matrix(points)) npts <- 1
      if(npts>=2) {
        covpoints <- var(points)
        meanpoints <- c(mean(points[,1]), mean(points[,2]))
        rconf  <- sqrt(2 * (npts-1) * qf(alpha, 2, npts-2)/(npts*(npts-2)))
        conf.elip <- ellipse::ellipse(covpoints, centre=meanpoints, level=alpha)
        if(!shade) lines(conf.elip, col=groups.col[index], lty=lty, lwd=lwd)
        if(shade) polygon(conf.elip, col=adjustcolor(groups.col[index], alpha.f=0.2), border=NA)
        if(shownames) text(meanpoints[1], meanpoints[2], labels=groupnames[index], col=groups.col[index], font=2, cex=cex)
        if(showcentr) points(meanpoints[1], meanpoints[2], col="white", bg=groups.col[index], pch=23, cex=cex+0.5)
      }
    }
  }
  # --------------------------------------------
  # ellipse = 0 normal theory confidence regions

  if(ellipse==0) {
    for(k in groups) {
      index <- index+1
      points <- cbind(x,y)[group==k,]
      npts <- nrow(points)
      if(!is.matrix(points)) npts <- 1
      if(npts>=2) {
        covpoints <- var(points)
        meanpoints <- c(mean(points[,1]), mean(points[,2]))
        rconf  <- sqrt(2 * (npts-1) * qf(alpha, 2, npts-2)/(npts*(npts-2)))
        conf.elip <- ellipse::ellipse(covpoints/npts, centre=meanpoints, level=alpha)
        if(!shade) lines(conf.elip, col=groups.col[index], lty=lty, lwd=lwd)
        if(shade) polygon(conf.elip, col=adjustcolor(groups.col[index], alpha.f=0.2), border=NA)
        if(shownames) text(meanpoints[1], meanpoints[2], labels=groupnames[index], col=groups.col[index], font=2, cex=cex)
        if(showcentr) points(meanpoints[1], meanpoints[2], col="white", bg=groups.col[index], pch=23, cex=cex+0.5)
      }
    }
  }

  # -----------------------------------------
  # ellipse = 1: bootstrap confidence regions

  if(ellipse==1) {
    for(k in groups) {
      index <- index+1
      points <- cbind(x,y)[group==k,]
      points.wt <- wt[group==k]
      npts <- nrow(points)
      if(!is.matrix(points)) npts <- 1
      if(npts>=2) {
        boots <- matrix(0, nrow=1000, ncol=2)
        for(iboot in 1:1000) {
          perm <- sample(1:npts, replace=T)
          foo <- points[perm,]
          foo.wt <- points.wt[perm]
          boots[iboot,] <- apply(foo * foo.wt, 2, sum) / sum(foo.wt)
        }
        covpoints <- var(boots)
        meanpoints <- c(sum(points[,1]*points.wt)/sum(points.wt),
                        sum(points[,2]*points.wt)/sum(points.wt))
        print(meanpoints)
        conf.elip <- ellipse::ellipse(covpoints, centre=meanpoints, level=alpha)
        if(!shade) lines(conf.elip, col=groups.col[index], lty=lty, lwd=lwd)
        if(shade) polygon(conf.elip, col=adjustcolor(groups.col[index], alpha.f=0.2), border=NA)
        if(shownames) text(meanpoints[1], meanpoints[2], labels=groupnames[index], col=groups.col[index], font=2, cex=cex)
        if(showcentr) points(meanpoints[1], meanpoints[2], col="white", bg=groups.col[index], pch=23, cex=cex+0.5)
      }
    }
  }

  # --------------------------------------------------
  # ellipse = 2: delta method, non-operative at moment

  # -----------------------------------------------------------------------
  # ellipse = 3: confidence intervals for separate variables, normal theory

  if(ellipse==3) {
    xrange <- xlim[2]-xlim[1]
    yrange <- ylim[2]-ylim[1]
    for(k in groups) {
      index <- index+1
      points <- cbind(x,y)[group==k,]
      npts <- nrow(points)
      if(!is.matrix(points)) npts <- 1
      if(npts>=2) {
        sdpoints <- apply(points, 2, sd)
        meanpoints <- apply(points, 2, mean)
        lines(c(meanpoints[1]-qt(alpha, npts-1)*sdpoints[1]/sqrt(npts),meanpoints[1]+qt(alpha, npts-1)*sdpoints[1]/sqrt(npts)), c(meanpoints[2],meanpoints[2]), col="gray", lwd=lwd, lty=lty)
        lines(c(meanpoints[1],meanpoints[1]), c(meanpoints[2]-qt(alpha, npts-1)*sdpoints[2]/sqrt(npts),meanpoints[2]+qt(alpha, npts-1)*sdpoints[2]/sqrt(npts)), col="gray", lwd=lwd, lty=lty)
        lines(c(meanpoints[1]-qt(alpha, npts-1)*sdpoints[1]/sqrt(npts),meanpoints[1]-qt(alpha, npts-1)*sdpoints[1]/sqrt(npts)), c(meanpoints[2]-frac*yrange,meanpoints[2]+frac*yrange), col="gray", lwd=lwd, lty=lty)
        lines(c(meanpoints[1]+qt(alpha, npts-1)*sdpoints[1]/sqrt(npts),meanpoints[1]+qt(alpha, npts-1)*sdpoints[1]/sqrt(npts)), c(meanpoints[2]-frac*yrange,meanpoints[2]+frac*yrange), col="gray", lwd=lwd, lty=lty)
        lines(c(meanpoints[1]-frac*xrange,meanpoints[1]+frac*xrange), c(meanpoints[2]-qt(alpha, npts-1)*sdpoints[2]/sqrt(npts),meanpoints[2]-qt(alpha, npts-1)*sdpoints[2]/sqrt(npts)), col="gray", lwd=lwd, lty=lty)
        lines(c(meanpoints[1]-frac*xrange,meanpoints[1]+frac*xrange), c(meanpoints[2]+qt(alpha, npts-1)*sdpoints[2]/sqrt(npts),meanpoints[2]+qt(alpha, npts-1)*sdpoints[2]/sqrt(npts)), col="gray", lwd=lwd, lty=lty)
        if(shownames) text(meanpoints[1], meanpoints[2], labels=groupnames[index], col=groups.col[index], font=2, cex=cex)
        if(showcentr) points(meanpoints[1], meanpoints[2], col="white", bg=groups.col[index], pch=23, cex=cex+0.5)
      }
    }
  }

  # --------------------------------------------------------------------------------
  # ellipse = 4: bootstrap confidence intervals for separate variables, by bootstrap

  if(ellipse==4) {
    xrange <- xlim[2]-xlim[1]
    yrange <- ylim[2]-ylim[1]
    for(k in groups) {
      index <- index+1
      points <- cbind(x,y)[group==k,]
      npts <- nrow(points)
      if(!is.matrix(points)) npts <- 1
      if(npts>=2) {
        boots <- matrix(0, nrow=1000, ncol=2)
        for(iboot in 1:1000) boots[iboot,] <- apply(points[sample(1:npts, replace=T),],2,mean)
        xquant <- quantile(boots[,1], c((1-alpha)/2,(1+alpha)/2))
        yquant <- quantile(boots[,2], c((1-alpha)/2,(1+alpha)/2))
        meanpoints <- apply(points, 2, mean)
        lines(c(xquant[1],xquant[2]), c(meanpoints[2],meanpoints[2]), col="gray", lwd=lwd, lty=lty)
        lines(c(meanpoints[1],meanpoints[1]), c(yquant[1], yquant[2]), col="gray", lwd=lwd, lty=lty)
        lines(c(xquant[1],xquant[1]), c(meanpoints[2]-frac*yrange,meanpoints[2]+frac*yrange), col="gray", lwd=lwd, lty=lty)
        lines(c(xquant[2],xquant[2]), c(meanpoints[2]-frac*yrange,meanpoints[2]+frac*yrange), col="gray", lwd=lwd, lty=lty)
        lines(c(meanpoints[1]-frac*xrange,meanpoints[1]+frac*xrange), c(yquant[1],yquant[1]), col="gray", lwd=lwd, lty=lty)
        lines(c(meanpoints[1]-frac*xrange,meanpoints[1]+frac*xrange), c(yquant[2],yquant[2]), col="gray", lwd=lwd, lty=lty)
        if(shownames) text(meanpoints[1], meanpoints[2], labels=groupnames[index], col=groups.col[index], font=2, cex=cex)
        if(showcentr) points(meanpoints[1], meanpoints[2], col="white", bg=groups.col[index], pch=23, cex=cex+0.5)
      }
    }
  }
}


#' PCA plot with ellipses:
#'
#' \code{microDA} function from the \code{funs_ggplot} group (expanding \code{baseplot} functionality).\cr
#' Builds a PCA plot adding centroids and ellipses using \code{\link[microDA]{CIplot_biv}} function.
#' @param pca PCA object from \code{\link{prcomp}}.
#' @param phenot Phenotypic data frame for factor correspondence.
#' @param factor Name of grouping variable.
#' @param fcolors Vector with color codes for groups.
#' @param comps Length=2 vector with PCA components to represent. Defaults to c(1,2).
#' @param title Whether to add plot title or not. Defaults TRUE.
#' @param level If title = TRUE, add a word with taxonomic (or other) level done by PCA. (must be changed to generalise).
#' @keywords PCA ellipse plot
#' @export

Ellipse_PCA <- function(pca, phenot, factor, fcolors, comps = c(1,2), title = TRUE, level = NULL){
  # For saving file: Change plot type "n" to "p", add pch and col plot line; Remove text() line
  phenot.sort <- phenot[rownames(pca$x),]
  classes <- phenot.sort[,factor]
  pca.coords <- pca$x[,comps]
  pinvec <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  #source("C:/Users/Adrian/Dropbox/Trabajo/0. INIA/5_R_Scripts/CIplot_biv.R")
  par(mar = c(4.5,4,1.5,1), mgp = c(2,0.7,0), font.lab = 2)
  plot(pca.coords, type = "p", asp = 1, cex.main=1.5,
       xlab = paste0("PCA", comps[1], "(", round(pinvec[comps[1]], 2), ")"),
       ylab = paste0("PCA", comps[2], "(", round(pinvec[comps[2]], 2), ")"),
       pch = 19, cex.lab = 1.5, col = fcolors[as.numeric(as.factor(classes))])
  if(isTRUE(title)){title(paste("PCA for LRA:", level, "level"))}
  abline(v = 0, h = 0, col = "gray", lty = 2)
  #text(pca.coords, labels = rownames(pca.coords), col = fcolors[as.numeric(as.factor(classes))], cex = 0.7)
  MicroDA::CIplot_biv(pca$x[,comps[1]], pca$x[,comps[2]], group = as.factor(classes), groupcols = fcolors, add = T,
                      shade = T, shownames = F, showcentr = T)
}


#' Venn diagram plot:
#'
#' \code{microDA} function from the \code{funs_ggplot} group (expanding \code{VennDiagram} functionality).\cr
#' Builds an expanded Venn Diagram using a list of string vectors (up to 4). Requires \code{VennDiagram} package.
#' @param vlist List of character strings (max list length = 4).
#' @param ... Other arguments from \code{VennDiagram} functions, such as \code{category} for headers of the diagram or \code{fill} for colors.
#' @keywords PCA ellipse plot
#' @export
#' @examples
#' # Create vectors with common and uncommon elements and group them in a named list:
#' a <- c("one", "two", "four", "ten")
#' b <- c("one", "two", "five", "seven", "eight", "ten")
#' c <- c("one", "three", "seven", "nine", "eleven")
#' vlist <- list(avec = a, bvec = b, cvec = c)
#' # Plot Venn diagram:
#' plotVenn(vlist, category = names(vlist), fill = c("red", "blue", "green"))

plotVenn <- function(vlist, ...) {
  # Function: sum every ovlp_sort element with a pattern
  ovlpsum <- function(pattern){
    out <- Reduce("+", lapply(ovlp_sort[grep(pattern, names(ovlp_sort))], length))
    return(out)
  }
  # Overlap calculation
  ovlp <- calculate.overlap(vlist)
  # Venn diagrams according to elements in the list (1 to 4)
  grid.newpage()
  if (length(vlist) == 1) {
    ovlp_sort <- list("A-U"=ovlp[[1]])
    out <- draw.single.venn(ovlpsum("A"), ...)
  }
  if (length(vlist) == 2) { # doesnt work as 3 or 4 diagram, only needs intersection sizes
    ovlp_sort <- list("A-U"=ovlp$a1, "B-U"=ovlp$a2, "A-B"=ovlp$a3)
    out <- draw.pairwise.venn(ovlpsum("A"), ovlpsum("B"), ovlpsum("A.+B"), ...)
  }
  if (length(vlist) == 3) {
    ovlp_sort <- list("A-U"=ovlp$a1, "B-U"=ovlp$a3, "C-U"=ovlp$a7, "A-B"=ovlp$a2, "A-C"=ovlp$a4,
                      "B-C"=ovlp$a6, "A-B-C"=ovlp$a5)
    out <- draw.triple.venn(ovlpsum("A"), ovlpsum("B"), ovlpsum("C"), ovlpsum("A.+B"), ovlpsum("B.+C"),
                            ovlpsum("A.+C"), ovlpsum("A.+B.+C"), ...)
  }
  if (length(vlist) == 4) {
    ovlp_sort <- list("A-U"=ovlp$a9, "B-U"=ovlp$a14, "C-U"=ovlp$a1, "D-U"=ovlp$a3, "A-B"=ovlp$a15,
                      "A-C"=ovlp$a4, "A-D"=ovlp$a10, "B-C"=ovlp$a13, "B-D"=ovlp$a8, "C-D"=ovlp$a2,
                      "A-B-C"=ovlp$a12, "A-B-D"=ovlp$a11, "A-C-D"=ovlp$a5, "B-C-D"=ovlp$a7, "A-B-C-D"=ovlp$a6)
    out <- draw.quad.venn(ovlpsum("A"), ovlpsum("B"), ovlpsum("C"), ovlpsum("D"), ovlpsum("A.+B"),
                          ovlpsum("A.+C"), ovlpsum("A.+D"), ovlpsum("B.+C"), ovlpsum("B.+D"),
                          ovlpsum("C.+D"), ovlpsum("A.+B.+C"), ovlpsum("A.+B.+D"), ovlpsum("A.+C.+D"),
                          ovlpsum("B.+C.+D"), ovlpsum("A.+B.+C.+D"), ...)
  }
  if (!exists("out"))
    out <- "Oops"
  return(out)
}
