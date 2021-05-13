#' Standard Error of the Mean
#'
#' Calculates standard error of the mean (sem) in one vector.
#' @param x Numeric data vector.
#' @keywords sem
#' @export
#' @examples
#' a <- seq(1:10)
#' sem(a)

sem <- function(x) sqrt(var(x)/length(x))

## SEM (no se q es)
#sem <- function(data, x, y){
#  sem <- list()
#  for(l in seq(x)){
#    fsem <- list()
#    for(i in levels(data[,x[l]])){
#      lvl <- which(levels(data[,x[l]]) %in% i)
#      dat <- data[data[,x[l]]==i, y]
#      fsem[[lvl]] <- sd(dat)/sqrt(length(dat))
#      names(fsem)[lvl] <- i
#    }
#    sem[[l]] <- fsem
#    names(sem)[l] <- x[l]
#  }
#  return(sem)
#}


#' Normality and homocedasticity test:
#'
#' Performs a normality and homocedasticity test for a dataset accounting each level of one factor. It also performs a
#' statistical test according to the normality obtained results. The \strong{\code{nortest}} package is required.
#' @param data Data frame with both data and factor as columns, and samples as rows.
#' @param factor String with the factor name.
#' @param vars String vector with names numeric variables to measure.
#' @details This function selects different test according to dataset properties. Being \strong{n} the number of samples and
#' \strong{k} the levels of the factor:\cr
#' @details Normality test: H0 = Normal distribution. With p < 0.05, H0 is rejected.\cr
#' \verb{   }n <= 50 -> \code{\link[nortest]{lillie.test}}\cr
#' \verb{   }n > 50 -> \code{\link{shapiro.test}}\cr
#' @details Variance equality test: H0 = Equal variances. With p < 0.05, H0 is rejected.\cr
#' \verb{   }k = 2 -> \code{\link{var.test}}\cr
#' \verb{   }k > 2 -> \code{\link{bartlett.test}}\cr
#' @details Statistical test:\cr
#' \verb{   }Normal distribution & k = 2: \code{\link{t.test}}\cr
#' \verb{   }Normal distribution & k > 2: \code{\link{aov}}\cr
#' \verb{   }Not normal distribution & k = 2: \code{\link{wilcox.test}}\cr
#' \verb{   }Not normal distribution & k > 2: \code{\link{kruskal.test}}
#' @keywords normality homocedasticity test
#' @export

Normtest_oneway <- function(data, factor, vars){
  ## Number of levels in factor
  k = length(levels(data[,factor]))

  ## Normality test: H0 = Normal distribution. With p < 0.05, H0 is rejected.
  nort <- NULL
  for (i in sort(levels(data[,factor]))){
    if(nrow(data) <= 50){
      nort[[i]] <- apply(data[data[,factor]==i,vars], 2, lillie.test)}
    else if(nrow(data) > 50){
      nort[[i]] <- apply(data[data[,factor]==i,vars], 2, shapiro.test)}
  }
  ### Normality test table:
  norm <- data.frame(matrix(NA, 0, 3))
  for (i in seq(vars)){
    tn <- data.frame(do.call("rbind", sapply(nort, "[", i))[,-4])
    tn$statistic <- as.double(tn$statistic)
    tn$p.value <- as.double(tn$p.value)
    tn$method <- as.character(tn$method)
    norm <- rbind(norm, tn)
  }
  notnorm <- norm[norm$p.value<=0.05,]
  notnorm.names <- unique(gsub("^.+\\.", "", row.names(notnorm)))
  norm.names <- vars[!(vars %in% notnorm.names)]

  ## Test for Variance equality: H0 = Equal variances. With p < 0.05, H0 is rejected
  if(k == 2){
    varT <- apply(data[,norm.names,drop = F], 2, function(x){var.test(x ~ data[,factor])})
    vartst <- data.frame(do.call("rbind", lapply(varT, '[', c('statistic', 'p.value', 'conf.int'))))
  } else if(k > 2){
    varT <- apply(data[,norm.names,drop = F], 2, function(x){bartlett.test(x ~ data[,factor])})
    vartst <- data.frame(do.call("rbind", lapply(varT, '[', c('statistic', 'p.value'))))
  }
  ### Variance equality test table:
  diffvar <- vartst[vartst$p.value<=0.05,]
  eqvar <- norm.names[!(norm.names %in% rownames(diffvar))]

  ## T-STUDENT/ONE-WAY ANOVA (k=2/>2, NORM)
  t.list <- NULL
  if(k == 2){
    for(i in seq(norm.names)){
      if(varT[[i]]$p.value <= 0.05){veq = F} else {veq = T}
      t.list[[i]] <- t.test(data[,norm.names[i]] ~ data[,factor], var.equal=veq)
      names(t.list)[[i]] <- norm.names[i]
      t.list[[i]]$data.name <- paste(norm.names[i], "by", factor)}
  } else if (k > 2){
    for (i in seq(norm.names)){
      t.list[[i]] <- aov(data[,norm.names[i]] ~ data[,factor])
      names(t.list)[[i]] <- norm.names[i]
      t.list[[i]]$data.name <- norm.names[i]}
  }

  ## MANN-WHITNEY/KRUSKAL-WALLIS (k=2/>2, NOT NORM)
  w.list <- NULL
  if(k==2){
    for (w in seq(notnorm.names)){
      w.list[[w]] <- wilcox.test(data[,notnorm.names[w]] ~ data[,factor])
      names(w.list)[[w]] <- notnorm.names[w]
      w.list[[w]]$data.name <- notnorm.names[w]}
  } else if (k > 2){
    for (w in seq(notnorm.names)){
      w.list[[w]] <- kruskal.test(data[,notnorm.names[w]] ~ data[,factor])
      names(w.list)[[w]] <- notnorm.names[w]
      w.list[[w]]$data.name <- notnorm.names[w]}
  }

  ## OUTPUT:
  out <- list("Norm.test"=nort, "Norm.table"=norm, "Var.test"=varT, "Var.table"=vartst, "Parametric"=t.list,
              "Non_Parametric"=w.list)
  cat("Normality hypothesis (H0) cannot be accepted in variables:", paste(notnorm.names, collapse=", "), "\n")
  cat("Variances are significantly different (H0 rejected) in variables:",
      paste(rownames(diffvar), collapse=", "), "\n")
  invisible(out)
}


#' Relative abundance calculation:
#'
#' Calculates relative abundance (RA) in percentage (sum = 100) from a vector and filters by a minimum value. After filtering,
#' it recalculates RA to sum 100.
#' @param vec Numeric vector to transform in abundances.
#' @param minRA Minimum RA value to filter. Defaults to 0 (no filtering).
#' @keywords abundance filter
#' @export
#' @examples
#' a <- seq(1:20)
#' abundance(a)
#' # check if abundance vector sums 100:
#' sum(abundance(a))

abundance <- function(vec, minRA = 0){
  ab <- vec/sum(vec)*100
  ignore <- which(ab < minRA)
  if (minRA == 0){
    return(ab)
  }
  else {
    ab <- vec/sum(vec[-ignore])*100
    ab[ignore] <- 0
    return(ab)
  }
}


#' relAbundance:
#'
#' Does something related to `abundance` function.
#' @param abundance ???
#' @param taxonomy_tree ???
#' @param level ???
#' @keywords abundance
#' @export

relAbundance <- function(abundance, taxonomy_tree, level){
  tabla <- table(taxonomy_tree[,level])
  #tabla<-tabla[tabla!=0]
  relAb <- vector(length=dim(tabla))
  for (i in 1:length(relAb)){
    relAb[i]<-sum(abundance[ taxonomy_tree[level]==names(tabla)[i] ] )
  }
  out <- list(taxonomy=names(tabla), relative_abundance=relAb)
  return(out)
}


#' Number of zeros:
#'
#' Counts the number of zeros in a vector.
#' @param arg Numeric vector.
#' @keywords zero count
#' @export
#' @examples
#' a <- c(1,1,0,0,1,0)
#' b <- seq(1:10)
#' countZero(a)
#' countZero(b)

countZero <- function(arg){
  return (length(which(arg==0)))
}
