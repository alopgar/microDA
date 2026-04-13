#' Centered log-ratio (CLR) transformation of a compositional dataset:
#'
#' \code{microDA} function from the \code{funs_analysis} group.\cr
#' Includes multiple options to obtain CLR transformation of a matrix, at several steps:\cr
#'   1. Zero replacement: no replacement, cmultRepl imputation, constant value replacement or detection limit-based replacement.\cr
#'   2. CLR weights: equal weights or custom.\cr
#'   3. CLR transformation: easyCODA or compositions packages.\cr
#' \code{easyCODA}, \code{zCompositions} and \code{compositions} packages required.\cr
#' @param data Data frame to transform.
#' @param zerorep Zero replacement method. Can take \strong{"cmultrepl"} (default, bayesian imputation), \strong{"const"} 
#' (replacement by 0.5\*DL), \strong{"unif"} (replacement by a uniform distribution between 0.1\*DL and DL), \strong{0} 
#' (no replacement) or any other number (constant replacement by specified value).
#' @param mode CLR transform method. Can take \strong{"eCODA"} (default, easyCODA::CLR function) and \strong{"compositions"} 
#' (compositions::clr function) values.
#' @param closure Relative abundances will be calculated to sum up to this number per sample. Defaults to 100. 
#' @param transpose Whether the matrix should be transposed or not, since samples must be in rows and features in columns.
#' Defaults TRUE.
#' @param clrweights Custom weights for each feature. Either a vector with weights or FALSE (default, no weights).
#' @keywords CLR easyCODA compositions
#' 
#' @export

CLR_calc <- function(data, zerorep="cmultrepl", mode="eCODA", closure=100, transpose=TRUE, clrweights=FALSE){
  # Transpose the matrix so that samples are in rows & features in columns:
  if(transpose) data = t(data)
  # Impute zeros:
  if(!(0 %in% unlist(data)) || zerorep == 0) {
    codata <- data
  } else if (is.numeric(zerorep)) {
    codata <- data + zerorep
  } else if(zerorep == "cmultrepl"){
    codata <- zCompositions::cmultRepl(data, output="p-counts", z.warning=0.9)
  } else if(zerorep %in% c("const", "unif")){
    codata <- t(apply(data, 1, REPLzeros, zerorep))
  } else stop("Error: zerorep argument must have one of these values: 'cmultrepl', 'const', 'unif', 0 or any number")
  # Compute RA:
  codata.pro <- as.matrix(apply(codata, 2, function(x){closure*x/sum(x)}))
  #codata.pro <- as.matrix(codata) / (rowSums(codata)/closure)
  # CLR weights:
  if(!clrweights){
    clrw = FALSE
  } else if(length(clrweights) > 1){
    clrw = clrweights
  } else if(isTRUE(clrweights)){
    clrw = TRUE
  }
  # CLR transformation:
  if(mode=="eCODA"){
    codata.clr <- easyCODA::CLR(codata.pro, weight = clrw)
  } else if(mode=="compositions"){
    codata.clr <- compositions::clr(codata.pro)
  }
  return(codata.clr)
}

#' Additive log-ratio (ALR) transformation of a compositional dataset:
#'
#' \code{microDA} function from the \code{funs_analysis} group.\cr
#' Includes multiple options to obtain ALR transformation of a matrix, at several steps:\cr
#'   1. Zero replacement: no replacement, cmultRepl imputation or constant value replacement.\cr
#'   2. ALR weights: equal weights or custom.\cr
#'   3. CLR transformation: easyCODA package.\cr
#' \code{easyCODA} and \code{zCompositions} packages required.\cr
#' @param data Data frame to transform.
#' @param reference Name of the feature to be used as reference.
#' @param zerorep Zero replacement method. Can take \strong{"cmultrepl"} (default, bayesian imputation), \strong{0} 
#' (no replacement) or any other number (constant replacement by specified value).
#' @param closure Relative abundances will be calculated to sum up to this number per sample. Defaults to 100. 
#' @param transpose Whether the matrix should be transposed or not, since samples must be in rows and features in columns.
#' Defaults TRUE.
#' @param alrweights Custom weights for each feature. Either a vector with weights or FALSE (default, no weights).
#' @keywords ALR easyCODA
#' @export

ALR_calc <- function(data, reference, zerorep="cmultrepl", closure=100, transpose=TRUE, alrweights=FALSE){
  # Transpose the matrix so that samples are in rows & features in columns:
  if(transpose) data = t(data)
  # Impute zeros:
  if(!(0 %in% unlist(data)) || zerorep == 0) codata <- data
  else if (is.numeric(zerorep)) codata <- data + zerorep
  else if(zerorep == "cmultrepl") codata <- zCompositions::cmultRepl(data, output="p-counts")
  else stop("Error: zerorep argument must have one of these values: 'cmultrepl', 0 or any number")
  # Compute RA:
  codata.pro <- as.matrix(apply(codata, 2, function(x){closure*x/sum(x)}))
  #codata.pro <- as.matrix(codata) / (rowSums(codata)/closure)
  # Get column index for reference:
  refindex <- which(colnames(codata.pro) == reference)
  # ALR transformation:
  if(!alrweights){
    alrw = FALSE
  } else if(length(alrweights) > 1){
    alrw = alrweights
  } else if(isTRUE(alrweights)){
    alrw = TRUE
  }
  codata.alr <- easyCODA::ALR(codata.pro, denom = refindex, weight = alrw)
  return(codata.alr)
}

#' Detection limit-based zero replacement on compositional data:
#'
#' \code{microDA} function from the \code{funs_analysis} group that applies detection limit (DL)-based zero replacement 
#' methods on individual compositions.\cr
#' @param comp Input composition (i.e., one sample with compositional features).
#' @param method Zero replacement method. Can take \strong{"const"} (replacement by 0.5\*DL) and \strong{"unif"} (replacement 
#' by a uniform distribution between 0.1\*DL and DL) values.
#' @details The two methods computed here are described in Lubbe et al (2021) (\url{https://doi.org/10.1016/j.chemolab.2021.104248}).
#' @keywords CLR easyCODA compositions
#' @export

REPLzeros <- function(comp, method){
  # See https://github.com/thomazbastiaanssen/Tjazi/blob/master/R/clr_lite.R
  DL = min(comp[comp != 0])
  if(DL==Inf){DL=0}
  if(method=="unif"){
    comp[comp == 0] <- apply(replicate(1000, runif(n=sum(comp == 0), min=0.1*DL, max=DL)), 1, median)
  } else if(method=="const"){
    comp[comp == 0] <- 0.5*DL
  }
  return(comp)
}
