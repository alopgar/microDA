CLR_calc <- function(data, zeroimp="cmultrepl", mode="eCODA", closure=100, transpose=TRUE, clrweights=FALSE){
  # Transpose the matrix so that samples are in rows & features in columns:
  if(transpose) data = t(data)
  # Impute zeros:
  if(!(0 %in% unlist(data)) || zeroimp == 0) {
    codata <- data
  } else if (is.numeric(zeroimp)) {
    codata <- data + zeroimp
  } else if(zeroimp == "cmultrepl"){
    codata <- zCompositions::cmultRepl(data, output="p-counts", z.warning=0.9)
  } else if(zeroimp %in% c("const", "unif")){
    codata <- t(apply(data, 1, REPLzeros, zeroimp))
  } else stop("Error: zeroimp argument must have one of these values: 'cmultrepl', 0 or any number")
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

ALR_calc <- function(data, reference, zeroimp="cmultrepl", closure=100, transpose=TRUE, alrweights=FALSE){
  # Transpose the matrix so that samples are in rows & features in columns:
  if(transpose) data = t(data)
  # Impute zeros:
  if(!(0 %in% unlist(data)) || zeroimp == 0) codata <- data
  else if (is.numeric(zeroimp)) codata <- data + zeroimp
  else if(zeroimp == "cmultrepl") codata <- zCompositions::cmultRepl(data, output="p-counts")
  else stop("Error: zeroimp argument must have one of these values: 'cmultrepl', 0 or any number")
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

FINDALR <- function(data, weight = FALSE) {
  require(easyCODA); require(vegan)
  ### -------------------------------------------------------------------
  ### function to identify the best reference for a set of ALRs
  ### various statistics are computed for each reference to assist
  ### the choice of best reference data is a normalized data matrix
  ### equal weighting is default for the logratio geometry here
  ### row (sample) weighting not implemented in this version
  
  if(sum(data[1,]!=1)) data <- data/rowSums(data)
  
  ### first compute the exact logratio geometry
  data.lra <- LRA(data, weight=weight)
  data.lra.rpc <- data.lra$rowpcoord
  tot.var <- sum(data.lra$sv^2)
  
  ### loop on all the potential references, computing Procrustes correlation
  ### of each set of ALRs with the exact geometry
  procrust.cor <- rep(0, ncol(data))
  dim <- min(nrow(data), ncol(data)) - 1
  progress <- txtProgressBar(min = 0, max = (ncol(data)), style = 3)
  it_count = 0
  for(j in 1:ncol(data)) {
    # ALR transformation
    alr <- ALR(data, denom=j, weight=weight)
    # ALR geometry using PCA 'by hand' using SVD, without or with weighting
    if(!weight) {
      alr.svd <- svd(sqrt(1/nrow(alr$LR)) * sweep(alr$LR, 2, colMeans(alr$LR)) * sqrt(1/ncol(alr$LR)))
      alr.rpc <- sqrt(nrow(alr$LR)) * alr.svd$u %*% diag(alr.svd$d)
    }
    if(weight) {
      c <- colMeans(data)
      cc <- c*c[j]
      cc <- cc[-j]
      alr.svd <- svd(sqrt(1/nrow(alr$LR)) * sweep(alr$LR, 2, colMeans(alr$LR)) %*%  diag(sqrt(cc)))
      alr.rpc <- sqrt(nrow(alr$LR)) * alr.svd$u %*% diag(alr.svd$d)
    }
    procrust.cor[j] <- protest(alr.rpc[,1:dim],data.lra.rpc, permutations=0)$t0
    it_count <- it_count + 1
    setTxtProgressBar(progress, it_count)
  }
  close(progress)
  
  ### the variances of the log-transformed parts
  var.log <- as.numeric(apply(log(data), 2, var))
  
  ### highest Procrustes correlation
  procrust.max <- max(procrust.cor)
  
  ### which reference gives maximum Procrustes
  procrust.ref <- which(procrust.cor==procrust.max)
  
  ### lowest log variance
  var.min <- min(var.log)
  
  ### which reference gives lowest log variance
  var.ref <- which(var.log==var.min)
  
  return(list(tot.var=tot.var, procrust.cor=procrust.cor, 
              procrust.max=procrust.max, procrust.ref=procrust.ref,
              var.log=var.log, var.min=var.min, var.ref=var.ref))
}

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

chiPower <- function(data, close=TRUE, power=1, chi=TRUE, BoxCox=TRUE, center=TRUE) {
  # See https://github.com/michaelgreenacre/CODAinPractice/blob/master/chiPower_miniscript.R
  # function to compute chiPower transformation with various options
  # close: close data
  # power: transformation (i.e., value of lambda)
  # chi: chi-square standardization
  # BoxCox: apply 1/power rescaling (without the subtraction of 1)
  # center: center the final result by column means
  foo <- data^power
  if(close)   foo <- foo / rowSums(foo)
  if(chi)     foo <- sweep(foo, 2, sqrt(colMeans(foo)), FUN="/")
  if(BoxCox)  foo <- (1/power) * foo
  if(center)  foo <- sweep(foo, 2, colMeans(foo))
  if(!center) foo <- foo - 1/sqrt(ncol(data))
  data.chiPower <- foo
  return(data.chiPower)
}

chiPower_test <- function(data, powers=seq(0.01,1,0.01), closure=100){
  CLR_tab <- CLR_calc(data, zeroimp=0.5, mode="compositions", closure=closure, transpose=F)
  CLR_PCA <- prcomp(CLR_tab, center=F, scale=F)
  CLR_coords <- CLR_PCA$x[,-(ncol(data)-1)]
  CLR_DIST <- dist(CLR_tab)/sqrt(ncol(data))
  
  clr_chip_proc <- rep(0,length(powers))
  clr_chip_corr <- rep(0,length(powers))
  for(pow in powers){
    print(paste("Testing power lambda:", pow))
    ChiP_tab <- chiPower(data, close=F, power=pow)
    ChiP_PCA <- prcomp(ChiP_tab, center=F, scale=F)
    ChiP_coords <- ChiP_PCA$x[,-(ncol(data)-1)]
    ChiP_DIST <- dist(ChiP_tab)
    clr_chip_proc[pow*100] <- protest(CLR_coords, ChiP_coords, permutations=0)$t0
    clr_chip_corr[pow*100] <- cor(CLR_DIST, ChiP_DIST, method="spearman")
  }
  return(list("Procrustes"=clr_chip_proc, "Spearman"=clr_chip_corr))
}

