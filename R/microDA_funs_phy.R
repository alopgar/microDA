#' Get max taxonomic rank:
#'
#' \code{microDA} function from the \code{funs_phy} group (expanding \code{phyloseq} functionality).\cr
#' Determines the lowest level of taxonomic classification, finds the last non-NA column in the taxonomy table and returns it.\cr
#' Function adapted from metagMisc repository (https://github.com/vmikk/metagMisc).
#' @param x Either a phyloseq object, or a data frame with columns as taxonomic ranks and rows as entries (e.g., OTUs). Columns in the data frame should be ordered from the highest level of classification (e.g., Kingdom) to the lowest level (e.g., Species), missing data are coded as NA.
#' @param return_rank_only Logical, if TRUE only name of the taxonomic rank will be returned
#' @keywords phyloseq taxonomy
#' @export

get_max_taxonomic_rank <- function(x, return_rank_only = FALSE){
  ## Input data
  inp_class <- class(x)
  ## If input is of class 'phyloseq'
  if("phyloseq" %in% inp_class || "taxonomyTable" %in% inp_class){

    if(is.null(phyloseq::tax_table(x, errorIfNULL = F))){
      stop("Error: taxonomy table slot is empty in the input data.\n")
    }
    otu_names <- phyloseq::taxa_names(x)
    x <- as.data.frame(phyloseq::tax_table(x), stringsAsFactors = F)
  }
  ## If there are factors in the input data frame this can cause an "Error: cannot allocate vector of size ... Gb"
  ## Convert all factors to character
  if("factor" %in% plyr::laply(.data = x, .fun = class)){
    x[] <- lapply(x, as.character)
  }
  ## Display progress bar?
  if(nrow(x) < 200){
    progr <- "none"
  } else {
    progr <- "text"
  }
  ## Find the indices of the last non-NA column for each row
  res <- plyr::adply(.data = x, .margins = 1, .fun = function(z){
    ## Test which tax ranks are not NAs
    rnk <- plyr::aaply(.data = z, .margins = 1, .fun = function(y) which(!is.na(y)) )
    ## Return last non-NA column number
    if(length(rnk) > 0){
      rez <- max(rnk)
    } else {
      rez <- 0   # if all tax ranks are NA
    }
    rez <- data.frame(RankColumn = rez)
    return(rez)
  }, .progress = progr)
  ## Table with correspondence of taxonomic ranks and their indices
  rnks <- data.frame(
    RankColumn = c(0, 1:ncol(x)),
    RankName = c(NA, colnames(x)),
    stringsAsFactors = F)
  ## Substitute column index with tax rank name
  res$RankName <- rnks$RankName[match(x = res$RankColumn, table = rnks$RankColumn)]
  ## Remove rank column
  res$RankColumn <- NULL
  ## Reorder taxonomic ranks
  res$RankName <- factor(res$RankName, levels = colnames(x))
  ## Add OTU name
  if("phyloseq" %in% inp_class || "taxonomyTable" %in% inp_class){
    res <- data.frame(TaxaName = otu_names, res, stringsAsFactors = F)
  }
  if(return_rank_only == TRUE){
    res <- res$RankName
  }
  return(res)
}


#' Phyloseq to data frame:
#'
#' \code{microDA} function from the \code{funs_phy} group (expanding \code{phyloseq} functionality).\cr
#' Converts a \code{phyloseq} object into a \code{data frame}.\cr
#' Function adapted from metagMisc repository (https://github.com/vmikk/metagMisc).
#' @param physeq \code{phyloseq} input object.
#' @param addtax Whether adding \code{tax_table(physeq)} to new data frame or not. Defaults TRUE.
#' @param addtot Whether adding total OTU abundance as a final column or not. Defaults FALSE.
#' @param addmaxrank Whether adding a max taxonomic rank column or not. Defaults FALSE.
#' @param sorting Reorder OTUs. Can take \strong{"abundance"} (default), \strong{"taxonomy"} or \strong{NULL} values.
#' @keywords phyloseq dataframe
#' @export

phyloseq_to_df <- function(physeq, addtax = T, addtot = F, addmaxrank = F, sorting = "abundance"){
  ## Data validation
  if(any(addtax == TRUE || sorting == "taxonomy")){
    if(is.null(phyloseq::tax_table(physeq, errorIfNULL = F))){
      stop("Error: taxonomy table slot is empty in the input data.\n")
    }
  }

  ## Prepare data frame
  if(taxa_are_rows(physeq) == TRUE){
    res <- data.frame(OTU = phyloseq::taxa_names(physeq), phyloseq::otu_table(physeq), stringsAsFactors = F)
  } else {
    res <- data.frame(OTU = phyloseq::taxa_names(physeq), t(phyloseq::otu_table(physeq)), stringsAsFactors = F)
  }

  ## Check if the sample names were silently corrected in the data.frame
  if(any(!phyloseq::sample_names(physeq) %in% colnames(res)[-1])){
    if(addtax == FALSE){
      warning("Warning: Sample names were converted to the syntactically valid column names
              in data.frame. See 'make.names'.\n")
    }
    if(addtax == TRUE){
      stop("Error: Sample names in 'physeq' could not be automatically converted to the
           syntactically valid column names in data.frame (see 'make.names'). Consider renaming
           with 'sample_names'.\n")
    }
  }

  ## Add taxonomy
  if(addtax == TRUE){
    ## Extract taxonomy table
    taxx <- as.data.frame(phyloseq::tax_table(physeq), stringsAsFactors = F)
    ## Reorder taxonomy table
    taxx <- taxx[match(x = res$OTU, table = rownames(taxx)), ]
    ## Add taxonomy table to the data
    res <- cbind(res, taxx)
    ## Add max tax rank column
    if(addmaxrank == TRUE){
      ## Determine the lowest level of taxonomic classification
      res$LowestTaxRank <- get_max_taxonomic_rank(taxx, return_rank_only = TRUE)
      ## Reorder columns (OTU name - Taxonomy - Max Rank - Sample Abundance)
      res <- res[, c("OTU", phyloseq::rank_names(physeq), "LowestTaxRank", phyloseq::sample_names(physeq))]
    } else {
      ## Reorder columns (OTU name - Taxonomy - Sample Abundance)
      res <- res[, c("OTU", phyloseq::rank_names(physeq), phyloseq::sample_names(physeq))]
    } # end of addmaxrank
  }   # end of addtax

  ## Reorder OTUs
  if(!is.null(sorting)){
    ## Sort by OTU abundance
    if(sorting == "abundance"){
      otus <- res[, which(colnames(res) %in% phyloseq::sample_names(physeq))]
      res <- res[order(rowSums(otus, na.rm = T), decreasing = T), ]
    }

    ## Sort by OTU taxonomy
    if(sorting == "taxonomy"){
      taxtbl <- as.data.frame( phyloseq::tax_table(physeq), stringsAsFactors = F )
      ## Reorder by all columns
      taxtbl <- taxtbl[do.call(order, taxtbl), ]
      # taxtbl <- data.table::setorderv(taxtbl, cols = colnames(taxtbl), na.last = T)
      res <- res[match(x = rownames(taxtbl), table = res$OTU), ]
    }
  }

  ## Add OTU total abundance
  if(addtot == TRUE){
    res$Total <- rowSums(res[, which(colnames(res) %in% phyloseq::sample_names(physeq))])
  }

  rownames(res) <- NULL
  return(res)
}


#' Phyloseq subsetting to core:
#'
#' \code{microDA} function from the \code{funs_phy} group (expanding \code{phyloseq} functionality).\cr
#' Subsets a \code{phyloseq} object given \strong{detection} and \strong{prevalence} thresholds.
#' @param physeq \code{phyloseq} input object.
#' @param transf Transformation to apply previous to the subsetting. See \code{\link[microbiome]{transform}} for more details.
#' Can also take NULL value (no transformation previous to subsetting).
#' @param detection Detection threshold for absence/presence. See \code{\link[microbiome]{core}} for more details.
#' @param prevalence Prevalence threshold (in \\[ 0, 1 \\]). See \code{\link[microbiome]{core}} for more details.
#' @details This function uses \code{\link[microbiome]{transform}} and \code{\link[microbiome]{core}} functions to get a
#' compositional-transformed phyloseq including only prevalent taxa.
#' @return \item{Reads}{Input phyloseq object.}
#' @return \item{RA}{Transformed initial phyloseq object. If \code{transf = NULL} it will be the same as \code{Reads} component.}
#' @return \item{Core}{Subsetted phyloseq. Transformed to compositional.}
#' @keywords phyloseq transform core
#' @export

phy_transform <- function(physeq, transf = "compositional", detection, prevalence){
  phyTR <- phyloseq

  ## Convert to compositional data
  if(!is.null(transf)){
    phyTR.rel <- microbiome::transform(phyTR, transf)
  } else {
    phyTR.rel <- phyTR
  }

  ## Pick core taxa with with the given prevalence and detection limits
  phyTR.core <- core(phyTR.rel, detection = detection, prevalence = prevalence, include.lowest = FALSE)

  ## Use relative abundances for the core
  phyTR.core <- microbiome::transform(phyTR.core, "compositional")
  phyTR.list <- list(phyTR, phyTR.rel, phyTR.core)
  names(phyTR.list) <- c("Reads", "RA", "Core")

  return(phyTR.list)
}


#' Filling blank taxonomy:
#'
#' \code{microDA} function from the \code{funs_phy} group (expanding \code{phyloseq} functionality).\cr
#' Fills blank or NA cells of a taxonomy table using the information of the superior level. Input data can be either a data frame
#' or a \code{tax_table} from a phyloseq object.
#' @param data Input dataset. Either a data frame or a \code{phyloseq} object.
#' @param physeq Whether input is a phyloseq object or not. Defaults TRUE.
#' @param unclassified Whether to use 'Unclassified' as prefix for filling (TRUE) or
#' \code{c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")} as prefix (FALSE). Defaults FALSE.
#' @details Time and again taxonomy tables from microbiome analysis software do not have a complete taxonomy, and one OTU may be
#' classified to Phylum level while other may be classified to Species level, leaving blank cells in the table. This function
#' fills the gaps using a prefix plus the information from the former cell.
#' @details As an example:\cr
#' OTU1: Kingdom: Bacteria; Phylum: NA
#' @details With \code{unclassified = FALSE}: \cr
#' OTU1: Kingdom: Bacteria; Phylum: Unclassified_Bacteria; Class: Unclassified_Bacteria.
#' @details With \code{unclassified = TRUE}: \cr
#' OTU1: Kingdom: Bacteria; Phylum: Kingdom_Bacteria; Class: Phylum_Bacteria.
#' @details Note that input data must have \strong{7 taxonomy columns}.
#' @keywords phyloseq taxonomy fill
#' @export

phy_filltax <- function(data, physeq = TRUE, unclassified = FALSE){
  if(isFALSE(unclassified)){
    prefxs <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
  }
  else if (isTRUE(unclassified)){
    prefxs <- rep("Unclassified", 6)
  }
  if(isTRUE(physeq)){tax.clean <- data.frame(tax_table(data))}
  else {tax.clean <- data}
  for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}
  tax.clean[is.na(tax.clean)] <- ""
  for (i in 1:nrow(tax.clean)){
    if (tax.clean[i,2] == ""){
      kingdom <- paste(prefxs[1], tax.clean[i,1], sep = "_")
      for (j in 2:7){
        tax.clean[i, j] <- kingdom
      }
    } else if (tax.clean[i,3] == ""){
      phylum <- paste(prefxs[2], tax.clean[i,2], sep = "_")
      for (j in 3:7){
        tax.clean[i, j] <- phylum
      }
    } else if (tax.clean[i,4] == ""){
      class <- paste(prefxs[3], tax.clean[i,3], sep = "_")
      for (j in 4:7){
        tax.clean[i, j] <- class
      }
    } else if (tax.clean[i,5] == ""){
      order <- paste(prefxs[4], tax.clean[i,4], sep = "_")
      for (j in 5:7){
        tax.clean[i, j] <- order
      }
    } else if (stringr::str_detect(tax.clean[i,5], "uncultured")){
      order <- paste(prefxs[4], tax.clean[i,4], "uncult", sep = "_")
      for (j in 5:7){
        tax.clean[i, j] <- order
      }
    } else if (tax.clean[i,6] == ""){
      family <- paste(prefxs[5], tax.clean[i,5], sep = "_")
      for (j in 6:7){
        tax.clean[i, j] <- family
      }
    } else if (stringr::str_detect(tax.clean[i,6], "uncultured")){
      family <- paste(prefxs[5], tax.clean[i,5], "uncult", sep = "_")
      for (j in 6:7){
        tax.clean[i, j] <- family
      }
    } else if (tax.clean[i,7] == ""){
      tax.clean[i,7] <- paste(prefxs[6], tax.clean[i,6], sep = "_")
    }
  }
  if(isTRUE(physeq)){tax.return <- tax_table(as.matrix(tax.clean))}
  else {tax.return <- tax.clean}
  return(tax.return)
}


#' Remove NA taxonomy columns:
#'
#' \code{microDA} function from the \code{funs_phy} group (expanding \code{phyloseq} functionality).\cr
#' This function is useful after applying \code{\link[phyloseq]{tax_glom}} to a phyloseq object. It removes all the taxonomy
#' columns from \code{tax_table(physeq)} composed exclusively by NAs.
#' @param physeq Input phyloseq object.
#' @keywords phyloseq taxonomy glom NA
#' @export

phy_rm_na_tax <- function(physeq){
  # rm_all <- function(x) { Filter(function(x)!all(is.na(x)), df) }
  rm_all <- function(df) { df[, !apply(is.na(df), 2, all)] }
  tax_table(physeq) <- rm_all( tax_table(physeq) )
  return(physeq)
}
