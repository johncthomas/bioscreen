# NOTE this is not a python file, but I've never found a consistent way to get data files
#  to actually get installed except by putting ".py" on the end

# These paths are the linux box system R, otherwise it uses the local conda
#  env, and I don't want to install a bunch of R packages into the general env
libs = c("/usr/local/lib/R/site-library", "/usr/lib/R/site-library", "/usr/lib/R/library")
# add them to the path
.libPaths(libs)

library(limma)
library(edgeR)
library(glue)

verbose = FALSE
vprint = function(thing) {
  if (verbose) {
    print(thing, quote=FALSE)
  }
}


do_voom = function(counts, design, block=NULL){
    v = voom(counts, design)
    corfit = NULL
    if (!is.null(block)) {
        # run this a couple of times
        corfit <- duplicateCorrelation(v, design, block = block)
        v = voom(counts, design, block=block, correlation = corfit$consensus.correlation)
        corfit <- duplicateCorrelation(v, design, block = block)
        v = voom(counts, design, block=block, correlation=corfit$consensus.correlation)
    }
    res = list(v, corfit$consensus.correlation)
    names(res) = c("vcounts", "correlation")
    return(res)
}



get_design = function(counts, details, test_groups) {
    # make sure the sample names match
  namesokay = all(rownames(details) == colnames(counts))
  stopifnot(namesokay)

  test_groups = as.factor(test_groups)

  design = model.matrix(~0+test_groups)

  colnames(design) = levels(test_groups)
  rownames(design) = colnames(counts)
  #as.data.frame(design)
  return(design)

}

prep_rnaseq_counts = function(counts, sample_groups, filter_exp=TRUE){

  # remove low expressing
  if (filter_exp) {
    groups = as.factor(sample_groups)
    keepers = filterByExpr(counts, group = groups)

    ngenes = dim(counts)[1]
    counts = counts[keepers,]
    remgenes = dim(counts)[1]
    remd = ngenes - remgenes
    s = paste(c("Removed", remd, "genes for low expression.", remgenes, "genes kept"),
              collapse = " ")
    vprint(s)
  }

  # Convert counts to a DGEList object
  dge <- DGEList(counts)

  # Calculate normalization factors using the TMM method
  dge <- calcNormFactors(dge)

  return(dge)

}



get_fit = function(dge, design,
                   block=NULL, correlation=NULL){

  if (!is.null(block)) {
    vprint("Removing covar from factors:")
    vprint(block)
    fit = lmFit(dge,
                design,
                block=block,
                correlation=correlation$consensus.correlation)
  } else {
    vprint("Not removing any covariance")
    fit = lmFit(dge,
                design,)
  }
  return(fit)
}

# contrasts is a character vector giving formulas, e.g. , "Treat1 - Ctrl1"
# adds $contrastsFit and $contrastsMat
fit_contrasts = function(fit, design, contrasts, contrast_names) {

  cm = makeContrasts(contrasts=contrasts, levels=design)
  colnames(cm) = contrast_names
  fit2 = contrasts.fit(fit, cm)
  fit2 = eBayes(fit2)
  res = list(fit2, cm)
  names(res) = c("contrastFit", "contrastMatrix")
  return(res)
}

get_toptables = function(contrastRes) {
  tables = list()
  contrast_names = colnames(contrastRes$contrastMatrix)

  for (c in contrast_names) {

    toptab = topTable(contrastRes$contrastFit, coef = c, number = Inf)
    toptab = as.data.frame(toptab)

    tables = append(tables, toptab)
  }

  names(tables) = contrast_names
  return(tables)
}

get_toptable = function(contrastRes, contrastName) {

  toptab = topTable(contrastRes$contrastFit, coef = contrastName, number = Inf)
  toptab = as.data.frame(toptab)

  return(toptab)

}

write_toptables = function(contrastRes, outprefix) {
  for (c in colnames(contrastRes$contrastMatrix)) {
    outfn = paste(outprefix, c, 'csv', sep = '.')
    print(outfn)
    write.csv(topTable(contrastRes$contrastFit, coef = c, number = Inf), outfn)
  }
}
