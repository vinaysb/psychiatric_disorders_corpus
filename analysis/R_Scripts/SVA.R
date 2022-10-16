#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)
suppressMessages(library("ArrayExpress"))
suppressMessages(library("GEOquery"))
suppressMessages(library("dplyr"))
suppressMessages(library("data.table"))
suppressMessages(library("optparse"))
suppressMessages(library("affycoretools"))
suppressMessages(library("sva"))
suppressMessages(library("pamr"))
suppressMessages(library("Biobase"))
suppressMessages(library("limma"))
suppressMessages(library("oligo"))

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="Name of the file containing dataset accession number (separated by a newline).", metavar="character"),
  make_option(c("-d", "--dir"), type="character", default=NULL,
              help="Directory of the files.", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("Both arguments must be supplied (--file missing).", call.=FALSE)
}

if (is.null(opt$dir)){
  print_help(opt_parser)
  stop("Both arguments must be supplied (--dir missing).", call.=FALSE)
}

dataset = opt$file
datadir = paste0(opt$dir, "/")

if (!(all(file.exists(paste0(datadir, dataset, "_diffex.tsv")), file.exists(paste0(datadir, dataset, "_annotated.tsv"))))) {
  print(paste("Starting processing ", dataset, " dataset."))

  tmp = read.table(paste0(datadir, dataset, "_data.tsv"))
  exprs = t(as.matrix(tmp))

  pdata = read.table(paste0(datadir, dataset, "_design.tsv"), row.names=1, header=TRUE)

  normData = new("ExpressionSet", exprs = exprs, phenoData = AnnotatedDataFrame(pdata))

  mod = model.matrix(~as.factor(Target), data=pdata)
  mod0 = model.matrix(~1, data=pdata)

  n.sv = num.sv(exprs,mod,method="leek")
  print(paste0("Number of SV's: ", n.sv))

  if (n.sv > 0){
    # SVA
    print(paste("Carrying out SVA on ", dataset))
    svobj = sva(exprs, mod, mod0, n.sv=n.sv)

    pvals = f.pvalue(exprs, mod, mod0)
    qvals = p.adjust(pvals, method="BH")

    print(paste0("Number of significant genes (pvals): ", length(pvals[pvals < 0.05])))
    print(paste0("Number of significant genes (qvals): ", length(qvals[qvals < 0.05])))

    modSv = cbind(mod, svobj$sv)
    mod0Sv = cbind(mod0, svobj$sv)
    pvals.sv = f.pvalue(exprs, modSv, mod0Sv)
    qvals.sv = p.adjust(pvals.sv, method="BH")

    print(paste0("Number of significant genes (pvals): ", length(pvals.sv[pvals.sv < 0.05])))
    print(paste0("Number of significant genes (qvals): ", length(qvals.sv[qvals.sv < 0.05])))
  }

  # Annotation
  annotatedData = paste0(datadir, dataset,"_annotated.tsv")

  if (startsWith(dataset, "E-")){
    print(paste("Re-downloading ", dataset, " dataset."))
    sink("/dev/null")
    if (dataset == "E-GEOD-53987"){
      path = paste(dirname(getwd()), "E-GEOD-53987", sep="/")
      expFiles = getAE("E-GEOD-53987", path, type = "raw", local=TRUE)
      rawset = ae2bioc(expFiles)
    } else {
      rawset = ArrayExpress(dataset)
    }
    sink()
    if (!is.null(rawset)) {
      print(paste("Re-normalising ", dataset, " dataset."))
      sink("/dev/null")
      tmp = oligo::rma(rawset)
      sink()
      featureNames(normData) = featureNames(tmp)
      annotation(normData) = annotation(tmp)

      print(paste("Annotating ", dataset))

      if(annotation(rawset) == "pd.hg.u133a"){
        annotated_eset = annotateEset(
          normData,
          "hgu133a.db",
          columns = c("PROBEID", "ENTREZID", "SYMBOL", "GENENAME"),
          multivals = "first"
        )
      }

      else if(annotation(rawset) == "pd.hg.u133.plus.2"){
        annotated_eset = annotateEset(
          normData,
          "hgu133plus2.db",
          columns = c("PROBEID", "ENTREZID", "SYMBOL", "GENENAME"),
          multivals = "first"
        )
      }

      else if(annotation(rawset) == "pd.hg.u95av2"){
        annotated_eset = annotateEset(
          normData,
          "hgu95av2.db",
          columns = c("PROBEID", "ENTREZID", "SYMBOL", "GENENAME"),
          multivals = "first"
        )
      }

      else{
        annotated_eset = annotateEset(
          normData,
          annotation(rawset),
          columns = c("PROBEID", "ENTREZID", "SYMBOL", "GENENAME"),
          multivals = "first"
        )
      }

      df = pData(featureData(annotated_eset))

      symbols = df[, "SYMBOL"]
      normData_annotated = exprs(normData)[!is.na(symbols),]
      symbols_filtered = symbols[!is.na(symbols)]

      normData_annotated = transpose(as.data.table(normData_annotated))
      rownames(normData_annotated) = colnames(normData)
      colnames(normData_annotated) = symbols_filtered

      fwrite(normData_annotated, file = annotatedData, row.names=TRUE, sep="\t")

      lm.exprs = t(as.matrix(normData_annotated))
      colnames(lm.exprs) = colnames(normData)
    }
  }

  else if (startsWith(dataset, "GSE")){
    print(paste("Re-downloading ", dataset, " dataset."))
    sink("/dev/null")
    rawset = getGEO(dataset, GSEMatrix=TRUE)[[1]]
    sink()

    print(paste("Annotating ", dataset))

    symbols = fData(rawset)[,'Gene Symbol']
    normData_annotated = exprs(normData)[!is.na(symbols),]
    symbols_filtered = symbols[!is.na(symbols)]

    normData_annotated = transpose(as.data.table(normData_annotated))
    rownames(normData_annotated) = colnames(normData)
    colnames(normData_annotated) = symbols_filtered

    fwrite(normData_annotated, file = annotatedData, row.names=TRUE, sep="\t")

    lm.exprs = t(as.matrix(normData_annotated))
    colnames(lm.exprs) = colnames(normData)
  }

  else{
    print(paste0("Error:", datasets[i],"is a bad accession number"))
    quit(save="no", status = 0)
  }

  op_file = paste0(datadir, dataset, "_diffex.tsv")

  # Limma
  if (n.sv > 0){
    print(paste("Calculating the differential expression on ", dataset, ", using limma (with SVA)."))

    fit = lmFit(lm.exprs, modSv)
    contrasts.matrix = cbind("C1"=c(-1,1,rep(0, svobj$n.sv)))
    fitContrasts = contrasts.fit(fit, contrasts.matrix)

    eb = eBayes(fitContrasts)
    limma_output = topTable(eb, adjust="BH")
    fwrite(limma_output, file = op_file, row.names=TRUE, sep="\t")
  }

  else {
    print(paste("Calculating the differential expression on ", dataset, ", using limma (without SVA)."))
    f = factor(pdata$Target, levels=unique(pdata$Target))
    design = model.matrix(~0+f)
    fit = lmFit(lm.exprs, design)
    print(fit$genes)

    contrasts.matrix = makeContrasts(paste(colnames(design)[1], "-", colnames(design)[2]), levels=design)
    fitContrasts = contrasts.fit(fit, contrasts.matrix)

    eb = eBayes(fitContrasts)
    limma_output = topTable(eb, adjust="BH", coef=1, number=Inf, confint=0.95)
    fwrite(limma_output, file = op_file, row.names=TRUE, sep="\t")
  }
}
