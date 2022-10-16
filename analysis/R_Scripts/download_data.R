#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)
suppressMessages(library("ArrayExpress"))
suppressMessages(library("GEOquery"))
suppressMessages(library("oligo"))
suppressMessages(library("dplyr"))
suppressMessages(library("data.table"))
suppressMessages(library("optparse"))
suppressMessages(library("affycoretools"))


option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="Name of the file containing dataset accession number (separated by a newline).", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).", call.=FALSE)
}

dataset = opt$file
datasetFile = paste0(dataset, ".tsv")

if (!(all(file.exists(paste0(dataset, ".tsv")), file.exists(paste0(dataset, "_annotated.tsv"))))) {
  print(paste0("Started processing dataset - ", dataset, "!"))

  fname_exprsSet = datasetFile

  if (file.exists(fname_exprsSet)) {
  	quit(save="no", status = 0)
  }

  # Process ArrayExpress Data
  if (startsWith(dataset, "E-")){
    rawset = ArrayExpress(dataset)
    if (!is.null(rawset)) {
			if (length(rawset)>1) {
				indx = which(names(rawset)=="A-AFFY-44")
				rawData = rawset[[indx]]
				normData = oligo::rma(rawData)
			}
			else {
				normData = oligo::rma(rawset)
			}
      write.table(normData, file = fname_exprsSet, row.names=TRUE, sep="\t")
    }
  }

  # Process GEO Data
  else if (startsWith(dataset, "GSE")){
		normData = getGEO(dataset, GSEMatrix=TRUE)[[1]]

    write.table(normData, file = fname_exprsSet, row.names=TRUE, sep="\t")
  }

  else{
  	print(paste0("Error:", dataset,"is a bad accession number"))
      	break
  }
  print("Finished processing dataset!")
}
