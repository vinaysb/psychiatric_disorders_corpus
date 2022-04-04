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

datasets = readLines(opt$file)
datadir = paste0(dirname(normalizePath(opt$file)), "/")
datasetFiles = paste0(datadir, datasets, ".tsv")


if (!(all(file.exists(paste0(datadir, datasets, ".tsv")), file.exists(paste0(datadir, datasets, "_annotated.tsv"))))) {
  print("Started processing datasets!")

  pb <- txtProgressBar(min = 0, max = length(datasets), style = 3)

  suppressMessages(
    for (i in 1:length(datasets)) {
      fname_exprsSet = datasetFiles[i]
      
      if (file.exists(fname_exprsSet)) {
      	next
      }

      # Process ArrayExpress Data
      if (startsWith(datasets[i], "E-")){
        rawset = ArrayExpress(datasets[i])
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

          annotatedData = paste0(datadir, datasets[i],"_annotated.tsv")

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
        }
      }
	    
      # Process GEO Data
	    else if (startsWith(datasets[i], "GSE")){
				normData = getGEO(datasets[i], GSEMatrix=TRUE)[[1]]

        write.table(normData, file = fname_exprsSet, row.names=TRUE, sep="\t")

        annotatedData = paste0(datadir, datasets[i],"_annotated.tsv")

        symbols = fData(normData)[,'Gene Symbol']
        normData_annotated = exprs(normData)[!is.na(symbols),]
        symbols_filtered = symbols[!is.na(symbols)]
        
        normData_annotated = transpose(as.data.table(normData_annotated))
        rownames(normData_annotated) = colnames(normData)
        colnames(normData_annotated) = symbols_filtered

        fwrite(normData_annotated, file = annotatedData, row.names=TRUE, sep="\t")

	    }
	    else{
	    	print(paste0("Error:", datasets[i],"is a bad accession number"))
          	break
	    }
      setTxtProgressBar(pb, i-1)
	  }
  )
  close(pb)
  print("Finished processing datasets!")
}