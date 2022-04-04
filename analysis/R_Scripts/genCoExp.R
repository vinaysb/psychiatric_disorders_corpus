#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)
suppressMessages(library("WGCNA"))
suppressMessages(library("data.table"))
suppressMessages(library("dplyr"))
suppressMessages(library("doParallel"))
suppressMessages(library("optparse"))
suppressMessages(library("tools"))

# Script options
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="Name of the file containing dataset gene expression which will be used to create a co-expression network.", metavar="character"),
  make_option(c("-d", "--dataset"), type="character", default=NULL,
        	  help="Name of the dataset/disease")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).", call.=FALSE)
}

if (is.null(opt$dataset)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (dataset name).", call.=FALSE)
}

fname_exprsSet = opt$file
dataset = tools::file_path_sans_ext(basename(opt$file))

datadir = paste0(dirname(normalizePath(opt$file)), "/")
resultdir = paste0("./results/", opt$dataset, "/")

# WGCNA network generation
header = unlist(strsplit(readLines(fname_exprsSet ,n=1), "\t"))
data = fread(file=fname_exprsSet, sep='\t', header=FALSE, quote="", data.table=FALSE)

column_names = data$V1

data = subset(data, select=-V1)
data = mutate_all(data, function(x) as.numeric(as.character(x)))
# data = transpose(data)

colnames(data) = header[header != "index"]
rownames(data) = column_names
data = data[-1,, drop=F]

threads = 30
net = blockwiseModules(t(data), power = 30, maxBlockSize = nrow(data),
                       minModuleSize = 20, verbose = 4,
                       saveTOMs = TRUE, loadTOM = FALSE,
                       saveTOMFileBase = paste0(datadir, dataset, "_coexp_network"), nThreads=threads)

save(net, file = paste0(datadir, dataset, "_coexp_network.rda"))

# Thresholding
load(paste0(datadir, dataset, "_coexp_network-block.1.RData"))
load(paste0(paste0(datadir, dataset, "_coexp_network.rda")))
print(1)

vis = exportNetworkToVisANT(TOM, weighted = TRUE, threshold = 0, probeToGene=data.frame(c(seq(from = 1, to=length(rownames(data)))), rownames(data)))

print(2)

threshold = quantile(vis$weight, 0.99)

print(threshold)

vis2 = subset(vis, weight>=threshold)

write.table(net$colors, file=paste0(resultdir, dataset, "_coexp_network_modules.tsv"), sep="\t", quote=F, row.names=T, col.names=F)
write.table(vis2, file=paste0(resultdir, dataset, "_coexp_network_edges.tsv"), sep="\t", quote=F, row.names=F, col.names=T)