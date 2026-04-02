###################################################
### GO terms Tree creation
###################################################
library(RamiGO)
library(DiagrammeR)
library(DiagrammeRsvg)
library(magrittr)
library(rsvg)
library("optparse")

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input path", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="Output path", metavar="character"),
  make_option(c("-w", "--work"), type="character", default=NULL, 
              help="Working directory path", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (is.null(opt$input) | is.null(opt$output)){
  print_help(opt_parser)
  stop("Provide all parameters.", call.=FALSE)
}

setwd(opt$work)

if (file.exists(opt$input)){
  rsvg_pdf(charToRaw(export_svg(grViz(diagram=opt$input))), file=opt$output)
}
