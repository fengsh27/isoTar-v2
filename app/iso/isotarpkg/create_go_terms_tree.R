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
              help="Directory path", metavar="character"),
  make_option(c("-g", "--go"), type="character", default=NULL, 
              help="GO terms list", metavar="character"),
  make_option(c("-w", "--work"), type="character", default=NULL, 
              help="Working directory path", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (is.null(opt$input) | is.null(opt$go) | is.null(opt$work)){
  print_help(opt_parser)
  stop("Provide all parameters.", call.=FALSE)
}

setwd(opt$work)

# Create the GO terms list
go_list <- unique(strsplit(opt$go, ',')[[1]])
  
# Create the GO terms dot file
dot_file <- getAmigoTree(goIDs=go_list, filename=opt$input, picType="dot", saveResult=TRUE)

  