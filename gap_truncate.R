#!/usr/bin/Rscript
#Matti Ruuskanen, Jul 2016
#gap truncate: truncate alignment from both ends to the first column with a base character in a set fraction of the sequences

ipak <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, repos = "http://cran.us.r-project.org", dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c("optparse", "seqinr")
ipak(packages)

option_list <- list(
  make_option(
    c("-a", "--alignment"),
    type = "character",
    default = NULL,
    help = "alignment file name",
    metavar = "character"
  ),
  make_option(
    c("-o", "--out"),
    type = "character",
    default = getwd(),
    help = "output folder [default = %default]",
    metavar = "character"
  ),
  make_option(
    c("-g", "--gaps"),
    type = "double",
    default = 0.5,
    help = "threshold for fraction of sequences with a base in the column [default = %default]",
    metavar = "double"
  )
)

opt_parser <- OptionParser(option_list = option_list)

opt <- parse_args(opt_parser)

if (is.null(opt$clusters)) {
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).", call. = FALSE)
}

setwd(opt$out)

file <- read.fasta(opt$fasta, forceDNAtolower = F)

i <- 1
x <- length(file)
while (x > (length(file) * opt$gaps)) {
  x <- table(getFrag(file, i, i) == "-")["TRUE"]
  i <- i + 1
}

file2 <- getFrag(file, i - 1, min(getLength(file)))

i <- min(getLength(file2))
x <- length(file2)
while (x > (length(file2) * opt$gaps)) {
  x <- table(getFrag(file2, i, i) == "-")["TRUE"]
  names(x) <- NULL
  i <- i - 1
}

file3 <- getFrag(file2, 1, i + 1)
fileName <- sub(".+/+(.+)\\.f.*", "\\1", opt$alignment)
write.fasta(
  sequences = file3,
  names = names(file),
  file.out = paste(fileName, "truncated.fasta", sep = "_")
)
