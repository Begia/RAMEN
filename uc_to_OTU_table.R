#!/usr/bin/Rscript
#Matti Ruuskanen, Jul 2016
#Makes an OTU table from a clustered *.uc file.
#Possibility to remove OTUs from the OTU table with less than a specified number of reads.
#If given a matching fasta file of cluster seeds, renames the sequences and removes
#OTUs with less than the user specified minimum number of reads also from the fasta file.

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("Biostrings"))

ipak <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, repos = "http://cran.us.r-project.org", dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c("optparse", "Biostrings")
ipak(packages)

option_list <- list(
  make_option(
    c("-c", "--clusters"),
    type = "character",
    default = NULL,
    help = "cluster file name",
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
    c("-f", "--fasta"),
    type = "character",
    default = NULL,
    help = "fasta file name containing cluster seeds",
    metavar = "character"
  ),
  make_option(
    c("-m", "--min"),
    type = "integer",
    default = 0,
    help = "minimum number of reads in OTU to be included [default = disabled]",
    metavar = "integer"
  )
)


opt_parser <- OptionParser(option_list = option_list)

opt <- parse_args(opt_parser)


if (is.null(opt$clusters)) {
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).", call. = FALSE)
}

setwd(opt$out)

df <- read.table(opt$clusters)
df <- df[,-c(6, 7)]

cSizes <- df[df$V1 == "C",]
cSizes <- cSizes[, 2:3]
cSizes[, 1] <- paste("OTU", cSizes[, 1], sep = "_")
colnames(cSizes) <- c("#OTU_ID", "count")
removedOTUs <- cSizes[cSizes$count <= opt$min,]
cSizes <- cSizes[!cSizes$`#OTU_ID` %in% removedOTUs$`#OTU_ID`,]

if (opt$min > 0) {
  print(paste(
    "Removed ",
    nrow(cSizes),
    " OTUs with size under the threshold of ",
    opt$min,
    sep = ""
  ))
  print("")
}

df <- df[!df$V1 == "C", ]
df <- df[, c(1, 2, 7)]
df <- df[order(df$V2), ]
colnames(df) <- c("rType", "OTU", "Sample")
df$Sample <- gsub("(.*)_.*", "\\1", df$Sample)
df$Dummy <- rep(1, nrow(df))
df2 <- aggregate(Dummy ~ Sample * OTU, FUN = length, data = df)

otuTable <- (xtabs(Dummy ~ OTU + Sample, data = df2))
names(attributes(otuTable)$dimnames) <- NULL
otuTable = as.data.frame.matrix(otuTable)

otuTable$OTU_ID = paste("OTU", row.names(otuTable), sep = "_")

otuTable = otuTable[, c(length(otuTable), c(1:length(otuTable) - 1))]
colnames(otuTable)[1] <- "#OTU_ID"

otuTable <- otuTable[otuTable$`#OTU_ID` %in% cSizes$`#OTU_ID`, ]

print(paste(
  "A total of ",
  colSums(cSizes[2]),
  " reads in ",
  nrow(cSizes),
  " OTUs, with:",
  sep = ""
))
print("")
for (i in 2:length(colnames(otuTable))) {
  print(paste(
    colnames(otuTable)[i],
    ": ",
    colSums(otuTable[i]),
    " reads in ",
    length(which(otuTable[, i] != 0)),
    " OTUs",
    sep = ""
  ))
}

write.table(
  otuTable,
  file = "otu_table.txt",
  quote = F,
  sep = "\t",
  row.names = F
)

if (!is.null(opt$fasta)) {
  cSeeds <- readDNAStringSet(opt$fasta)
  names(cSeeds) <-
    paste("OTU", seq(1:length(names(cSeeds))), sep = "_")
  cSeeds <- cSeeds[which(names(cSeeds) %in% otuTable$`#OTU_ID`)]
  fileName <- sub("(.+)\\.f.*", "\\1", opt$fasta)
  writeXStringSet(
    cSeeds,
    filepath = paste(fileName, "renamed.fasta", sep = "_"),
    format = "fasta"
  )
}
