#!/usr/bin/Rscript
#Matti Ruuskanen, Jul 2016

ipak <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, repos = "http://cran.us.r-project.org", dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- "optparse"
ipak(packages)

option_list <- list(
  make_option(
    c("-f", "--file"),
    type = "character",
    default = NULL,
    help = "dataset file name",
    metavar = "character"
  ),
  make_option(
    c("-o", "--out"),
    type = "character",
    default = getwd(),
    help = "output folder [default = %default]",
    metavar = "character"
  ),
)


opt_parser <- OptionParser(option_list = option_list)

opt <- parse_args(opt_parser)


if (is.null(opt$file)) {
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).", call. = FALSE)
}

setwd(opt$out)

df <- read.table(opt$file)
df <- df[,-c(6, 7)]
cSizes <- df[df$V1 == "C",]
df <- df[!df$V1 == "C",]
df <- df[, c(1, 2, 7)]
df <- df[order(df$V2),]
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
write.table(
  otuTable,
  file = "otu_table.txt",
  quote = F,
  sep = "\t",
  row.names = F
)
