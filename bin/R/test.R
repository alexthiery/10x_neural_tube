#!/usr/bin/env Rscript
args <- commandArgs(TRUE)

# read all files from folder and keep only those from chr_edit
files <- Filter(function(x) grepl("chr_edit", x), list.files(args[2], recursive = T, full.names = T)
)

# remove file suffix
file.path <- dirname(files)[!duplicated(dirname(files))]

# make dataframe with tissue matching directory
tissue = c("hh4", "hh6", "4ss", "8ss")
matches <- sapply(tissue, function(x) file.path[grep(pattern = x, x = file.path)])
sample.paths <- data.frame(tissue = names(matches), path = matches, row.names = NULL)
sample.paths$tissue = c("hh4", "hh6", "ss4", "ss8")

sapply(list.files(args[1], full.names = T), source)

plot_path <- "plots/"
dir.create(plot_path, recursive = T)
processed_data <- "processed_data/"
dir.create(processed_data, recursive = T)

sample.paths