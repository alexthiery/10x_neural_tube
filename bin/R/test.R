#!/usr/bin/env Rscript
args <- commandArgs(TRUE)

cat("R script is running")


plot_path <- "plots/"
dir.create(plot_path, recursive = T)
processed_data <- "processed_data/"
dir.create(processed_data, recursive = T)

sapply(list.files(args[1], full.names = T), source)
