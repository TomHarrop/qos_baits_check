#!/usr/bin/env Rscript

library(jsonlite)
library(data.table)

unclustered_kmer_count_files <- list.files(
    "output/020_kmer-counts",
    full.names = TRUE,
    pattern = "kmer_counts.json",
    recursive = TRUE
)

names(unclustered_kmer_count_files) <- unlist(
    strsplit(
        unclustered_kmer_count_files, "/",
        fixed = TRUE
    )
)[[3]]

unclustered_kmer_count_list <- sapply(
    unclustered_kmer_count_files, function(x) {
        fromJSON(x)$Stats$`#Unique_k-mers`
    }
)

unclustered_kmer_counts <- data.table(
    file = names(unclustered_kmer_count_list),
    identity = 100,
    count = unclustered_kmer_count_list)

