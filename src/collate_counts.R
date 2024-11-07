#!/usr/bin/env Rscript

if (exists("snakemake")) {
    log <- file(snakemake@log[[1]], open = "wt")
    sink(log, type = "message")
    sink(log, append = TRUE, type = "output")

    # inputs
    unclustered_kmer_count_files <- snakemake@input[["unclustered_kmer_counts"]]
    cluster_count_files <- snakemake@input[["clustered_kmer_counts"]]
    # outputs
    csv_file <- snakemake@output[["csv"]]
} else {
    unclustered_kmer_count_files <- list.files(
        "output/020_kmer-counts",
        full.names = TRUE,
        pattern = "kmer_counts.json",
        recursive = TRUE
    )
    cluster_count_files <- list.files(
        "output/030_usearch-clusters",
        full.names = TRUE,
        pattern = "count.csv",
        recursive = TRUE
    )
    csv_file <- "counts.csv"
}


library(jsonlite)
library(data.table)


# Collate unclustered counts
names(unclustered_kmer_count_files) <- sapply(unclustered_kmer_count_files, function(x) {
    unlist(strsplit(
        x,
        "/",
        fixed = TRUE
    ))[[3]]
})

unclustered_kmer_count_list <- sapply(
    unclustered_kmer_count_files, function(x) {
        fromJSON(x)$Stats$`#Unique_k-mers`
    }
)

unclustered_kmer_counts <- data.table(
    file = names(unclustered_kmer_count_list),
    identity = 100,
    count = unclustered_kmer_count_list
)

# Collate clustered counts
cluster_counts <- rbindlist(
    lapply(
        cluster_count_files,
        fread,
        col.names = c("file", "identity", "count")
    )
)

all_counts <- rbind(unclustered_kmer_counts, cluster_counts)
setkey(all_counts, file, identity)

fwrite(all_counts, csv_file)

sessionInfo()
