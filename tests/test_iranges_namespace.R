#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(GenomicRanges)
    library(IRanges)
})

exports <- getNamespaceExports("IRanges")
stopifnot(all(c("overlapsAny", "findOverlaps") %in% exports))

query <- GRanges(seqnames = c("chr1", "chr1"),
                 ranges = IRanges(start = c(10L, 100L), width = 10L),
                 strand = "*")
subject <- GRanges(seqnames = "chr1",
                   ranges = IRanges(start = 15L, width = 10L),
                   strand = "*")

stopifnot(identical(
    IRanges::overlapsAny(query, subject, ignore.strand = TRUE),
    c(TRUE, FALSE)
))
hits <- IRanges::findOverlaps(query, subject, ignore.strand = TRUE)
stopifnot(length(hits) == 1L, queryHits(hits)[[1L]] == 1L,
          subjectHits(hits)[[1L]] == 1L)

cat("OK   IRanges overlap generics dispatch correctly for GRanges\n")
