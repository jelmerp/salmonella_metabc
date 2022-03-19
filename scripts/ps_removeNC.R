#!/usr/bin/env Rscript

# SBATCH --account=PAS0471
# SBATCH --time=5

## Report
message("\n## Starting script ps_removeNC.R")
Sys.time()
message()

## Load packages
library(phyloseq)

## Process command-line args
args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
outfile <- args[2]

## Remove NC sample
ps_in <- readRDS(infile)
ps_out <- subset_samples(ps_in, treatment != "NC")
saveRDS(ps_out, outfile)

## Report
message("\n## Input object:")
print(ps_in)
message("\n## Output object:")
print(ps_out)

message("\n## Listing output file:")
system(paste("ls -lh", outfile))
message("## Done with script ps_removeNC.R")
Sys.time()
message()