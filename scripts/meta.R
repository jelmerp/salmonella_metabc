#!/usr/bin/env Rscript

#SBATCH --account=PAS0471
#SBATCH --output=slurm-meta-%j.out

## Report
message("## Starting script meta.R")
Sys.time()
message()

## Load packages
if (!"pacman" %in% installed.packages()) install.packages("pacman")
pkgs <- c("tidyverse")
pacman::p_load(pkgs, character.only = TRUE)

## Process command-line arguments
args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]      # infile <- "metadata/meta_raw.tsv"
outfile <- args[2]     # outfile <- "metadata/meta.tsv"

## Read and process raw metadata
meta <- read_tsv(infile, col_names = c("s_nr", "well", "sampleID")) %>%
  mutate(sampleID = gsub(" ", "-", sampleID),  # Replace spaces and + by underscores
         # There is a weird mismatch in the water samples S-numbers and filename nrs,
         # so we add 2:
         sampleID = ifelse(sampleID %in% c("water", "Water"),
                           paste0(s_nr + 2, "-Water"), sampleID),
         treatment_org = gsub("\\d", "", sampleID),
         treatment_org = gsub("-$|^-", "", treatment_org),
         treatment = case_when(treatment_org == "LGG" ~ "L",
                               treatment_org == "PC+S" ~ "S",
                               TRUE ~ treatment_org),
         lacto = ifelse(treatment %in% c("L", "L+S"), TRUE, FALSE),
         salmo = ifelse(treatment %in% c("S", "L+S"), TRUE, FALSE),
         neg_control = ifelse(treatment == "Water", TRUE, FALSE)) %>%
  arrange(sampleID) %>%
  select(sampleID, treatment, lacto, salmo, treatment_org, neg_control, s_nr)

## Print number of samples
message("\n## Total nr of samples: ", nrow(meta))
message("## Nr of samples by treatment:")
print(meta %>% count(treatment))
message("\n## Showing head of meta df:")
print(head(meta))

## Write output file
write_tsv(meta, outfile)

## Report
message("\n## Done with script meta.R")
Sys.time()
message()