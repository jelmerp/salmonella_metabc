contam_method <- "prevalence"    # Method for contam detection: "prevalence" (neg. control), "frequency" (DNA conc.), "either", "both", or "NA"
contam_thres <- 0.1              # P-value for an ASV to be considered a contaminant 

conc_column <- NA                # Name of the column in the metadata containing DNA concentrations (use 'NA' if none)
batch_column <- NA               # Name of the column in the metadata containing batch IDs (use 'NA' if none)
neg_column <- "neg_control"      # Name of the column in the metadata indicating neg. control status with TRUE/FALSE; specify either `neg_column` or `neg_ids` to identify negative controls (use 'NA' if none)
neg_ids <- NA                    # IDs of samples that are neg. controls; specify either `neg_column` or `neg_ids` to identify negative controls (use 'NA' if none) 

rm_offtarget <- TRUE             # Whether to remove off-target taxa: Chloroplasts, Mitochondria & Eukaryotes

min_ASV <- 1000                  # Min. total ASV count for a sample; sample will be excluded if it has a lower value 
