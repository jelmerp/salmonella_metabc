## Running DESeq at the ASV level returns an error ('every gene contains at least one zero')
## Therefore, we first filter the data to make sure this is no longer the case
prep_taxrank_asv <- function() {
  counts <- t(otu_table(ps)) %>% as(., "matrix") %>% round()
  counts <- counts[apply(counts, 1, prev_fun), ]
  
  asv_w_least_miss <- names(sort(apply(counts, 1, function(x) sum(x == 0))))[1]
  samp_miss <- which(counts[asv_w_least_miss, ] == 0)
  counts2 <- counts[, -samp_miss]
  meta2 <- meta %>% filter(! sampleID %in% names(samp_miss))
  dds <- mk_deseq(counts2, meta2)
  dds <- estimateSizeFactors(dds)
    
  counts_norm <- sweep(counts2, 2, sizeFactors(dds), "/")
    
  tax <- as.data.frame(tax_table(ps)) %>%
    rownames_to_column("ID") %>%
    dplyr::select(-domain) %>%
    filter(ID %in% row.names(counts))
  
  taxdat <- list(ps = ps, dds = dds, counts_norm = counts_norm,
                 tax = tax, meta = meta2)
  return(taxdat)
}

if (taxrank == "ASV") {
  taxdat <- prep_taxrank_asv()
} else {
  taxdat <- prep_taxrank(taxrank)
}