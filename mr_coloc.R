#!/usr/bin/env Rscript
# Colocalization: TNFRSF1A region + TNF/MHC region
# Uses eQTLGen VCF + IEU OpenGWAS API
suppressPackageStartupMessages({
  library(coloc)
  library(data.table)
  library(ieugwasr)
})

OUT <- "coloc_results"
dir.create(OUT, showWarnings = FALSE)

# ---------------------------------------------------------------
# Helper: parse eQTL VCF
# ---------------------------------------------------------------
parse_eqtl_vcf <- function(vcf_path) {
  con <- gzfile(vcf_path, "rt")
  lines <- readLines(con)
  close(con)
  data_start <- grep("^#CHROM", lines) + 1
  data_lines <- lines[data_start:length(lines)]
  
  n <- length(data_lines)
  chr <- character(n); pos <- integer(n); id <- character(n)
  es <- numeric(n); se <- numeric(n); lp <- numeric(n); af <- numeric(n)
  
  for (i in seq_len(n)) {
    parts <- strsplit(data_lines[i], "\t")[[1]]
    chr[i] <- parts[1]
    pos[i] <- as.integer(parts[2])
    id[i] <- parts[3]
    gt <- strsplit(parts[10], ":")[[1]]
    es[i] <- as.numeric(gt[1])
    se[i] <- as.numeric(gt[2])
    lp[i] <- as.numeric(gt[3])
    af[i] <- as.numeric(gt[4])
  }
  data.frame(chr, pos, id, es, se, lp, af, stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------
# Helper: run coloc for one region
# ---------------------------------------------------------------
run_coloc_region <- function(vcf_path, out_id, n_eqtl, out_type, out_n, out_s, label) {
  cat(sprintf("\n========== %s ==========\n", label))
  
  # Parse eQTL VCF
  d <- parse_eqtl_vcf(vcf_path)
  d <- d[!is.na(d$lp) & is.finite(d$lp), ]
  cat(sprintf("eQTL SNPs: %d\n", nrow(d)))
  
  # Get top SNPs and query outcome data
  d <- d[order(-d$lp), ]
  top_n <- min(200, nrow(d))
  top_snps <- d$id[1:top_n]
  
  aau_list <- list()
  for (i in seq(1, top_n, by = 50)) {
    end <- min(i + 49, top_n)
    r <- try(associations(top_snps[i:end], out_id), silent = TRUE)
    if (!inherits(r, "try-error") && !is.null(r) && nrow(r) > 0) {
      aau_list[[length(aau_list) + 1]] <- r
    }
  }
  
  if (length(aau_list) == 0) { cat("No outcome data\n"); return(NULL) }
  aau <- rbindlist(aau_list)
  
  # Merge
  m <- merge(d, aau, by.x = "id", by.y = "rsid")
  m <- m[complete.cases(m$es, m$se, m$beta, m$se), ]
  cat(sprintf("Overlapping SNPs: %d\n", nrow(m)))
  
  if (nrow(m) < 10) { cat("Too few SNPs for coloc\n"); return(NULL) }
  
  maf <- pmin(m$af, 1 - m$af)
  D1 <- list(beta = m$es, varbeta = m$se^2, snp = m$id,
             type = "quant", N = n_eqtl, MAF = maf, position = m$pos)
  D2 <- list(beta = m$beta, varbeta = m$se^2, snp = m$id,
             type = out_type, s = out_s, N = out_n, MAF = maf,
             position = m$position)
  
  res <- coloc.abf(D1, D2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  cat("Coloc results:\n")
  print(res$summary)
  
  saveRDS(list(summary = res$summary, results = res$results, n_snps = nrow(m)),
          file.path(OUT, sprintf("%s_coloc.rds", gsub(" ", "_", label))))
  
  # Prior sensitivity
  for (p12 in c(1e-4, 1e-6)) {
    res2 <- coloc.abf(D1, D2, p1 = 1e-4, p2 = 1e-4, p12 = p12)
    cat(sprintf("  Prior p12=%s: PP.H4=%.3f\n", format(p12, scientific = TRUE),
                res2$summary[6]))
  }
  
  return(res)
}

# ---------------------------------------------------------------
# Region 1: TNFRSF1A (chr12:6.2-6.7 Mb)
# ---------------------------------------------------------------
cat("\n========================================")
cat("\nColocalization analysis")
cat("\n========================================\n")

# TNFRSF1A: chr12:6.2-6.7 Mb
# eQTL VCF from eQTLGen
vcf_tnfrsf1a <- "data/tnfrsf1a_eqtl.vcf.gz"
if (file.exists(vcf_tnfrsf1a)) {
  run_coloc_region(vcf_tnfrsf1a, "finn-b-H7_IRIDOCYCLITIS",
                   n_eqtl = 31684, out_type = "cc",
                   out_n = 212909, out_s = 3622 / 212909,
                   label = "TNFRSF1A chr12")
} else {
  cat("\nTNFRSF1A VCF not found at", vcf_tnfrsf1a, "\n")
  cat("Download eQTLGen data for chr12:6.2-6.7 Mb first.\n")
}

# ---------------------------------------------------------------
# Region 2: TNF/MHC (chr6:31.3-31.8 Mb)
# ---------------------------------------------------------------
vcf_tnf <- "data/tnf_mhc_eqtl.vcf.gz"
if (file.exists(vcf_tnf)) {
  run_coloc_region(vcf_tnf, "finn-b-H7_IRIDOCYCLITIS",
                   n_eqtl = 14263, out_type = "cc",
                   out_n = 212909, out_s = 3622 / 212909,
                   label = "TNF MHC chr6")
} else {
  cat("\nTNF/MHC VCF not found at", vcf_tnf, "\n")
  cat("Download eQTLGen data for chr6:31.3-31.8 Mb first.\n")
}

cat("\nDone. Results saved to", OUT, "/\n")
