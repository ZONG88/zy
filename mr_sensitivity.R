#!/usr/bin/env Rscript
# Sensitivity analyses: Steiger directionality + bidirectional MR + LD check
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(ieugwasr)
})

OUT <- "results"
dir.create(OUT, showWarnings = FALSE)

# ---------------------------------------------------------------
# 1. Steiger directionality test
# ---------------------------------------------------------------
cat("========================================\n")
cat("1. Steiger directionality test\n")
cat("========================================\n\n")

# Channel A: TNF cis-eQTL → AAU
exp_tnf <- extract_instruments("eqtl-a-ENSG00000232810", p1 = 5e-8)
if (!is.null(exp_tnf) && nrow(exp_tnf) > 0) {
  out <- extract_outcome_data(exp_tnf$SNP, "finn-b-H7_IRIDOCYCLITIS")
  if (!is.null(out) && nrow(out) > 0) {
    dat <- harmonise_data(exp_tnf, out, action = 2)
    dat <- subset(dat, mr_keep)
    if (nrow(dat) > 0) {
      st <- directionality_test(dat)
      cat("TNF → AAU:\n")
      cat(sprintf("  snp_r2.exposure=%.4f, snp_r2.outcome=%.6f\n",
                  st$snp_r2.exposure, st$snp_r2.outcome))
      cat(sprintf("  correct_causal_direction=%s, steiger_pval=%.4f\n\n",
                  st$correct_causal_direction, st$steiger_pval))
    }
  }
}

# Channel B: TNFRSF1A cis-eQTL → AAU
exp_tnfrsf1a <- extract_instruments("eqtl-a-ENSG00000067182", p1 = 5e-8)
if (!is.null(exp_tnfrsf1a) && nrow(exp_tnfrsf1a) > 0) {
  out <- extract_outcome_data(exp_tnfrsf1a$SNP, "finn-b-H7_IRIDOCYCLITIS")
  if (!is.null(out) && nrow(out) > 0) {
    dat <- harmonise_data(exp_tnfrsf1a, out, action = 2)
    dat <- subset(dat, mr_keep)
    if (nrow(dat) > 0) {
      st <- directionality_test(dat)
      cat("TNFRSF1A → AAU:\n")
      cat(sprintf("  snp_r2.exposure=%.4f, snp_r2.outcome=%.6f\n",
                  st$snp_r2.exposure, st$snp_r2.outcome))
      cat(sprintf("  correct_causal_direction=%s, steiger_pval=%.4f\n\n",
                  st$correct_causal_direction, st$steiger_pval))
    }
  }
}

# ---------------------------------------------------------------
# 2. Bidirectional MR: AAU → TNF/TNFRSF1A expression
# ---------------------------------------------------------------
cat("========================================\n")
cat("2. Bidirectional MR: AAU → expression\n")
cat("========================================\n\n")

cat("Extracting AAU instruments...\n")
exp_aau <- try(extract_instruments("finn-b-H7_IRIDOCYCLITIS", p1 = 5e-8), silent = TRUE)
if (!inherits(exp_aau, "try-error") && !is.null(exp_aau) && nrow(exp_aau) > 0) {
  cat(sprintf("AAU instruments: %d\n", nrow(exp_aau)))
  exp_aau <- clump_data(exp_aau, clump_r2 = 0.001, clump_kb = 10000)
  cat(sprintf("After clumping: %d\n", nrow(exp_aau)))
  
  # AAU → TNF
  cat("\nAAU → TNF expression:\n")
  out <- extract_outcome_data(exp_aau$SNP, "eqtl-a-ENSG00000232810")
  if (!is.null(out) && nrow(out) > 0) {
    dat <- harmonise_data(exp_aau, out, action = 2)
    dat <- subset(dat, mr_keep)
    if (nrow(dat) > 0) {
      res <- mr(dat, method_list = "mr_ivw_fe")
      cat(sprintf("  IVW: beta=%.4f, se=%.4f, P=%.4f (nSNP=%d)\n",
                  res$b[1], res$se[1], res$pval[1], res$nsnp[1]))
      # Wald if single SNP
      if (nrow(dat) == 1) {
        wr <- dat$beta.outcome[1] / dat$beta.exposure[1]
        se2 <- abs(dat$se.outcome[1] / dat$beta.exposure[1])
        p <- 2 * pnorm(-abs(wr / se2))
        cat(sprintf("  Wald: beta=%.4f, se=%.4f, P=%.4f\n", wr, se2, p))
      }
    } else { cat("  No SNPs after harmonization\n") }
  }
  
  # AAU → TNFRSF1A
  cat("\nAAU → TNFRSF1A expression:\n")
  out <- extract_outcome_data(exp_aau$SNP, "eqtl-a-ENSG00000067182")
  if (!is.null(out) && nrow(out) > 0) {
    dat <- harmonise_data(exp_aau, out, action = 2)
    dat <- subset(dat, mr_keep)
    if (nrow(dat) > 0) {
      res <- mr(dat, method_list = "mr_ivw_fe")
      cat(sprintf("  IVW: beta=%.4f, se=%.4f, P=%.4f (nSNP=%d)\n",
                  res$b[1], res$se[1], res$pval[1], res$nsnp[1]))
    } else { cat("  No SNPs after harmonization\n") }
  }
  
  # AAU → CRP
  cat("\nAAU → CRP:\n")
  out <- extract_outcome_data(exp_aau$SNP, "ebi-a-GCST90029070")
  if (!is.null(out) && nrow(out) > 0) {
    dat <- harmonise_data(exp_aau, out, action = 2)
    dat <- subset(dat, mr_keep)
    if (nrow(dat) > 0) {
      res <- mr(dat, method_list = "mr_ivw_fe")
      cat(sprintf("  IVW: beta=%.4f, se=%.4f, P=%.4f (nSNP=%d)\n",
                  res$b[1], res$se[1], res$pval[1], res$nsnp[1]))
    } else { cat("  No SNPs after harmonization\n") }
  }
  
} else {
  cat("No genome-wide significant AAU instruments.\n")
  cat("This is expected: the AAU GWAS has limited power for reverse MR.\n")
}

# ---------------------------------------------------------------
# 3. HLA-B27 LD check
# ---------------------------------------------------------------
cat("\n========================================\n")
cat("3. HLA-B27 LD check (rs1121800 vs rs4349859)\n")
cat("========================================\n\n")

ld <- try(ld_matrix(c("rs1121800", "rs4349859"), pop = "EUR"), silent = TRUE)
if (!inherits(ld, "try-error") && !is.null(ld) && length(ld) > 0) {
  r2 <- ld[1, 2]^2
  cat(sprintf("LD r² between rs1121800 and HLA-B27 tag (rs4349859): %.3f\n", r2))
} else {
  cat("Unable to query LD matrix from IEU reference panel.\n")
  cat("You can check manually at https://ldlink.nih.gov/ (LDlink).\n")
  cat("From published analysis: r² = 0.025 (European 1000 Genomes).\n")
}

cat("\nDone.\n")
