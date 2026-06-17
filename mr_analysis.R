#!/usr/bin/env Rscript
#
# TNF/TNFRSF1A drug-target cis-MR → AAU
# Single cis-SNP Wald ratio per Karhunen et al. (2026)
# Data from IEU OpenGWAS API
#
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
})

set.seed(20260602)

# -----------------------------------------------------------------
# Instrument info (all within ±10kb of target gene)
# -----------------------------------------------------------------
# TNF eQTL (eQTLGen): rs1121800, chr6:31535074, beta=-0.2915, F=599
# TNFRSF1A eQTL (eQTLGen): rs1800692, chr12:6442346, beta=-0.370, F=974
# CRP proxy (within TNFRSF1A): rs1800693, chr12:6440009, beta=-0.0218

IVS <- list(
  tnf = list(
    snp = "rs1121800",
    data_id = "eqtl-a-ENSG00000232810",
    beta = -0.2915, se = 0.0119,
    effect_allele = "T", other_allele = "A",
    eaf = 0.601,
    label = "TNF cis-eQTL"
  ),
  tnfrsf1a = list(
    snp = "rs1800692",
    data_id = "eqtl-a-ENSG00000067182",
    beta = -0.3700, se = 0.0119,
    effect_allele = "C", other_allele = "T",
    eaf = 0.619,
    label = "TNFRSF1A cis-eQTL"
  )
)

OUTCOMES <- list(
  c("AAU (iridocyclitis)", "finn-b-H7_IRIDOCYCLITIS", "Primary"),
  c("Ankylosing spondylitis", "finn-b-M13_ANKYLOSPON", "Positive control (protective)"),
  c("Crohn's disease", "ieu-a-12", "Positive control (protective)"),
  c("Rheumatoid arthritis", "finn-b-M13_RHEUMA", "Positive control (protective)"),
  c("Multiple sclerosis", "finn-b-G6_MS", "Positive control (risk)"),
  c("Type 2 diabetes", "ebi-a-GCST006867", "Specificity"),
  c("Hypertension", "ukb-b-12493", "Specificity")
)

# -----------------------------------------------------------------
# Wald ratio for single cis-SNP
# -----------------------------------------------------------------
wald_ratio <- function(bx, sx, by, sy) {
  w <- by / bx
  se_w <- abs(sqrt((sy^2 / bx^2) + (by^2 * sx^2 / bx^4)))
  or <- exp(w)
  ci_l <- exp(w - 1.96 * se_w)
  ci_u <- exp(w + 1.96 * se_w)
  p <- 2 * pnorm(-abs(w / se_w))
  f <- (bx / sx)^2
  list(wald = w, se = se_w, OR = or, CI_low = ci_l, CI_high = ci_u, pval = p, F = f)
}

# -----------------------------------------------------------------
# Helper: run one channel
# -----------------------------------------------------------------
run_channel <- function(iv, outcomes, label) {
  cat(sprintf("\n========== %s ==========\n", label))
  cat(sprintf("Instrument: %s\n", iv$snp))
  cat(sprintf("Exposure: %s\n\n", iv$data_id))
  
  # Extract exposure
  cat("Extracting exposure...")
  exp <- extract_instruments(iv$data_id, p1 = 5e-8)
  if (is.null(exp) || nrow(exp) == 0) {
    # Manual fallback
    cat(" manual construction\n")
    exp <- data.frame(
      SNP = iv$snp,
      chr.exposure = NA, pos.exposure = NA,
      effect_allele.exposure = iv$effect_allele,
      other_allele.exposure = iv$other_allele,
      beta.exposure = iv$beta,
      se.exposure = iv$se,
      pval.exposure = 2 * pnorm(-abs(iv$beta / iv$se)),
      samplesize.exposure = ifelse(iv$data_id == "eqtl-a-ENSG00000232810", 14263, 31684),
      exposure = iv$label,
      id.exposure = iv$data_id,
      mr_keep.exposure = TRUE,
      stringsAsFactors = FALSE
    )
  } else {
    cat(sprintf(" %d SNPs found\n", nrow(exp)))
  }
  
  results <- data.frame()
  
  for (oc in outcomes) {
    outc_name <- oc[1]
    out_id <- oc[2]
    role <- oc[3]
    
    cat(sprintf("  %s (%s) ... ", outc_name, out_id))
    
    out <- tryCatch(
      extract_outcome_data(snps = unique(exp$SNP), outcomes = out_id),
      error = function(e) NULL
    )
    if (is.null(out) || nrow(out) == 0) {
      cat("no outcome data\n")
      next
    }
    
    dat <- harmonise_data(exp, out, action = 2)
    dat <- subset(dat, mr_keep)
    
    if (nrow(dat) == 0) {
      cat("harmonization failed\n")
      next
    }
    
    wr <- wald_ratio(
      dat$beta.exposure[1], dat$se.exposure[1],
      dat$beta.outcome[1], dat$se.outcome[1]
    )
    
    cat(sprintf("OR=%.2f [%.2f, %.2f] P=%.4f\n",
                wr$OR, wr$CI_low, wr$CI_high, wr$pval))
    cat(sprintf("    Wald=%.4f, SE=%.4f, F=%.1f\n",
                wr$wald, wr$se, wr$F))
    
    # Harmonization info
    cat(sprintf("    EA_exposure=%s EA_outcome=%s palindromic=%s\n",
                dat$effect_allele.exposure[1],
                dat$effect_allele.outcome[1],
                dat$palindromic[1]))
    
    results <- rbind(results, data.frame(
      Channel = label,
      SNP = iv$snp,
      Outcome = outc_name,
      Role = role,
      Wald = wr$wald,
      OR = wr$OR,
      CI_low = wr$CI_low,
      CI_high = wr$CI_high,
      P = wr$pval,
      F = wr$F,
      Harmonized = nrow(dat),
      stringsAsFactors = FALSE
    ))
  }
  
  return(results)
}

# -----------------------------------------------------------------
# Run
# -----------------------------------------------------------------
cat("============================================================\n")
cat("TNF/TNFRSF1A drug-target cis-MR analysis\n")
cat("============================================================\n")
cat(sprintf("R %s | TwoSampleMR %s\n\n",
            R.version.string,
            as.character(packageVersion("TwoSampleMR"))))

dir.create("results", showWarnings = FALSE)

r1 <- run_channel(IVS$tnf, OUTCOMES, "TNF")
r2 <- run_channel(IVS$tnfrsf1a, OUTCOMES, "TNFRSF1A")

all_res <- rbind(r1, r2)

# -----------------------------------------------------------------
# Print final table
# -----------------------------------------------------------------
cat("\n\n============================================================\n")
cat("FINAL RESULTS\n")
cat("============================================================\n\n")

cat(sprintf("%-12s %-30s %-8s %-6s %-20s %-10s\n",
            "Channel", "Outcome", "Role", "OR", "95% CI", "P"))
cat(paste(rep("-", 90), collapse = ""), "\n")

for (i in 1:nrow(all_res)) {
  r <- all_res[i, ]
  p_str <- ifelse(r$P < 0.0001, "<0.0001", sprintf("%.4f", r$P))
  cat(sprintf("%-12s %-30s %-8s %-6.2f [%-6.2f, %-6.2f] %s\n",
              r$Channel, r$Outcome, r$Role, r$OR, r$CI_low, r$CI_high, p_str))
}

# -----------------------------------------------------------------
# Save
# -----------------------------------------------------------------
write.csv(all_res, "results/final_cis_mr_results.csv", row.names = FALSE)
cat(sprintf("\nSaved to results/final_cis_mr_results.csv\n"))
cat("Done.\n")
