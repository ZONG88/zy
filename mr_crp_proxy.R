#!/usr/bin/env Rscript
# CRP proxy channel: rs1800693 (within TNFRSF1A) → AAU
# Qualitative sensitivity only (small CRP beta → inflated CI)
suppressPackageStartupMessages({
  library(TwoSampleMR)
})

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

# rs1800693: CRP-lowering allele (C) proxies TNF/TNFR1 inhibition
# Beta = -0.0218, SE = 0.002, F = 118.8
crp_exp <- data.frame(
  SNP = "rs1800693",
  chr.exposure = 12,
  pos.exposure = 6440009,
  effect_allele.exposure = "C",
  other_allele.exposure = "A",
  beta.exposure = -0.0218,
  se.exposure = 0.002,
  pval.exposure = 9.6e-27,
  samplesize.exposure = 575531,
  exposure = "CRP proxy (TNF/TNFR1 inhibition)",
  id.exposure = "crp_proxy",
  mr_keep.exposure = TRUE,
  stringsAsFactors = FALSE
)

outcomes <- list(
  c("AAU", "finn-b-H7_IRIDOCYCLITIS"),
  c("Ankylosing spondylitis", "finn-b-M13_ANKYLOSPON"),
  c("Crohn's disease", "ieu-a-12"),
  c("Rheumatoid arthritis", "finn-b-M13_RHEUMA"),
  c("Multiple sclerosis", "finn-b-G6_MS")
)

cat("CRP proxy (rs1800693) → outcomes\n")
cat("================================\n\n")

for (oc in outcomes) {
  name <- oc[1]; id <- oc[2]
  cat(sprintf("  %s (%s) ... ", name, id))
  
  out <- tryCatch(
    extract_outcome_data("rs1800693", id),
    error = function(e) NULL
  )
  if (is.null(out) || nrow(out) == 0) { cat("no data\n"); next }
  
  dat <- harmonise_data(crp_exp, out, action = 2)
  dat <- subset(dat, mr_keep)
  if (nrow(dat) == 0) { cat("harmonization failed\n"); next }
  
  wr <- wald_ratio(dat$beta.exposure[1], dat$se.exposure[1],
                   dat$beta.outcome[1], dat$se.outcome[1])
  cat(sprintf("OR=%.2f [%.2f, %.2f] P=%.4f (F=%.1f)\n",
              wr$OR, wr$CI_low, wr$CI_high, wr$pval, wr$F))
}

cat("\nNote: CRP proxy estimates are directionally informative only.\n")
cat("Small exposure beta leads to inflated OR and wide CI.\n")
