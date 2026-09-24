# GDR-3627: does the native backend register through the existing fitting layer,
# with no change to the pipeline?
#
# Registers a native (drc-free) fit function as a fit_fn and writes it into the
# same Metrics assay under its own fit_source, next to the native gDR rows.

suppressMessages({
  library(data.table)
  library(gDRutils)
  library(MultiAssayExperiment)
})

# --- native backend, same model / bounds as .prepareFitModel(), ec50 on log10 --

.native_one_curve <- function(conc, y, x_0 = 1, cap = 0.1, pcutoff = 0.05, n_point_cutoff = 4L) {
  keep <- !is.na(y) & !is.na(conc) & conc > 0
  conc <- conc[keep]; y <- pmin(y[keep], x_0 + cap)
  if (uniqueN(conc) < n_point_cutoff) {
    return(list(fit_type = "DRCTooFewPointsToFit", x_mean = mean(y), r2 = NA_real_, xc50 = NA_real_))
  }
  if (uniqueN(y) == 1L) {
    # same guard logisticFit() applies before it ever calls a solver
    return(list(fit_type = "DRCConstantFitResult", x_mean = y[1], r2 = 0, ec50 = 0, h = 1e-04,
                xc50 = Inf, x_inf = y[1], x_0 = y[1]))
  }
  med <- stats::median(conc); mn <- min(conc)
  start <- c(2, 0.4, log10(med))
  lowerl <- c(0.1, 0, log10(mn / 10))
  upperl <- c(5, min(x_0 + cap, 1), log10(max(conc) * 10))
  start <- pmin(pmax(start, lowerl), upperl)
  sse <- function(p) sum((y - (p[2] + (x_0 - p[2]) / (1 + exp(p[1] * (log(conc) - log(10 ^ p[3])))))) ^ 2)
  res <- tryCatch(stats::optim(start, sse, method = "L-BFGS-B", lower = lowerl, upper = upperl,
                               control = list(factr = 1e7, maxit = 500)), error = function(e) NULL)
  if (is.null(res)) {
    # the state gDRutils documents but never reaches: solver tried and failed
    return(list(fit_type = "DRCInvalidFitResult", x_mean = mean(y), r2 = NA_real_, xc50 = NA_real_))
  }
  h <- res$par[1]; x_inf <- res$par[2]; ec50 <- 10 ^ res$par[3]
  rss <- res$value
  rss1 <- sum((y - mean(y)) ^ 2)
  r2 <- 1 - rss / rss1
  p <- gDRutils:::.calculate_f_pval(2L, length(y) - 2L, rss1, rss)
  if (is.na(p) || p >= pcutoff) {
    return(list(fit_type = "DRCConstantFitResult", x_mean = mean(y), r2 = 0, ec50 = 0,
                h = 1e-04, xc50 = Inf, x_inf = mean(y), x_0 = mean(y)))
  }
  list(fit_type = "DRC3pHillFitModelFixS0", x_mean = mean(y), r2 = r2, p_value = p,
       h = h, x_inf = x_inf, x_0 = x_0, ec50 = ec50, rss = rss,
       xc50 = gDRutils::predict_conc_from_efficacy(0.5, x_inf = x_inf, x_0 = x_0, ec50 = ec50, h = h))
}

native_fit_fn <- function(dt) {
  .native_one_curve(dt$Concentration, dt$x)
}

# --- run it through the existing layer, no pipeline change ------------------

mae <- gDRutils::get_synthetic_data("finalMAE_small.qs2")
se <- mae[["single-agent"]]

cat("Metrics before:", paste(SummarizedExperiment::assayNames(se), collapse = ", "), "\n")
before <- as.data.table(gDRutils::convert_se_assay_to_dt(se, "Metrics"))
cat("native gDR rows:", nrow(before), "| fit_source column present:",
    "fit_source" %in% names(before), "\n")

se2 <- gDRutils::apply_fit(
  se, native_fit_fn, "single-agent",
  output_assay = "Metrics", fit_source = "native_optim", merge = "merge", on_error = "warn"
)

after <- as.data.table(gDRutils::convert_se_assay_to_dt(se2, "Metrics"))
cat("\nrows after apply_fit:", nrow(after), "\n")
print(after[, .N, by = fit_source])
print(after[, .N, by = .(fit_source, fit_type)][order(fit_source, -N)])

# --- compare the two sources on the same curves -----------------------------

keys <- intersect(c("rId", "cId", "normalization_type"), names(after))
a <- after[is.na(fit_source) | fit_source != "native_optim"]
b <- after[fit_source == "native_optim"]
cmp <- merge(a[, c(keys, "fit_type", "xc50", "ec50", "h", "x_inf"), with = FALSE],
             b[, c(keys, "fit_type", "xc50", "ec50", "h", "x_inf"), with = FALSE],
             by = keys, suffixes = c("_drc", "_native"))
cat("\nmatched curves:", nrow(cmp), "\n")
cat("fit_type agreement:", sum(cmp$fit_type_drc == cmp$fit_type_native), "of", nrow(cmp), "\n")
sig <- cmp[fit_type_drc == "DRC3pHillFitModelFixS0" & fit_type_native == "DRC3pHillFitModelFixS0"]
if (nrow(sig)) {
  rel <- function(x, y) abs(x - y) / pmax(abs(y), 1e-12)
  cat(sprintf("xc50 median rel.diff %.2e | max %.2e (n = %d)\n",
              stats::median(rel(sig$xc50_native, sig$xc50_drc), na.rm = TRUE),
              max(rel(sig$xc50_native, sig$xc50_drc), na.rm = TRUE), nrow(sig)))
}
cat("\ndisagreements by normalization_type:\n")
print(cmp[fit_type_drc != fit_type_native, .N, by = c("normalization_type", "fit_type_drc", "fit_type_native")])
cat("\nassays after:", paste(SummarizedExperiment::assayNames(se2), collapse = ", "), "\n")
