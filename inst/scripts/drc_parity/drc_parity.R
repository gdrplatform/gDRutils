# Parity probe: a native bounded fit against drc::drm(), compared on what reaches
# the user rather than on the coefficients (GDR-3627).
#
# For each curve and each backend (drc / own / own_log) the script reproduces the
# rest of gDRutils::logisticFit(): F-test against the flat model, the constant-fit
# rule (p >= pcutoff or ec50 NA), r2 and xc50. It then compares
#   (a) fit_type   - sigmoidal vs constant, which is what the report shows
#   (b) xc50       - on curves both backends call sigmoidal
#   (c) parameters - restricted to those same curves
# and times both backends.

suppressMessages({
  library(data.table)
  library(gDRutils)
  library(MultiAssayExperiment)
})

OUT_DIR <- Sys.getenv("PARITY_OUT", file.path(tempdir(), "drc-parity"))
N_POINT_CUTOFF <- 4L
CAP <- 0.1
X_0 <- 1
PCUTOFF <- 0.05
# Path to the examples/ directory of a gDRworkshops checkout:
# git clone https://github.com/gdrplatform/gDRworkshops
WS <- Sys.getenv("GDRWORKSHOPS_EXAMPLES", "../gDRworkshops/examples")

DATASETS <- c(
  PRISM = file.path(WS, "PRISMBroadScreen_Hagenbeek_NatComm_2026/gDR_data/gDR_mae.qs2"),
  Zhou = file.path(WS, "SmallDrugCombo_Zhou_CellChemBio_2026/gDR_data/gDR_mae.qs2"),
  Goetz = file.path(WS, "LargeDrugCombo_Goetz_Cancers_2024/gDR_data/gDR_mae.qs2")
)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

ll_mean <- function(conc, h, x_inf, x_0, ec50) {
  x_inf + (x_0 - x_inf) / (1 + exp(h * (log(conc) - log(ec50))))
}

sse <- function(par, conc, y, x_0, log_ec50) {
  ec50 <- if (log_ec50) 10 ^ par[3] else par[3]
  sum((y - ll_mean(conc, par[1], par[2], x_0, ec50)) ^ 2)
}

fit_drc <- function(conc, y, x_0, start, lowerl, upperl) {
  controls <- drc::drmc(relTol = 1e-04, errorm = FALSE, noMessage = TRUE, rmNA = TRUE)
  fct <- drc::LL.3u(upper = x_0, names = c("h", "x_inf", "ec50"))
  fit <- tryCatch(
    drc::drm(y ~ conc, data = data.table(y = y, conc = conc), logDose = NULL, fct = fct,
             start = start, lowerl = lowerl, upperl = upperl,
             control = controls, na.action = stats::na.omit),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  co <- stats::coef(fit)
  list(h = co[[1]], x_inf = co[[2]], ec50 = co[[3]],
       rss = sum(stats::residuals(fit) ^ 2, na.rm = TRUE))
}

fit_own <- function(conc, y, x_0, start, lowerl, upperl, log_ec50 = FALSE) {
  if (log_ec50) {
    start[3] <- log10(start[3]); lowerl[3] <- log10(lowerl[3]); upperl[3] <- log10(upperl[3])
  }
  res <- tryCatch(
    stats::optim(par = start, fn = sse, method = "L-BFGS-B", lower = lowerl, upper = upperl,
                 conc = conc, y = y, x_0 = x_0, log_ec50 = log_ec50,
                 control = list(factr = 1e7, maxit = 500)),
    error = function(e) NULL
  )
  if (is.null(res)) return(NULL)
  list(h = res$par[1], x_inf = res$par[2],
       ec50 = if (log_ec50) 10 ^ res$par[3] else res$par[3], rss = res$value)
}

# rest of logisticFit(): F-test, constant-fit rule, r2, xc50
finish <- function(fit, y, x_0) {
  if (is.null(fit)) {
    return(list(fit_type = "DRCInvalidFitResult", r2 = NA_real_, p = NA_real_, xc50 = NA_real_))
  }
  rss1 <- sum((y - mean(y, na.rm = TRUE)) ^ 2, na.rm = TRUE)
  r2 <- 1 - fit$rss / rss1
  nparam <- 3L                       # 3-param fixed-x_0 path, as in production
  df1 <- nparam - 1L
  df2 <- length(stats::na.omit(y)) - nparam + 1L
  p <- gDRutils:::.calculate_f_pval(df1, df2, rss1, fit$rss)
  if (is.na(p) || p >= PCUTOFF || is.na(fit$ec50)) {
    return(list(fit_type = "DRCConstantFitResult", r2 = r2, p = p, xc50 = NA_real_))
  }
  xc50 <- gDRutils::predict_conc_from_efficacy(efficacy = 0.5, x_inf = fit$x_inf,
                                               x_0 = x_0, ec50 = fit$ec50, h = fit$h)
  list(fit_type = "DRC3pHillFitModelFixS0", r2 = r2, p = p, xc50 = xc50)
}

run_dataset <- function(label, path) {
  mae <- qs2::qs_read(path)
  exps <- intersect(names(mae), c("single-agent", "combination"))
  out <- list()
  for (e in exps) {
    dt <- tryCatch(gDRutils::convert_se_assay_to_dt(mae[[e]], "Averaged"), error = function(err) NULL)
    if (is.null(dt)) next
    setDT(dt)
    if (!all(c("Concentration", "x", "normalization_type") %in% names(dt))) next
    dt <- dt[!is.na(x) & !is.nan(x) & Concentration > 0]
    if (!nrow(dt)) next
    keys <- intersect(c("rId", "cId", "normalization_type"), names(dt))
    dt[, .exp := e]
    out[[e]] <- list(dt = dt, keys = keys)
  }
  rows <- list()
  t_drc <- t_own <- 0
  for (e in names(out)) {
    dt <- out[[e]]$dt
    keys <- out[[e]]$keys
    for (nm in names(split(dt, by = keys, drop = TRUE))) {
      cur <- split(dt, by = keys, drop = TRUE)[[nm]]
      norm_type <- cur$normalization_type[1]
      cur <- cur[, .(x = mean(x, na.rm = TRUE)), by = Concentration]
      if (uniqueN(cur$Concentration) < N_POINT_CUTOFF) next
      conc <- cur$Concentration
      y <- pmin(cur$x, X_0 + CAP)
      med_conc <- stats::median(conc); min_conc <- min(conc)
      if (norm_type == "GR") {
        start <- c(2, 0.1, med_conc); lowerl <- c(0.1, -1, min_conc / 10)
      } else {
        start <- c(2, 0.4, med_conc); lowerl <- c(0.1, 0, min_conc / 10)
      }
      upperl <- c(5, min(X_0 + CAP, 1), max(conc) * 10)
      start <- pmin(pmax(start, lowerl), upperl)

      t0 <- proc.time()[["elapsed"]]; d <- fit_drc(conc, y, X_0, start, lowerl, upperl)
      t_drc <- t_drc + proc.time()[["elapsed"]] - t0
      t0 <- proc.time()[["elapsed"]]; o <- fit_own(conc, y, X_0, start, lowerl, upperl, FALSE)
      t_own <- t_own + proc.time()[["elapsed"]] - t0
      ol <- fit_own(conc, y, X_0, start, lowerl, upperl, TRUE)

      fd <- finish(d, y, X_0); fo <- finish(o, y, X_0); fol <- finish(ol, y, X_0)
      rows[[length(rows) + 1L]] <- data.table(
        dataset = label, experiment = e, curve = nm, norm_type = norm_type,
        n_conc = uniqueN(conc),
        drc_type = fd$fit_type, own_type = fo$fit_type, own_log_type = fol$fit_type,
        drc_p = fd$p, own_p = fo$p,
        drc_r2 = fd$r2, own_r2 = fo$r2,
        drc_xc50 = fd$xc50, own_xc50 = fo$xc50, own_log_xc50 = fol$xc50,
        drc_h = if (is.null(d)) NA_real_ else d$h, own_h = if (is.null(o)) NA_real_ else o$h,
        drc_x_inf = if (is.null(d)) NA_real_ else d$x_inf,
        own_x_inf = if (is.null(o)) NA_real_ else o$x_inf,
        drc_ec50 = if (is.null(d)) NA_real_ else d$ec50,
        own_ec50 = if (is.null(o)) NA_real_ else o$ec50,
        drc_rss = if (is.null(d)) NA_real_ else d$rss,
        own_rss = if (is.null(o)) NA_real_ else o$rss,
        own_log_rss = if (is.null(ol)) NA_real_ else ol$rss
      )
    }
  }
  list(res = rbindlist(rows), t_drc = t_drc, t_own = t_own)
}

all_res <- list(); timing <- list()
for (lbl in names(DATASETS)) {
  if (!file.exists(DATASETS[[lbl]])) next
  message("== ", lbl)
  r <- run_dataset(lbl, DATASETS[[lbl]])
  if (!nrow(r$res)) next
  all_res[[lbl]] <- r$res
  timing[[lbl]] <- data.table(dataset = lbl, curves = nrow(r$res),
                              drc_s = round(r$t_drc, 1), own_s = round(r$t_own, 1))
}
res <- rbindlist(all_res)
fwrite(res, file.path(OUT_DIR, "poc_fit_comparison_v2.csv"))

rel <- function(a, b) abs(a - b) / pmax(abs(b), 1e-12)

cat("\n==== native backend vs drc - outcome-level comparison ====\n\n")
print(rbindlist(timing))

cat("\n-- fit_type agreement (what the report shows)\n")
print(res[, .N, by = .(dataset, drc_type, own_type)][order(dataset, -N)])

cat("\n-- own_log vs drc fit_type\n")
print(res[, .N, by = .(drc_type, own_log_type)][order(-N)])

sig <- res[drc_type == "DRC3pHillFitModelFixS0" & own_type == "DRC3pHillFitModelFixS0"]
cat(sprintf("\n-- curves both call sigmoidal: %d of %d\n", nrow(sig), nrow(res)))
for (p in c("h", "x_inf", "ec50")) {
  r <- rel(sig[[paste0("own_", p)]], sig[[paste0("drc_", p)]])
  cat(sprintf("   %-6s median %.2e | p95 %.2e | max %.2e | within 1e-6: %5.1f%%\n",
              p, stats::median(r, na.rm = TRUE), stats::quantile(r, 0.95, na.rm = TRUE),
              max(r, na.rm = TRUE), 100 * mean(r < 1e-6, na.rm = TRUE)))
}
rx <- rel(sig$own_xc50, sig$drc_xc50)
cat(sprintf("   %-6s median %.2e | p95 %.2e | max %.2e | within 1e-6: %5.1f%% | within 1%%: %5.1f%%\n",
            "xc50", stats::median(rx, na.rm = TRUE), stats::quantile(rx, 0.95, na.rm = TRUE),
            max(rx, na.rm = TRUE), 100 * mean(rx < 1e-6, na.rm = TRUE), 100 * mean(rx < 1e-2, na.rm = TRUE)))
d_rss <- sig$own_rss - sig$drc_rss
cat(sprintf("   RSS    own better %d | tie %d | drc better %d | worst excess %.2e\n",
            sum(d_rss < -1e-12), sum(abs(d_rss) <= 1e-12), sum(d_rss > 1e-12), max(d_rss, na.rm = TRUE)))

sigl <- res[drc_type == "DRC3pHillFitModelFixS0" & own_log_type == "DRC3pHillFitModelFixS0"]
cat(sprintf("\n-- own_log: curves both call sigmoidal: %d\n", nrow(sigl)))
rxl <- rel(sigl$own_log_xc50, sigl$drc_xc50)
cat(sprintf("   xc50   median %.2e | p95 %.2e | within 0.1%%: %5.1f%% | within 1%%: %5.1f%%\n",
            stats::median(rxl, na.rm = TRUE), stats::quantile(rxl, 0.95, na.rm = TRUE),
            100 * mean(rxl < 1e-3, na.rm = TRUE), 100 * mean(rxl < 1e-2, na.rm = TRUE)))
dl <- sigl$own_log_rss - sigl$drc_rss
cat(sprintf("   RSS    own_log better %d | tie %d | drc better %d | worst excess %.2e | median excess %.2e\n",
            sum(dl < -1e-12), sum(abs(dl) <= 1e-12), sum(dl > 1e-12),
            max(dl, na.rm = TRUE), stats::median(dl, na.rm = TRUE)))

cat("\n-- curves where the backends disagree on fit_type\n")
dis <- res[drc_type != own_type]
print(dis[, .(dataset, norm_type, n_conc, drc_type, own_type,
              drc_p = signif(drc_p, 3), own_p = signif(own_p, 3),
              drc_r2 = signif(drc_r2, 3), own_r2 = signif(own_r2, 3))][1:min(10, nrow(dis))])

cat("\nwrote:", file.path(OUT_DIR, "poc_fit_comparison_v2.csv"), "\n")
