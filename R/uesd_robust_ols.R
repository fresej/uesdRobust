#'  OLS UESD with robust or cluster‐robust SE and optional placebo inference
#'
#' @param df A data.frame containing your data.
#' @param outcome String name of the numeric outcome variable.
#' @param date_num String name of the numeric running date variable.
#' @param bw Numeric half‐width of the window around `cutoff`.
#' @param cutoff Numeric cutoff on `date_num` defining the treated (≥cutoff) vs control (<cutoff) groups. Defaults to 0.
#' @param cluster NULL (default) for heteroskedasticity‐robust SE, or a string naming a column in `df` on which to cluster.
#' @param controls NULL (default) or a character vector of other column-names in `df` to add as linear controls.
#' @param placebos Logical; if TRUE, runs regressions at every feasible cutoff (skipping `excluded` after the main), and computes a placebo‐robust p‐value and CI.
#' @param excluded Integer; number of subsequent cutoffs after the main one to skip in the placebo loop (default = 3).
#' @return Invisibly returns a list with elements:
#'  * `estimate`, `se`, `t_value`, `ci` (length‐2 numeric), `p_value`, `n_obs`
#'  * if `placebos = TRUE`: also `placebo_t_values`, `placebo_p_value`, `ci_placebo`
#' @importFrom stats lm coef qnorm pnorm quantile
#' @importFrom sandwich vcovHC vcovCL
#' @importFrom lmtest coeftest
#' @export

uesd_robust_ols <- function(df, outcome, date_num, bw = 15, cutoff = 0,
                                                        cluster = NULL,
                                                        controls    = NULL,
                                                        placebos    = FALSE,
                                                        excluded    = 3) {
  # 0) Drop any NA in the running variable
  df_clean <- df[!is.na(df[[date_num]]), , drop = FALSE]

  # 1) Main window subset and ITT variable
  dat_main <- df_clean[abs(df_clean[[date_num]] - cutoff) <= bw, , drop = FALSE]
  if (nrow(dat_main) == 0) {
    stop("No observations within the main window. Increase `bw` or adjust `cutoff`.")
  }
  dat_main$ITT <- ifelse(dat_main[[date_num]] >= cutoff, 1L, 0L)

  # 2) Build regression formula, including controls if requested
  if (is.null(controls) || length(controls) == 0) {
    fml <- as.formula(paste(outcome, "~ ITT"))
  } else {
    ctrl_terms <- paste(controls, collapse = " + ")
    fml <- as.formula(paste(outcome, "~ ITT +", ctrl_terms))
  }


  fit_lm <- stats::lm(fml, data = dat_main)

  # 3) Covariance for main
  if (is.null(cluster)) {
    V_main  <- sandwich::vcovHC(fit_lm, type = "HC1")
    se_type <- "Heteroskedasticity‐robust"
  } else {
    V_main  <- sandwich::vcovCL(fit_lm, cluster = dat_main[[cluster]], type = "HC1")
    se_type <- paste0("Cluster‐robust (by ", cluster, ")")
  }

  # 4) Extract main stats
  est_main <- stats::coef(fit_lm)["ITT"]
  se_main  <- sqrt(V_main["ITT","ITT"])
  t_main   <- est_main / se_main
  z        <- stats::qnorm(0.975)
  ci_main  <- c(est_main - z * se_main, est_main + z * se_main)
  p_main   <- 2 * (1 - stats::pnorm(abs(t_main)))
  n_main   <- nrow(dat_main)

  # 5) Print main results
  cat("\nOLS UESD for outcome:", outcome, "\n")
  cat("  Window: [", cutoff - bw, ",", cutoff + bw,
      "]  (n =", n_main, ")\n")
  cat("  SE type:", se_type, "\n\n")
  print(lmtest::coeftest(fit_lm, vcov. = V_main))
  cat(sprintf("\nUnadjusted 95%% CI for ITT: [%.4f, %.4f]\n", ci_main[1], ci_main[2]))

  # 6) Placebo loop if requested
  placebo_p      <- NULL
  placebo_t_vals <- NULL
  ci_placebo     <- NULL
  if (placebos) {
    cuts <- sort(unique(df_clean[[date_num]]))
    # determine which cutoffs to exclude
    pos_main <- match(cutoff, cuts)
    if (!is.na(pos_main) && excluded > 0) {
      exclude_inds <- seq(pos_main + 1, pos_main + excluded)
      exclude_inds <- exclude_inds[exclude_inds <= length(cuts)]
      exclude_cuts <- cuts[exclude_inds]
    } else {
      exclude_cuts <- numeric(0)
    }

    t_vals <- numeric(0)
    for (c0 in cuts) {
      if (c0 %in% exclude_cuts) next
      n_lower <- sum(df_clean[[date_num]] <  c0, na.rm = TRUE)
      n_upper <- sum(df_clean[[date_num]] >= c0, na.rm = TRUE)
      if (n_lower < 10 || n_upper < 10) next

      dat_p <- df_clean[abs(df_clean[[date_num]] - c0) <= bw, , drop = FALSE]
      dat_p$ITT <- ifelse(dat_p[[date_num]] >= c0, 1L, 0L)

      fit_p <- tryCatch(stats::lm(fml, data = dat_p), error = function(e) NULL)
      if (is.null(fit_p)) next
      if (!"ITT" %in% names(stats::coef(fit_p))) next

      Vp <- if (is.null(cluster)) {
        sandwich::vcovHC(fit_p, type = "HC1")
      } else {
        sandwich::vcovCL(fit_p, cluster = dat_p[[cluster]], type = "HC1")
      }
      if (!("ITT" %in% rownames(Vp) && "ITT" %in% colnames(Vp))) next

      se_p  <- sqrt(Vp["ITT","ITT"])
      est_p <- stats::coef(fit_p)["ITT"]
      t_vals <- c(t_vals, est_p / se_p)
    }

    placebo_t_vals <- t_vals
    placebo_p      <- mean(abs(t_vals) >= abs(t_main))

    # Symmetric placebo‐robust CI using the 97.5% quantile of |t|
    q_abs <- stats::quantile(abs(t_vals), 0.975, na.rm = TRUE)
    ci_placebo_lower <- est_main - q_abs * se_main
    ci_placebo_upper <- est_main + q_abs * se_main
    ci_placebo       <- c(lower = ci_placebo_lower, upper = ci_placebo_upper)

    cat(sprintf("Placebo‐robust 95%% CI for ITT: [%.4f, %.4f]\n\n",
                ci_placebo_lower, ci_placebo_upper))
    cat(sprintf("Unadjusted p-value for ITT: %.5f\n", p_main))
    cat(sprintf("Placebo‐robust p-value for ITT: %.5f\n", placebo_p))
  }

  # 7) Return invisibly, now including t_value and placebo CI
  out <- list(
    estimate = est_main,
    se        = se_main,
    t_value   = t_main,
    ci        = setNames(ci_main, c("lower","upper")),
    p_value   = p_main,
    n_obs     = n_main
  )
  if (placebos) {
    out$placebo_t_values <- placebo_t_vals
    out$placebo_p_value  <- placebo_p
    out$ci_placebo       <- ci_placebo
  }
  invisible(out)
}
