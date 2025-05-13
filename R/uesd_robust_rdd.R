#'  RDD UESD with conventional, bias-corrected, or robust inference and optional placebo tests
#'
#' @param df A data.frame containing your data.
#' @param outcome String name of the numeric outcome variable.
#' @param running String name of the numeric running (forcing) variable.
#' @param cutoff Numeric discontinuity point on `running`. Defaults to 0.
#' @param bw Numeric half-width of the RDD window (h). If NULL (default), uses rdrobust’s optimal bandwidth.
#' @param placebos Logical; if TRUE, runs placebo RDDs at every feasible cutoff (skipping `excluded` after the main) and computes a placebo-robust p-value and CI.
#' @param excluded Integer; number of subsequent cutoffs after the main one to skip in the placebo loop (default = 3).
#' @param se_type Character; which rdrobust row to use for placebo inference. One of
#'    "conventional" (row 1), "bias-corrected" (row 2), or "robust" (row 3).  Default: "conventional".
#' @return Invisibly returns a list with:
#'    - h_used, n_obs, se_type
#'    - **All three** rdrobust estimates in coef_all, se_all, z_all, p_all, ci_all  (named by row)
#'     - **Selected** row statistics in estimate, se, z_value, p_value, ci
#'    - if placebos = TRUE, also z_values, placebo_p_value, ci_placebo
#' @importFrom rdrobust rdrobust
#' @importFrom stats pnorm qnorm quantile
#' @export

uesd_robust_rdd <- function(df, outcome, running, cutoff = 0,
                            bw = NULL,
                            placebos = FALSE,
                            excluded = 3,
                            se_type = c("conventional","bias-corrected","robust")) {

  se_type <- match.arg(se_type)
  # map to row index
  row_idx <- switch(se_type,
                    conventional   = 1,
                    "bias-corrected" = 2,
                    robust          = 3)

  # 1) Drop NAs in running variable
  df_clean <- df[!is.na(df[[running]]), , drop = FALSE]

  # 2) Main rdrobust call
  rdout <- if (is.null(bw)) {
    rdrobust::rdrobust(y = df_clean[[outcome]],
                       x = df_clean[[running]],
                       c = cutoff)
  } else {
    rdrobust::rdrobust(y = df_clean[[outcome]],
                       x = df_clean[[running]],
                       c = cutoff,
                       h = bw)
  }

  # 3) Report bandwidth used (one‐sided)
  h_used <- rdout$bws[1]

  # 4) Print full summary (includes all three rows)

  cat("\n RDD UESD for outcome:", outcome, "\n")
  summary(rdout)
  cat("\n")

  # 5) Extract all three main estimates
  coefs   <- rdout$coef        # length-3
  ses     <- rdout$se
  z_all   <- coefs / ses
  p_all   <- 2*(1 - stats::pnorm(abs(z_all)))
  ci_mat  <- rdout$ci          # 3×2 matrix, columns = lower,upper

  # total obs in window
  n_main  <- sum(abs(df_clean[[running]] - cutoff) <= h_used, na.rm = TRUE)

  # 6) Now pick the user‐requested row for inference
  est_sel     <- coefs[row_idx]
  se_sel      <- ses[row_idx]
  z_sel       <- z_all[row_idx]
  p_sel       <- p_all[row_idx]
  ci_sel      <- ci_mat[row_idx, , drop = TRUE]  # named length-2

  cat(sprintf("Using [%s] row for placebos:   estimate=%.4f   se=%.4f   z=%.4f   p=%.4f\n\n",
              se_type, est_sel, se_sel, z_sel, p_sel))

  # 7) Show conventional 95% CI (z = 1.96) for the selected row
  zcrit <- stats::qnorm(0.975)
  ci_conv <- c(est_sel - zcrit*se_sel,
               est_sel + zcrit*se_sel)
  cat(sprintf("Unadjusted 95%% CI for ITT: [%.4f, %.4f]\n\n",
              ci_conv[1], ci_conv[2]))

  # 8) Placebo loop if requested
  z_vals     <- NULL
  placebo_p  <- NULL
  ci_placebo <- NULL

  if (placebos) {
    cuts     <- sort(unique(df_clean[[running]]))
    pos_main <- match(cutoff, cuts)

    # compute exclusion set
    if (!is.na(pos_main) && excluded > 0) {
      exc_idx    <- seq(pos_main + 1, pos_main + excluded)
      exc_idx    <- exc_idx[exc_idx <= length(cuts)]
      exclude_cuts <- cuts[exc_idx]
    } else {
      exclude_cuts <- numeric(0)
    }

    for (c0 in cuts) {
      if (c0 %in% exclude_cuts) next
      n_lo <- sum(df_clean[[running]] <  c0, na.rm = TRUE)
      n_hi <- sum(df_clean[[running]] >= c0, na.rm = TRUE)
      if (n_lo < 10 || n_hi < 10) next

      rd_p <- tryCatch({
        suppressWarnings(
          rdrobust::rdrobust(
            y = df_clean[[outcome]],
            x = df_clean[[running]],
            c = c0,
            h = h_used
          )
        )
      }, error = function(e) {
        NULL
      })
      if (is.null(rd_p)) next

      # extract same row
      coefs_p <- rd_p$coef
      ses_p   <- rd_p$se
      if (length(coefs_p) < row_idx) next  # guard
      z_vals  <- c(z_vals, coefs_p[row_idx] / ses_p[row_idx])
    }

    # placebo‐robust p‐value
    placebo_p <- mean(abs(z_vals) >= abs(z_sel))

    # symmetric placebo CI: quantile of |z|
    qz <- stats::quantile(abs(z_vals), 0.975, na.rm = TRUE)
    ci_placebo <- c(est_sel - qz*se_sel,
                    est_sel + qz*se_sel)
    cat(sprintf("Placebo‐robust 95%% CI for ITT: [%.4f, %.4f]\n\n",
                ci_placebo[1], ci_placebo[2]))
    cat(sprintf("Unadjusted p-value for ITT: %.5f\n",   p_sel))
    cat(sprintf("Placebo‐robust p-value for ITT: %.5f\n",   placebo_p))

  }

  # 9) Return invisibly
  out <- list(
    h_used     = h_used,
    n_obs      = n_main,
    se_type    = se_type,
    coef_all   = setNames(coefs,   c("conventional","bias-corrected","robust")),
    se_all     = setNames(ses,     c("conventional","bias-corrected","robust")),
    z_all      = setNames(z_all,   c("conventional","bias-corrected","robust")),
    p_all      = setNames(p_all,   c("conventional","bias-corrected","robust")),
    ci_all     = structure(ci_mat,
                           dimnames = list(c("conventional","bias-corrected","robust"),
                                           c("lower","upper"))),
    # selected
    estimate   = est_sel,
    se         = se_sel,
    z_value    = z_sel,
    p_value    = p_sel,
    ci         = setNames(ci_sel, c("lower","upper"))
  )
  if (placebos) {
    out$z_values        <- z_vals
    out$placebo_p_value <- placebo_p
    out$ci_placebo      <- setNames(ci_placebo, c("lower","upper"))
  }

  invisible(out)
}
