#' RDD UESD with conventional, bias‐corrected, or robust inference,
#' optional weights, and flexible placebo exclusions
#'
#' @param df A data.frame containing your data.
#' @param outcome String name of the numeric outcome variable.
#' @param running String name of the numeric running (forcing) variable.
#' @param cutoff Numeric discontinuity point on `running`. Defaults to 0.
#' @param bw Numeric half‐width of the RDD window (h). If `NULL` (default), uses rdrobust’s optimal bandwidth.
#' @param weight NULL (default) or a string naming a column in `df` of observation weights to pass to `rdrobust()`.
#' @param placebos Logical; if `TRUE`, runs placebo RDDs at every feasible cutoff (skipping those in `excluded`) and computes a placebo‐robust p‐value and CI.
#' @param excluded Numeric vector of length 2, `c(post_excl, pre_excl)`.
#'   Defaults to `c(Inf, 0)`, i.e. exclude **all** post‐cutoff placebo dates, none pre‐cutoff.
#' @param se_type Character; which rdrobust row to use for placebo inference. One of
#'   `"conventional"` (row 1), `"bias-corrected"` (row 2), or `"robust"` (row 3). Default: `"conventional"`.
#' @return Invisibly returns a list with:
#'   - `h_used`, `n_obs`, `se_type`
#'   - **All three** rdrobust estimates in `coef_all`, `se_all`, `z_all`, `p_all`, `ci_all`
#'   - **Selected** row statistics in `estimate`, `se`, `z_value`, `p_value`, `ci`
#'   - if `placebos=TRUE`: also `z_values`, `placebo_p_value`, `ci_placebo`
#' @importFrom rdrobust rdrobust
#' @importFrom stats pnorm qnorm quantile
#' @export

uesd_robust_rdd <- function(df,
                            outcome,
                            running,
                            cutoff    = 0,
                            bw        = NULL,
                            weight    = NULL,
                            placebos  = FALSE,
                            excluded  = c(Inf, 0),
                            se_type   = c("conventional","bias-corrected","robust")) {
  se_type <- match.arg(se_type)
  row_idx <- switch(se_type,
                    conventional    = 1,
                    "bias-corrected"= 2,
                    robust          = 3)

  # 1) Clean
  df_clean <- df[!is.na(df[[running]]), , drop = FALSE]

  # 2) Main rdrobust call
  rd_args <- list(
    y = df_clean[[outcome]],
    x = df_clean[[running]],
    c = cutoff
  )
  if (!is.null(bw))    rd_args$h       <- bw
  if (!is.null(weight)) rd_args$weights <- df_clean[[weight]]
  rdout <- do.call(rdrobust::rdrobust, rd_args)

  # 3) Bandwidth used (one‐sided)
  h_used <- rdout$bws[1]

  # 4) Print summary
  cat("\nRDD UESD for outcome:", outcome, "\n")
  if (!is.null(weight)) cat("Weights:", weight, "\n")
  summary(rdout); cat("\n")

  # 5) Extract all three rows
  coefs  <- rdout$coef
  ses    <- rdout$se
  z_all  <- coefs / ses
  p_all  <- 2*(1 - stats::pnorm(abs(z_all)))
  ci_mat <- rdout$ci
  n_main <- sum(abs(df_clean[[running]] - cutoff) <= h_used, na.rm = TRUE)

  # 6) Select the user‐requested row
  est_sel <- coefs[row_idx]
  se_sel  <- ses[row_idx]
  z_sel   <- z_all[row_idx]
  p_sel   <- p_all[row_idx]
  ci_sel  <- ci_mat[row_idx, , drop = TRUE]

  cat(sprintf("Using [%s] row for placebos: estimate=%.4f  se=%.4f  z=%.4f  p=%.4f\n\n",
              se_type, est_sel, se_sel, z_sel, p_sel))
  zcrit   <- stats::qnorm(0.975)
  ci_conv <- c(est_sel - zcrit*se_sel, est_sel + zcrit*se_sel)
  cat(sprintf("Unadjusted 95%% CI for ITT: [%.4f, %.4f]\n\n",
              ci_conv[1], ci_conv[2]))

  # 7) Placebo loop if requested
  z_vals <- numeric(0)
  if (placebos) {
    cuts <- sort(unique(c(df_clean[[date_num]], cutoff)))
    pos_main <- match(cutoff, cuts)
    post_excl <- excluded[1]
    pre_excl  <- excluded[2]

    # build post‐cutoff exclusion indices
    post_idx <- if (post_excl > 0) {
      if (is.infinite(post_excl)) {
        seq(pos_main+1, length(cuts))
      } else {
        seq(pos_main+1, pos_main+post_excl)
      }
    } else integer(0)

    # build pre‐cutoff exclusion indices
    pre_idx <- if (pre_excl > 0) {
      if (is.infinite(pre_excl)) {
        seq(1, pos_main-1)
      } else {
        seq(pos_main-pre_excl, pos_main-1)
      }
    } else integer(0)

    all_idx <- unique(c(post_idx, pre_idx))
    all_idx <- all_idx[all_idx >= 1 & all_idx <= length(cuts)]
    exclude_cuts <- cuts[all_idx]

    for (c0 in cuts) {
      if (c0 %in% exclude_cuts) next
      n_lo <- sum(df_clean[[running]] <  c0, na.rm = TRUE)
      n_hi <- sum(df_clean[[running]] >= c0, na.rm = TRUE)
      if (n_lo < 10 || n_hi < 10) next

      args_p <- list(
        y = df_clean[[outcome]],
        x = df_clean[[running]],
        c = c0,
        h = h_used
      )
      if (!is.null(weight)) args_p$weights <- df_clean[[weight]]

      rd_p <- tryCatch(
        suppressWarnings(do.call(rdrobust::rdrobust, args_p)),
        error = function(e) NULL
      )
      if (is.null(rd_p)) next
      co_p <- rd_p$coef; se_p <- rd_p$se
      if (length(co_p) < row_idx) next
      z_vals <- c(z_vals, co_p[row_idx] / se_p[row_idx])
    }

    placebo_p  <- mean(abs(z_vals) >= abs(z_sel))
    qz         <- stats::quantile(abs(z_vals), 0.975, na.rm = TRUE)
    ci_placebo <- c(est_sel - qz*se_sel, est_sel + qz*se_sel)

    cat(sprintf("Placebo‐robust p-value: %.4f\n",   placebo_p))
    cat(sprintf("Placebo‐robust 95%% CI:     [%.4f, %.4f]\n\n",
                ci_placebo[1], ci_placebo[2]))
  }

  # 8) Return invisibly
  out <- list(
    h_used    = h_used,
    n_obs     = n_main,
    se_type   = se_type,
    coef_all  = setNames(coefs,  c("conventional","bias-corrected","robust")),
    se_all    = setNames(ses,    c("conventional","bias-corrected","robust")),
    z_all     = setNames(z_all,  c("conventional","bias-corrected","robust")),
    p_all     = setNames(p_all,  c("conventional","bias-corrected","robust")),
    ci_all    = structure(ci_mat,
                          dimnames = list(
                            c("conventional","bias-corrected","robust"),
                            c("lower","upper")
                          )),
    estimate  = est_sel,
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
