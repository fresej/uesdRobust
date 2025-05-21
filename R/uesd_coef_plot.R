#' Plot original vs. placebo‐robust UESD results
#'
#' Given the list returned (invisibly) by `uesd_robust_ols(..., placebos=TRUE)` or
#' `uesd_robust_rdd(..., placebos=TRUE)`, draws a side‐by‐side coefficient plot
#' showing the original 95% CI and the placebo‐robust 95% CI, plus a horizontal
#' reference line at zero.
#'
#' @param res A list returned by one of the UESD functions with components:
#'   * `estimate`: numeric, the point estimate
#'   * `ci`: named length‐2 numeric for the conventional 95% CI
#'     (names must include “lower” and “upper” somewhere)
#'   * `ci_placebo`: named length‐2 numeric for the placebo‐robust 95% CI
#'     (names must include “lower” and “upper” somewhere)
#' @return Invisibly returns the ggplot object.
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar geom_hline scale_x_discrete labs theme_classic
#' @export

uesd_coef_plot <- function(res) {
  if (is.null(res$estimate) ||
      is.null(res$ci) ||
      is.null(res$ci_placebo)) {
    stop("`res` must contain `estimate`, `ci`, and `ci_placebo` elements.")
  }
  # helper to pull out lower/upper regardless of exact name
  get_bound <- function(x, which = c("lower","upper")) {
    nm <- names(x)
    idx <- grep(which, nm, ignore.case = TRUE)
    if (length(idx)==0) stop("Cannot find '", which, "' in names: ", paste(nm,collapse=","))
    x[idx[1]]
  }
  est    <- res$estimate
  main_lo <- get_bound(res$ci,         "lower")
  main_hi <- get_bound(res$ci,         "upper")
  pl_lo   <- get_bound(res$ci_placebo, "lower")
  pl_hi   <- get_bound(res$ci_placebo, "upper")

  df_plot <- data.frame(
    Type     = factor(c("Original","Placebo-robust"),
                      levels = c("Original","Placebo-robust")),
    Estimate = est,
    CI_lower = c(main_lo, pl_lo),
    CI_upper = c(main_hi, pl_hi),
    row.names = NULL
  )

  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Type, y = Estimate)) +
    ggplot2::geom_hline(yintercept = 0,
                        linetype   = "dashed",
                        color      = "red",
                        size       = 1) +
    ggplot2::geom_point(size = 6) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = CI_lower, ymax = CI_upper),
      width = 0,
      size  = 2
    ) +
    ggplot2::scale_x_discrete(expand = c(0.2, 0)) +
    ggplot2::labs(
      x     = "",
      y     = "Estimate",
      title = "Original vs. Placebo-Robust 95% CI"
    ) +
    ggplot2::theme_classic(base_size = 20)

  print(p)
  invisible(p)
}
