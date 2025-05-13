#' Create a LaTeX table of UESD estimates (original vs. placebo‐robust)
#'
#' Takes one or more result lists from `uesd_robust_ols(...,placebos=TRUE)` or
#' `uesd_robust_rdd(...,placebos=TRUE)` and produces a LaTeX table with columns
#' for each test and rows for:
#'   1. Coefficient (SE)
#'   2. Unadjusted 95\\% CI
#'   3. Unadjusted p‐value
#'   4. Placebo‐robust 95\\% CI
#'   5. Placebo‐robust p‐value
#'   6. Sample size
#'   7. Estimator (“OLS” or “RDD”)
#'   8. SE type (e.g. “Heteroskedasticity‐robust”, “clustered”, or rdrobust row)
#'
#' @param res One of:
#'   * a single result‐list (with `estimate`,`se`,`ci`,`p_value`,
#'     `ci_placebo`,`placebo_p_value`,`n_obs`, and `se_type`),
#'   * a list of such result‐lists, or
#'   * a character vector of object names in your calling environment.
#' @param digits Integer: number of decimal places for non‐p‐value entries (default: 3).
#' @param digits_p Integer: number of decimal places for p‐value entries (default: 5).
#' @return Invisibly returns the LaTeX code (character vector).
#' @importFrom knitr kable
#' @export

uesd_table <- function(res, digits = 3, digits_p = 5) {
  # helper to pull out lower/upper regardless of exact name
  get_bound <- function(x, which) {
    nm  <- names(x)
    idx <- grep(which, nm, ignore.case = TRUE)
    if (length(idx) == 0) {
      stop("Cannot find '", which, "' in names: ", paste(nm, collapse = ","))
    }
    unname(x[idx[1]])
  }

  # normalize res into a named list of result‐lists
  if (is.character(res)) {
    res_list <- mget(res, envir = parent.frame())
  } else if (
    is.list(res) &&
    all(c("estimate","se","ci","p_value",
          "ci_placebo","placebo_p_value","n_obs") %in% names(res))
  ) {
    res_list <- list(Test1 = res)
  } else if (
    is.list(res) &&
    all(sapply(res, function(r)
      is.list(r) &&
      all(c("estimate","se","ci","p_value",
            "ci_placebo","placebo_p_value","n_obs") %in% names(r))
    ))
  ) {
    res_list <- res
  } else {
    stop("`res` must be a single result‐list, a list of them, or a character vector.")
  }

  # ensure names
  if (is.null(names(res_list)) || any(names(res_list) == "")) {
    names(res_list) <- paste0("Test", seq_along(res_list))
  }

  tests   <- names(res_list)
  metrics <- c(
    "Coefficient (SE)",
    "Unadjusted 95\\% CI",
    "Unadjusted p‐value",
    "Placebo‐robust 95\\% CI",
    "Placebo‐robust p‐value",
    "Sample size",
    "Estimator",
    "SE type"
  )

  # build empty matrix
  m <- matrix("", nrow = length(metrics), ncol = length(tests),
              dimnames = list(metrics, tests))

  # formatters
  fmt_gen <- function(x) {
    format(round(x, digits), nsmall = digits)
  }
  fmt_p <- function(x) {
    tiny <- paste0("0.", strrep("0", digits_p-1), "1")
    if (x < 10^(-digits_p)) {
      paste0("<", tiny)
    } else {
      format(round(x, digits_p), nsmall = digits_p)
    }
  }
  # significance stars
  star <- function(p) {
    if      (p < 0.001) "***"
    else if (p < 0.01)  "**"
    else if (p < 0.05)  "*"
    else                ""
  }

  for (j in seq_along(tests)) {
    r         <- res_list[[j]]
    est       <- r$estimate
    se        <- r$se
    ci_lo     <- get_bound(r$ci,         "lower")
    ci_hi     <- get_bound(r$ci,         "upper")
    p_unadj   <- r$p_value
    ci_pl_lo  <- get_bound(r$ci_placebo, "lower")
    ci_pl_hi  <- get_bound(r$ci_placebo, "upper")
    p_pl      <- r$placebo_p_value
    n_obs     <- r$n_obs

    estimator <- if ("t_value" %in% names(r)) "OLS" else "RDD"
    se_type   <- if (!is.null(r$se_type)) {
      r$se_type
    } else if (estimator == "OLS") {
      "Heteroskedasticity-robust"
    } else {
      "conventional"
    }

    m[1, j] <- paste0(fmt_gen(est), " (", fmt_gen(se), ")")
    m[2, j] <- paste0("[", fmt_gen(ci_lo), ", ", fmt_gen(ci_hi), "]")
    m[3, j] <- paste0(fmt_p(p_unadj), star(p_unadj))
    m[4, j] <- paste0("[", fmt_gen(ci_pl_lo), ", ", fmt_gen(ci_pl_hi), "]")
    m[5, j] <- paste0(fmt_p(p_pl), star(p_pl))
    m[6, j] <- as.character(n_obs)
    m[7, j] <- estimator
    m[8, j] <- se_type
  }

  # render LaTeX table
  tab <- knitr::kable(
    m,
    format    = "latex",
    booktabs  = TRUE,
    escape    = FALSE,
    caption   = "UESD Estimates: original vs. placebo-robust",
    label     = "tab:uesd_estimates"
  )
  cat(tab, "\n")
  invisible(tab)
}
