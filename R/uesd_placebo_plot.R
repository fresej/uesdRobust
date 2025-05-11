#'  Plot the placebo distribution of t‐ or z‐statistics
#'
#'  Given the output from `uesd_robust_ols(..., placebos = TRUE)` or
#'  `uesd_robust_rdd(..., placebos = TRUE)`, draws a histogram of the
#'  placebo t- or z-values, highlights in red those exceeding the main
#'  statistic, and marks the main value with a vertical line.
#'
#' @param res A list returned (invisibly) by `uesd_robust_ols(..., placebos = TRUE)`
#'             or `uesd_robust_rdd(..., placebos = TRUE)`.
#' @param binwidth Numeric. Width of the histogram bins (default = 0.1).
#' @return Invisibly returns the ggplot object.
#' @importFrom ggplot2 ggplot aes geom_histogram scale_fill_manual geom_vline labs theme_classic
#' @export

uesd_placebo_plot <- function(res, binwidth = 0.1) {
  # Decide whether we have t‐values or z‐values
  if (!is.null(res$z_values)) {
    stat_vec <- res$z_values
    main_stat <- res$z_value
    stat_name <- "Z-values"
    # read out the se_type if present
    subtitle_prefix <- if (!is.null(res$se_type)) {
      paste0("RDD (", res$se_type, ") Main Z-value = ")
    } else {
      "Main z = "
    }
  } else if (!is.null(res$placebo_t_values)) {
    stat_vec <- res$placebo_t_values
    main_stat <- res$t_value
    stat_name <- "T-values"
    subtitle_prefix <- "OLS Main T-value = "
  } else {
    stop("`res` must contain either `$placebo_t_values` or `$z_values`.")
  }

  if (length(stat_vec) == 0) {
    stop("No placebo statistics to plot.")
  }

  df <- data.frame(stat = stat_vec)
  df$exceed <- abs(df$stat) >= abs(main_stat)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = stat, fill = exceed)) +
    ggplot2::geom_histogram(binwidth = binwidth, color = "black") +
    ggplot2::scale_fill_manual(
      values = c("FALSE" = "grey", "TRUE" = "red"),
      guide  = FALSE
    ) +
    ggplot2::geom_vline(
      xintercept = main_stat,
      color      = "red",
      size       = 1
    ) +
    ggplot2::labs(
      x        = paste("Placebo", stat_name),
      y        = "Count",
      title    = paste("Placebo Distribution of", stat_name),
      subtitle = paste0(
        subtitle_prefix,
        round(main_stat, 3),
        if (!is.null(res$placebo_p_value)) {
          paste0("\nPlacebo-robust p-value = ", round(res$placebo_p_value, 3))
        } else ""
      )
    ) +
    ggplot2::theme_classic(base_size = 20)

  print(p)
  invisible(p)
}
