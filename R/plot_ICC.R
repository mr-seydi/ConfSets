# Suppress global variable warnings for CRAN checks
utils::globalVariables(c("x", "Lower", "Upper", "ICC"))
#' @title Plot ICC with Simultaneous Confidence Bands
#'
#' @description This function generates a plot of the estimated ICC with simultaneous confidence bands based on Rank method.
#'
#' @param conf_sets_output The output of Conf_sets function when it the type = "Rank".
#' @param estimated_ICC_color The color of the estimated ICC line. The default value is "black".
#' @param bands_color The color of the confidence bands shaded area. The default value is "darkgrey".
#' @param plot_title The title of the plot. The default value is "Estimated ICC with Confidence Bands", it could be "" for having no title.
#' @param xlable The x-axis label. The default value is "Domain".
#' @param ylable The y-axis label. The default value is "ICC".
#' @param plot_legend_position The position of the legend. The default value is "top". it could be "none" for having no legend.
#'
#' @return A ggplot object.
#'
#' @examples
#' # The package contains the datasets, TEST_data, and RETEST_data which are the
#' #test and retest data, respectively.
#'
#'
#' set.seed(123)
#' confsets_outputs <- Conf_sets(test = TEST_data, retest = RETEST_data,
#'                              type = c("Rank", "Normal"), n_boot=1000,
#'                              designation = "ICC(A,1)", alpha = 0.05)
#'
#' plot_ICC(confsets_outputs, estimated_ICC_color = "black",
#' bands_color = "darkgrey",
#' plot_title = "Estimated ICC with Confidence Bands",
#' xlable = "Domain", ylable = "ICC",
#' plot_legend_position = "top")
#'
#' @export
plot_ICC <- function(conf_sets_output, estimated_ICC_color = "black",
                     bands_color = "darkgrey",
                     plot_title = "Estimated ICC with Confidence Bands",
                     xlable = "Domain", ylable = "ICC",
                     plot_legend_position = "top") {


  # Check if necessary elements exist in the input
  if (!all(c("pointwise_ICC", "CSets") %in% names(conf_sets_output)) ||
      !all(c("Rank") %in% names(conf_sets_output$CSets)) ||
      !all(c("upper", "lower") %in% names(conf_sets_output$CSets$Rank))) {
    stop("Invalid input structure: Ensure the output contains pointwise_ICC and CSets$Rank$upper/lower.")
  }

  # Create a data frame for plotting
  df <- data.frame(
    x = 0:(length(conf_sets_output$pointwise_ICC)-1),
    ICC = conf_sets_output$pointwise_ICC,
    Upper = conf_sets_output$CSets$Rank$upper,
    Lower = conf_sets_output$CSets$Rank$lower
  )

  # Generate the plot
  # Generate the plot
  ggplot2::ggplot(df, ggplot2::aes(x = x)) +
    # Shaded area for confidence bands
    ggplot2::geom_ribbon(ggplot2::aes(ymin = Lower, ymax = Upper, fill = "Confidence bands"), alpha = 0.3) +
    # Line for estimated ICC
    ggplot2::geom_line(ggplot2::aes(y = ICC, color = "Estimated ICC"), linewidth = 1) +
    # Custom legend and colors
    ggplot2::scale_color_manual(values = c("Estimated ICC" = estimated_ICC_color)) +
    ggplot2::scale_fill_manual(values = c("Confidence bands" = bands_color)) +
    # ylim (-1,1)
    ggplot2::coord_cartesian(ylim = c(-1, 1)) +
    ggplot2::labs(
      title = plot_title,
      x = xlable,
      y = ylable,
      color = "Legend",
      fill = "Legend"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = plot_legend_position,
      plot.title = ggplot2::element_text(hjust = 0.5),
      legend.title = ggplot2::element_blank()
    )
}
