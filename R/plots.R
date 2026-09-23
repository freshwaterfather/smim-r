utils::globalVariables(c("t_s", "C", "site"))

#' Axis labels for breakthrough-curve figures (plotmath for static, HTML for plotly)
#'
#' Units are taken from the data (the package's example data are in µg L⁻¹ and s).
#' @param html logical; return HTML markup (for `plotly`, which ignores plotmath) instead
#'   of plotmath expressions.
#' @return list `conc`, `time`.
#' @export
btc_labels <- function(html = FALSE) {
  if (html) {
    list(conc = "Rhodamine WT concentration, <i>C</i> (µg L<sup>−1</sup>)",
         time = "Time since release, <i>t</i> (s)")
  } else {
    list(conc = expression(paste("Rhodamine WT concentration, ", italic(C), " (", mu, "g L"^-1, ")")),
         time = expression(paste("Time since release, ", italic(t), " (s)")))
  }
}

#' Observed breakthrough curves and SMIM fits, linear and log panels
#'
#' Observed data are drawn as open points and the fitted model as a black line, one
#' column per site, in two figures (linear and log10 concentration) so the tail fit is
#' visible. Legible in grayscale.
#'
#' @param obs data frame with columns `site`, `t_s`, `C`: the corrected, background-
#'   subtracted observations inside the fit window (µg L⁻¹, s).
#' @param fit data frame with columns `site`, `t_s`, `C`: the fitted model on a fine time
#'   grid, scaled to the observation units.
#' @param labels list from [btc_labels()].
#' @param log_floor concentrations at or below this are omitted from the log panel.
#' @return list of two `ggplot` objects, `linear` and `log`.
#' @export
plot_btc_fits <- function(obs, fit, labels = btc_labels(), log_floor = 1e-3) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required")
  base <- function(o, f) ggplot2::ggplot() +
    ggplot2::geom_point(data = o, ggplot2::aes(x = t_s, y = C), shape = 21, size = 1.4,
                        colour = "grey25", fill = "white", stroke = 0.5) +
    ggplot2::geom_line(data = f, ggplot2::aes(x = t_s, y = C), colour = "black", linewidth = 0.7) +
    ggplot2::facet_wrap(~site, nrow = 1, labeller = ggplot2::labeller(site = function(s) paste("Stream", s))) +
    ggplot2::labs(x = labels$time, y = labels$conc) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(strip.background = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(face = "bold"))
  list(linear = base(obs, fit),
       log = base(obs[obs$C > log_floor, ], fit[fit$C > log_floor, ]) +
         ggplot2::scale_y_log10(labels = function(x) format(x, scientific = FALSE, drop0trailing = TRUE)) +
         ggplot2::annotation_logticks(sides = "l"))
}

#' Save a figure as PDF, 600-dpi PNG and TIFF, and an interactive HTML (plotly)
#'
#' The static files use the plotmath labels of `p`; the HTML version is rebuilt with
#' `plotly::ggplotly()` and the HTML labels from [btc_labels()], because plotly does not
#' render plotmath.
#' @param p a ggplot object. @param path_base file path without extension.
#' @param width,height inches. @param html_labels list from `btc_labels(html = TRUE)`.
#' @return invisibly, the vector of files written.
#' @export
save_btc_figure <- function(p, path_base, width = 7.2, height = 3.2, html_labels = btc_labels(html = TRUE)) {
  files <- character()
  ggplot2::ggsave(paste0(path_base, ".pdf"), p, width = width, height = height, device = grDevices::cairo_pdf)
  ggplot2::ggsave(paste0(path_base, ".png"), p, width = width, height = height, dpi = 600, bg = "white")
  ggplot2::ggsave(paste0(path_base, ".tiff"), p, width = width, height = height, dpi = 600, bg = "white", compression = "lzw")
  files <- paste0(path_base, c(".pdf", ".png", ".tiff"))
  if (requireNamespace("plotly", quietly = TRUE) && requireNamespace("htmlwidgets", quietly = TRUE)) {
    ph <- p + ggplot2::labs(x = html_labels$time, y = html_labels$conc)
    w <- plotly::ggplotly(ph)
    w <- plotly::layout(w, xaxis = list(title = list(text = html_labels$time)),
                        yaxis = list(title = list(text = html_labels$conc)))
    htmlwidgets::saveWidget(w, paste0(path_base, ".html"), selfcontained = TRUE)
    files <- c(files, paste0(path_base, ".html"))
  }
  invisible(files)
}
