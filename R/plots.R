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

#' Observed breakthrough curves and SMIM fits: linear–linear and log–log panels
#'
#' Follows the presentation of Volponi et al. (2025, Figs 2 and 4): a linear–linear panel
#' with the time axis starting at the release (t = 0), and a log–log panel that makes the
#' power-law tail visible. Observed data are drawn as open points and the fitted model as a
#' black line, one column per site. Legible in grayscale.
#'
#' @param obs data frame with columns `site`, `t_s`, `C` and logical `in_window`: all
#'   background-corrected observations from the release onward (µg L⁻¹, s); rows with
#'   `in_window = TRUE` are the fitted points. The linear panel shows every row; the
#'   log–log panel shows every positive row (`t > 0`, `C > log_floor`), so the low
#'   pre-arrival values fill the decades before the rising limb, as in the original
#'   figures.
#' @param fit data frame with columns `site`, `t_s`, `C`: the fitted model on a fine time
#'   grid, scaled to the observation units.
#' @param labels list from [btc_labels()].
#' @param log_floor concentrations at or below this are omitted from the log–log panel.
#' @param x_breaks,y_breaks major tick positions of the log–log panel.
#' @param x_logticks logical; draw log minor ticks on the x axis of the log–log panel.
#' @return list of two `ggplot` objects, `linear` and `loglog`.
#' @export
plot_btc_fits <- function(obs, fit, labels = btc_labels(), log_floor = 1e-3,
                          x_breaks = c(10, 100, 1000), y_breaks = 10^(-2:3), x_logticks = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required")
  if (is.null(obs$in_window)) obs$in_window <- TRUE
  base <- function(o, f) ggplot2::ggplot() +
    ggplot2::geom_point(data = o, ggplot2::aes(x = t_s, y = C), shape = 21, size = 1.4,
                        colour = "grey25", fill = "white", stroke = 0.5) +
    ggplot2::geom_line(data = f, ggplot2::aes(x = t_s, y = C), colour = "black", linewidth = 0.7) +
    ggplot2::facet_wrap(~site, nrow = 1, labeller = ggplot2::labeller(site = function(s) paste("Stream", s))) +
    ggplot2::labs(x = labels$time, y = labels$conc) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(strip.background = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(face = "bold"))
  fmt <- function(x) format(x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE)
  list(linear = base(obs[obs$t_s >= 0, ], fit) +
         ggplot2::scale_x_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0, 0.03))),
       loglog = base(obs[obs$t_s > 0 & obs$C > log_floor, ], fit[fit$C > log_floor & fit$t_s > 0, ]) +
         ggplot2::scale_x_log10(breaks = x_breaks, labels = fmt, limits = c(min(x_breaks), NA)) +
         ggplot2::scale_y_log10(breaks = y_breaks, labels = fmt) +
         (if (x_logticks) ggplot2::annotation_logticks(sides = "b") else NULL))
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
