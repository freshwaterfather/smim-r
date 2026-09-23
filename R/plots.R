utils::globalVariables(c("t_s", "C", "site"))

#' Axis labels for breakthrough-curve figures (plotmath)
#'
#' Units are taken from the data (the package's example data are in µg L⁻¹ and s).
#' @return list `conc`, `time` of plotmath expressions.
#' @export
btc_labels <- function() {
  list(conc = expression(paste("Rhodamine WT concentration (", mu, "g L"^-1, ")")),
       time = "Time since release (s)")
}

#' Figure theme following the Methods in Stream Ecology figure guidelines
#'
#' Arial lettering, 14-pt axis numbers, 16-pt axis labels and bold axis lines, on top of
#' [ggplot2::theme_classic()].
#' @param family font family (default "Arial"; Helvetica and Calibri are also accepted by
#'   the publisher).
#' @return a `ggplot2` theme.
#' @export
theme_mise <- function(family = "Arial") {
  ggplot2::theme_classic(base_size = 14, base_family = family) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 14, colour = "black"),
                   axis.title = ggplot2::element_text(size = 16),
                   axis.line = ggplot2::element_line(linewidth = 0.9, colour = "black"),
                   axis.ticks = ggplot2::element_line(linewidth = 0.7, colour = "black"),
                   axis.ticks.length = ggplot2::unit(4, "pt"),
                   strip.background = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(size = 14, face = "bold"))
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
    theme_mise() +
    ggplot2::theme(panel.spacing = ggplot2::unit(0.8, "lines"))
  fmt <- function(x) format(x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE)
  list(linear = base(obs[obs$t_s >= 0, ], fit) +
         ggplot2::scale_x_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0, 0.03))),
       loglog = base(obs[obs$t_s > 0 & obs$C > log_floor, ], fit[fit$C > log_floor & fit$t_s > 0, ]) +
         ggplot2::scale_x_log10(breaks = x_breaks, labels = fmt, limits = c(min(x_breaks), NA)) +
         ggplot2::scale_y_log10(breaks = y_breaks, labels = fmt) +
         (if (x_logticks) ggplot2::annotation_logticks(sides = "b") else NULL))
}

#' Save a figure as PDF, 600-dpi PNG and 600-dpi TIFF
#'
#' @param p a ggplot object. @param path_base file path without extension.
#' @param width,height inches.
#' @return invisibly, the vector of files written.
#' @export
save_btc_figure <- function(p, path_base, width = 7.2, height = 3.2) {
  ggplot2::ggsave(paste0(path_base, ".pdf"), p, width = width, height = height, device = grDevices::cairo_pdf)
  ggplot2::ggsave(paste0(path_base, ".png"), p, width = width, height = height, dpi = 600, bg = "white")
  ggplot2::ggsave(paste0(path_base, ".tiff"), p, width = width, height = height, dpi = 600, bg = "white", compression = "lzw")
  invisible(paste0(path_base, c(".pdf", ".png", ".tiff")))
}

#' Fitted SMIM parameters per stream as a four-panel bar figure
#'
#' One bar per stream for each of the four fitted parameters, in the layout and with the
#' axis labels of Volponi et al. (2025, Fig. 5), velocity written as V: velocity (top left), dispersion (top
#' right), exchange rate (bottom left) and power-law slope (bottom right). Thin error bars
#' show ± one standard error when the `se_*` columns are present.
#'
#' @param results data frame with columns `site`, `v_m_s`, `D_m2_s`, `Lambda_1_s`, `beta`
#'   and optionally `se_v`, `se_D`, `se_Lambda`, `se_beta` (as written by
#'   `analysis/example_four_streams.R`).
#' @param fill bar fill colour (grayscale-safe default).
#' @return a `ggplot` object (2 × 2 facets with the parameter label on the y axis).
#' @export
plot_smim_params <- function(results, fill = "grey55") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required")
  pars <- c(v_m_s = "Velocity~'['*italic(V)*', m/s'*']'",
            D_m2_s = "Dispersion~'['*italic(D)*', m'^2*'/s'*']'",
            Lambda_1_s = "Exchange~Rate~'['*Lambda*', 1/s'*']'",
            beta = "Power~Law~Slope~'['*beta*']'")
  ses <- c(v_m_s = "se_v", D_m2_s = "se_D", Lambda_1_s = "se_Lambda", beta = "se_beta")
  long <- do.call(rbind, lapply(names(pars), function(k) data.frame(
    site = results$site, parameter = pars[[k]], value = results[[k]],
    se = if (ses[[k]] %in% names(results)) results[[ses[[k]]]] else NA_real_)))
  long$parameter <- factor(long$parameter, levels = unname(pars))
  p <- ggplot2::ggplot(long, ggplot2::aes(x = site, y = value)) +
    ggplot2::geom_col(fill = fill, colour = "black", linewidth = 0.3, width = 0.65) +
    ggplot2::facet_wrap(~parameter, nrow = 2, scales = "free_y", strip.position = "left",
                        labeller = ggplot2::label_parsed) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.08))) +
    ggplot2::labs(x = "Stream", y = NULL) +
    theme_mise() +
    ggplot2::theme(strip.placement = "outside",
                   strip.text.y.left = ggplot2::element_text(angle = 90, size = 16, face = "plain"),
                   panel.spacing = ggplot2::unit(1, "lines"))
  if (any(!is.na(long$se))) {
    p <- p + ggplot2::geom_errorbar(ggplot2::aes(ymin = value - se, ymax = value + se), width = 0.2, linewidth = 0.4)
  }
  p
}
