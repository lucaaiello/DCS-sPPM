#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_argument)) {
  normalizePath(sub("^--file=", "", script_argument[[1L]]), mustWork = TRUE)
} else {
  normalizePath("Manuscript_figures/scripts/generate_conceptual_scheme.R", mustWork = TRUE)
}

project_root <- dirname(script_path)
while (!file.exists(file.path(project_root, "R", "project_paths.R"))) {
  parent <- dirname(project_root)
  if (identical(parent, project_root)) stop("Could not locate the manuscript project.")
  project_root <- parent
}
figure_dir <- file.path(project_root, "Luc_i_a", "fig")
asset_dir <- file.path(project_root, "Manuscript_figures", "results")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(asset_dir, recursive = TRUE, showWarnings = FALSE)

pdf_file <- file.path(figure_dir, "dcs_sppm_conceptual_scheme.pdf")
png_file <- file.path(figure_dir, "dcs_sppm_conceptual_scheme.png")
caption_file <- file.path(asset_dir, "dcs_sppm_conceptual_scheme_caption.txt")

colour <- list(
  navy = "#173F62", ink = "#18324A", muted = "#566674",
  border = "#C9D5DE", blue = "#2A83B8", blue_light = "#E4F2FA",
  green = "#3D8B57", green_light = "#E9F5E9",
  red = "#C9574D", red_light = "#FCE9E7",
  purple = "#72509A", purple_light = "#F1EAF8",
  grey_light = "#F4F7F9"
)

transform_xy <- function(data, xmin, xmax, ymin, ymax) {
  transform(data, x = xmin + x * (xmax - xmin), y = ymin + y * (ymax - ymin))
}

region <- rbind(
  data.frame(group = "Regime 1", x = c(.02, .10, .32, .45, .38, .21),
             y = c(.27, .76, .96, .63, .35, .05), fill = "#1895C3"),
  data.frame(group = "Regime 2", x = c(.32, .64, .68, .61, .45),
             y = c(.96, .91, .70, .57, .63), fill = "#74B66A"),
  data.frame(group = "Regime 3", x = c(.21, .38, .45, .61, .73, .50),
             y = c(.05, .35, .63, .57, .18, .01), fill = "#F1BE3E"),
  data.frame(group = "Regime 4", x = c(.64, .91, .99, .92, .73, .61, .68),
             y = c(.91, .75, .48, .25, .18, .57, .70), fill = "#DE574B")
)
outline <- data.frame(
  x = c(.02, .10, .32, .64, .91, .99, .92, .73, .50, .21),
  y = c(.27, .76, .96, .91, .75, .48, .25, .18, .01, .05)
)
outline <- rbind(outline, outline[1L, ])
station_points <- data.frame(
  x = c(.08, .15, .22, .25, .31, .35, .40, .44, .49, .55, .60,
        .69, .70, .76, .82, .87, .91, .20, .33, .52, .68, .79),
  y = c(.39, .66, .24, .82, .54, .12, .73, .42, .88, .25, .67,
        .45, .82, .30, .61, .34, .47, .51, .33, .58, .25, .72)
)

map_a <- transform_xy(region, .45, 5.25, 4.88, 7.18)
outline_a <- transform_xy(outline, .45, 5.25, 4.88, 7.18)
stations_a <- transform_xy(station_points, .45, 5.25, 4.88, 7.18)
map_b <- transform_xy(region, 8.42, 11.78, 5.25, 7.18)
outline_b <- transform_xy(outline, 8.42, 11.78, 5.25, 7.18)
stations_b <- transform_xy(station_points, 8.42, 11.78, 5.25, 7.18)
star <- transform_xy(data.frame(x = .73, y = .66), 8.42, 11.78, 5.25, 7.18)
neighbor_index <- c(9, 11, 13, 14, 16, 22)
neighbor_lines <- data.frame(
  x = stations_b$x[neighbor_index], y = stations_b$y[neighbor_index],
  xend = star$x, yend = star$y
)

surface_grid <- expand.grid(ix = 1:10, iy = 1:5)
surface_grid$value <- with(surface_grid, sin(ix / 2.2) + cos(iy / 1.4))
surface_palette <- grDevices::hcl.colors(100, "Spectral", rev = TRUE)
surface_grid$fill <- surface_palette[
  cut(surface_grid$value, breaks = 100, labels = FALSE, include.lowest = TRUE)
]
surface_grid$xmin <- 8.62 + (surface_grid$ix - 1) * .29
surface_grid$xmax <- surface_grid$xmin + .30
surface_grid$ymin <- 2.74 + (surface_grid$iy - 1) * .20
surface_grid$ymax <- surface_grid$ymin + .21
small_surface <- surface_grid
small_surface$xmin <- 2.78 + (small_surface$ix - 1) * .23
small_surface$xmax <- small_surface$xmin + .24
small_surface$ymin <- 1.35 + (small_surface$iy - 1) * .14
small_surface$ymax <- small_surface$ymin + .15

wave_x <- seq(0, 1, length.out = 100)
local_waves <- do.call(rbind, lapply(seq_len(4), function(series) {
  data.frame(
    x = 2.85 + 2.15 * wave_x,
    y = .62 + series * .16 + .055 * sin(5 * pi * wave_x + series),
    group = series,
    colour = c("#2693C1", "#6DB36A", "#E5AE2E", "#D95046")[series]
  )
}))
error_wave <- data.frame(
  x = 5.80 + 1.42 * wave_x,
  y = 1.42 + .12 * sin(22 * pi * wave_x) + .055 * sin(51 * pi * wave_x)
)
latent_wave <- data.frame(
  x = 12.42 + 2.78 * wave_x,
  y = 2.74 + .10 * sin(10 * pi * wave_x) + .045 * sin(27 * pi * wave_x)
)
response_wave <- data.frame(
  x = 9.42 + 5.20 * wave_x,
  y = .80 + .10 * sin(11 * pi * wave_x) + .055 * sin(31 * pi * wave_x)
)
arrow_style <- grid::arrow(length = grid::unit(2.2, "mm"), type = "closed")

scheme <- ggplot() +
  annotate("rect", xmin = .08, xmax = 7.92, ymin = .08, ymax = 7.92,
           fill = "white", colour = colour$border, linewidth = .65) +
  annotate("rect", xmin = 8.08, xmax = 15.92, ymin = .08, ymax = 7.92,
           fill = "white", colour = colour$border, linewidth = .65) +
  annotate("rect", xmin = .08, xmax = 7.92, ymin = 7.48, ymax = 7.92,
           fill = colour$navy, colour = NA) +
  annotate("rect", xmin = 8.08, xmax = 15.92, ymin = 7.48, ymax = 7.92,
           fill = colour$navy, colour = NA) +
  annotate("text", x = .28, y = 7.70,
           label = "Observed monitoring network: model decomposition",
           hjust = 0, colour = "white", family = "sans", fontface = "bold", size = 4.15) +
  annotate("text", x = 8.28, y = 7.70,
           label = "Unmonitored location: posterior prediction",
           hjust = 0, colour = "white", family = "sans", fontface = "bold", size = 4.15) +

  geom_polygon(data = map_a, aes(x, y, group = group, fill = fill),
               colour = "white", linewidth = .45, alpha = .88) +
  geom_path(data = outline_a, aes(x, y), colour = "#111820", linewidth = .8) +
  geom_point(data = stations_a, aes(x, y), shape = 21, size = 2.4,
             fill = "#111820", colour = "white", stroke = .45) +
  annotate("rect", xmin = 5.42, xmax = 7.55, ymin = 4.90, ymax = 7.17,
           fill = colour$grey_light, colour = NA) +
  annotate("text", x = 5.66, y = 6.55, label = "Spatial partition (sPPM)",
           hjust = 0, family = "sans", fontface = "bold", size = 3.05, colour = colour$ink) +
  annotate("rect", xmin = rep(5.67, 4), xmax = rep(5.93, 4),
           ymin = c(6.15, 5.83, 5.51, 5.19), ymax = c(6.38, 6.06, 5.74, 5.42),
           fill = c("#1895C3", "#74B66A", "#F1BE3E", "#DE574B"), colour = NA) +
  annotate("text", x = 6.04, y = c(6.265, 5.945, 5.625, 5.305),
           label = paste("Regime", 1:4), hjust = 0, family = "sans", size = 2.9,
           colour = colour$ink) +

  annotate("rect", xmin = .42, xmax = 7.58, ymin = 4.13, ymax = 4.72,
           fill = "white", colour = colour$border, linewidth = .45) +
  annotate("text", x = c(2.30, 3.65, 4.90, 6.25), y = 4.43,
           label = c("y[t](bold(s)[i])", "bold(x)[ti]^T*bold(beta)",
                     "w[t](bold(s)[i])", "epsilon[t](bold(s)[i])"),
           parse = TRUE, family = "sans", size = 4.15, colour = colour$ink) +
  annotate("text", x = c(2.96, 4.28, 5.56), y = 4.43,
           label = c("=", "+", "+"), family = "sans", size = 4.15,
           colour = colour$ink) +

  annotate("rect", xmin = .30, xmax = 2.58, ymin = .36, ymax = 3.92,
           fill = colour$blue_light, colour = "#A9D4EB", linewidth = .6) +
  annotate("rect", xmin = 2.62, xmax = 5.42, ymin = .36, ymax = 3.92,
           fill = colour$green_light, colour = "#B8DDBB", linewidth = .6) +
  annotate("rect", xmin = 5.57, xmax = 7.58, ymin = .36, ymax = 3.92,
           fill = colour$red_light, colour = "#EDB9B4", linewidth = .6) +
  annotate("text", x = 1.445, y = 3.62, label = "Measured systematic variation",
           family = "sans", fontface = "bold", size = 3.05, colour = "#17669B") +
  annotate("text", x = 4.02, y = 3.62, label = "Latent residual process",
           family = "sans", fontface = "bold", size = 3.05, colour = "#24713D") +
  annotate("text", x = 6.575, y = 3.62, label = "Measurement error",
           family = "sans", fontface = "bold", size = 3.05, colour = "#B83D34") +
  annotate("text", x = 1.445, y = 3.28, label = "bold(x)[ti]^T*bold(beta)",
           parse = TRUE, family = "sans", size = 4.0, colour = colour$ink) +
  annotate("text", x = 4.02, y = 3.28,
           label = "w[t](bold(s)[i]) == phi[i]*w[t-1](bold(s)[i]) + u[t](bold(s)[i]) + eta[t](bold(s)[i])",
           parse = TRUE, family = "sans", size = 2.65, colour = colour$ink) +
  annotate("text", x = 6.575, y = 3.28, label = "epsilon[t](bold(s)[i])",
           parse = TRUE, family = "sans", size = 4.0, colour = colour$ink) +
  annotate("rect",
           xmin = c(.70, .84, .98, 1.12), xmax = c(1.80, 1.94, 2.08, 2.22),
           ymin = c(1.18, 1.34, 1.50, 1.66), ymax = c(2.38, 2.54, 2.70, 2.86),
           fill = c("#D95F45", "#E9B949", "#4DA6C8", "#4472A7"),
           colour = "white", linewidth = .3, alpha = .86) +
  geom_rect(data = small_surface,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
            colour = NA, alpha = .82) +
  geom_line(data = local_waves, aes(x, y, group = group, colour = colour),
            linewidth = .55) +
  geom_line(data = error_wave, aes(x, y), colour = "#202830", linewidth = .55) +
  annotate("segment", x = 3.65, y = 4.25, xend = 1.45, yend = 3.92,
           colour = colour$blue, linewidth = .55, arrow = arrow_style) +
  annotate("segment", x = 4.90, y = 4.25, xend = 4.02, yend = 3.92,
           colour = colour$green, linewidth = .55, arrow = arrow_style) +
  annotate("segment", x = 6.25, y = 4.25, xend = 6.58, yend = 3.92,
           colour = colour$red, linewidth = .55, arrow = arrow_style) +

  geom_polygon(data = map_b, aes(x, y, group = group, fill = fill),
               colour = "white", linewidth = .35, alpha = .38) +
  geom_path(data = outline_b, aes(x, y), colour = "#111820", linewidth = .8) +
  geom_segment(data = neighbor_lines, aes(x, y, xend = xend, yend = yend),
               colour = colour$navy, linewidth = .45, linetype = "22") +
  geom_point(data = stations_b, aes(x, y), shape = 21, size = 2.0,
             fill = "#111820", colour = "white", stroke = .4) +
  geom_point(data = star, aes(x, y), shape = 8, size = 5.5,
             colour = colour$purple, fill = "white", stroke = 1.15) +
  annotate("text", x = star$x + .25, y = star$y - .03,
           label = "s[0]", parse = TRUE, family = "sans", size = 3.4,
           colour = colour$ink) +

  annotate("rect", xmin = 12.02, xmax = 15.58, ymin = 5.98, ymax = 7.18,
           fill = colour$blue_light, colour = "#B8D8E9", linewidth = .55) +
  annotate("point", x = 12.28, y = 6.92, shape = 21, size = 6.5,
           fill = colour$navy, colour = colour$navy, stroke = 0) +
  annotate("text", x = 12.28, y = 6.92, label = "1", family = "sans",
           fontface = "bold", size = 3.0, colour = "white") +
  annotate("text", x = 12.58, y = 6.94, label = "Predictive allocation (spatial partition)",
           hjust = 0, family = "sans", fontface = "bold", size = 3.0,
           colour = colour$ink) +
  annotate("text", x = 13.80, y = 6.46,
           label = "Pr(c[0]==k*' | '*rho[n]*','*bold(s)[0]) %prop% Delta[C*g](k)",
           parse = TRUE, family = "sans", size = 3.15, colour = colour$ink) +
  annotate("rect", xmin = 12.12, xmax = 13.76, ymin = 5.08, ymax = 5.82,
           fill = colour$green_light, colour = "#B9DCBE", linewidth = .5) +
  annotate("rect", xmin = 13.94, xmax = 15.58, ymin = 5.08, ymax = 5.82,
           fill = "#FFF2D7", colour = "#E7C77D", linewidth = .5) +
  annotate("text", x = 12.94, y = 5.64, label = "Existing regime",
           family = "sans", fontface = "bold", size = 2.65, colour = colour$green) +
  annotate("text", x = 12.94, y = 5.34,
           label = "plain('inherit')~(phi[k]*','*tau[k]^2)",
           parse = TRUE, family = "sans", size = 2.35, colour = colour$muted) +
  annotate("text", x = 14.76, y = 5.64, label = "New regime",
           family = "sans", fontface = "bold", size = 2.65, colour = "#9A6814") +
  annotate("text", x = 14.76, y = 5.34, label = "draw from the base measure",
           family = "sans", size = 2.30, colour = colour$muted) +

  annotate("rect", xmin = 8.42, xmax = 11.82, ymin = 2.33, ymax = 5.04,
           fill = colour$blue_light, colour = "#ADD5EA", linewidth = .6) +
  annotate("point", x = 8.68, y = 4.77, shape = 21, size = 6.5,
           fill = colour$navy, colour = colour$navy, stroke = 0) +
  annotate("text", x = 8.68, y = 4.77, label = "2", family = "sans",
           fontface = "bold", size = 3.0, colour = "white") +
  annotate("text", x = 9.02, y = 4.79, label = "Spatial interpolation (common shock)",
           hjust = 0, family = "sans", fontface = "bold", size = 3.0,
           colour = colour$ink) +
  annotate("text", x = 10.12, y = 4.34,
           label = "u[t](bold(s)[0])~'|'~bold(u)[t](N(0)) %~% N(a[0]^T*bold(u)[t](N(0))*','*delta[0])",
           parse = TRUE, family = "sans", size = 2.75, colour = colour$ink) +
  geom_rect(data = surface_grid,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
            colour = NA, alpha = .86) +
  annotate("point", x = 10.78, y = 3.25, shape = 8, size = 4.2,
           colour = colour$purple, stroke = 1.0) +

  annotate("rect", xmin = 12.12, xmax = 15.58, ymin = 2.33, ymax = 4.98,
           fill = colour$green_light, colour = "#B9DCBE", linewidth = .6) +
  annotate("point", x = 12.38, y = 4.71, shape = 21, size = 6.5,
           fill = colour$green, colour = colour$green, stroke = 0) +
  annotate("text", x = 12.38, y = 4.71, label = "3", family = "sans",
           fontface = "bold", size = 3.0, colour = "white") +
  annotate("text", x = 12.72, y = 4.74, label = "Local dynamics (predicted trajectory)",
           hjust = 0, family = "sans", fontface = "bold", size = 3.0,
           colour = colour$ink) +
  annotate("text", x = 13.85, y = 4.20,
           label = "w[1](bold(s)[0]) %~% N(0*','*tau[0]^2/(1-phi[0]^2))",
           parse = TRUE, family = "sans", size = 2.75, colour = colour$ink) +
  annotate("text", x = 13.85, y = 3.75,
           label = "w[t](bold(s)[0]) == phi[0]*w[t-1](bold(s)[0]) + u[t](bold(s)[0]) + eta[t](bold(s)[0])",
           parse = TRUE, family = "sans", size = 2.52, colour = colour$ink) +
  geom_ribbon(data = latent_wave, aes(x = x, ymin = y - .14, ymax = y + .14),
              fill = "#9DD4AA", alpha = .35) +
  geom_line(data = latent_wave, aes(x, y), colour = "#187A45", linewidth = .65) +

  annotate("rect", xmin = 9.12, xmax = 15.30, ymin = .38, ymax = 2.02,
           fill = colour$purple_light, colour = "#CDBCE0", linewidth = .6) +
  annotate("point", x = 9.38, y = 1.76, shape = 21, size = 6.5,
           fill = colour$purple, colour = colour$purple, stroke = 0) +
  annotate("text", x = 9.38, y = 1.76, label = "4", family = "sans",
           fontface = "bold", size = 3.0, colour = "white") +
  annotate("text", x = 9.72, y = 1.78, label = "Predicted response",
           hjust = 0, family = "sans", fontface = "bold", size = 3.0,
           colour = colour$ink) +
  annotate("text", x = 12.20, y = 1.38,
           label = "hat(y)[t](bold(s)[0]) == bold(x)[t0]^T*bold(beta) + w[t](bold(s)[0]) + epsilon[t](bold(s)[0])",
           parse = TRUE, family = "sans", size = 2.9, colour = colour$ink) +
  geom_ribbon(data = response_wave, aes(x = x, ymin = y - .14, ymax = y + .14),
              fill = "#BDA4D5", alpha = .35) +
  geom_line(data = response_wave, aes(x, y), colour = colour$purple, linewidth = .65) +

  annotate("segment", x = 10.10, y = 5.25, xend = 10.10, yend = 5.04,
           colour = colour$blue, linewidth = .65, arrow = arrow_style) +
  annotate("segment", x = 13.80, y = 5.98, xend = 12.94, yend = 5.82,
           colour = colour$green, linewidth = .55, arrow = arrow_style) +
  annotate("segment", x = 13.80, y = 5.98, xend = 14.76, yend = 5.82,
           colour = "#B77B18", linewidth = .55, arrow = arrow_style) +
  annotate("segment", x = 12.94, y = 5.08, xend = 13.45, yend = 4.98,
           colour = colour$green, linewidth = .55, linetype = "22", arrow = arrow_style) +
  annotate("segment", x = 14.76, y = 5.08, xend = 14.25, yend = 4.98,
           colour = "#B77B18", linewidth = .55, linetype = "22", arrow = arrow_style) +
  annotate("segment", x = 11.82, y = 3.15, xend = 12.12, yend = 3.15,
           colour = colour$blue, linewidth = .65, arrow = arrow_style) +
  annotate("curve", x = 13.85, y = 2.33, xend = 13.10, yend = 2.02,
           curvature = .12, colour = colour$purple, linewidth = .65,
           arrow = arrow_style) +

  scale_fill_identity() +
  scale_colour_identity() +
  coord_cartesian(xlim = c(0, 16), ylim = c(0, 8), clip = "off") +
  theme_void(base_size = 11, base_family = "sans") +
  theme(
    plot.margin = margin(8, 8, 8, 8),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

ggsave(pdf_file, scheme, device = grDevices::cairo_pdf,
       width = 14, height = 7.2, units = "in", bg = "white")
ggsave(png_file, scheme, width = 14, height = 7.2,
       units = "in", dpi = 320, bg = "white")

caption <- paste(
  "Conceptual decomposition underlying the DCS-sPPM. Spatially shared",
  "variation is represented through a continuous common-shock process, while",
  "residual local temporal behavior is characterized by regime-specific",
  "persistence and innovation variability. At an unmonitored location, these",
  "components are extended through spatial interpolation and predictive regime",
  "allocation, respectively."
)
writeLines(caption, caption_file)

message("Created: ", pdf_file)
message("Created: ", png_file)
message("Created: ", caption_file)
