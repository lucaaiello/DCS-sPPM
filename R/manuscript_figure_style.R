manuscript_palette <- list(
  ink = "#172A3A",
  muted = "#5C6770",
  land = "grey92",
  border = "white",
  background = "white",
  continuous = "#2C7FB8",
  discrete = "#238B45",
  accent = "#D95F0E"
)

manuscript_map_theme <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(
        fill = manuscript_palette$background,
        colour = NA
      ),
      plot.background = ggplot2::element_rect(
        fill = manuscript_palette$background,
        colour = NA
      ),
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = base_size - 2, colour = "grey25"),
      axis.ticks = ggplot2::element_line(colour = "grey65", linewidth = 0.25),
      legend.position = "right",
      legend.direction = "vertical",
      legend.box = "vertical",
      legend.justification = "center",
      legend.title = ggplot2::element_text(
        size = base_size - 1,
        family = "sans",
        face = "plain"
      ),
      legend.text = ggplot2::element_text(
        size = base_size - 2,
        family = "sans"
      ),
      legend.key.height = grid::unit(0.48, "cm"),
      legend.key.width = grid::unit(0.42, "cm"),
      legend.spacing.y = grid::unit(0.12, "cm"),
      plot.margin = ggplot2::margin(3, 3, 3, 3)
    )
}

manuscript_colourbar <- function(title.position = "top") {
  ggplot2::guide_colourbar(
    title.position = title.position,
    title.hjust = 0.5,
    direction = "vertical",
    barheight = grid::unit(3.1, "cm"),
    barwidth = grid::unit(0.38, "cm"),
    ticks = TRUE,
    frame.colour = "grey70"
  )
}

manuscript_country_layer <- function(world, fill = manuscript_palette$land) {
  ggplot2::geom_polygon(
    data = world,
    ggplot2::aes(long, lat, group = group),
    inherit.aes = FALSE,
    fill = fill,
    colour = manuscript_palette$border,
    linewidth = 0.25
  )
}

save_manuscript_pair <- function(
    plot,
    stem,
    directories,
    width,
    height,
    dpi = 320
) {
  directories <- unique(normalizePath(
    directories,
    winslash = "/",
    mustWork = FALSE
  ))
  for (directory in directories) {
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(
      file.path(directory, paste0(stem, ".pdf")),
      plot,
      width = width,
      height = height,
      units = "in",
      device = grDevices::cairo_pdf,
      bg = "white"
    )
    ggplot2::ggsave(
      file.path(directory, paste0(stem, ".png")),
      plot,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      bg = "white"
    )
  }
}

english_month_labels <- function(dates) {
  month.abb[as.integer(format(dates, "%m"))]
}
