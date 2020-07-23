sem <- function (x) {
    sd(x) / sqrt(length(x))
}

qqplot.data <- function (vec) # argument: vector of numbers
{
  # following four lines from base R's qqline()
  y <- quantile(vec[!is.na(vec)], c(0.25, 0.75))
  x <- qnorm(c(0.25, 0.75))
  slope <- diff(y)/diff(x)
  int <- y[1L] - slope * x[1L]

  d <- data.frame(resids = vec)

  residplot <- ggplot(d, aes(sample = resids)) +
    stat_qq() + geom_abline(slope = slope, intercept = int)

  return(residplot)
}

left_join_NA <- function(x, y, ...) {
  left_join(x = x, y = y, by = ...) %>% 
    mutate_each(funs(replace(., which(is.na(.)), 0)))
}

# To create multi-panel plots (option 2)
grid_arrange_shared_legend <- function(..., ncol = length(list(...)),
                                       nrow = 1, position = c("bottom", "right")) {
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
#  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)

  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))

  # grid.newpage()
  # grid.draw(combined)
  return(combined)

  # return gtable invisibly
  # invisible(combined)

}

# Basic plotting theme settings
basic.theme <- theme(
	panel.background = element_rect(
		fill = "transparent",colour = NA),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	plot.background = element_rect(
		fill = "transparent",colour = NA),
	legend.background = element_rect(
		fill="transparent"),
	legend.text = element_text(size=30),
	legend.title = element_text(size=30),
	legend.key = element_rect(colour = NA, fill = NA),
	legend.key.height = unit(2, "lines"),
	axis.text.x = element_text(size=30),
	axis.title.x = element_text(size=30),
	axis.text.y = element_text(size=30),
	axis.title.y = element_text(size=30),
	strip.text = element_text(size=30),
	panel.spacing = unit(2, "lines"),
	plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
