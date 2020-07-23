# rm(list = ls())
library(tidyverse)
library(gtools)

basic.theme <- theme(
	panel.background = element_rect(
		fill = "transparent",colour = NA),
	panel.grid.major = element_line(colour = "gray50"),
#	panel.grid.minor = element_blank(),
	plot.background = element_rect(
		fill = "transparent",colour = NA),
	legend.background = element_rect(
		fill="transparent"),
	legend.key = element_rect(colour = NA, fill = NA),
	legend.key.height = unit(2, "lines"),
	panel.spacing = unit(2, "lines"))

retrieve.summary <- function(sample, version, measures, model) {
  ################################################################################
  # Set up
  ################################################################################
  
  # sample <- "Random"
  # version <- "Casillas, Brown, & Levinson (submitted March 2019)"
  # measures <- "TCDS"
  # model <- "Yes"

  all.sum.stats <- read_csv("shiny_input/quantity-scores_rand-and-tt.csv")
  sum.stat.tbl.a <- filter(all.sum.stats, Sample == "Turn taking") %>%
    mutate(age_months = substr(paste0("0", as.character(age_mo_round)),
                               nchar(age_mo_round), nchar(age_mo_round)+1)) %>%
    select(-age_mo_round) %>%
    group_by(age_months) %>%
    summarise(mean = mean(tds_mph),
              median = median(tds_mph),
              min = min(tds_mph),
              max = max(tds_mph),
              n_clips = n())
  sum.stat.tbl.b <- filter(all.sum.stats, Sample == "Turn taking") %>%
    summarise(age_months = "All",
              mean = mean(tds_mph),
              median = median(tds_mph),
              min = min(tds_mph),
              max = max(tds_mph),
              n_clips = n())
  sum.stat.tbl <- bind_rows(sum.stat.tbl.a, sum.stat.tbl.b)
  
  sum.stat.fig <- ggplot(data=all.sum.stats,
                         aes(x = age_mo_round, y = tds_mph)) +
    geom_boxplot(aes(group = age_mo_round), fill = "forestgreen",
                 alpha = 0.4, outlier.shape = NA) +
    geom_jitter(aes(group = age_mo_round), color = "forestgreen",
                alpha = 0.4) +
    xlab("Age (months)") +
    ylab("TCDS min/hr") +
    basic.theme
    
  
  all.models <- read_csv("shiny_input/all_model_tables.csv")
  model.output.tbl <- filter(all.models, model == "TCDS_turntaking_nb")
  
  model.res.fig <- "c_o.tpm_random_log_gaus.res.plot.png"
  

    return(list(
      sum.stat.tbl = sum.stat.tbl,
      sum.stat.fig = sum.stat.fig,
      model.output.tbl = model.output.tbl,
      model.res.fig = model.res.fig
    ))
}

