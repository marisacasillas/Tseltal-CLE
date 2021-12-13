# Create summaries, adding in zeroes for annotated clips with no speech

# ... at the clip level
summarize.over.clips <- function(utt.info,
                                 clip.dur.min, all.clips.CHIOTR.tiers,
                                 rec.info.stds) {
  clip.summary <- utt.info %>%
    group_by(recording, tier.name, clip.onset, spkr.type.corr) %>%
    summarize(
      spch.min = sum(dur.sec)/60,
      min.p.hr = spch.min*(60/clip.dur.min)
    ) %>%
    pivot_wider(names_from = spkr.type.corr,
                values_from = c(spch.min, min.p.hr)) %>%
    mutate_at(vars(starts_with(c("spch", "min"))), 
              funs(if_else(is.na(.), 0, .))) %>%
    mutate(tot.spch.min = rowSums(across(starts_with(c("spch")))),
           tot.min.p.hr = tot.spch.min*(60/clip.dur.min)) %>%
    mutate(`spch.min_niños` = `spch.min_niña` + `spch.min_niño`,
           `min.p.hr_niños` = `spch.min_niños`*(60/clip.dur.min)) %>%
    full_join(mutate(all.clips.CHIOTR.tiers,
                     clip.onset.uncorrected = as.numeric(clip.onset.uncorrected),
                     clip.onset.corrected = as.numeric(clip.onset.corrected)),
              by = c("recording", "tier.name",
                     "clip.onset" = "clip.onset.uncorrected")) %>%
    ungroup() %>%
    select(-clip.onset) %>%
    rename(clip.onset = clip.onset.corrected) %>%
    full_join(rec.info.stds, by = "recording") %>%
    mutate(across(everything(), ~replace_na(.x, 0))) %>%
    # add non-zero marker
    mutate(
      nz.spch = ifelse(tot.min.p.hr > 0, 1, 0)
    ) %>%
    # add standardized/recoded predictors
    mutate(
      mated.std = recode_factor(mat_ed,
                                "none" = "none",
                                "primary" = "primary",
                                "secondary" = "secondary",
                                "preparatory" = "preparatory"),
      mated.bin = recode_factor(mat_ed,
                                "none" = "0-5", "primary" = "0-5",
                                "secondary" = "6+", "preparatory" = "6+"),
      motyr.std = ((mother_dob - motyr.m)/motyr.sd),
      nsb.std = ((number_older_sibs - nsb.m)/nsb.sd),
      hsz.std = ((household_size - hsz.m)/hsz.sd),
      start.hr = as.numeric(substr(clip.onset, 1, 2)),
      stthr.std = (start.hr - 12)/12,
      stthr.tri = ifelse(start.hr < 11, "morning",
                         ifelse(start.hr > 13, "afternoon", "midday"))
    )
  CHI.vocs.temp <- clip.summary %>%
    filter(tier.name == "CHI") %>%
    # add NAs where needed to avoid human tibble-reviewing confusion
    mutate_at(vars(starts_with(c("spch", "min")) & !(ends_with("CHI"))), 
              funs(if_else(. >= 0, NA, NA)))
  OTR.vocs.temp <- clip.summary %>%
    filter(tier.name == "OTR") %>%
    # add NAs where needed to avoid human tibble-reviewing confusion
    mutate_at(vars(starts_with(c("spch", "min")) & (ends_with("CHI"))), 
              funs(if_else(. >= 0, NA, NA)))
  clip.summary <- CHI.vocs.temp %>% #
    bind_rows(OTR.vocs.temp) %>%
    mutate(stthr.tri = as_factor(stthr.tri),
           stthr.tri = fct_relevel(stthr.tri, c("midday",
                                                "morning", "afternoon")),
           stthr.tri.a = fct_relevel(stthr.tri, c("afternoon",
                                                  "midday", "morning")),
           stthr.tri.o = fct_relevel(stthr.tri, c("morning",
                                                  "afternoon", "midday"))) %>%
    arrange(aclew_id, tier.name, segment.num)
  return(clip.summary)
}


# ... and then with averages across clips at the recording level
summarize.over.recs <- function(clip.summary, rec.info.stds) {
  rec.summary <- clip.summary %>%
    group_by(aclew_id, tier.name,
             mated.std, mated.bin, motyr.std, nsb.std, hsz.std) %>%
    summarize(
      total.clips = n(),
      tot.spchmin = sum(tot.spch.min),
      tot.fa.spchmin = sum(spch.min_mujer),
      tot.ma.spchmin = sum(spch.min_hombre),
      tot.fc.spchmin = sum(`spch.min_niña`),
      tot.mc.spchmin = sum(`spch.min_niño`),
      tot.oc.spchmin = sum(`spch.min_niños`),
      tot.un.spchmin = sum(`spch.min_other/unknown`),
      mean.tot.spchmin = mean(tot.spch.min),
      median.tot.spchmin = median(tot.spch.min),
      min.tot.spchmin = min(tot.spch.min),
      max.tot.spchmin = max(tot.spch.min),
      sd.tot.spchmin = sd(tot.spch.min),
      mean.fa.spchmin = mean(spch.min_mujer),
      median.fa.spchmin = median(spch.min_mujer),
      min.fa.spchmin = min(spch.min_mujer),
      max.fa.spchmin = max(spch.min_mujer),
      sd.fa.spchmin = sd(spch.min_mujer),
      mean.ma.spchmin = mean(spch.min_hombre),
      median.ma.spchmin = median(spch.min_hombre),
      min.ma.spchmin = min(spch.min_hombre),
      max.ma.spchmin = max(spch.min_hombre),
      sd.ma.spchmin = sd(spch.min_hombre),
      mean.fc.spchmin = mean(`spch.min_niña`),
      median.fc.spchmin = median(`spch.min_niña`),
      min.fc.spchmin = min(`spch.min_niña`),
      max.fc.spchmin = max(`spch.min_niña`),
      sd.fc.spchmin = sd(`spch.min_niña`),
      mean.mc.spchmin = mean(`spch.min_niño`),
      median.mc.spchmin = median(`spch.min_niño`),
      min.mc.spchmin = min(`spch.min_niño`),
      max.mc.spchmin = max(`spch.min_niño`),
      sd.mc.spchmin = sd(`spch.min_niño`),
      mean.oc.spchmin = mean(`spch.min_niños`),
      median.oc.spchmin = median(`spch.min_niños`),
      min.oc.spchmin = min(`spch.min_niños`),
      max.oc.spchmin = max(`spch.min_niños`),
      sd.oc.spchmin = sd(`spch.min_niños`),
      mean.un.spchmin = mean(`spch.min_other/unknown`),
      median.un.spchmin = median(`spch.min_other/unknown`),
      min.un.spchmin = min(`spch.min_other/unknown`),
      max.un.spchmin = max(`spch.min_other/unknown`),
      sd.un.spchmin = sd(`spch.min_other/unknown`),
      mean.tot.mph = mean(tot.min.p.hr),
      median.tot.mph = median(tot.min.p.hr),
      min.tot.mph = min(tot.min.p.hr),
      max.tot.mph = max(tot.min.p.hr),
      sd.tot.mph = sd(tot.min.p.hr),
      mean.fa.mph = mean(min.p.hr_mujer),
      median.fa.mph = median(min.p.hr_mujer),
      min.fa.mph = min(min.p.hr_mujer),
      max.fa.mph = max(min.p.hr_mujer),
      sd.fa.mph = sd(min.p.hr_mujer),
      mean.ma.mph = mean(min.p.hr_hombre),
      median.ma.mph = median(min.p.hr_hombre),
      min.ma.mph = min(min.p.hr_hombre),
      max.ma.mph = max(min.p.hr_hombre),
      sd.ma.mph = sd(min.p.hr_hombre),
      mean.fc.mph = mean(`min.p.hr_niña`),
      median.fc.mph = median(`min.p.hr_niña`),
      min.fc.mph = min(`min.p.hr_niña`),
      max.fc.mph = max(`min.p.hr_niña`),
      sd.fc.mph = sd(`min.p.hr_niña`),
      mean.mc.mph = mean(`min.p.hr_niño`),
      median.mc.mph = median(`min.p.hr_niño`),
      min.mc.mph = min(`min.p.hr_niño`),
      max.mc.mph = max(`min.p.hr_niño`),
      sd.mc.mph = sd(`min.p.hr_niño`),
      mean.oc.mph = mean(`min.p.hr_niños`),
      median.oc.mph = median(`min.p.hr_niños`),
      min.oc.mph = min(`min.p.hr_niños`),
      max.oc.mph = max(`min.p.hr_niños`),
      sd.oc.mph = sd(`min.p.hr_niños`),
      mean.un.mph = mean(`min.p.hr_other/unknown`),
      median.un.mph = median(`min.p.hr_other/unknown`),
      min.un.mph = min(`min.p.hr_other/unknown`),
      max.un.mph = max(`min.p.hr_other/unknown`),
      sd.un.mph = sd(`min.p.hr_other/unknown`),
    ) %>%
    left_join(rec.info.stds, by = "aclew_id")
  return(rec.summary)
}

# # corr btwn OTR and CHI at clip and recording levels
# clip.CHI.OTR.rates <- clip.summary %>%
#   select(-tot.spch.min, -nz.spch) %>%
#   pivot_wider(names_from = "tier.name", values_from = "min.p.hr")
#     
# rec.CHI.OTR.rates <- rec.summary %>%
#   select(-median.mph, -min.mph, -max.mph, -sd.mph, -total.clips) %>%
#   pivot_wider(names_from = "tier.name", values_from = "mean.mph")
# 
# OTR.mean <- mean(rec.CHI.OTR.rates$OTR)
# OTR.sd <- sd(rec.CHI.OTR.rates$OTR)
# 
# rec.CHI.OTR.rates <- rec.CHI.OTR.rates %>%
#   mutate(
#     Tertile = factor(case_when(
#       OTR <= OTR.mean - OTR.sd ~ "low TCDS",
#       OTR >= OTR.mean + OTR.sd ~ "high TCDS",
#       TRUE ~ "middle TCDS"
#     ), levels = c("low TCDS", "middle TCDS", "high TCDS"))
#   ) 
