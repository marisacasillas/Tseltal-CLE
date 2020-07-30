library(papaja)
library(ggplot2)
library(tidyverse)
library(viridis)
library(grid)
library(gridExtra)
library(ggpirate)
library(lme4)
library(glmmTMB)
library(DHARMa)
library(bbmle)
library(broom.mixed)
source("0-Helper.R")


options(scipen=999)
data.path <- "transcripts/anon/" # text files exported from ELAN
metadata.path <- "metadata/"
plot.path <- "plots/" # output plots
shiny.input.path <- "shiny_input/"
shiny.input.img.path <- "www/"
seg.index.file <- paste0(metadata.path, "Segment-order-inventory.csv")
ptcp.info.file <- paste0(metadata.path, "recording-info.csv")
comparison.file <- paste0(metadata.path, "comparison_studies.csv")
samplelabels <- c("High activity  ", "Random  ")
col.sample.bu <- list(
  scale_fill_manual(labels=samplelabels, values=viridis(2)),
  scale_color_manual(labels=samplelabels, values=viridis(2)))
col.sample.bu3 <- list(
  scale_fill_manual(labels=samplelabels, values=viridis(3)),
  scale_color_manual(labels=samplelabels, values=viridis(3)))
allowed.overlap <- 1000 #ms
allowed.gap <- 2000 #ms
waking.hours <- 14

# Read in annotation files
files <- list.files(path=data.path,pattern="*.txt")
all.data <- data.frame()
for (i in 1:length(files)) {
#  print(files[i])
  newfile <- read_csv(paste0(data.path, files[i]),
                      col_types = cols(val = col_character()))
  newfile$aclew_child_id <- unlist(strsplit(files[i], '\\.'))[1]
  all.data <- rbind(all.data, newfile)
}
all.data$row <- c(1:nrow(all.data))

# Read in supplementary data
ptcp.info <- read_csv(ptcp.info.file, col_types = cols()) %>%
  dplyr::select(-row)
seg.info <- read_csv(seg.index.file)

# Extract and convert start time of each sample
seg.info$start.hhmmss <- regmatches(seg.info$Media,
                                    regexpr("[[:digit:]]{6}", seg.info$Media))
seg.info$start.sec <- as.numeric(substr(seg.info$start.hhmmss,1,2))*3600 +
  as.numeric(substr(seg.info$start.hhmmss,3,4))*60 +
  as.numeric(substr(seg.info$start.hhmmss,5,6))
seg.info$start.hr <- round(seg.info$start.sec/3600, 3)

seg.info$clipoffset.hhmmss <- regmatches(seg.info$Media,
                                    regexpr("(?<=[[:digit:]]{6}_)[[:digit:]]{6}",
                                            seg.info$Media, perl = TRUE))
seg.info$clipoffset.sec <- as.numeric(substr(seg.info$clipoffset.hhmmss,1,2))*3600 +
  as.numeric(substr(seg.info$clipoffset.hhmmss,3,4))*60 +
  as.numeric(substr(seg.info$clipoffset.hhmmss,5,6))
seg.info$clipoffset.hr <- round(seg.info$clipoffset.sec/3600, 3)

# Add mean and sd values for participant-level predictors to ptcp.info
ptcp.info <- ptcp.info %>%
  mutate(
    tchiyr.m = mean(age_mo_round),
    motyr.m = mean(mother_age),
    nsb.m = mean(number_older_sibs),
    hsz.m = mean(household_size),
    tchiyr.sd = sd(age_mo_round),
    motyr.sd = sd(mother_age),
    nsb.sd = sd(number_older_sibs),
    hsz.sd = sd(household_size)
    )

# Merge in participant and segment info to the main data table
codes <- all.data %>% filter(tier == "code")

all.data <- all.data %>%
  filter(speaker != "") %>%
  left_join(ptcp.info, by = "aclew_child_id") %>%
  mutate(segment = "", sample = "",
         sample_type = "", segment_dur = 0)

for (i in 1:nrow(codes)) {
  rec <- codes$aclew_child_id[i]
  seg <- as.character(codes$val[i])
  seg.on <- codes$start[i]
  seg.off <- codes$stop[i]
  seg.idx <- which(all.data$aclew_child_id == rec &
                     all.data$start < seg.off &
                     all.data$stop > seg.on)
  all.data$segment[seg.idx] <- seg
}

# Label samples
all.data$sample[which(
  grepl('^random', all.data$segment))] <- "random"
all.data$sample[which(
  grepl('tt', all.data$segment))] <- "turn-taking"
all.data$sample[which(
  grepl('va', all.data$segment))] <- "high-activity"

# Label sample types and durations
random.samples <- which(grepl('^random', all.data$segment))
all.data$sample_type[random.samples] <- "random"
all.data$segment_dur[random.samples] <- 5

ext.samples <- which(grepl('^extension', all.data$segment))
all.data$sample_type[ext.samples] <- "extension"
all.data$segment_dur[ext.samples] <- 5

tt.samples <- which(grepl('^tt', all.data$segment))
all.data$sample_type[tt.samples] <- "turn-taking"
all.data$segment_dur[tt.samples] <- 1

va.samples <- which(grepl('^va', all.data$segment))
all.data$sample_type[va.samples] <- "turn-taking"
all.data$segment_dur[va.samples] <- 1

# Add in segment start time
all.data <- all.data %>%
  left_join(dplyr::select(seg.info, c("aclew_id", "CodeName", "start.hhmmss",
                                      "start.sec", "start.hr")),
            by = c("aclew_child_id" = "aclew_id", "segment" = "CodeName"))

avg.utt.len.tseltal <- all.data %>%
  filter(speaker != "CHI") %>%
  dplyr::select(dur) %>%
  summarise(mean.utt.len = mean(dur),
            median.utt.len = median(dur))

sample.demog <- ptcp.info %>%
  mutate(mat_age_rd = as.integer(round(mother_age, 0)),
         fat_age_rd = as.integer(round(father_age, 0)),
         hszrd = as.integer(round(household_size + 1, 0))) %>%
  select(age_childes, child_sex, mat_age_rd, mat_ed, hszrd) %>%
  rename("Age" = age_childes, "Sex" = child_sex,
         "Mot age" = mat_age_rd, "Mot edu" = mat_ed,
         "People in house" = hszrd)

ptcp.info <- mutate(ptcp.info,
                    rec.start.hr = lubridate::hour(start_of_recording) +
                      lubridate::minute(start_of_recording)/60 +
                      lubridate::second(start_of_recording)/3600,
                    rec.stop.hr = rec.start.hr + length_of_recording/3600) %>%
  arrange(age_mo_round) %>%
  mutate(order = seq(1:10))

used.clips <- seg.info %>%
  filter(Include == 1) %>%
  left_join(ptcp.info, by = c("aclew_id" = "aclew_child_id")) %>%
  mutate(clip.dur = ifelse(grepl(('random'), CodeName), 5,
                           ifelse(grepl('extension', CodeName), 6, 1)),
         sample.type = ifelse(grepl(('random'), CodeName), "Random",
                           ifelse(grepl('tt', CodeName),
                                  "Turn taking", "Vocal activity")))
used.clips$sample.type = factor(used.clips$sample.type,
                                levels = c("Random", "Turn taking", "Vocal activity"))
  
clip.distribution <- ggplot() +
  geom_segment(data = ptcp.info,
               aes(x = rec.start.hr, y = order, xend = rec.stop.hr, yend = order)) +
  geom_segment(data = used.clips,
               aes(x = start.hr, y = order, xend = start.hr + clip.dur/60, yend = order,
                   color = sample.type), size = 3) +
  theme_apa() +
  theme(legend.position = "bottom") +
  scale_x_continuous(breaks = 7:21) +
  scale_y_continuous(breaks = 1:10, labels = ptcp.info$age_mo_round) +
  ylab("Child age (mo)") + xlab("Time of day (hr)") + labs(color = "Sample type") +
  scale_color_manual(values = c(viridis(3)[1], viridis(3)[2], "black"))
clip.distribution

# RANDOM
# Get min/hr speech measures
n.unique.rand.segs <- length(unique(all.data$segment[grepl("random", all.data$segment)]))
n.unique.recs <- length(unique(all.data$aclew_child_id))
all.rand.segments <- tibble(
  aclew_child_id = rep(unique(all.data$aclew_child_id),
                n.unique.rand.segs),
  segment = rep(unique(all.data$segment[grepl("random", all.data$segment)]),
                n.unique.recs))
# XDS
xds.per.seg.rand <- all.data %>%
  filter(sample == "random" & speaker != "CHI" &
           grepl("xds@", tier)) %>%
  group_by(aclew_child_id, segment, segment_dur) %>%
  summarise(xds_min = round(sum(dur)/60000,3)) %>%
  full_join(all.rand.segments, by = c("aclew_child_id", "segment")) %>%
  replace_na(list(segment_dur = 5, xds_min = 0)) %>%
  mutate(xds_mph = (xds_min/segment_dur)*60) %>%
  arrange(aclew_child_id, segment)
# ODS
ods.per.seg.rand <- all.data %>%
  filter(sample == "random" & speaker != "CHI" &
           grepl("xds@", tier) & val != "T") %>%
  group_by(aclew_child_id, segment, segment_dur) %>%
  summarise(ods_min = round(sum(dur)/60000,3)) %>%
  full_join(all.rand.segments, by = c("aclew_child_id", "segment")) %>%
  replace_na(list(segment_dur = 5, ods_min = 0)) %>%
  mutate(ods_mph = (ods_min/segment_dur)*60) %>%
  arrange(aclew_child_id, segment)
# TDS
tds.per.seg.rand <- all.data %>%
  filter(sample == "random" & speaker != "CHI" &
           grepl("xds@", tier) & val == "T") %>%
  group_by(aclew_child_id, segment, segment_dur) %>%
  summarise(tds_min = round(sum(dur)/60000,3)) %>%
  full_join(all.rand.segments, by = c("aclew_child_id", "segment")) %>%
  replace_na(list(segment_dur = 5, tds_min = 0)) %>%
  mutate(tds_mph = (tds_min/segment_dur)*60) %>%
  arrange(aclew_child_id, segment)
# All CDS
cds.per.seg.rand <- all.data %>%
  filter(sample == "random" & speaker != "CHI" &
           grepl("xds@", tier) & (val == "T" | val == "C")) %>%
  group_by(aclew_child_id, segment, segment_dur) %>%
  summarise(cds_min = round(sum(dur)/60000,3)) %>%
  full_join(all.rand.segments, by = c("aclew_child_id", "segment")) %>%
  replace_na(list(segment_dur = 5, cds_min = 0)) %>%
  mutate(cds_mph = (cds_min/segment_dur)*60) %>%
  arrange(aclew_child_id, segment)
# Number of speakers per clip
spkrs.per.seg.rand <- all.data %>%
  filter(sample == "random" & speaker != "CHI" &
           !(grepl("@", tier))) %>%
  group_by(aclew_child_id, segment) %>%
  summarise(n_spkrs_clip = length(unique(speaker)))
# All together
quantity.rand <- xds.per.seg.rand %>%
  full_join(ods.per.seg.rand, by = c("aclew_child_id", "segment", "segment_dur")) %>%
  full_join(tds.per.seg.rand, by = c("aclew_child_id", "segment", "segment_dur")) %>%
  full_join(cds.per.seg.rand, by = c("aclew_child_id", "segment", "segment_dur")) %>%
  full_join(spkrs.per.seg.rand, by = c("aclew_child_id", "segment")) %>%
  left_join(dplyr::select(seg.info, c("aclew_id", "CodeName", "start.hr")),
            by = c("aclew_child_id" = "aclew_id", "segment" = "CodeName")) %>%
  full_join(ptcp.info, by = "aclew_child_id") %>%
  replace_na(list(xds_min = 0, xds_mph = 0,
                  tds_min = 0, tds_mph = 0,
                  cds_min = 0, cds_mph = 0,
                  n_spkrs_clip = 0)) %>%
  mutate(prop_tds = tds_min/xds_min)
  # Don't replace NAs with 0s in this case; proportion is not meaningful w/o any speech
quantity.rand.bychild <- quantity.rand %>%
  group_by(aclew_child_id) %>%
  summarise(
    xds_min = mean(xds_min),
    xds_mph = mean(xds_mph),
    ods_min = mean(ods_min),
    ods_mph = mean(ods_mph),
    tds_min = mean(tds_min),
    tds_mph = mean(tds_mph),
    cds_min = mean(cds_min),
    cds_mph = mean(cds_mph),
    prop_tds = mean(prop_tds, na.rm = TRUE),
    m_n_spkrs = mean(n_spkrs_clip)) %>%
  full_join(ptcp.info, by = "aclew_child_id")

# Get xds and tds min/hr by speaker type
all.data$SpkrAge <- "Not known"
all.data$SpkrAge[grepl("FA|MA|UA", all.data$speaker)] <- "Adult"
all.data$SpkrAge[grepl("FC|MC|UC", all.data$speaker)] <- "Child"
all.rand.segments.sa <- tibble(
  aclew_child_id = rep(unique(all.data$aclew_child_id),
                2*n.unique.rand.segs),
  segment = rep(unique(all.data$segment[grepl("random", all.data$segment)]),
                2*n.unique.recs),
  SpkrAge = c(rep("Adult", (n.unique.rand.segs * n.unique.recs)),
              rep("Child", (n.unique.rand.segs * n.unique.recs))))
# XDS
xds.per.seg.rand.sa <- all.data %>%
  filter(sample == "random" & speaker != "CHI" & SpkrAge != "Not known" &
           grepl("xds@", tier)) %>%
  group_by(aclew_child_id, SpkrAge, segment, segment_dur) %>%
  summarise(xds_min.sa = round(sum(dur)/60000,3)) %>%
  full_join(all.rand.segments.sa, by = c("aclew_child_id", "segment", "SpkrAge")) %>%
  replace_na(list(segment_dur = 5, xds_min.sa = 0)) %>%
  mutate(xds_mph.sa = (xds_min.sa/segment_dur)*60) %>%
  arrange(aclew_child_id, segment, SpkrAge)
# ODS
ods.per.seg.rand.sa <- all.data %>%
  filter(sample == "random" & speaker != "CHI" & SpkrAge != "Not known" &
           grepl("xds@", tier) & val != "T") %>%
  group_by(aclew_child_id, SpkrAge, segment, segment_dur) %>%
  summarise(ods_min.sa = round(sum(dur)/60000,3)) %>%
  full_join(all.rand.segments.sa, by = c("aclew_child_id", "segment", "SpkrAge")) %>%
  replace_na(list(segment_dur = 5, ods_min.sa = 0)) %>%
  mutate(ods_mph.sa = (ods_min.sa/segment_dur)*60) %>%
  arrange(aclew_child_id, segment, SpkrAge)
# TDS
tds.per.seg.rand.sa <- all.data %>%
  filter(sample == "random" & speaker != "CHI" & SpkrAge != "Not known" &
           grepl("xds@", tier) & val == "T") %>%
  group_by(aclew_child_id, SpkrAge, segment, segment_dur) %>%
  summarise(tds_min.sa = round(sum(dur)/60000,3)) %>%
  full_join(all.rand.segments.sa, by = c("aclew_child_id", "segment", "SpkrAge")) %>%
  replace_na(list(segment_dur = 5, tds_min.sa = 0)) %>%
  mutate(tds_mph.sa = (tds_min.sa/segment_dur)*60) %>%
  arrange(aclew_child_id, segment, SpkrAge)
# All CDS
cds.per.seg.rand.sa <- all.data %>%
  filter(sample == "random" & speaker != "CHI" & SpkrAge != "Not known" &
           grepl("xds@", tier) & (val == "T" | val == "C")) %>%
  group_by(aclew_child_id, SpkrAge, segment, segment_dur) %>%
  summarise(cds_min.sa = round(sum(dur)/60000,3)) %>%
  full_join(all.rand.segments.sa, by = c("aclew_child_id", "segment", "SpkrAge")) %>%
  replace_na(list(segment_dur = 5, cds_min.sa = 0)) %>%
  mutate(cds_mph.sa = (cds_min.sa/segment_dur)*60) %>%
  arrange(aclew_child_id, segment, SpkrAge)
# Number of speakers per clip
spkrs.per.seg.rand.sa <- all.data %>%
  filter(sample == "random" & speaker != "CHI" & SpkrAge != "Not known" &
           !(grepl("@", tier))) %>%
  group_by(aclew_child_id, SpkrAge, segment) %>%
  summarise(n_spkrs_clip = length(unique(speaker)))
# All together
quantity.rand.sa <- xds.per.seg.rand.sa %>%
  full_join(ods.per.seg.rand.sa, by = c("aclew_child_id", "SpkrAge",
                                        "segment", "segment_dur")) %>%
  full_join(tds.per.seg.rand.sa, by = c("aclew_child_id", "SpkrAge",
                                        "segment", "segment_dur")) %>%
  full_join(cds.per.seg.rand.sa, by = c("aclew_child_id", "SpkrAge",
                                        "segment", "segment_dur")) %>%
  full_join(dplyr::select(quantity.rand, c("aclew_child_id", "segment", "tds_min")),
            by = c("aclew_child_id", "segment")) %>%
  full_join(spkrs.per.seg.rand.sa, by = c("aclew_child_id", "SpkrAge", "segment")) %>%
  left_join(dplyr::select(seg.info, c("aclew_id", "CodeName", "start.hr")),
            by = c("aclew_child_id" = "aclew_id", "segment" = "CodeName")) %>%
  full_join(ptcp.info, by = "aclew_child_id") %>%
  replace_na(list(xds_min.sa = 0, xds_mph.sa = 0,
                  ods_min.sa = 0, ods_mph.sa = 0,
                  tds_min.sa = 0, tds_mph.sa = 0,
                  cds_min.sa = 0, cds_mph.sa = 0,
                  n_spkrs_clip = 0)) %>%
  mutate(prop_tds.sa = tds_min.sa/xds_min.sa,
         prop_sa.tds = tds_min.sa/tds_min)
  # Don't replace NAs with 0s in this case; proportion is not meaningful w/o any speech


# NON-RANDOM
# Get min/hr speech measures
all.nonrand.segments <- seg.info %>%
  filter(!(grepl("random", CodeName))) %>%
  select(aclew_id, CodeName) %>%
  rename(aclew_child_id = aclew_id, segment = CodeName)
all.nonrand.segments$sample <- ifelse(grepl("va", all.nonrand.segments$segment),
                                      "high-activity","turn-taking")
all.nonrand.segments$sample_type <- ifelse(
  grepl("^va", all.nonrand.segments$segment), "high-activity",
  ifelse(grepl("^tt", all.nonrand.segments$segment),"turn-taking",
                                        "extension"))
# XDS
xds.per.seg.nonrand <- all.data %>%
  filter(sample != "random" & speaker != "CHI" &
           grepl("xds@", tier)) %>%
  group_by(aclew_child_id, segment, segment_dur) %>%
  summarise(xds_min = round(sum(dur)/60000,3)) %>%
  full_join(all.nonrand.segments, by = c("aclew_child_id", "segment")) %>%
  replace_na(list(segment_dur = 1, xds_min = 0)) %>%
  mutate(segment_dur = ifelse(grepl("ext", segment), 5, 1),
         xds_mph = (xds_min/segment_dur)*60) %>%
  arrange(aclew_child_id, segment)
# ODS
ods.per.seg.nonrand <- all.data %>%
  filter(sample != "random" & speaker != "CHI" &
           grepl("xds@", tier) & val != "T") %>%
  group_by(aclew_child_id, segment, segment_dur) %>%
  summarise(ods_min = round(sum(dur)/60000,3)) %>%
  full_join(all.nonrand.segments, by = c("aclew_child_id", "segment")) %>%
  replace_na(list(segment_dur = 1, ods_min = 0)) %>%
  mutate(segment_dur = ifelse(grepl("ext", segment), 5, 1),
         ods_mph = (ods_min/segment_dur)*60) %>%
  arrange(aclew_child_id, segment)
# TDS
tds.per.seg.nonrand <- all.data %>%
  filter(sample != "random" & speaker != "CHI" &
           grepl("xds@", tier) & val == "T") %>%
  group_by(aclew_child_id, segment, segment_dur) %>%
  summarise(tds_min = round(sum(dur)/60000,3)) %>%
  full_join(all.nonrand.segments, by = c("aclew_child_id", "segment")) %>%
  replace_na(list(segment_dur = 1, tds_min = 0)) %>%
  mutate(segment_dur = ifelse(grepl("ext", segment), 5, 1),
         tds_mph = (tds_min/segment_dur)*60) %>%
  arrange(aclew_child_id, segment)
# All CDS
cds.per.seg.nonrand <- all.data %>%
  filter(sample != "random" & speaker != "CHI" &
           grepl("xds@", tier) & (val == "T" | val == "C")) %>%
  group_by(aclew_child_id, segment, segment_dur) %>%
  summarise(cds_min = round(sum(dur)/60000,3)) %>%
  full_join(all.nonrand.segments, by = c("aclew_child_id", "segment")) %>%
  replace_na(list(segment_dur = 1, cds_min = 0)) %>%
  mutate(segment_dur = ifelse(grepl("ext", segment), 5, 1),
         cds_mph = (cds_min/segment_dur)*60) %>%
  arrange(aclew_child_id, segment)
# Number of speakers per clip
spkrs.per.seg.nonrand <- all.data %>%
  filter(sample != "random" & speaker != "CHI" & SpkrAge != "Not known" &
           !(grepl("@", tier))) %>%
  group_by(aclew_child_id, segment) %>%
  summarise(n_spkrs_clip = length(unique(speaker)))
# All together
quantity.nonrand <- xds.per.seg.nonrand %>%
  full_join(ods.per.seg.nonrand, by = c("aclew_child_id", "segment",
                                        "segment_dur", "sample", "sample_type")) %>%
  full_join(tds.per.seg.nonrand, by = c("aclew_child_id", "segment",
                                        "segment_dur", "sample", "sample_type")) %>%
  full_join(cds.per.seg.nonrand, by = c("aclew_child_id", "segment",
                                        "segment_dur", "sample", "sample_type")) %>%
  full_join(ptcp.info, by = "aclew_child_id") %>%
  full_join(spkrs.per.seg.nonrand, by = c("aclew_child_id", "segment")) %>%
  left_join(dplyr::select(seg.info, c("aclew_id", "CodeName", "start.hr")),
            by = c("aclew_child_id" = "aclew_id", "segment" = "CodeName")) %>%
  replace_na(list(xds_min = 0, xds_mph = 0,
                  ods_min = 0, ods_mph = 0,
                  tds_min = 0, tds_mph = 0,
                  cds_min = 0, cds_mph = 0,
                  n_spkrs_clip = 0)) %>%
  mutate(prop_tds = tds_min/xds_min)
  # Don't replace NAs with 0s in this case; proportion is not meaningful w/o any speech
quantity.nonrand.bychild <- quantity.nonrand %>%
  group_by(aclew_child_id, sample) %>%
  summarise(
    xds_min = mean(xds_min),
    xds_mph = mean(xds_mph),
    ods_min = mean(ods_min),
    ods_mph = mean(ods_mph),
    tds_min = mean(tds_min),
    tds_mph = mean(tds_mph),
    cds_min = mean(cds_min),
    cds_mph = mean(cds_mph),
    prop_tds = mean(prop_tds, na.rm = TRUE),
    m_n_spkrs = mean(n_spkrs_clip)) %>%
  full_join(ptcp.info, by = "aclew_child_id")

# Get xds and tds min/hr by speaker type
all.nonrand.segments.sa <- tibble(
  aclew_child_id = rep(all.nonrand.segments$aclew_child_id, 2),
  segment = rep(all.nonrand.segments$segment, 2),
  sample = rep(all.nonrand.segments$sample, 2),
  sample_type = rep(all.nonrand.segments$sample_type, 2),
  SpkrAge = c(rep("Adult", nrow(all.nonrand.segments)),
              rep("Child", nrow(all.nonrand.segments))))
# XDS
xds.per.seg.nonrand.sa <- all.data %>%
  filter(sample != "random" & speaker != "CHI" & SpkrAge != "Not known" &
           grepl("xds@", tier)) %>%
  group_by(aclew_child_id, SpkrAge, segment, segment_dur) %>%
  summarise(xds_min.sa = round(sum(dur)/60000,3)) %>%
  full_join(all.nonrand.segments.sa, by = c("aclew_child_id", "segment", "SpkrAge")) %>%
  replace_na(list(segment_dur = 1, xds_min.sa = 0)) %>%
  mutate(segment_dur = ifelse(grepl("ext", segment), 5, 1),
         xds_mph.sa = (xds_min.sa/segment_dur)*60) %>%
  arrange(aclew_child_id, segment, SpkrAge)
# ODS
ods.per.seg.nonrand.sa <- all.data %>%
  filter(sample != "random" & speaker != "CHI" & SpkrAge != "Not known" &
           grepl("xds@", tier) & val != "T") %>%
  group_by(aclew_child_id, SpkrAge, segment, segment_dur) %>%
  summarise(ods_min.sa = round(sum(dur)/60000,3)) %>%
  full_join(all.nonrand.segments.sa, by = c("aclew_child_id", "segment", "SpkrAge")) %>%
  replace_na(list(segment_dur = 1, ods_min.sa = 0)) %>%
  mutate(segment_dur = ifelse(grepl("ext", segment), 5, 1),
         ods_mph.sa = (ods_min.sa/segment_dur)*60) %>%
  arrange(aclew_child_id, segment, SpkrAge)
# TDS
tds.per.seg.nonrand.sa <- all.data %>%
  filter(sample != "random" & speaker != "CHI" & SpkrAge != "Not known" &
           grepl("xds@", tier) & val == "T") %>%
  group_by(aclew_child_id, SpkrAge, segment, segment_dur) %>%
  summarise(tds_min.sa = round(sum(dur)/60000,3)) %>%
  full_join(all.nonrand.segments.sa, by = c("aclew_child_id", "segment", "SpkrAge")) %>%
  replace_na(list(segment_dur = 1, tds_min.sa = 0)) %>%
  mutate(segment_dur = ifelse(grepl("ext", segment), 5, 1),
         tds_mph.sa = (tds_min.sa/segment_dur)*60) %>%
  arrange(aclew_child_id, segment, SpkrAge)
# All CDS
cds.per.seg.nonrand.sa <- all.data %>%
  filter(sample != "random" & speaker != "CHI" & SpkrAge != "Not known" &
           grepl("xds@", tier) & (val == "T" | val == "C")) %>%
  group_by(aclew_child_id, SpkrAge, segment, segment_dur) %>%
  summarise(cds_min.sa = round(sum(dur)/60000,3)) %>%
  full_join(all.nonrand.segments.sa, by = c("aclew_child_id", "segment", "SpkrAge")) %>%
  replace_na(list(segment_dur = 1, cds_min.sa = 0)) %>%
  mutate(segment_dur = ifelse(grepl("ext", segment), 5, 1),
         cds_mph.sa = (cds_min.sa/segment_dur)*60) %>%
  arrange(aclew_child_id, segment, SpkrAge)
# Number of speakers per clip
spkrs.per.seg.nonrand.sa <- all.data %>%
  filter(sample != "random" & speaker != "CHI" & SpkrAge != "Not known" &
           !(grepl("@", tier))) %>%
  group_by(aclew_child_id, SpkrAge, segment) %>%
  summarise(n_spkrs_clip = length(unique(speaker)))
# All together
quantity.nonrand.sa <- xds.per.seg.nonrand.sa %>%
  full_join(ods.per.seg.nonrand.sa, by = c("aclew_child_id", "SpkrAge",
                                           "segment", "segment_dur",
                                           "sample", "sample_type")) %>%
  full_join(tds.per.seg.nonrand.sa, by = c("aclew_child_id", "SpkrAge",
                                           "segment", "segment_dur",
                                           "sample", "sample_type")) %>%
  full_join(cds.per.seg.nonrand.sa, by = c("aclew_child_id", "SpkrAge",
                                           "segment", "segment_dur",
                                           "sample", "sample_type")) %>%
  full_join(select(quantity.nonrand, c("aclew_child_id", "segment", "tds_min")),
            by = c("aclew_child_id", "segment")) %>%
  full_join(ptcp.info, by = "aclew_child_id") %>%
  full_join(spkrs.per.seg.nonrand.sa, by = c("aclew_child_id", "SpkrAge", "segment")) %>%
  left_join(dplyr::select(seg.info, c("aclew_id", "CodeName", "start.hr")),
            by = c("aclew_child_id" = "aclew_id", "segment" = "CodeName")) %>%
  replace_na(list(xds_min.sa = 0, xds_mph.sa = 0,
                  ods_min.sa = 0, ods_mph.sa = 0,
                  tds_min.sa = 0, tds_mph.sa = 0,
                  cds_min.sa = 0, cds_mph.sa = 0,
                  n_spkrs_clip = 0)) %>%
  mutate(prop_tds.sa = tds_min.sa/xds_min.sa,
         prop_sa.tds = tds_min.sa/tds_min)
  # Don't replace NAs with 0s in this case; proportion is not meaningful w/o any speech

# Subset the non-random samples (used for differnt purposes)
quantity.nonrand.tt <- filter(quantity.nonrand, sample == "turn-taking")
quantity.nonrand.tt.sa <- filter(quantity.nonrand.sa, sample == "turn-taking")
quantity.nonrand.va <- filter(quantity.nonrand, sample != "turn-taking")
quantity.nonrand.va.sa <- filter(quantity.nonrand.sa, sample != "turn-taking")


## Get variables ready for modeling
# random sample
quantity.rand$child_sex <- as.factor(quantity.rand$child_sex)
quantity.rand$mat_ed <- as.factor(quantity.rand$mat_ed)
nspkrs.m <- mean(quantity.rand$n_spkrs_clip)
nspkrs.sd <- sd(quantity.rand$n_spkrs_clip)
quantity.rand <- quantity.rand %>%
  mutate(
    xds_mph.nz = ifelse(xds_mph > 0, 1, 0),
    ods_mph.nz = ifelse(ods_mph > 0, 1, 0),
    tds_mph.nz = ifelse(tds_mph > 0, 1, 0),
    cds_mph.nz = ifelse(cds_mph > 0, 1, 0),
    tchiyr.std = ((age_mo_round - tchiyr.m)/tchiyr.sd),
    chisx.std = recode_factor(child_sex,
                              "M" = "M", "F" = "F"),
    mated.std = recode_factor(mat_ed,
                              "none" = "none", "primary" = "primary",
                              "secondary" = "secondary", "preparatory" = "preparatory"),
    mated.bin = recode_factor(mat_ed,
                              "none" = "0-5", "primary" = "0-5",
                              "secondary" = "6+", "preparatory" = "6+"),
    motyr.std = ((mother_age - motyr.m)/motyr.sd),
    nsb.std = ((number_older_sibs - nsb.m)/nsb.sd),
    hsz.std = ((household_size - hsz.m)/hsz.sd),
    nsk.std = ((n_spkrs_clip - nspkrs.m)/nspkrs.sd),
    stthr.std = (start.hr - 12)/12,
    stthr.tri = ifelse(start.hr < 11, "morning",
                       ifelse(start.hr > 13, "afternoon", "midday")))
quantity.rand$stthr.tri <- factor(quantity.rand$stthr.tri,
                                  levels = c("midday", "morning", "afternoon"))
quantity.rand$stthr.tri.a <- factor(quantity.rand$stthr.tri,
                                  levels = c("afternoon", "midday", "morning"))

quantity.rand.sa$child_sex <- as.factor(quantity.rand.sa$child_sex)
quantity.rand.sa$mat_ed <- as.factor(quantity.rand.sa$mat_ed)
nspkrs.sa.m <- mean(quantity.rand.sa$n_spkrs_clip)
nspkrs.sa.sd <- sd(quantity.rand.sa$n_spkrs_clip)
quantity.rand.sa <- quantity.rand.sa %>%
  mutate(
    xds_mph.sa.nz = ifelse(xds_mph.sa > 0, 1, 0),
    ods_mph.sa.nz = ifelse(ods_mph.sa > 0, 1, 0),
    tds_mph.sa.nz = ifelse(tds_mph.sa > 0, 1, 0),
    cds_mph.sa.nz = ifelse(cds_mph.sa > 0, 1, 0),
    tchiyr.std = ((age_mo_round - tchiyr.m)/tchiyr.sd),
    chisx.std = recode_factor(child_sex,
                              "M" = "M", "F" = "F"),
    mated.std = recode_factor(mat_ed,
                              "none" = "none", "primary" = "primary",
                              "secondary" = "secondary", "preparatory" = "preparatory"),
    mated.bin = recode_factor(mat_ed,
                              "none" = "0-5", "primary" = "0-5",
                              "secondary" = "6+", "preparatory" = "6+"),
    motyr.std = ((mother_age - motyr.m)/motyr.sd),
    nsb.std = ((number_older_sibs - nsb.m)/nsb.sd),
    hsz.std = ((household_size - hsz.m)/hsz.sd),
    nsk.std = ((n_spkrs_clip - nspkrs.sa.m)/nspkrs.sa.sd),
    stthr.std = (start.hr - 12)/12,
    stthr.tri = ifelse(start.hr < 11, "morning",
                       ifelse(start.hr > 13, "afternoon", "midday")))
quantity.rand.sa$stthr.tri <- factor(quantity.rand.sa$stthr.tri,
                                     levels = c("midday", "morning", "afternoon"))
quantity.rand.sa$stthr.tri.a <- factor(quantity.rand.sa$stthr.tri,
                                       levels = c("afternoon", "midday", "morning"))


# tt sample
quantity.nonrand.tt$child_sex <- as.factor(quantity.nonrand.tt$child_sex)
quantity.nonrand.tt$mat_ed <- as.factor(quantity.nonrand.tt$mat_ed)
nspkrs.m <- mean(quantity.nonrand.tt$n_spkrs_clip)
nspkrs.sd <- sd(quantity.nonrand.tt$n_spkrs_clip)
quantity.nonrand.tt <- quantity.nonrand.tt %>%
  mutate(
    xds_mph.nz = ifelse(xds_mph > 0, 1, 0),
    ods_mph.nz = ifelse(ods_mph > 0, 1, 0),
    tds_mph.nz = ifelse(tds_mph > 0, 1, 0),
    cds_mph.nz = ifelse(cds_mph > 0, 1, 0),
    tchiyr.std = ((age_mo_round - tchiyr.m)/tchiyr.sd),
    chisx.std = recode_factor(child_sex,
                              "M" = "M", "F" = "F"),
    mated.std = recode_factor(mat_ed,
                              "none" = "none", "primary" = "primary",
                              "secondary" = "secondary", "preparatory" = "preparatory"),
    mated.bin = recode_factor(mat_ed,
                              "none" = "0-5", "primary" = "0-5",
                              "secondary" = "6+", "preparatory" = "6+"),
    motyr.std = ((mother_age - motyr.m)/motyr.sd),
    nsb.std = ((number_older_sibs - nsb.m)/nsb.sd),
    hsz.std = ((household_size - hsz.m)/hsz.sd),
    nsk.std = ((n_spkrs_clip - nspkrs.m)/nspkrs.sd),
    stthr.std = (start.hr - 12)/12,
    stthr.tri = ifelse(start.hr < 11, "morning",
                       ifelse(start.hr > 13, "afternoon", "midday")))
quantity.nonrand.tt$stthr.tri <- factor(quantity.nonrand.tt$stthr.tri,
                                        levels = c("midday", "morning", "afternoon"))
quantity.nonrand.tt$stthr.tri.a <- factor(quantity.nonrand.tt$stthr.tri,
                                          levels = c("afternoon", "midday", "morning"))

quantity.nonrand.tt.sa$child_sex <- as.factor(quantity.nonrand.tt.sa$child_sex)
quantity.nonrand.tt.sa$mat_ed <- as.factor(quantity.nonrand.tt.sa$mat_ed)
nspkrs.sa.m <- mean(quantity.nonrand.tt.sa$n_spkrs_clip)
nspkrs.sa.sd <- sd(quantity.nonrand.tt.sa$n_spkrs_clip)
quantity.nonrand.tt.sa <- quantity.nonrand.tt.sa %>%
  mutate(
    xds_mph.sa.nz = ifelse(xds_mph.sa > 0, 1, 0),
    ods_mph.sa.nz = ifelse(ods_mph.sa > 0, 1, 0),
    tds_mph.sa.nz = ifelse(tds_mph.sa > 0, 1, 0),
    cds_mph.sa.nz = ifelse(cds_mph.sa > 0, 1, 0),
    tchiyr.std = ((age_mo_round - tchiyr.m)/tchiyr.sd),
    chisx.std = recode_factor(child_sex,
                              "M" = "M", "F" = "F"),
    mated.std = recode_factor(mat_ed,
                              "none" = "none", "primary" = "primary",
                              "secondary" = "secondary", "preparatory" = "preparatory"),
    mated.bin = recode_factor(mat_ed,
                              "none" = "0-5", "primary" = "0-5",
                              "secondary" = "6+", "preparatory" = "6+"),
    motyr.std = ((mother_age - motyr.m)/motyr.sd),
    nsb.std = ((number_older_sibs - nsb.m)/nsb.sd),
    hsz.std = ((household_size - hsz.m)/hsz.sd),
    nsk.std = ((n_spkrs_clip - nspkrs.sa.m)/nspkrs.sa.sd),
    stthr.std = (start.hr - 12)/12,
    stthr.tri = ifelse(start.hr < 11, "morning",
                       ifelse(start.hr > 13, "afternoon", "midday")))
quantity.nonrand.tt.sa$stthr.tri <- factor(quantity.nonrand.tt.sa$stthr.tri,
                                           levels = c("midday", "morning", "afternoon"))
quantity.nonrand.tt.sa$stthr.tri.a <- factor(quantity.nonrand.tt.sa$stthr.tri,
                                             levels = c("afternoon", "midday", "morning"))


quantity.nonrand.tt.minimum <- dplyr::select(quantity.nonrand.tt,
                                             age_mo_round, xds_mph, ods_mph, tds_mph,
                                             prop_tds, n_spkrs_clip) %>%
                                             mutate(Sample = "Turn taking")
quantity.rand.minimum <- dplyr::select(quantity.rand,
                                       age_mo_round, xds_mph, ods_mph, tds_mph,
                                       prop_tds, n_spkrs_clip) %>%
                                       mutate(Sample = "Random")
quantity.rand_and_tt <- bind_rows(quantity.nonrand.tt.minimum, quantity.rand.minimum)

quantity.nonrand.sa.tt.minimum <- dplyr::select(quantity.nonrand.tt.sa,
                                             age_mo_round, prop_sa.tds, SpkrAge) %>%
                                             mutate(Sample = "Turn taking")
quantity.rand.sa.minimum <- dplyr::select(quantity.rand.sa,
                                       age_mo_round, prop_sa.tds, SpkrAge) %>%
                                       mutate(Sample = "Random")
quantity.sa.rand_and_tt <- bind_rows(
  quantity.nonrand.sa.tt.minimum, quantity.rand.sa.minimum)

write_csv(quantity.rand_and_tt,
          paste0(shiny.input.path,
                 "quantity-scores_rand-and-tt.csv"))

# ODS min/hr
odsmph.segments.rand_and_tt <- ggplot(quantity.rand_and_tt,
                          aes(x = age_mo_round, y = ods_mph, lty = Sample)) +
  geom_boxplot(aes(group = interaction(age_mo_round, Sample),
                   color = Sample), fill = "white", outlier.shape = NA,
               lty = "solid", alpha = 0.4) +
  geom_smooth(aes(fill = Sample, color = Sample), method = "lm") +
  ylab("ODS (min/hr)") + xlab("Child age (mo)")	+
  scale_y_continuous(limits=c(-10,80),
                     breaks=seq(0,80,20)) +
  scale_x_continuous(limits=c(0,38),
                     breaks=seq(0,38,6)) +
  coord_cartesian(ylim=c(0,80),xlim=c(0,38)) +
  scale_color_manual(values = viridis(3)) +
  scale_fill_manual(values = viridis(3)) +
  theme_apa() +
  theme(legend.position="none",
        axis.line = element_line(color="black", size = 0.4))

# TDS min/hr
tdsmph.segments.rand_and_tt <- ggplot(quantity.rand_and_tt,
                          aes(x = age_mo_round, y = tds_mph, lty = Sample)) +
  geom_boxplot(aes(group = interaction(age_mo_round, Sample),
                   color = Sample), fill = "white", outlier.shape = NA,
               lty = "solid", alpha = 0.4) +
  geom_smooth(aes(fill = Sample, color = Sample), method = "lm") +
  ylab("TCDS (min/hr)") + xlab("Child age (mo)")	+
  scale_y_continuous(limits=c(0,80),
                     breaks=seq(0,80,20)) +
  scale_x_continuous(limits=c(0,38),
                     breaks=seq(0,38,6)) +
  coord_cartesian(ylim=c(0,80),xlim=c(0,38)) +
  scale_color_manual(values = viridis(3)) +
  scale_fill_manual(values = viridis(3)) +
  theme_apa() +
  theme(legend.position="none",
        axis.line = element_line(color="black", size = 0.4))

# TDS min/hr - zoomed in
tdsmph.segments.rand_and_tt.zoomedin <- ggplot(quantity.rand_and_tt,
                          aes(x = age_mo_round, y = tds_mph, lty = Sample)) +
  geom_boxplot(aes(group = interaction(age_mo_round, Sample),
                   color = Sample), fill = "white", outlier.shape = NA,
               lty = "solid", alpha = 0.4) +
  geom_smooth(aes(fill = Sample, color = Sample), method = "lm") +
  ylab("TCDS (min/hr)") + xlab("Child age (mo)")	+
  scale_y_continuous(limits=c(0,40),
                     breaks=seq(0,40,10)) +
  scale_x_continuous(limits=c(0,38),
                     breaks=seq(0,38,6)) +
  coord_cartesian(ylim=c(0,40),xlim=c(0,38)) +
  scale_color_manual(values = viridis(3)) +
  scale_fill_manual(values = viridis(3)) +
  theme_apa() +
  theme(legend.position="none",
        axis.line = element_line(color="black", size = 0.4))

# TDS prop
tdsprp.segments.rand_and_tt <- ggplot(quantity.rand_and_tt,
                          aes(x = age_mo_round, y = prop_tds, lty = Sample)) +
  geom_boxplot(aes(group = interaction(age_mo_round, Sample),
                   color = Sample), fill = "white", outlier.shape = NA,
               lty = "solid", alpha = 0.4) +
  geom_smooth(aes(fill = Sample, color = Sample), method = "lm") +
  ylab("TCDS/All spch") + xlab("Child age (mo)")	+
  scale_y_continuous(limits=c(-.2,1.2),
                     breaks=seq(0,1,0.2)) +
  scale_x_continuous(limits=c(0,38),
                     breaks=seq(0,38,6)) +
  coord_cartesian(ylim=c(0,1),xlim=c(0,38)) +
  scale_color_manual(values = viridis(3)) +
  scale_fill_manual(values = viridis(3)) +
  theme_apa() +
  theme(legend.position="none",
        axis.line = element_line(color="black", size = 0.4))

# prop TDS from children
tdsprp.segments.rand_and_tt.sa <- ggplot(subset(quantity.sa.rand_and_tt,
                                                SpkrAge == "Child"),
                             aes(x = age_mo_round, y = prop_sa.tds, lty = Sample)) +
  geom_boxplot(aes(group = interaction(age_mo_round, Sample),
                   color = Sample), fill = "white", outlier.shape = NA,
               lty = "solid", alpha = 0.4) +
  geom_smooth(aes(fill = Sample, color = Sample), method = "lm") +
  ylab("Prop of TCDS") + xlab("Child age (mo)")	+
  scale_y_continuous(limits=c(-.2,1.2),
                     breaks=seq(0,1,0.2)) +
  scale_x_continuous(limits=c(0,38),
                     breaks=seq(0,38,6)) +
  coord_cartesian(ylim=c(0,1), xlim=c(0,38)) +
  scale_color_manual(values = viridis(3)) +
  scale_fill_manual(values = viridis(3)) +
  theme_apa() +
  theme(legend.position="none",
        axis.line = element_line(color="black", size = 0.4))

## TOD plots
quantity.nonrand.tt$Sample <- "Turn taking"
quantity.rand$Sample <- "Random"
quantity.rand_and_tt.all <- bind_rows(quantity.nonrand.tt, quantity.rand)
quantity.rand_and_tt.all$stthr.tri.centered <- factor(
  quantity.rand_and_tt.all$stthr.tri,
  levels = c("morning", "midday", "afternoon"))
quantity.rand_and_tt.all$SplitAge <- ifelse(
  quantity.rand_and_tt.all$age_mo_round < 13, "Under 1;0", "1;0+")
quantity.rand_and_tt.all$SplitAge <- factor(quantity.rand_and_tt.all$SplitAge,
  levels = c("Under 1;0", "1;0+"))

write_csv(quantity.rand_and_tt.all,
          paste0(shiny.input.path,
                 "quantity-scores_rand-and-tt-TOD.csv"))

tod.tcds <- ggplot(data = quantity.rand_and_tt.all, aes(
  y = round(tds_mph, 0), x = stthr.tri.centered,
  group = interaction(SplitAge, stthr.tri.centered),
  color = Sample,
  fill = Sample,
  alpha = SplitAge)) +
  geom_jitter() +
  geom_boxplot(color = "black") +
  facet_grid(~ Sample) +
  ylab("TCDS (min/hr)") + xlab("Time of day")	+
  labs(alpha = "Age") +
  scale_y_continuous(limits=c(-10,40),
                     breaks=c(0,10,20, 30, 40)) +
  scale_color_manual(guide = FALSE, values = viridis(3)) +
  scale_fill_manual(guide = FALSE, values = viridis(3)) +
  scale_alpha_discrete(range = c(0.2, 0.8)) +
  coord_cartesian(ylim=c(0,40)) +
  theme_apa() +
  theme(legend.position="right",
        axis.line = element_line(color="black", size = 0.4)) +
  guides(alpha = guide_legend(override.aes = list(
    fill = c("gray80", "gray30"))))

tod.ods <- ggplot(data = quantity.rand_and_tt.all, aes(
  y = round(ods_mph, 0), x = stthr.tri.centered,
  group = interaction(SplitAge, stthr.tri.centered),
  color = Sample,
  fill = Sample,
  alpha = SplitAge)) +
  geom_jitter() +
  geom_boxplot(color = "black") +
  facet_grid(~ Sample) +
  ylab("ODS (min/hr)") + xlab("Time of day")	+
  labs(alpha = "Age") +
  scale_y_continuous(limits=c(-10,80),
                     breaks=seq(0,80,20)) +
  scale_color_manual(guide = FALSE, values = viridis(3)) +
  scale_fill_manual(guide = FALSE, values = viridis(3)) +
  scale_alpha_discrete(range = c(0.2, 0.8)) +
  coord_cartesian(ylim=c(0,80)) +
  theme_apa() +
  theme(legend.position="right",
        axis.line = element_line(color="black", size = 0.4)) +
  guides(alpha = guide_legend(override.aes = list(
    fill = c("gray80", "gray30"))))

# individual TOD plots
tod.tcds.rand <- ggplot(
  data = subset(quantity.rand_and_tt.all, Sample == "Random"), aes(
  y = round(tds_mph, 0), x = stthr.tri.centered,
  group = interaction(SplitAge, stthr.tri.centered),
  color = Sample,
  fill = Sample,
  alpha = SplitAge)) +
  geom_jitter() +
  geom_boxplot(color = "black") +
  ylab("TCDS (min/hr)") + xlab("Time of day")	+
  labs(alpha = "Age") +
  scale_y_continuous(limits=c(-10,40),
                     breaks=c(0,10,20, 30, 40)) +
  scale_color_manual(guide = FALSE, values = viridis(3)) +
  scale_fill_manual(guide = FALSE, values = viridis(3)) +
  scale_alpha_discrete(range = c(0.2, 0.8)) +
  coord_cartesian(ylim=c(0,40)) +
  theme_apa() +
  theme(legend.position="none",
        axis.line = element_line(color="black", size = 0.4))

tod.ods.rand <- ggplot(
  data = subset(quantity.rand_and_tt.all, Sample == "Random"), aes(
  y = round(ods_mph, 0), x = stthr.tri.centered,
  group = interaction(SplitAge, stthr.tri.centered),
  color = Sample,
  fill = Sample,
  alpha = SplitAge)) +
  geom_jitter() +
  geom_boxplot(color = "black") +
  ylab("ODS (min/hr)") + xlab("Time of day")	+
  labs(alpha = "Age") +
  scale_y_continuous(limits=c(-10,80),
                     breaks=seq(0,80,20)) +
  scale_color_manual(guide = FALSE, values = viridis(3)) +
  scale_fill_manual(guide = FALSE, values = viridis(3)) +
  scale_alpha_discrete(range = c(0.2, 0.8)) +
  coord_cartesian(ylim=c(0,80)) +
  theme_apa() +
  theme(legend.position="none",
        axis.line = element_line(color="black", size = 0.4))

tod.tcds.tt <- ggplot(
  data = subset(quantity.rand_and_tt.all, Sample == "Turn taking"), aes(
  y = round(tds_mph, 0), x = stthr.tri.centered,
  group = interaction(SplitAge, stthr.tri.centered),
  color = Sample,
  fill = Sample,
  alpha = SplitAge)) +
  geom_jitter() +
  geom_boxplot(color = "black") +
  ylab("TCDS (min/hr)") + xlab("Time of day")	+
  labs(alpha = "Age") +
  scale_y_continuous(limits=c(-10,40),
                     breaks=c(0,10,20, 30, 40)) +
  scale_color_manual(guide = FALSE,
                     values = c("#21908CFF", "#21908CFF", "#FDE725FF")) +
  scale_fill_manual(guide = FALSE,
                    values = c("#21908CFF", "#21908CFF", "#FDE725FF")) +
  scale_alpha_discrete(range = c(0.2, 0.8)) +
  coord_cartesian(ylim=c(0,40)) +
  theme_apa() +
  theme(legend.position="none",
        axis.line = element_line(color="black", size = 0.4))

tod.ods.tt <- ggplot(
  data = subset(quantity.rand_and_tt.all, Sample == "Turn taking"), aes(
  y = round(ods_mph, 0), x = stthr.tri.centered,
  group = interaction(SplitAge, stthr.tri.centered),
  color = Sample,
  fill = Sample,
  alpha = SplitAge)) +
  geom_jitter() +
  geom_boxplot(color = "black") +
  ylab("ODS (min/hr)") + xlab("Time of day")	+
  labs(alpha = "Age") +
  scale_y_continuous(limits=c(-10,80),
                     breaks=seq(0,80,20)) +
  scale_color_manual(guide = FALSE,
                     values = c("#21908CFF", "#21908CFF", "#FDE725FF")) +
  scale_fill_manual(guide = FALSE,
                    values = c("#21908CFF", "#21908CFF", "#FDE725FF")) +
  scale_alpha_discrete(range = c(0.2, 0.8)) +
  coord_cartesian(ylim=c(0,80)) +
  theme_apa() +
  theme(legend.position="none",
        axis.line = element_line(color="black", size = 0.4))

grid.arrange(tdsmph.segments.rand_and_tt.zoomedin,
             odsmph.segments.rand_and_tt, nrow=1, ncol=2)

# Prep comparison data for matching figures
comparison.data <- read_csv(comparison.file)
tds.mph.tofill <- which(is.na(comparison.data$tds_mph) &
                          !(is.na(comparison.data$tds_uph)))
comparison.data$tds_mph[tds.mph.tofill] <- ((comparison.data$tds_uph[tds.mph.tofill] *
                                               avg.utt.len.tseltal$median.utt.len)/60000)
ods.mph.tofill <- which(is.na(comparison.data$ods_mph) &
                          !(is.na(comparison.data$ods_uph)))
comparison.data$ods_mph[ods.mph.tofill] <- ((comparison.data$ods_uph[ods.mph.tofill] *
                                               avg.utt.len.tseltal$median.utt.len)/60000)
xds.mph.tofill <- which(is.na(comparison.data$xds_mph) &
                          !(is.na(comparison.data$ods_mph)) &
                          !(is.na(comparison.data$tds_mph)))
comparison.data$xds_mph[xds.mph.tofill] <- comparison.data$ods_mph[xds.mph.tofill] +
  comparison.data$tds_mph[xds.mph.tofill]
# Add in current data
by.chi.rand.avgs <- quantity.rand %>%
  ungroup() %>%
  select(age_mo_round, ods_mph, cds_mph, tds_mph, xds_mph, prop_tds) %>%
  group_by(age_mo_round) %>%
  summarise(ods_mph = mean(ods_mph),
            cds_mph = mean(cds_mph),
            tds_mph = mean(tds_mph),
            xds_mph = mean(xds_mph),
            prp_tds = mean(prop_tds, na.rm = TRUE)) %>%
  rename(AgeMonths = age_mo_round) %>%
  mutate(ods_uph = NA,
         tds_uph = NA,
         xds_uph = NA,
         prp_tds_adu = NA,
         prp_tds_chi = NA,
         shape = 8,
         N = 1,
         Type = "Rural, non-WEIRD",
         Site = "Tseltal",
         Reference = "The current study",
         Notes = "")
comparison.data <- bind_rows(comparison.data, by.chi.rand.avgs)

comparison.data$Site <- factor(comparison.data$Site,
                               levels = unique(comparison.data$Site))
comparison.data$Site <- factor(comparison.data$Site, levels =
                                 c("US/Canada", "US", "Dutch", "MozambiqueUrb",
                                   "MozambiqueRur", "Tsimane", "Yukatek", "Tseltal"))
comparison.data.shapes <- unique(comparison.data$shape)

# TDS min/hr comparison
tdsmph.comparison <- ggplot(comparison.data,
                          aes(x = AgeMonths, y = tds_mph, stroke = 1,
                              color = Type, shape = Site, size = N)) +
  geom_point(fill = "gray50", alpha = 0.7) +
  scale_shape_manual(values=comparison.data.shapes) +
  ylab("TCDS (min/hr)") + xlab("Child age (mo)")	+
  scale_y_continuous(limits=c(0,20),
                     breaks=seq(0,20,5)) +
  scale_x_continuous(limits=c(0,38),
                     breaks=seq(0,38,6)) +
  coord_cartesian(ylim=c(0,20),xlim=c(0,38)) +
  scale_color_manual(values = c("black", "gray40")) +
  theme_apa() +
  guides(size = FALSE, color = FALSE) +
  theme(axis.line = element_line(color="black", size = 0.4))
tdsmph.comparison

## TCDS random sample ####
tcds.random.distribution <- ggplot(quantity.rand,
                            aes(round(tds_mph,0))) +
  geom_histogram() +
  ylab("# of clips") +
  xlab ("TCDS min/hr") +
  basic.theme
ggsave(paste0(shiny.input.img.path, "TCDS_random_distribution.png"),
       plot = tcds.random.distribution,
       width = 8, height = 6, dpi = 300)
#sd(round(quantity.rand$tds_mph,0))^2
#mean(round(quantity.rand$tds_mph,0))
# mean is much smaller than variance
tds.rand.zinb <- glmmTMB(round(tds_mph,0) ~
                           tchiyr.std +
                           stthr.tri +
                           hsz.std +
                           nsk.std +
                           tchiyr.std:stthr.tri +
                           tchiyr.std:hsz.std +
                           tchiyr.std:nsk.std +
                           (1|aclew_child_id),
                         data=quantity.rand,
                         ziformula=~nsk.std,
                         family="nbinom1")
#tds.rand.zinb.res = simulateResiduals(tds.rand.zinb)
#plot(tds.rand.zinb.res, rank = T) # (manually saved)
#summary(tds.rand.zinb)
#Conditional model:
#                              Estimate Std. Error z value Pr(>|z|)   
#(Intercept)                    0.82318    0.38815   2.121  0.03394 * 
#tchiyr.std                     0.44447    0.42262   1.052  0.29294   
#stthr.trimorning               0.81894    0.39733   2.061  0.03929 * 
#stthr.triafternoon             0.49070    0.37401   1.312  0.18952   
#hsz.std                       -0.08792    0.26286  -0.334  0.73802   
#nsk.std                       -0.12782    0.16269  -0.786  0.43206   
#tchiyr.std:stthr.trimorning   -0.23599    0.39305  -0.600  0.54823   
#tchiyr.std:stthr.triafternoon -0.81223    0.37874  -2.145  0.03199 * 
#tchiyr.std:hsz.std            -0.21337    0.32104  -0.665  0.50629   
#tchiyr.std:nsk.std             0.60579    0.19780   3.063  0.00219 **
#---
#Zero-inflation model:
#            Estimate Std. Error z value Pr(>|z|)
#(Intercept)   -56.90   14003.31  -0.004    0.997
#nsk.std       -55.17   14243.76  -0.004    0.997
# save for reporting
tds.rand.zinb.COEF.morn <-
  coef(summary(tds.rand.zinb))[[1]]["stthr.trimorning",] 
tds.rand.zinb.COEF.aft <-
  coef(summary(tds.rand.zinb))[[1]]["stthr.triafternoon",] 
tds.rand.zinb.COEF.agetime <-
  coef(summary(tds.rand.zinb))[[1]]["tchiyr.std:stthr.triafternoon",]
tds.rand.zinb.COEF.agensk <-
  coef(summary(tds.rand.zinb))[[1]]["tchiyr.std:nsk.std",] 

# test the other two-way effects of time of day
tds.rand.zinb.v2 <- glmmTMB(round(tds_mph,0) ~
                           tchiyr.std +
                           stthr.tri.a +
                           hsz.std +
                           nsk.std +
                           tchiyr.std:stthr.tri.a +
                           tchiyr.std:hsz.std +
                           tchiyr.std:nsk.std +
                           (1|aclew_child_id),
                         data=quantity.rand,
                         ziformula=~nsk.std,
                         family="nbinom1")
#summary(tds.rand.zinb.v2)
#stthr.tri.amidday              -0.4886     0.3787  -1.290       0.19700    
#stthr.tri.amorning              0.3044     0.2879   1.057       0.29048    
#tchiyr.std:stthr.tri.amidday    0.7326     0.3588   2.042       0.04114 *  
#tchiyr.std:stthr.tri.amorning   0.4552     0.2766   1.646       0.09984 .  
# save for reporting
tds.rand.zinb.v2.COEF.aft.vs.morn <-
  coef(summary(tds.rand.zinb.v2))[[1]]["stthr.tri.amorning",] 
tds.rand.zinb.v2.COEF.age.aft.vs.morn <-
  coef(summary(tds.rand.zinb.v2))[[1]]["tchiyr.std:stthr.tri.amorning",] 
tds.rand.zinb.v2.COEF.age.aft.vs.midd <-
  coef(summary(tds.rand.zinb.v2))[[1]]["tchiyr.std:stthr.tri.amidday",] 

## TCDS tt sample ####
tcds.tt.distribution <- ggplot(quantity.nonrand.tt,
                            aes(round(tds_mph,0))) +
  geom_histogram() +
  ylab("# of clips") +
  xlab ("TCDS min/hr") +
  basic.theme
ggsave(paste0(shiny.input.img.path, "TCDS_turntaking_distribution.png"),
       plot = tcds.tt.distribution,
       width = 8, height = 6, dpi = 300)
#sd(round(quantity.nonrand.tt$tds_mph,0))^2
#mean(round(quantity.nonrand.tt$tds_mph,0))
# mean is much smaller than variance
# not zero-inflated (nature of this sample)
tds.tt.nb <- glmmTMB(round(tds_mph,0) ~
                       tchiyr.std +
                       stthr.tri +
                       hsz.std +
                       nsk.std +
                       tchiyr.std:stthr.tri +
                       tchiyr.std:hsz.std +
                       tchiyr.std:nsk.std +
                       (1|aclew_child_id),
                     data=quantity.nonrand.tt,
                     family="nbinom1")
#tds.tt.nb.res = simulateResiduals(tds.tt.nb)
#plot(tds.tt.nb.res, rank = T) # (manually saved)
#summary(tds.tt.nb)
#Conditional model:
#                              Estimate Std. Error z value  Pr(>|z|)    
#(Intercept)                    2.18313    0.21298  10.250  <0.0001 ***
#tchiyr.std                    -0.22751    0.19781  -1.150  0.2501    
#stthr.trimorning               0.35533    0.23475   1.514  0.1301    
#stthr.triafternoon             0.14484    0.21939   0.660  0.5091    
#hsz.std                       -0.25282    0.15437  -1.638  0.1015    
#nsk.std                       -0.06740    0.09198  -0.733  0.4637    
#tchiyr.std:stthr.trimorning   -0.24161    0.22835  -1.058  0.2900    
#tchiyr.std:stthr.triafternoon -0.15237    0.20725  -0.735  0.4622    
#tchiyr.std:hsz.std            -0.46538    0.19241  -2.419  0.0156 *  
#tchiyr.std:nsk.std             0.04160    0.11620   0.358  0.7203    
#--
#Zero-inflation model:
#            Estimate Std. Error z value Pr(>|z|)    
#(Intercept)  -3.7900     1.0050  -3.771 0.000162 ***
#nsk.std       0.8721     0.7178   1.215 0.224397
# save for reporting
tds.tt.nb.COEF.age.hsz <-
  coef(summary(tds.tt.nb))[[1]]["tchiyr.std:hsz.std",]

# test the other two-way effects of time of day
tds.tt.nb.v2 <- glmmTMB(round(tds_mph,0) ~
                       tchiyr.std +
                       stthr.tri.a +
                       hsz.std +
                       nsk.std +
                       tchiyr.std:stthr.tri.a +
                       tchiyr.std:hsz.std +
                       tchiyr.std:nsk.std +
                       (1|aclew_child_id),
                     data=quantity.nonrand.tt,
                     family="nbinom1")
#summary(tds.tt.nb.v2)
#stthr.tri.amidday             -0.06471    0.23241  -0.278              0.7807    
#stthr.tri.amorning             0.26769    0.21593   1.240              0.2151    
#tchiyr.std:stthr.tri.amidday   0.03344    0.21400   0.156              0.8758    
#tchiyr.std:stthr.tri.amorning -0.23644    0.21540  -1.098              0.2723    
# no results to save for reporting

## TCDS random sample ####
tds.rand.gaus <- glmmTMB(log(tds_mph+1) ~
                           tchiyr.std +
                           stthr.tri +
                           hsz.std +
                           nsk.std +
                           tchiyr.std:stthr.tri +
                           tchiyr.std:hsz.std +
                           tchiyr.std:nsk.std +
                           (1|aclew_child_id),
                         data=quantity.rand)
#tds.rand.gaus.res = simulateResiduals(tds.rand.gaus)
#plot(tds.rand.gaus.res, rank = T) # (manually saved)
#summary(tds.rand.gaus)
#Conditional model:
#                              Estimate Std. Error z value Pr(>|z|)    
#(Intercept)                    0.77507    0.22523   3.441 0.000579 ***
#tchiyr.std                     0.49352    0.25950   1.902 0.057198 .  
#stthr.trimorning               0.50824    0.25050   2.029 0.042465 *  
#stthr.triafternoon             0.29315    0.22126   1.325 0.185203    
#hsz.std                       -0.19656    0.19595  -1.003 0.315809    
#nsk.std                        0.23224    0.11838   1.962 0.049776 *  
#tchiyr.std:stthr.trimorning   -0.16002    0.26992  -0.593 0.553276    
#tchiyr.std:stthr.triafternoon -0.68051    0.23874  -2.850 0.004366 ** 
#tchiyr.std:hsz.std            -0.08496    0.23915  -0.355 0.722381    
#tchiyr.std:nsk.std             0.24923    0.14810   1.683 0.092402 .  
# DIFFERENCES?
# broadly similar to z-i nb model

# test the other two-way effects of time of day
tds.rand.gaus.v2 <- glmmTMB(log(tds_mph+1) ~
                           tchiyr.std +
                           stthr.tri.a +
                           hsz.std +
                           nsk.std +
                           tchiyr.std:stthr.tri.a +
                           tchiyr.std:hsz.std +
                           tchiyr.std:nsk.std +
                           (1|aclew_child_id),
                         data=quantity.rand)
#summary(tds.rand.gaus.v2)
#tthr.tri.amidday             -0.29314    0.22126  -1.325       0.18521    
#stthr.tri.amorning             0.21510    0.21916   0.981       0.32635    
#tchiyr.std:stthr.tri.amidday   0.68050    0.23874   2.850       0.00437 ** 
#tchiyr.std:stthr.tri.amorning  0.52049    0.23242   2.239       0.02512 *  
# DIFFERENCES?
# broadly similar to z-i nb model

## TCDS tt sample ####
tds.tt.gaus <- glmmTMB(log(tds_mph+1) ~
                       tchiyr.std +
                       stthr.tri +
                       hsz.std +
                       nsk.std +
                       tchiyr.std:stthr.tri +
                       tchiyr.std:hsz.std +
                       tchiyr.std:nsk.std +
                       (1|aclew_child_id),
                     data=quantity.nonrand.tt)
#tds.tt.gaus.res = simulateResiduals(tds.tt.gaus)
#plot(tds.tt.gaus.res, rank = T) # (manually saved)
#summary(tds.tt.gaus)
#Conditional model:
#                               Estimate Std. Error z value  Pr(>|z|)    
#(Intercept)                    2.081105   0.243378   8.551  <0.0001 ***
#tchiyr.std                    -0.128497   0.233182  -0.551  0.5816    
#stthr.trimorning               0.380764   0.298155   1.277  0.2016    
#stthr.triafternoon             0.105890   0.267468   0.396  0.6922    
#hsz.std                       -0.147307   0.172597  -0.853  0.3934    
#nsk.std                       -0.077111   0.114714  -0.672  0.5015    
#tchiyr.std:stthr.trimorning   -0.344770   0.296548  -1.163  0.2450    
#tchiyr.std:stthr.triafternoon -0.004403   0.262779  -0.017  0.9866    
#tchiyr.std:hsz.std            -0.491258   0.219527  -2.238  0.0252 *  
#tchiyr.std:nsk.std             0.169108   0.149388   1.132  0.2576 
# DIFFERENCES?
# very similar to z-i nb model

# test the other two-way effects of time of day
tds.tt.gaus.v2 <- glmmTMB(log(tds_mph+1) ~
                       tchiyr.std +
                       stthr.tri.a +
                       hsz.std +
                       nsk.std +
                       tchiyr.std:stthr.tri.a +
                       tchiyr.std:hsz.std +
                       tchiyr.std:nsk.std +
                       (1|aclew_child_id),
                     data=quantity.nonrand.tt)
#summary(tds.tt.gaus.v2)
#stthr.tri.amidday             -0.105890   0.267468  -0.396              0.6922    
#stthr.tri.amorning             0.274873   0.265438   1.036              0.3004    
#tchiyr.std:stthr.tri.amidday   0.004403   0.262779   0.017              0.9866    
#tchiyr.std:stthr.tri.amorning -0.340366   0.277150  -1.228              0.2194
# DIFFERENCES?
# very similar to z-i nb model

# Write model results out for input to shiny
TCDS.models <- bind_rows(
  broom.mixed::tidy(tds.rand.zinb) %>%
    mutate(model = "TCDS_random_z-inb"),
  broom.mixed::tidy(tds.rand.zinb.v2) %>%
    mutate(model = "TCDS_random_z-inb.v2"),
  broom.mixed::tidy(tds.rand.gaus) %>%
    mutate(model = "TCDS_random_log_gaus"),
  broom.mixed::tidy(tds.rand.gaus.v2) %>%
    mutate(model = "TCDS_random_log_gaus.v2"),
  broom.mixed::tidy(tds.tt.nb) %>%
    mutate(model = "TCDS_turntaking_nb"),
  broom.mixed::tidy(tds.tt.nb.v2) %>%
    mutate(model = "TCDS_turntaking_nb.v2"),
  broom.mixed::tidy(tds.tt.gaus) %>%
    mutate(model = "TCDS_turntaking_log_gaus"),
  broom.mixed::tidy(tds.tt.gaus.v2) %>%
    mutate(model = "TCDS_turntaking_log_gaus.v2"))


grid.arrange(tod.tcds.rand, tod.ods.rand,
             tod.tcds.tt, tod.ods.tt, nrow=2, ncol=2)

# Aggregate to one point per child and then test a correlation with age
# random
propchitcds.rand <- subset(quantity.rand.sa[
  which(!is.na(quantity.rand.sa$prop_sa.tds)),],
  SpkrAge == "Child")
propchitcds.rand.corr_age <- propchitcds.rand %>%
  group_by(aclew_child_id, tchiyr.std, age_mo_round) %>%
  summarise(avg_prpchitcds = mean(prop_sa.tds))
propchitcds.rand.corr_age.test <- cor.test(
  ~ age_mo_round + avg_prpchitcds,
  data = propchitcds.rand.corr_age, method = "spearman")
# turn taking
propchitcds.tt <- subset(quantity.nonrand.tt.sa[
  which(!is.na(quantity.nonrand.tt.sa$prop_sa.tds)),],
  SpkrAge == "Child")
propchitcds.tt.corr_age <- propchitcds.tt %>%
  group_by(aclew_child_id, tchiyr.std, age_mo_round) %>%
  summarise(avg_prpchitcds = mean(prop_sa.tds))
propchitcds.tt.corr_age.test <- cor.test(
  ~ age_mo_round + avg_prpchitcds,
  data = propchitcds.tt.corr_age, method = "spearman")


## ODS random sample ####
ods.random.distribution <- ggplot(quantity.rand,
                            aes(round(ods_mph,0))) +
  geom_histogram() +
  ylab("# of clips") +
  xlab ("ODS min/hr") +
  basic.theme
ggsave(paste0(shiny.input.img.path, "ODS_random_distribution.png"),
       plot = ods.random.distribution,
       width = 8, height = 6, dpi = 300)
#sd(round(quantity.rand$ods_mph,0))^2
#mean(round(quantity.rand$ods_mph,0))
# mean is much smaller than variance
ods.rand.zinb <- glmmTMB(round(ods_mph,0) ~
                           tchiyr.std +
                           stthr.tri +
                           hsz.std +
                           nsk.std +
                           tchiyr.std:stthr.tri +
                           tchiyr.std:hsz.std +
                           tchiyr.std:nsk.std +
                           (1|aclew_child_id),
                         data=quantity.rand,
                         ziformula=~nsk.std,
                         family="nbinom1")
#ods.rand.zinb.res = simulateResiduals(ods.rand.zinb)
#plot(ods.rand.zinb.res, rank = T) # (manually saved)
#summary(ods.rand.zinb)
#Conditional model:
#                              Estimate Std. Error z value  Pr(>|z|)    
#(Intercept)                    2.86750    0.15974  17.951  < 0.001 ***
#tchiyr.std                    -0.12805    0.18369  -0.697  0.4857    
#stthr.trimorning               0.36202    0.17348   2.087  0.0369 *  
#stthr.triafternoon             0.29431    0.15573   1.890  0.0588 .  
#hsz.std                        0.04324    0.09867   0.438  0.6612    
#nsk.std                        0.64982    0.08871   7.325  < 0.001 ***
#tchiyr.std:stthr.trimorning    0.09915    0.20802   0.477  0.6336    
#tchiyr.std:stthr.triafternoon  0.37800    0.17103   2.210  0.0271 *  
#tchiyr.std:hsz.std             0.31660    0.13159   2.406  0.0161 *  
#tchiyr.std:nsk.std            -0.01992    0.13093  -0.152  0.8791    
#---
#Zero-inflation model:
#            Estimate Std. Error z value Pr(>|z|)
#(Intercept)   -50.25   10421.88  -0.005    0.996
#nsk.std       -53.75   10600.83  -0.005    0.996
# save for reporting
ods.rand.zinb.COEF.morn <-
  coef(summary(ods.rand.zinb))[[1]]["stthr.trimorning",] 
ods.rand.zinb.COEF.aft <-
  coef(summary(ods.rand.zinb))[[1]]["stthr.triafternoon",] 
ods.rand.zinb.COEF.nsk <-
  coef(summary(ods.rand.zinb))[[1]]["nsk.std",] 
ods.rand.zinb.COEF.age.morn <-
  coef(summary(ods.rand.zinb))[[1]]["tchiyr.std:stthr.trimorning",]
ods.rand.zinb.COEF.age.aft <-
  coef(summary(ods.rand.zinb))[[1]]["tchiyr.std:stthr.triafternoon",]
ods.rand.zinb.COEF.age.hsz <-
  coef(summary(ods.rand.zinb))[[1]]["tchiyr.std:hsz.std",] 

# test the other two-way effects of time of day
ods.rand.zinb.v2 <- glmmTMB(round(ods_mph,0) ~
                              tchiyr.std +
                              stthr.tri.a +
                              hsz.std +
                              nsk.std +
                              tchiyr.std:stthr.tri.a +
                              tchiyr.std:hsz.std +
                              tchiyr.std:nsk.std +
                              (1|aclew_child_id),
                            data=quantity.rand,
                            ziformula=~nsk.std,
                            family="nbinom1")
# save for reporting
#summary(ods.rand.zinb.v2)
#stthr.tri.amidday             -0.29431    0.15573  -1.890               0.0588 .  
#stthr.tri.amorning             0.06770    0.14238   0.476               0.6344    
#tchiyr.std:stthr.tri.amidday  -0.37800    0.17103  -2.210               0.0271 *  
#tchiyr.std:stthr.tri.amorning -0.27885    0.17162  -1.625               0.1042    
ods.rand.zinb.v2.COEF.aft.vs.morn <-
  coef(summary(ods.rand.zinb.v2))[[1]]["stthr.tri.amorning",] 
ods.rand.zinb.v2.COEF.age.aft.vs.morn <-
  coef(summary(ods.rand.zinb.v2))[[1]]["tchiyr.std:stthr.tri.amorning",] 

## ODS tt sample ####
ods.tt.distribution <- ggplot(quantity.nonrand.tt,
                            aes(round(ods_mph,0))) +
  geom_histogram() +
  ylab("# of clips") +
  xlab ("ODS min/hr") +
  basic.theme
ggsave(paste0(shiny.input.img.path, "ODS_turntaking_distribution.png"),
       plot = ods.tt.distribution,
       width = 8, height = 6, dpi = 300)
#sd(round(quantity.nonrand.tt$ods_mph,0))^2
#mean(round(quantity.nonrand.tt$ods_mph,0))
# mean is much smaller than variance
# still zero-inflated
ods.tt.zinb <- glmmTMB(round(ods_mph,0) ~
                         tchiyr.std +
                         stthr.tri +
                         hsz.std +
                         nsk.std +
                         tchiyr.std:stthr.tri +
                         tchiyr.std:hsz.std +
                         tchiyr.std:nsk.std +
                         (1|aclew_child_id),
                       data=quantity.nonrand.tt,
                       ziformula=~nsk.std,
                       family="nbinom1")
#ods.tt.zinb.res = simulateResiduals(ods.tt.zinb)
#plot(ods.tt.zinb.res, rank = T) # (manually saved)
#summary(ods.tt.zinb)
#Conditional model:
#                              Estimate Std. Error z value  Pr(>|z|)    
#(Intercept)                    2.47226    0.21894  11.292  < 0.001 ***
#tchiyr.std                    -0.44969    0.20523  -2.191  0.02844 *  
#stthr.trimorning               0.02291    0.25858   0.089  0.92941    
#stthr.triafternoon            -0.69960    0.29240  -2.393  0.01673 *  
#hsz.std                       -0.43839    0.16870  -2.599  0.00936 ** 
#nsk.std                        0.70801    0.10686   6.626  < 0.001 ***
#tchiyr.std:stthr.trimorning   -0.56061    0.28206  -1.988  0.04686 *  
#tchiyr.std:stthr.triafternoon -0.14065    0.30202  -0.466  0.64142    
#tchiyr.std:hsz.std            -0.37513    0.21525  -1.743  0.08137 .  
#tchiyr.std:nsk.std             0.10420    0.14355   0.726  0.46792    
#---
#Zero-inflation model:
#            Estimate Std. Error z value Pr(>|z|)
#(Intercept)   -32.21   12233.79  -0.003    0.998
#nsk.std       -31.55   12037.74  -0.003    0.998
# save for reporting
ods.tt.zinb.COEF.age <-
  coef(summary(ods.tt.zinb))[[1]]["tchiyr.std",] 
ods.tt.zinb.COEF.morn <-
  coef(summary(ods.tt.zinb))[[1]]["stthr.trimorning",] 
ods.tt.zinb.COEF.aft <-
  coef(summary(ods.tt.zinb))[[1]]["stthr.triafternoon",] 
ods.tt.zinb.COEF.hsz <-
  coef(summary(ods.tt.zinb))[[1]]["hsz.std",] 
ods.tt.zinb.COEF.nsk <-
  coef(summary(ods.tt.zinb))[[1]]["nsk.std",] 
ods.tt.zinb.COEF.age.morn <-
  coef(summary(ods.tt.zinb))[[1]]["tchiyr.std:stthr.trimorning",] 
ods.tt.zinb.COEF.age.aft <-
  coef(summary(ods.tt.zinb))[[1]]["tchiyr.std:stthr.triafternoon",] 
ods.tt.zinb.COEF.age.hsz <-
  coef(summary(ods.tt.zinb))[[1]]["tchiyr.std:hsz.std",] 

# test the other two-way effects of time of day
ods.tt.zinb.v2 <- glmmTMB(round(ods_mph,0) ~
                            tchiyr.std +
                            stthr.tri.a +
                            hsz.std +
                            nsk.std +
                            tchiyr.std:stthr.tri.a +
                            tchiyr.std:hsz.std +
                            tchiyr.std:nsk.std +
                            (1|aclew_child_id),
                          data=quantity.nonrand.tt,
                          ziformula=~nsk.std,
                          family="nbinom1")
# save for reporting
#summary(ods.tt.zinb.v2)
#stthr.tri.amidday               0.6996     0.2924   2.393         0.01673 *  
#stthr.tri.amorning              0.7225     0.2485   2.907         0.00364 ** 
#tchiyr.std:stthr.tri.amidday    0.1407     0.3020   0.466         0.64142    
#tchiyr.std:stthr.tri.amorning  -0.4200     0.2720  -1.544         0.12266   
ods.tt.zinb.v2.COEF.aft.vs.morn <-
  coef(summary(ods.tt.zinb.v2))[[1]]["stthr.tri.amorning",] 
ods.tt.zinb.v2.COEF.aft.vs.midd <-
  coef(summary(ods.tt.zinb.v2))[[1]]["stthr.tri.amidday",] 
ods.tt.zinb.v2.COEF.age.aft.vs.morn <-
  coef(summary(ods.tt.zinb.v2))[[1]]["tchiyr.std:stthr.tri.amorning",] 

## ODS random sample ####
ods.rand.gaus <- glmmTMB(log(ods_mph+1) ~
                           tchiyr.std +
                           stthr.tri +
                           hsz.std +
                           nsk.std +
                           tchiyr.std:stthr.tri +
                           tchiyr.std:hsz.std +
                           tchiyr.std:nsk.std +
                           (1|aclew_child_id),
                         data=quantity.rand)
#ods.rand.gaus.res = simulateResiduals(ods.rand.gaus)
#plot(ods.rand.gaus.res, rank = T) # (manually saved)
#summary(ods.rand.gaus)
#Conditional model:
#                               Estimate Std. Error z value  Pr(>|z|)    
#(Intercept)                    2.209209   0.173268  12.750  <0.0001 ***
#tchiyr.std                    -0.082121   0.198366  -0.414  0.6789    
#stthr.trimorning               0.214457   0.210419   1.019  0.3081    
#stthr.triafternoon             0.336292   0.186481   1.803  0.0713 .  
#hsz.std                       -0.222994   0.137259  -1.625  0.1042    
#nsk.std                        1.534691   0.094447  16.249  <0.0001 ***
#tchiyr.std:stthr.trimorning   -0.006369   0.225841  -0.028  0.9775    
#tchiyr.std:stthr.triafternoon  0.422003   0.201088   2.099  0.0359 *  
#tchiyr.std:hsz.std             0.322522   0.169461   1.903  0.0570 .  
#tchiyr.std:nsk.std             0.082755   0.122301   0.677  0.4986    
# DIFFERENCES?
# broadly similar to z-i nb model


# test the other two-way effects of time of day
ods.rand.gaus.v2 <- glmmTMB(log(ods_mph+1) ~
                              tchiyr.std +
                              stthr.tri.a +
                              hsz.std +
                              nsk.std +
                              tchiyr.std:stthr.tri.a +
                              tchiyr.std:hsz.std +
                              tchiyr.std:nsk.std +
                              (1|aclew_child_id),
                            data=quantity.rand)
#summary(ods.rand.gaus.v2)
# save for reporting
#stthr.tri.amidday             -0.33629    0.18648  -1.803 0.0713 .  
#stthr.tri.amorning            -0.12183    0.18372  -0.663 0.5072    
#tchiyr.std:stthr.tri.amidday  -0.42200    0.20109  -2.099 0.0359 *  
#tchiyr.std:stthr.tri.amorning -0.42837    0.19584  -2.187 0.0287 *
# DIFFERENCES?
# broadly similar to z-i nb model

## ODS tt sample ####
ods.tt.gaus <- glmmTMB(log(ods_mph+1) ~
                         tchiyr.std +
                         stthr.tri +
                         hsz.std +
                         nsk.std +
                         tchiyr.std:stthr.tri +
                         tchiyr.std:hsz.std +
                         tchiyr.std:nsk.std +
                         (1|aclew_child_id),
                       data=quantity.nonrand.tt)
#ods.tt.gaus.res = simulateResiduals(ods.tt.gaus)
#plot(ods.tt.gaus.res, rank = T) # (manually saved)
#summary(ods.tt.gaus)
#Conditional model:
#                              Estimate Std. Error z value Pr(>|z|)    
#(Intercept)                    2.24127    0.23554   9.516 < 0.0001 ***
#tchiyr.std                    -0.25514    0.22567  -1.131 0.25822    
#stthr.trimorning              -0.18879    0.28855  -0.654 0.51294    
#stthr.triafternoon            -0.74887    0.25885  -2.893 0.00382 ** 
#hsz.std                       -0.38815    0.16704  -2.324 0.02014 *  
#nsk.std                        1.11687    0.11102  10.060 < 0.0001 ***
#tchiyr.std:stthr.trimorning   -0.28973    0.28699  -1.010 0.31271    
#tchiyr.std:stthr.triafternoon  0.04321    0.25431   0.170 0.86509    
#tchiyr.std:hsz.std            -0.21240    0.21245  -1.000 0.31742    
#tchiyr.std:nsk.std            -0.03232    0.14457  -0.224 0.82311    
# DIFFERENCES?
# Fairly similar, but some differences: no effect of child age or child age * morning

# test the other two-way effects of time of day
ods.tt.gaus.v2 <- glmmTMB(log(ods_mph+1) ~
                            tchiyr.std +
                            stthr.tri.a +
                            hsz.std +
                            nsk.std +
                            tchiyr.std:stthr.tri.a +
                            tchiyr.std:hsz.std +
                            tchiyr.std:nsk.std +
                            (1|aclew_child_id),
                          data=quantity.nonrand.tt)
#summary(ods.tt.gaus.v2)
# save for reporting
#stthr.tri.amidday              0.74887    0.25885   2.893  0.00382 ** 
#stthr.tri.amorning             0.56008    0.25688   2.180  0.02924 *  
#tchiyr.std:stthr.tri.amidday  -0.04321    0.25431  -0.170  0.86509    
#tchiyr.std:stthr.tri.amorning -0.33294    0.26822  -1.241  0.21450    
# DIFFERENCES?
# very similar to z-i nb model

# Write model results out for input to shiny
ODS.models <- bind_rows(
  broom.mixed::tidy(ods.rand.zinb) %>%
    mutate(model = "ODS_random_z-inb"),
  broom.mixed::tidy(ods.rand.zinb.v2) %>%
    mutate(model = "ODS_random_z-inb.v2"),
  broom.mixed::tidy(ods.rand.gaus) %>%
    mutate(model = "ODS_random_log_gaus"),
  broom.mixed::tidy(ods.rand.gaus.v2) %>%
    mutate(model = "ODS_random_log_gaus.v2"),
  broom.mixed::tidy(ods.tt.zinb) %>%
    mutate(model = "ODS_turntaking_z-inb"),
  broom.mixed::tidy(ods.tt.zinb.v2) %>%
    mutate(model = "ODS_turntaking_z-inb.v2"),
  broom.mixed::tidy(ods.tt.gaus) %>%
    mutate(model = "ODS_turntaking_log_gaus"),
  broom.mixed::tidy(ods.tt.gaus.v2) %>%
    mutate(model = "ODS_turntaking_log_gaus.v2"))


# All segment info
all.segments <- all.rand.segments %>%
  mutate(sample = "random", sample_type = "random") %>%
  bind_rows(all.nonrand.segments) %>%
  arrange(aclew_child_id, segment)

# Number of speakers per clip
spkrs.per.seg.all <- all.data %>%
  filter(speaker != "CHI" &
           !(grepl("@", tier))) %>%
  group_by(aclew_child_id, segment) %>%
  summarise(n_spkrs_clip = length(unique(speaker))) %>%
  full_join(all.segments, by = c("aclew_child_id", "segment")) %>%
  replace_na(list(n_spkrs_clip = 0)) %>%
  arrange(aclew_child_id, segment)

# Turn transitions:
# other-to-child and child-to-other
turn.transitions <- tibble()
for (i in 1:nrow(all.segments)) {
  subdata <- all.data %>%
    filter(aclew_child_id == all.segments$aclew_child_id[i] &
             segment == all.segments$segment[i])
  c_utt <- subdata %>%
    filter(tier == "CHI")
  used.tm1s <- c()
  used.tp1s <- c()
  child <- all.segments$aclew_child_id[i]
  # Save turn-by-turn info
  chi.turn.info <- c_utt %>%
    select(aclew_child_id, segment, tier, speaker, start, stop) %>%
    mutate(tm1.tier = NA, tm1.speaker = NA, tm1.start = NA, tm1.stop = NA, tm1.val = NA,
           tp1.tier = NA, tp1.speaker = NA, tp1.start = NA, tp1.stop = NA, tp1.val = NA)
  for (j in 1:nrow(c_utt)) {
    # Find CHI-OTH transitions
    # "T" responses that start:
    #    - earliest: when the child starts vocalizing, with a limit on vocal overlap
    #    - latest: before the maximum allowed gap after the child's voc ends
    tp1.start <- max((c_utt$stop[j] - allowed.overlap), c_utt$start[j])
    tp1.stop <- c_utt$stop[j] + allowed.gap
    t.plus1 <- which(subdata$speaker != "CHI" &
                    (subdata$val == "T") &
                    subdata$start >= tp1.start &
                    subdata$start <= tp1.stop)
    # Find OTH-CHI transitions
    # "T" prompts that start:
    #    - earliest: within the maximum gap allowed before the child begins vocalizing
    #    - latest: when the child stops vocalizing, with a limit on vocal overlap
    tm1.start <- c_utt$start[j] - allowed.gap
    tm1.stop <- min((c_utt$start[j] + allowed.overlap), c_utt$stop[j])
    t.minus1 <- which(subdata$speaker != "CHI" &
                    (subdata$val == "T") &
                    subdata$stop <= tm1.stop &
                    subdata$stop >= tm1.start)
    if(length(t.plus1) > 0) {
      tp1.match <- 0
      for (turn in t.plus1) {
        # Each OTH turn can only be a response once 
        if (!(turn %in% used.tp1s) & tp1.match == 0) {
          used.tp1s <- c(used.tp1s, turn)
          chi.turn.info$tp1.tier[j] <- subdata$tier[turn]
          chi.turn.info$tp1.speaker[j] <- subdata$speaker[turn]
          chi.turn.info$tp1.start[j] <- subdata$start[turn]
          chi.turn.info$tp1.stop[j] <- subdata$stop[turn]
          chi.turn.info$tp1.val[j] <- subdata$val[turn]
          tp1.match <- 1
        }
      }
    }
    if(length(t.minus1) > 0) {
      tm1.match <- 0
      for (turn in t.minus1) {
        # Each OTH turn can only be a prompt once 
        if (!(turn %in% used.tm1s) & tm1.match == 0) {
          used.tm1s <- c(used.tm1s, turn)
          chi.turn.info$tm1.tier[j] <- subdata$tier[turn]
          chi.turn.info$tm1.speaker[j] <- subdata$speaker[turn]
          chi.turn.info$tm1.start[j] <- subdata$start[turn]
          chi.turn.info$tm1.stop[j] <- subdata$stop[turn]
          chi.turn.info$tm1.val[j] <- subdata$val[turn]
          tm1.match <- 1
        }
      }
    }
  }
  turn.transitions <- bind_rows(turn.transitions, chi.turn.info)
}

turn.transitions.overview.o_c <- turn.transitions %>%
  filter(!(is.na(tm1.val))) %>%
  group_by(aclew_child_id, segment) %>%
  summarise(n.o_c.tts = n()) %>%
  full_join(all.segments, by = c("aclew_child_id", "segment")) %>%
  replace_na(list(n.o_c.tts = 0)) %>%
  arrange(aclew_child_id, segment) %>%
  mutate(duration = ifelse(grepl('extension', sample_type), 5,
                           ifelse(grepl('random', sample_type), 5, 1)))

turn.transitions.overview.c_o <- turn.transitions %>%
  filter(!(is.na(tp1.val))) %>%
  group_by(aclew_child_id, segment) %>%
  summarise(n.c_o.tts = n()) %>%
  full_join(all.segments, by = c("aclew_child_id", "segment")) %>%
  replace_na(list(n.c_o.tts = 0)) %>%
  arrange(aclew_child_id, segment) %>%
  mutate(duration = ifelse(grepl('extension', sample_type), 5,
                           ifelse(grepl('random', sample_type), 5, 1)))

# Combine the turn-taking info into one table
turn.transitions.overview <- turn.transitions.overview.o_c %>%
  full_join(dplyr::select(turn.transitions.overview.c_o,
                          c("aclew_child_id", "segment", "n.c_o.tts")),
            by = c("aclew_child_id", "segment")) %>%
  mutate(n.o_c.tpm = n.o_c.tts/duration,
         n.c_o.tpm = n.c_o.tts/duration) %>%
  left_join(ptcp.info, by = "aclew_child_id") %>%
  left_join(spkrs.per.seg.all, by = c("aclew_child_id", "segment",
                                      "sample", "sample_type")) %>%
  left_join(dplyr::select(seg.info, c("aclew_id", "CodeName", "start.hr")),
            by = c("aclew_child_id" = "aclew_id", "segment" = "CodeName"))

# Summaries of turn taking by sample
turn.taking.by.sample <- turn.transitions.overview %>%
  group_by(sample) %>%
  summarise(mean.n.o_c.tpm = mean(n.o_c.tpm),
         median.n.o_c.tpm = median(n.o_c.tpm),
         min.n.o_c.tpm = min(n.o_c.tpm),
         max.n.o_c.tpm = max(n.o_c.tpm),
         mean.n.c_o.tpm = mean(n.c_o.tpm),
         median.n.c_o.tpm = median(n.c_o.tpm),
         min.n.c_o.tpm = min(n.c_o.tpm),
         max.n.c_o.tpm = max(n.c_o.tpm))
turn.taking.by.child.bysample <- turn.transitions.overview %>%
  group_by(aclew_child_id, sample) %>%
  summarise(mean.n.o_c.tpm = mean(n.o_c.tpm),
         median.n.o_c.tpm = median(n.o_c.tpm),
         min.n.o_c.tpm = min(n.o_c.tpm),
         max.n.o_c.tpm = max(n.o_c.tpm),
         mean.n.c_o.tpm = mean(n.c_o.tpm),
         median.n.c_o.tpm = median(n.c_o.tpm),
         min.n.c_o.tpm = min(n.c_o.tpm),
         max.n.c_o.tpm = max(n.c_o.tpm))

# Interactional sequences:
# continuous sequences only including target-child vocalization and
# target-child-directed talk
turn.sequences <- tibble()
for (i in 1:nrow(all.segments)) {
  # retrieve all the target child vocs and xds == T vocalizations
  subdata.tcds <- all.data %>%
    filter(aclew_child_id == all.segments$aclew_child_id[i] &
             segment == all.segments$segment[i] &
             grepl('xds@', tier) &
             val == "T") %>%
    mutate(uttid = paste0(speaker, start))
  subdata.chi <- all.data %>%
    filter(aclew_child_id == all.segments$aclew_child_id[i] &
             segment == all.segments$segment[i] &
             tier == "CHI") %>%
    mutate(uttid = paste0(speaker, start))
  subdata <- bind_rows(subdata.tcds, subdata.chi) %>%
    arrange(start)
  # work through each target child voc in the clip
  if (nrow(subdata) > 0) {
    if (nrow(subdata.chi) > 0) {
      chi.vocs <- mutate(subdata.chi,
          seq.num = rep(NA,nrow(subdata.chi)),
          seq.start = rep(NA,nrow(subdata.chi)),
          seq.start.spkr = rep(NA,nrow(subdata.chi)),
          seq.stop = rep(NA,nrow(subdata.chi)),
          seq.stop.spkr = rep(NA,nrow(subdata.chi)))
      seq.num <- 1
      for (j in 1:nrow(chi.vocs)) {
        curr.start <- chi.vocs$start[j]
        curr.stop <- chi.vocs$stop[j]
        curr.spk <- "CHI"
        curr.utt <- chi.vocs$uttid[j]
        # First check if this chi voc is already in a sequence. If so, skip it
        if (seq.num > 1) {
          prev.seq <- which(chi.vocs$seq.num == seq.num - 1)[1]
          prev.seq.start <- chi.vocs$seq.start[prev.seq]
          prev.seq.stop <- chi.vocs$seq.stop[prev.seq]
          if (curr.start >= prev.seq.start & curr.stop <= prev.seq.stop) {
            chi.vocs$seq.start[j] <- chi.vocs$seq.start[prev.seq]
            chi.vocs$seq.start.spkr[j] <- chi.vocs$seq.start.spkr[prev.seq]
            chi.vocs$seq.stop[j] <- chi.vocs$seq.stop[prev.seq]
            chi.vocs$seq.stop.spkr[j] <- chi.vocs$seq.stop.spkr[prev.seq]
            chi.vocs$seq.num[j] <- chi.vocs$seq.num[prev.seq]
            next
          }
        }
        # We start a loop to look for related utterances to the LEFT
        stop.looking.left <- FALSE
        # Look first for prior turns from that STOP
        #    - earliest: within the maximum gap allowed before the curr utt begins
        #    - latest: up to the limit on vocal overlap or the end of
        #              the curr utt (whichever comes first)
        while (stop.looking.left == FALSE) {
          candidate.vocs <- subdata %>%
            filter(uttid != curr.utt &
                     start <= curr.start &
                     stop >= curr.start - allowed.gap &
                     stop <= min(curr.stop, (curr.start + allowed.overlap))) %>%
            arrange(start)
          if (nrow(candidate.vocs) > 0) {
            # continue with the earliest (first occurring) candidate
            curr.start <- candidate.vocs$start[1]
            curr.stop <- candidate.vocs$stop[1]
            curr.spk <- candidate.vocs$speaker[1]
            curr.utt <- candidate.vocs$uttid[1]
          } else {
            stop.looking.left <- TRUE
          }
        }
        seq.start <- curr.start
        seq.start.spkr <- curr.spk

        # We start a loop to look for related utterances to the RIGHT
        curr.start <- chi.vocs$start[j]
        curr.stop <- chi.vocs$stop[j]
        curr.spk <- "CHI"
        stop.looking.right <- FALSE
        # Look first for turn transitions that START
        #    - earliest: within the allowed vocal overlap over the current utt,
        #              up to its start time (whichever comes later)
        #    - latest: before the maximum allowed gap after the curr utt ends
        while(stop.looking.right == FALSE) {
          candidate.vocs <- subdata %>%
            filter(uttid != curr.utt &
                     stop >= curr.stop &
                     start >= max(curr.start, (curr.stop - allowed.overlap)) &
                     start <= curr.stop + allowed.gap) %>%
            arrange(-stop)
          if (nrow(candidate.vocs) > 0) {
            # continue with the earliest (first occurring) candidate
            curr.start <- candidate.vocs$start[1]
            curr.stop <- candidate.vocs$stop[1]
            curr.spk <- candidate.vocs$speaker[1]
            curr.utt <- candidate.vocs$uttid[1]
          } else {
            stop.looking.right <- TRUE
          }
        }
        seq.stop <- curr.stop
        seq.stop.spkr <- curr.spk

        # Check if the "sequence" is actually just 1+ CHI utterances
        if (seq.start == chi.vocs$start[j] & seq.stop == chi.vocs$stop[j]) {
          chi.vocs$seq.stop[j] <-  0
          chi.vocs$seq.stop.spkr[j] <- "NONE"
          chi.vocs$seq.start[j] <- 0
          chi.vocs$seq.start.spkr[j] <- "NONE"
        } else if (length(unique(subset(subdata, start >= seq.start &
                                        stop <= seq.stop)$speaker)) == 1) {
          chi.vocs$seq.stop[j] <-  0
          chi.vocs$seq.stop.spkr[j] <- "NONE"
          chi.vocs$seq.start[j] <- 0
          chi.vocs$seq.start.spkr[j] <- "NONE"
        } else {
          chi.vocs$seq.stop[j] <- seq.stop
          chi.vocs$seq.stop.spkr[j] <- seq.stop.spkr
          chi.vocs$seq.start[j] <- seq.start
          chi.vocs$seq.start.spkr[j] <- seq.start.spkr
          # Write the sequence number
          chi.vocs$seq.num[j] <- seq.num
          seq.num <- seq.num + 1
        }
      }
      turn.sequences <- bind_rows(turn.sequences, chi.vocs)
    }
  }
}

# Sequence duration summary
turn.sequences <- turn.sequences %>%
  filter(!(is.na(seq.num))) %>%
  mutate(seq.dur = (seq.stop - seq.start)/60000)
turn.sequences.overview <- turn.sequences %>%
  group_by(aclew_child_id, sample, segment, seq.num, seq.start, seq.stop, seq.dur) %>%
  summarise(n_cvcs_seq = n()) %>%
  left_join(ptcp.info, by = "aclew_child_id") %>%
  left_join(spkrs.per.seg.all, by = c("aclew_child_id", "segment", "sample")) %>%
  left_join(dplyr::select(seg.info, c("aclew_id", "CodeName", "start.hr")),
            by = c("aclew_child_id" = "aclew_id", "segment" = "CodeName"))
turn.seq.by.sample <- turn.sequences.overview %>%
  group_by(sample) %>%
  summarise(mean.tsq_dur = mean(seq.dur),
            median.tsq_dur = median(seq.dur),
            min.tsq_dur = min(seq.dur),
            max.tsq_dur = max(seq.dur),
            mean.cvcs = mean(n_cvcs_seq),
            median.cvcs = median(n_cvcs_seq),
            min.cvcs = min(n_cvcs_seq),
            max.cvcs = max(n_cvcs_seq))
turn.seq.by.sample.by.child <- turn.sequences.overview %>%
  group_by(aclew_child_id, sample) %>%
  summarise(mean.tsq_dur = mean(seq.dur),
            median.tsq_dur = median(seq.dur),
            min.tsq_dur = min(seq.dur),
            max.tsq_dur = max(seq.dur),
            mean.cvcs = mean(n_cvcs_seq),
            median.cvcs = median(n_cvcs_seq),
            min.cvcs = min(n_cvcs_seq),
            max.cvcs = max(n_cvcs_seq))

# Delta of chi voc and TCDS minutes
cvc.ovc.per.seg.rand <- all.data %>%
  filter(sample == "random" & speaker == "CHI" & speaker == tier) %>%
  group_by(aclew_child_id, segment, segment_dur) %>%
  summarise(cvc_min = round(sum(dur)/60000,3)) %>%
  full_join(all.rand.segments, by = c("aclew_child_id", "segment")) %>%
  replace_na(list(segment_dur = 5, cvc_min = 0)) %>%
  mutate(cvc_mph = (cvc_min/segment_dur)*60) %>%
  left_join(dplyr::select(tds.per.seg.rand, c("aclew_child_id", "segment", "tds_mph")),
            by = c("aclew_child_id", "segment")) %>%
  left_join(ptcp.info, by = "aclew_child_id") %>%
  mutate(cvc.ovc.delta = cvc_mph-tds_mph,
         cvc.ovc.delta.norm = cvc.ovc.delta/(cvc_mph+tds_mph)) %>%
  left_join(spkrs.per.seg.all, by = c("aclew_child_id", "segment")) %>%
  left_join(dplyr::select(seg.info, c("aclew_id", "CodeName", "start.hr")),
            by = c("aclew_child_id" = "aclew_id", "segment" = "CodeName")) %>%
  arrange(aclew_child_id, segment)

cvc.ovc.per.seg.nonrand <- all.data %>%
  filter(sample != "random" & speaker == "CHI" & speaker == tier) %>%
  group_by(aclew_child_id, segment, segment_dur) %>%
  summarise(cvc_min = round(sum(dur)/60000,3)) %>%
  full_join(all.nonrand.segments, by = c("aclew_child_id", "segment")) %>%
  replace_na(list(segment_dur = 5, cvc_min = 0)) %>%
  mutate(cvc_mph = (cvc_min/segment_dur)*60) %>%
  left_join(dplyr::select(tds.per.seg.nonrand, c("aclew_child_id", "segment", "tds_mph")),
            by = c("aclew_child_id", "segment")) %>%
  left_join(ptcp.info, by = "aclew_child_id") %>%
  mutate(cvc.ovc.delta = cvc_mph-tds_mph,
         cvc.ovc.delta.norm = cvc.ovc.delta/(cvc_mph+tds_mph)) %>%
  left_join(spkrs.per.seg.all,
            by = c("aclew_child_id", "segment", "sample", "sample_type")) %>%
  left_join(dplyr::select(seg.info, c("aclew_id", "CodeName", "start.hr")),
            by = c("aclew_child_id" = "aclew_id", "segment" = "CodeName")) %>%
  arrange(aclew_child_id, segment)

# Subset the samples for analysis
turn.transitions.rand_and_tt <- filter(turn.transitions.overview,
                                       sample == "random" | sample == "turn-taking")
turn.sequences.rand_and_tt <- filter(turn.sequences.overview,
                                     sample == "random" | sample == "turn-taking")
cvc.ovc.per.seg.rand_and_tt <- bind_rows(mutate(cvc.ovc.per.seg.rand, sample = "random"),
                                         filter(cvc.ovc.per.seg.nonrand,
                                                sample == "turn-taking"))
cvc.ovc.per.seg.rand_and_tt.nonas <- filter(cvc.ovc.per.seg.rand_and_tt,
                                            !(is.na(cvc.ovc.delta.norm)))

## Prepare variables for modeling
turn.transitions.rand <- filter(turn.transitions.rand_and_tt, sample == "random")
turn.transitions.rand$child_sex <- as.factor(turn.transitions.rand$child_sex)
turn.transitions.rand$mat_ed <- as.factor(turn.transitions.rand$mat_ed)
nspkrs.m <- mean(turn.transitions.rand$n_spkrs_clip)
nspkrs.sd <- sd(turn.transitions.rand$n_spkrs_clip)
turn.transitions.rand <- turn.transitions.rand %>%
  mutate(
    tchiyr.std = ((age_mo_round - tchiyr.m)/tchiyr.sd),
    chisx.std = recode_factor(child_sex,
                              "M" = "M", "F" = "F"),
    mated.std = recode_factor(mat_ed,
                              "none" = "none", "primary" = "primary",
                              "secondary" = "secondary", "preparatory" = "preparatory"),
    mated.bin = recode_factor(mat_ed,
                              "none" = "0-5", "primary" = "0-5",
                              "secondary" = "6+", "preparatory" = "6+"),
    motyr.std = ((mother_age - motyr.m)/motyr.sd),
    nsb.std = ((number_older_sibs - nsb.m)/nsb.sd),
    hsz.std = ((household_size - hsz.m)/hsz.sd),
    nsk.std = ((n_spkrs_clip - nspkrs.m)/nspkrs.sd),
    stthr.std = (start.hr - 12)/12,
    stthr.tri = ifelse(start.hr < 11, "morning",
                       ifelse(start.hr > 13, "afternoon", "midday")))
turn.transitions.rand$stthr.tri <- factor(turn.transitions.rand$stthr.tri,
                                          levels = c("midday", "morning", "afternoon"))
turn.transitions.rand$stthr.tri.a <- factor(turn.transitions.rand$stthr.tri,
                                            levels = c("afternoon", "midday", "morning"))

turn.sequences.rand <- filter(turn.sequences.rand_and_tt, sample == "random")
turn.sequences.rand$child_sex <- as.factor(turn.sequences.rand$child_sex)
turn.sequences.rand$mat_ed <- as.factor(turn.sequences.rand$mat_ed)
nspkrs.m <- mean(turn.sequences.rand$n_spkrs_clip)
nspkrs.sd <- sd(turn.sequences.rand$n_spkrs_clip)
turn.sequences.rand <- turn.sequences.rand %>%
  mutate(
    tchiyr.std = ((age_mo_round - tchiyr.m)/tchiyr.sd),
    chisx.std = recode_factor(child_sex,
                              "M" = "M", "F" = "F"),
    mated.std = recode_factor(mat_ed,
                              "none" = "none", "primary" = "primary",
                              "secondary" = "secondary", "preparatory" = "preparatory"),
    mated.bin = recode_factor(mat_ed,
                              "none" = "0-5", "primary" = "0-5",
                              "secondary" = "6+", "preparatory" = "6+"),
    motyr.std = ((mother_age - motyr.m)/motyr.sd),
    nsb.std = ((number_older_sibs - nsb.m)/nsb.sd),
    hsz.std = ((household_size - hsz.m)/hsz.sd),
    nsk.std = ((n_spkrs_clip - nspkrs.m)/nspkrs.sd),
    stthr.std = (start.hr - 12)/12,
    stthr.tri = ifelse(start.hr < 11, "morning",
                       ifelse(start.hr > 13, "afternoon", "midday")))
turn.sequences.rand$stthr.tri <- factor(turn.sequences.rand$stthr.tri,
                                        levels = c("midday", "morning", "afternoon"))
turn.sequences.rand$stthr.tri.a <- factor(turn.sequences.rand$stthr.tri,
                                          levels = c("afternoon", "midday", "morning"))

cvc.ovc.per.seg.randnonas <- filter(cvc.ovc.per.seg.rand_and_tt.nonas, sample == "random")
cvc.ovc.per.seg.randnonas$child_sex <- as.factor(cvc.ovc.per.seg.randnonas$child_sex)
cvc.ovc.per.seg.randnonas$mat_ed <- as.factor(cvc.ovc.per.seg.randnonas$mat_ed)
nspkrs.m <- mean(turn.sequences.rand$n_spkrs_clip)
nspkrs.sd <- sd(turn.sequences.rand$n_spkrs_clip)
cvc.ovc.per.seg.randnonas <- cvc.ovc.per.seg.randnonas %>%
  mutate(
    tchiyr.std = ((age_mo_round - tchiyr.m)/tchiyr.sd),
    chisx.std = recode_factor(child_sex,
                              "M" = "M", "F" = "F"),
    mated.std = recode_factor(mat_ed,
                              "none" = "none", "primary" = "primary",
                              "secondary" = "secondary", "preparatory" = "preparatory"),
    mated.bin = recode_factor(mat_ed,
                              "none" = "0-5", "primary" = "0-5",
                              "secondary" = "6+", "preparatory" = "6+"),
    motyr.std = ((mother_age - motyr.m)/motyr.sd),
    nsb.std = ((number_older_sibs - nsb.m)/nsb.sd),
    hsz.std = ((household_size - hsz.m)/hsz.sd),
    nsk.std = ((n_spkrs_clip - nspkrs.m)/nspkrs.sd),
    stthr.std = (start.hr - 12)/12)

cvc.ovc.per.seg.randnonas.bychild <- cvc.ovc.per.seg.randnonas %>%
  group_by(aclew_child_id) %>%
  summarise(mean.ratio = mean(cvc.ovc.delta.norm))

# tt clips
turn.transitions.tt <- filter(turn.transitions.rand_and_tt, sample == "turn-taking")
turn.transitions.tt$child_sex <- as.factor(turn.transitions.tt$child_sex)
turn.transitions.tt$mat_ed <- as.factor(turn.transitions.tt$mat_ed)
nspkrs.m <- mean(turn.transitions.tt$n_spkrs_clip)
nspkrs.sd <- sd(turn.transitions.tt$n_spkrs_clip)
turn.transitions.tt <- turn.transitions.tt %>%
  mutate(
    tchiyr.std = ((age_mo_round - tchiyr.m)/tchiyr.sd),
    chisx.std = recode_factor(child_sex,
                              "M" = "M", "F" = "F"),
    mated.std = recode_factor(mat_ed,
                              "none" = "none", "primary" = "primary",
                              "secondary" = "secondary", "preparatory" = "preparatory"),
    mated.bin = recode_factor(mat_ed,
                              "none" = "0-5", "primary" = "0-5",
                              "secondary" = "6+", "preparatory" = "6+"),
    motyr.std = ((mother_age - motyr.m)/motyr.sd),
    nsb.std = ((number_older_sibs - nsb.m)/nsb.sd),
    hsz.std = ((household_size - hsz.m)/hsz.sd),
    nsk.std = ((n_spkrs_clip - nspkrs.m)/nspkrs.sd),
    stthr.std = (start.hr - 12)/12,
    stthr.tri = ifelse(start.hr < 11, "morning",
                       ifelse(start.hr > 13, "afternoon", "midday")))
turn.transitions.tt$stthr.tri <- factor(turn.transitions.tt$stthr.tri,
                                        levels = c("midday", "morning", "afternoon"))
turn.transitions.tt$stthr.tri.a <- factor(turn.transitions.tt$stthr.tri,
                                          levels = c("afternoon", "midday", "morning"))

turn.sequences.tt <- filter(turn.sequences.rand_and_tt, sample == "turn-taking")
turn.sequences.tt$child_sex <- as.factor(turn.sequences.tt$child_sex)
turn.sequences.tt$mat_ed <- as.factor(turn.sequences.tt$mat_ed)
nspkrs.m <- mean(turn.sequences.tt$n_spkrs_clip)
nspkrs.sd <- sd(turn.sequences.tt$n_spkrs_clip)
turn.sequences.tt <- turn.sequences.tt %>%
  mutate(
    tchiyr.std = ((age_mo_round - tchiyr.m)/tchiyr.sd),
    chisx.std = recode_factor(child_sex,
                              "M" = "M", "F" = "F"),
    mated.std = recode_factor(mat_ed,
                              "none" = "none", "primary" = "primary",
                              "secondary" = "secondary", "preparatory" = "preparatory"),
    mated.bin = recode_factor(mat_ed,
                              "none" = "0-5", "primary" = "0-5",
                              "secondary" = "6+", "preparatory" = "6+"),
    motyr.std = ((mother_age - motyr.m)/motyr.sd),
    nsb.std = ((number_older_sibs - nsb.m)/nsb.sd),
    hsz.std = ((household_size - hsz.m)/hsz.sd),
    nsk.std = ((n_spkrs_clip - nspkrs.m)/nspkrs.sd),
    stthr.std = (start.hr - 12)/12,
    stthr.tri = ifelse(start.hr < 11, "morning",
                       ifelse(start.hr > 13, "afternoon", "midday")))
turn.sequences.tt$stthr.tri <- factor(turn.sequences.tt$stthr.tri,
                                      levels = c("midday", "morning", "afternoon"))
turn.sequences.tt$stthr.tri.a <- factor(turn.sequences.tt$stthr.tri,
                                        levels = c("afternoon", "midday", "morning"))

cvc.ovc.per.seg.ttnonas <- filter(cvc.ovc.per.seg.rand_and_tt.nonas,
                                  sample == "turn-taking")
cvc.ovc.per.seg.ttnonas$child_sex <- as.factor(cvc.ovc.per.seg.ttnonas$child_sex)
cvc.ovc.per.seg.ttnonas$mat_ed <- as.factor(cvc.ovc.per.seg.ttnonas$mat_ed)
nspkrs.m <- mean(cvc.ovc.per.seg.ttnonas$n_spkrs_clip)
nspkrs.sd <- sd(cvc.ovc.per.seg.ttnonas$n_spkrs_clip)
cvc.ovc.per.seg.ttnonas <- cvc.ovc.per.seg.ttnonas %>%
  mutate(
    tchiyr.std = ((age_mo_round - tchiyr.m)/tchiyr.sd),
    chisx.std = recode_factor(child_sex,
                              "M" = "M", "F" = "F"),
    mated.std = recode_factor(mat_ed,
                              "none" = "none", "primary" = "primary",
                              "secondary" = "secondary", "preparatory" = "preparatory"),
    mated.bin = recode_factor(mat_ed,
                              "none" = "0-5", "primary" = "0-5",
                              "secondary" = "6+", "preparatory" = "6+"),
    motyr.std = ((mother_age - motyr.m)/motyr.sd),
    nsb.std = ((number_older_sibs - nsb.m)/nsb.sd),
    hsz.std = ((household_size - hsz.m)/hsz.sd),
    nsk.std = ((n_spkrs_clip - nspkrs.m)/nspkrs.sd),
    stthr.std = (start.hr - 12)/12)

cvc.ovc.per.seg.ttnonas.bychild <- cvc.ovc.per.seg.ttnonas %>%
  group_by(aclew_child_id) %>%
  summarise(mean.ratio = mean(cvc.ovc.delta.norm))

# Graph the basic turn taking rate info
# CHI-OTH transitions per minute
chi.oth.tts.rand_and_tt <- ggplot(turn.transitions.rand_and_tt,
                          aes(x = age_mo_round, y = n.c_o.tpm, lty = sample)) +
  geom_boxplot(aes(group = interaction(age_mo_round, sample),
                   color = sample), fill = "white", outlier.shape = NA,
               lty = "solid", alpha = 0.4) +
  geom_smooth(aes(fill = sample, color = sample), method = "lm") +
  ylab("CHI-OTH tts/min") + xlab("")	+
  scale_y_continuous(limits=c(0,30),
                     breaks=seq(0,30,5)) +
  scale_x_continuous(limits=c(0,38),
                     breaks=seq(0,38,6)) +
  coord_cartesian(ylim=c(0,30),xlim=c(0,38)) +
  scale_color_manual(values = viridis(3)) +
  scale_fill_manual(values = viridis(3)) +
  theme_apa() +
  theme(legend.position="none",
        axis.line = element_line(color="black", size = 0.4))

# OTH-CHI transitions per minute
oth.chi.tts.rand_and_tt <- ggplot(turn.transitions.rand_and_tt,
                          aes(x = age_mo_round, y = n.o_c.tpm, lty = sample)) +
  geom_boxplot(aes(group = interaction(age_mo_round, sample),
                   color = sample), fill = "white", outlier.shape = NA,
               lty = "solid", alpha = 0.4) +
  geom_smooth(aes(fill = sample, color = sample), method = "lm") +
  ylab("OTH-CHI tts/min") + xlab("Child age (mo)")	+
  scale_y_continuous(limits=c(0,30),
                     breaks=seq(0,30,5)) +
  scale_x_continuous(limits=c(0,38),
                     breaks=seq(0,38,6)) +
  coord_cartesian(ylim=c(0,30),xlim=c(0,38)) +
  scale_color_manual(values = viridis(3)) +
  scale_fill_manual(values = viridis(3)) +
  theme_apa() +
  theme(legend.position="none",
        axis.line = element_line(color="black", size = 0.4))

# Graph the basic sequence duration info
# plot per-clip averages so it's consistent with the rest
turn.sequences.rand_and_tt.byclip <- turn.sequences.rand_and_tt %>%
  group_by(aclew_child_id, age_mo_round, sample, segment) %>%
  summarise(m.seqdur.sec = mean(seq.dur*60))
seq.dur.rand_and_tt <- ggplot(turn.sequences.rand_and_tt.byclip,
                          aes(x = age_mo_round, y = m.seqdur.sec, lty = sample)) +
  geom_boxplot(aes(group = interaction(age_mo_round, sample),
                   color = sample), fill = "white", outlier.shape = NA,
               lty = "solid", alpha = 0.4) +
  geom_smooth(aes(fill = sample, color = sample), method = "lm") +
  ylab("Seq. dur. (sec)") + xlab("")	+
  scale_y_continuous(limits=c(0,60),
                     breaks=seq(0,60,20)) +
  scale_x_continuous(limits=c(0,38),
                     breaks=seq(0,38,6)) +
  coord_cartesian(ylim=c(0,60),xlim=c(0,38)) +
  scale_color_manual(values = viridis(3)) +
  scale_fill_manual(values = viridis(3)) +
  theme_apa() +
  theme(legend.position="none",
        axis.line = element_line(color="black", size = 0.4))

grid.arrange(chi.oth.tts.rand_and_tt, oth.chi.tts.rand_and_tt,
             seq.dur.rand_and_tt, nrow=1, ncol=3)

## CHI-OTH tts/min random sample ####
c_o.tpm.random.distribution <- ggplot(turn.transitions.rand,
                            aes(round(n.c_o.tpm,0))) +
  geom_histogram() +
  ylab("# of clips") +
  xlab ("TC-O transitions/min") +
  basic.theme
ggsave(paste0(shiny.input.img.path, "c_o_tpm_random_distribution.png"),
       plot = c_o.tpm.random.distribution,
       width = 8, height = 6, dpi = 300)
#sd(round(turn.transitions.rand$n.c_o.tpm,0))^2
#mean(round(turn.transitions.rand$n.c_o.tpm,0))
# mean isn't much smaller than variance
# zero inflated
c_o.tpm.rand.zinb <- glmmTMB(round(n.c_o.tpm,0) ~
                               tchiyr.std +
                               stthr.tri +
                               hsz.std +
                               nsk.std +
                               tchiyr.std:stthr.tri +
                               tchiyr.std:hsz.std +
                               tchiyr.std:nsk.std +
                               (1|aclew_child_id),
                             data=turn.transitions.rand,
                             ziformula=~nsk.std,
                             family="nbinom1")
#c_o.tpm.rand.zinb.res = simulateResiduals(c_o.tpm.rand.zinb)
#plot(c_o.tpm.rand.zinb.res, rank = T) # (manually saved)
#summary(c_o.tpm.rand.zinb)
#Conditional model:
#                              Estimate Std. Error z value Pr(>|z|)  
#(Intercept)                    -0.1257     0.5037  -0.250   0.8029  
#tchiyr.std                      0.8867     0.6063   1.462   0.1436  
#stthr.trimorning                0.4825     0.4500   1.072   0.2836  
#stthr.triafternoon              0.3405     0.3993   0.853   0.3938  
#hsz.std                        -0.1728     0.4541  -0.381   0.7035  
#nsk.std                        -0.1769     0.1746  -1.014   0.3107  
#tchiyr.std:stthr.trimorning    -0.1398     0.4751  -0.294   0.7686  
#tchiyr.std:stthr.triafternoon  -1.0820     0.4443  -2.436   0.0149 *
#tchiyr.std:hsz.std              0.1109     0.5569   0.199   0.8421  
#tchiyr.std:nsk.std              0.5639     0.2270   2.484   0.0130 *
#---
#Zero-inflation model:
#            Estimate Std. Error z value Pr(>|z|)
#(Intercept)   -116.7    53056.2  -0.002    0.998
#nsk.std       -100.0    52343.8  -0.002    0.998
# save for reporting
c_o.tpm.rand.zinb.COEF.age.morn <-
  coef(summary(c_o.tpm.rand.zinb))[[1]]["tchiyr.std:stthr.trimorning",] 
c_o.tpm.rand.zinb.COEF.age.aft <-
  coef(summary(c_o.tpm.rand.zinb))[[1]]["tchiyr.std:stthr.triafternoon",] 
c_o.tpm.rand.zinb.COEF.age.nsk <-
  coef(summary(c_o.tpm.rand.zinb))[[1]]["tchiyr.std:nsk.std",]

# test the other two-way effects of time of day
c_o.tpm.rand.zinb.v2 <- glmmTMB(round(n.c_o.tpm,0) ~
                                  tchiyr.std +
                                  stthr.tri.a +
                                  hsz.std +
                                  nsk.std +
                                  tchiyr.std:stthr.tri.a +
                                  tchiyr.std:hsz.std +
                                  tchiyr.std:nsk.std +
                                  (1|aclew_child_id),
                                data=turn.transitions.rand,
                                ziformula=~nsk.std,
                                family="nbinom1")
#summary(c_o.tpm.rand.zinb.v2)
#stthr.tri.amidday              -0.3405     0.3993  -0.853   0.3938  
#stthr.tri.amorning              0.1420     0.3254   0.436   0.6626  
#tchiyr.std:stthr.tri.amidday    1.0820     0.4443   2.436   0.0149 *
#tchiyr.std:stthr.tri.amorning   0.9422     0.3746   2.515   0.0119 *
# save for reporting
c_o.tpm.rand.zinb.v2.COEF.aft.vs.morn <-
  coef(summary(c_o.tpm.rand.zinb.v2))[[1]]["stthr.tri.amorning",] 
c_o.tpm.rand.zinb.v2.COEF.age.aft.vs.morn <-
  coef(summary(c_o.tpm.rand.zinb.v2))[[1]]["tchiyr.std:stthr.tri.amorning",] 
c_o.tpm.rand.zinb.v2.COEF.age.aft.vs.midd <-
  coef(summary(c_o.tpm.rand.zinb.v2))[[1]]["tchiyr.std:stthr.tri.amidday",] 

## CHI-OTH tts/min tt sample ####
c_o.tpm.tt.distribution <- ggplot(turn.transitions.tt,
                            aes(round(n.c_o.tpm,0))) +
  geom_histogram() +
  ylab("# of clips") +
  xlab ("TC-O transitions/min") +
  basic.theme
ggsave(paste0(shiny.input.img.path, "c_o_tpm_turntaking_distribution.png"),
       plot = c_o.tpm.tt.distribution,
       width = 8, height = 6, dpi = 300)
#ggplot(turn.transitions.tt, aes(round(n.c_o.tpm,0))) + geom_histogram()
#sd(round(turn.transitions.tt$n.c_o.tpm,0))^2
#mean(round(turn.transitions.tt$n.c_o.tpm,0))
# mean isn't much smaller than variance
# not zero-inflated (nature of sample)
c_o.tpm.tt.nb <- glmmTMB(round(n.c_o.tpm,0) ~
                           tchiyr.std +
                           stthr.tri +
                           hsz.std +
                           nsk.std +
                           tchiyr.std:stthr.tri +
                           tchiyr.std:hsz.std +
                           tchiyr.std:nsk.std +
                           (1|aclew_child_id),
                         data=turn.transitions.tt,
                         family="nbinom1")
#c_o.tpm.tt.nb.res = simulateResiduals(c_o.tpm.tt.nb)
#plot(c_o.tpm.tt.nb.res, rank = T) # (manually saved)
#summary(c_o.tpm.tt.nb)
#Conditional model:
#                              Estimate Std. Error z value Pr(>|z|)    
#(Intercept)                    1.78448    0.28828   6.190 <0.0001 ***
#tchiyr.std                     0.11866    0.29273   0.405 0.6852    
#stthr.trimorning               0.05827    0.31932   0.182 0.8552    
#stthr.triafternoon             0.11001    0.26991   0.408 0.6836    
#hsz.std                       -0.03938    0.22998  -0.171 0.8640    
#nsk.std                       -0.21384    0.11634  -1.838 0.0661 .  
#tchiyr.std:stthr.trimorning   -0.03954    0.29729  -0.133 0.8942    
#tchiyr.std:stthr.triafternoon -0.18839    0.25152  -0.749 0.4539    
#tchiyr.std:hsz.std            -0.19196    0.28233  -0.680 0.4966    
#tchiyr.std:nsk.std             0.05446    0.13186   0.413 0.6796    
# save for reporting
c_o.tpm.tt.nb.COEF.nsk <-
  coef(summary(c_o.tpm.tt.nb))[[1]]["nsk.std",] 

# test the other two-way effects of time of day
c_o.tpm.tt.nb.v2 <- glmmTMB(round(n.c_o.tpm,0) ~
                              tchiyr.std +
                              stthr.tri.a +
                              hsz.std +
                              nsk.std +
                              tchiyr.std:stthr.tri.a +
                              tchiyr.std:hsz.std +
                              tchiyr.std:nsk.std +
                              (1|aclew_child_id),
                            data=turn.transitions.tt,
                            family="nbinom1")
#summary(c_o.tpm.tt.nb.v2)
#stthr.tri.amidday             -0.11002    0.26991  -0.408              0.6835    
#stthr.tri.amorning            -0.05176    0.24239  -0.214              0.8309    
#tchiyr.std:stthr.tri.amidday   0.18840    0.25152   0.749              0.4538    
#tchiyr.std:stthr.tri.amorning  0.14885    0.26612   0.559              0.5759    
# no results to save for reporting

## CHI-OTH tts/min random sample ####
c_o.tpm.rand.gaus <- glmmTMB(log(n.c_o.tpm+2) ~
                               tchiyr.std +
                               stthr.tri +
                               hsz.std +
                               nsk.std +
                               tchiyr.std:stthr.tri +
                               tchiyr.std:hsz.std +
                               tchiyr.std:nsk.std +
                               (1|aclew_child_id),
                             data=turn.transitions.rand)
#c_o.tpm.rand.gaus.res = simulateResiduals(c_o.tpm.rand.gaus)
#plot(c_o.tpm.rand.gaus.res, rank = T) # (manually saved)
#summary(c_o.tpm.rand.gaus)
#Conditional model:
#                               Estimate Std. Error z value           Pr(>|z|)    
#(Intercept)                    1.056534   0.135493   7.798 0.0000000000000063 ***
#tchiyr.std                     0.250150   0.157480   1.588            0.11218    
#stthr.trimorning               0.137954   0.116248   1.187            0.23534    
#stthr.triafternoon             0.013519   0.102405   0.132            0.89498    
#hsz.std                       -0.110842   0.135817  -0.816            0.41444    
#nsk.std                        0.087300   0.056058   1.557            0.11940    
#tchiyr.std:stthr.trimorning   -0.025362   0.125456  -0.202            0.83979    
#tchiyr.std:stthr.triafternoon -0.303130   0.111001  -2.731            0.00632 ** 
#tchiyr.std:hsz.std             0.002589   0.164390   0.016            0.98744    
#tchiyr.std:nsk.std             0.094462   0.070742   1.335            0.18178
# DIFFERENCES?
# broadly similar to z-i nb model; no age * nsk effect

# test the other two-way effects of time of day
c_o.tpm.rand.gaus.v2 <- glmmTMB(log(n.c_o.tpm+2) ~
                                  tchiyr.std +
                                  stthr.tri.a +
                                  hsz.std +
                                  nsk.std +
                                  tchiyr.std:stthr.tri.a +
                                  tchiyr.std:hsz.std +
                                  tchiyr.std:nsk.std +
                                  (1|aclew_child_id),
                                data=turn.transitions.rand)
#summary(c_o.tpm.rand.gaus.v2)
#stthr.tri.amidday             -0.013520   0.102406  -0.132 0.89496    
#stthr.tri.amorning             0.124433   0.100790   1.235 0.21699    
#tchiyr.std:stthr.tri.amidday   0.303131   0.111001   2.731 0.00632 ** 
#tchiyr.std:stthr.tri.amorning  0.277769   0.107167   2.592 0.00954 ** 
# DIFFERENCES
# very similar to z-i nb model


## CHI-OTH tts/min tt sample ####
c_o.tpm.tt.gaus <- glmmTMB(log(n.c_o.tpm+2) ~
                           tchiyr.std +
                           stthr.tri.a +
                           hsz.std +
                           nsk.std +
                           tchiyr.std:stthr.tri.a +
                           tchiyr.std:hsz.std +
                           tchiyr.std:nsk.std +
                           (1|aclew_child_id),
                         data=turn.transitions.tt)
#c_o.tpm.tt.gaus.res = simulateResiduals(c_o.tpm.tt.gaus)
#plot(c_o.tpm.tt.gaus.res, rank = T) # (manually saved)
#summary(c_o.tpm.tt.gaus)
#Conditional model:
#                              Estimate Std. Error z value Pr(>|z|)    
#(Intercept)                    2.06738    0.20284  10.192 <0.0001 ***
#tchiyr.std                    -0.01995    0.22908  -0.087 0.9306    
#stthr.tri.amidday             -0.10853    0.20987  -0.517 0.6050    
#stthr.tri.amorning            -0.08067    0.20332  -0.397 0.6915    
#hsz.std                       -0.02492    0.20253  -0.123 0.9021    
#nsk.std                       -0.17907    0.09640  -1.857 0.0632 .  
#tchiyr.std:stthr.tri.amidday   0.11148    0.20262   0.550 0.5822    
#tchiyr.std:stthr.tri.amorning  0.09163    0.21556   0.425 0.6708    
#tchiyr.std:hsz.std            -0.16794    0.24828  -0.676 0.4988    
#tchiyr.std:nsk.std             0.09030    0.11495   0.786 0.4321  
# DIFFERENCES
# very similar to nb model

# test the other two-way effects of time of day
c_o.tpm.tt.gaus.v2 <- glmmTMB(log(n.c_o.tpm+2) ~
                              tchiyr.std +
                              stthr.tri.a +
                              hsz.std +
                              nsk.std +
                              tchiyr.std:stthr.tri.a +
                              tchiyr.std:hsz.std +
                              tchiyr.std:nsk.std +
                              (1|aclew_child_id),
                            data=turn.transitions.tt)
#summary(c_o.tpm.tt.gaus.v2)
#stthr.tri.amidday             -0.10853    0.20987  -0.517  0.6050    
#stthr.tri.amorning            -0.08067    0.20332  -0.397  0.6915    
#tchiyr.std:stthr.tri.amidday   0.11148    0.20262   0.550  0.5822    
#tchiyr.std:stthr.tri.amorning  0.09163    0.21556   0.425  0.6708    
# DIFFERENCES?
# very similar to nb model

# Write model results out for input to shiny
c_o.tpm.models <- bind_rows(
  broom.mixed::tidy(c_o.tpm.rand.zinb) %>%
    mutate(model = "c_o.tpm_random_z-inb"),
  broom.mixed::tidy(c_o.tpm.rand.zinb.v2) %>%
    mutate(model = "c_o.tpm_random_z-inb.v2"),
  broom.mixed::tidy(c_o.tpm.rand.gaus) %>%
    mutate(model = "c_o.tpm_random_log_gaus"),
  broom.mixed::tidy(c_o.tpm.rand.gaus.v2) %>%
    mutate(model = "c_o.tpm_random_log_gaus.v2"),
  broom.mixed::tidy(c_o.tpm.tt.nb) %>%
    mutate(model = "c_o.tpm_turntaking_nb"),
  broom.mixed::tidy(c_o.tpm.tt.nb.v2) %>%
    mutate(model = "c_o.tpm_turntaking_nb.v2"),
  broom.mixed::tidy(c_o.tpm.tt.gaus) %>%
    mutate(model = "c_o.tpm_turntaking_log_gaus"),
  broom.mixed::tidy(c_o.tpm.tt.gaus.v2) %>%
    mutate(model = "c_o.tpm_turntaking_log_gaus.v2"))

## OTH-CHI tts/min random sample ####
o_c.tpm.random.distribution <- ggplot(turn.transitions.rand,
                            aes(round(n.o_c.tpm,0))) +
  geom_histogram() +
  ylab("# of clips") +
  xlab ("O-TC transitions/min") +
  basic.theme
ggsave(paste0(shiny.input.img.path, "o_c_tpm_random_distribution.png"),
       plot = o_c.tpm.random.distribution,
       width = 8, height = 6, dpi = 300)
#sd(round(turn.transitions.rand$n.o_c.tpm,0))^2
#mean(round(turn.transitions.rand$n.o_c.tpm,0))
# mean isn't much smaller than variance
# zero-inflated
o_c.tpm.rand.zinb <- glmmTMB(round(n.o_c.tpm,0) ~
                               tchiyr.std +
                               stthr.tri +
                               hsz.std +
                               nsk.std +
                               tchiyr.std:stthr.tri +
                               tchiyr.std:hsz.std +
                               tchiyr.std:nsk.std +
                               (1|aclew_child_id),
                             data=turn.transitions.rand,
                             ziformula=~nsk.std,
                             family="nbinom1")
#o_c.tpm.rand.zinb.res = simulateResiduals(o_c.tpm.rand.zinb)
#plot(o_c.tpm.rand.zinb.res, rank = T) # (saved manually)
#summary(o_c.tpm.rand.zinb)
#Conditional model:
#                              Estimate Std. Error z value Pr(>|z|)   
#(Intercept)                    -0.4586     0.5436  -0.844  0.39890   
#tchiyr.std                      1.1417     0.6570   1.738  0.08226 . 
#stthr.trimorning                0.3179     0.4878   0.652  0.51462   
#stthr.triafternoon              0.4948     0.4082   1.212  0.22542   
#hsz.std                        -0.2025     0.4947  -0.409  0.68226   
#nsk.std                        -0.1414     0.1784  -0.793  0.42796   
#tchiyr.std:stthr.trimorning    -0.1238     0.5059  -0.245  0.80671   
#tchiyr.std:stthr.triafternoon  -1.4569     0.4647  -3.135  0.00172 **
#tchiyr.std:hsz.std              0.1412     0.6070   0.233  0.81610   
#tchiyr.std:nsk.std              0.5175     0.2247   2.303  0.02130 * 
#---
#Zero-inflation model:
#            Estimate Std. Error z value Pr(>|z|)
#(Intercept)  -115.46   43943.60  -0.003    0.998
#nsk.std       -98.63   42142.48  -0.002    0.998
# save for reporting
o_c.tpm.rand.zinb.COEF.age <-
  coef(summary(o_c.tpm.rand.zinb))[[1]]["tchiyr.std",]
o_c.tpm.rand.zinb.COEF.age.morn <-
  coef(summary(o_c.tpm.rand.zinb))[[1]]["tchiyr.std:stthr.trimorning",]
o_c.tpm.rand.zinb.COEF.age.aft <-
  coef(summary(o_c.tpm.rand.zinb))[[1]]["tchiyr.std:stthr.triafternoon",]
o_c.tpm.rand.zinb.COEF.age.nsk <-
  coef(summary(o_c.tpm.rand.zinb))[[1]]["tchiyr.std:nsk.std",]

# test the other two-way effects of time of day
o_c.tpm.rand.zinb.v2 <- glmmTMB(round(n.o_c.tpm,0) ~
                                  tchiyr.std +
                                  stthr.tri.a +
                                  hsz.std +
                                  nsk.std +
                                  tchiyr.std:stthr.tri.a +
                                  tchiyr.std:hsz.std +
                                  tchiyr.std:nsk.std +
                                  (1|aclew_child_id),
                                data=turn.transitions.rand,
                                ziformula=~nsk.std,
                                family="nbinom1")
#summary(o_c.tpm.rand.zinb.v2)
#stthr.tri.amidday             -0.49479    0.40816  -1.212  0.22542   
#stthr.tri.amorning            -0.17693    0.36072  -0.490  0.62379   
#tchiyr.std:stthr.tri.amidday   1.45689    0.46472   3.135  0.00172 **
#tchiyr.std:stthr.tri.amorning  1.33312    0.41759   3.192  0.00141 **
# save for reporting
o_c.tpm.rand.zinb.v2.COEF.morn <-
  coef(summary(o_c.tpm.rand.zinb.v2))[[1]]["stthr.tri.amorning",]
o_c.tpm.rand.zinb.v2.COEF.age.morn <-
  coef(summary(o_c.tpm.rand.zinb.v2))[[1]]["tchiyr.std:stthr.tri.amorning",]
o_c.tpm.rand.zinb.v2.COEF.age.midd <-
  coef(summary(o_c.tpm.rand.zinb.v2))[[1]]["tchiyr.std:stthr.tri.amidday",]

## OTH-CHI tts/min tt sample ####
o_c.tpm.tt.distribution <- ggplot(turn.transitions.tt,
                            aes(round(n.o_c.tpm,0))) +
  geom_histogram() +
  ylab("# of clips") +
  xlab ("O-TC transitions/min") +
  basic.theme
ggsave(paste0(shiny.input.img.path, "o_c_tpm_turntaking_distribution.png"),
       plot = o_c.tpm.tt.distribution,
       width = 8, height = 6, dpi = 300)
#sd(round(turn.transitions.tt$n.o_c.tpm,0))^2
#mean(round(turn.transitions.tt$n.o_c.tpm,0))
# mean is much smaller than variance
# not really zero-inflated
o_c.tpm.tt.nb <- glmmTMB(round(n.o_c.tpm,0) ~
                           tchiyr.std +
                           stthr.tri +
                           hsz.std +
                           nsk.std +
                           tchiyr.std:stthr.tri +
                           tchiyr.std:hsz.std +
                           tchiyr.std:nsk.std +
                           (1|aclew_child_id),
                         data=turn.transitions.tt,
                         family="nbinom1")
#o_c.tpm.tt.nb.res = simulateResiduals(o_c.tpm.tt.nb)
#plot(o_c.tpm.tt.nb.res, rank = T) # (manually saved)
#summary(o_c.tpm.tt.nb)
#Conditional model:
#                              Estimate Std. Error z value Pr(>|z|)    
#(Intercept)                    1.84121    0.29135   6.319 < 0.0001 ***
#tchiyr.std                     0.05962    0.29556   0.202 0.8401    
#stthr.trimorning              -0.04921    0.33317  -0.148 0.8826    
#stthr.triafternoon             0.02709    0.28235   0.096 0.9236    
#hsz.std                       -0.03631    0.23216  -0.156 0.8757    
#nsk.std                       -0.22871    0.12170  -1.879 0.0602 .  
#tchiyr.std:stthr.trimorning    0.12037    0.30574   0.394 0.6938    
#tchiyr.std:stthr.triafternoon -0.07718    0.26007  -0.297 0.7667    
#tchiyr.std:hsz.std            -0.18813    0.28410  -0.662 0.5078    
#tchiyr.std:nsk.std             0.08522    0.13646   0.624 0.5323
# save for reporting
o_c.tpm.tt.nb.COEF.nsk <-
  coef(summary(o_c.tpm.tt.nb))[[1]]["nsk.std",]

# test the other two-way effects of time of day
o_c.tpm.tt.nb.v2 <- glmmTMB(round(n.o_c.tpm,0) ~
                              tchiyr.std +
                              stthr.tri.a +
                              hsz.std +
                              nsk.std +
                              tchiyr.std:stthr.tri.a +
                              tchiyr.std:hsz.std +
                              tchiyr.std:nsk.std +
                              (1|aclew_child_id),
                            data=turn.transitions.tt,
                            family="nbinom1")
#summary(o_c.tpm.tt.nb.v2)
#stthr.tri.amidday             -0.02712    0.28236  -0.096             0.9235    
#stthr.tri.amorning            -0.07632    0.25647  -0.298             0.7660    
#tchiyr.std:stthr.tri.amidday   0.07718    0.26007   0.297             0.7666    
#tchiyr.std:stthr.tri.amorning  0.19755    0.27015   0.731             0.4646    
# no results to save for reporting

## OTH-CHI tts/min random sample ####
o_c.tpm.rand.gaus <- glmmTMB(log(n.o_c.tpm+2) ~
                               tchiyr.std +
                               stthr.tri +
                               hsz.std +
                               nsk.std +
                               tchiyr.std:stthr.tri +
                               tchiyr.std:hsz.std +
                               tchiyr.std:nsk.std +
                               (1|aclew_child_id),
                             data=turn.transitions.rand)
#o_c.tpm.rand.gaus.res = simulateResiduals(o_c.tpm.rand.gaus)
#plot(o_c.tpm.rand.gaus.res, rank = T) # (saved manually)
#summary(o_c.tpm.rand.gaus)
#Conditional model:
#                                 Estimate  Std. Error z value Pr(>|z|)    
#(Intercept)                    1.00641393  0.12448510   8.085 <0.0001 ***
#tchiyr.std                     0.22740872  0.14462973   1.572 0.11587    
#stthr.trimorning               0.10010749  0.10843729   0.923 0.35591    
#stthr.triafternoon             0.00009137  0.09549393   0.001 0.99924    
#hsz.std                       -0.11953003  0.12408216  -0.963 0.33539    
#nsk.std                        0.06939389  0.05214713   1.331 0.18328    
#tchiyr.std:stthr.trimorning   -0.01939578  0.11702461  -0.166 0.86836    
#tchiyr.std:stthr.triafternoon -0.28596439  0.10350542  -2.763 0.00573 ** 
#tchiyr.std:hsz.std            -0.02916843  0.15025795  -0.194 0.84608    
#tchiyr.std:nsk.std             0.07454357  0.06594419   1.130 0.25831 
# DIFFERENCES?
# broadly similar to z-i nb model, but with no age * nsk interaction

# test the other two-way effects of time of day
o_c.tpm.rand.gaus.v2 <- glmmTMB(log(n.o_c.tpm+2) ~
                                  tchiyr.std +
                                  stthr.tri.a +
                                  hsz.std +
                                  nsk.std +
                                  tchiyr.std:stthr.tri.a +
                                  tchiyr.std:hsz.std +
                                  tchiyr.std:nsk.std +
                                  (1|aclew_child_id),
                                data=turn.transitions.rand)
#summary(o_c.tpm.rand.gaus.v2)
#stthr.tri.amidday             -0.00009122  0.09549413  -0.001 0.99924    
#stthr.tri.amorning             0.10001752  0.09406177   1.063 0.28764    
#tchiyr.std:stthr.tri.amidday   0.28596268  0.10350565   2.763 0.00573 ** 
#tchiyr.std:stthr.tri.amorning  0.26656607  0.10003252   2.665 0.00770 ** 
# DIFFERENCES?
# very similar to z-i nb model

## OTH-CHI tts/min tt sample ####
o_c.tpm.tt.gaus <- glmmTMB(log(n.o_c.tpm+2) ~
                           tchiyr.std +
                           stthr.tri +
                           hsz.std +
                           nsk.std +
                           tchiyr.std:stthr.tri +
                           tchiyr.std:hsz.std +
                           tchiyr.std:nsk.std +
                           (1|aclew_child_id),
                         data=turn.transitions.tt)
#o_c.tpm.tt.gaus.res = simulateResiduals(o_c.tpm.tt.gaus)
#plot(o_c.tpm.tt.gaus.res, rank = T) # (manually saved)
#summary(o_c.tpm.tt.gaus)
#Conditional model:
#                               Estimate Std. Error z value Pr(>|z|)    
#(Intercept)                    2.007746   0.247768   8.103 <0.0001 ***
#tchiyr.std                     0.039817   0.246349   0.162 0.8716    
#stthr.trimorning              -0.091549   0.270425  -0.339 0.7350    
#stthr.triafternoon             0.035551   0.223546   0.159 0.8736    
#hsz.std                       -0.002974   0.208196  -0.014 0.9886    
#nsk.std                       -0.202931   0.102724  -1.975 0.0482 *  
#tchiyr.std:stthr.trimorning    0.166181   0.253668   0.655 0.5124    
#tchiyr.std:stthr.triafternoon -0.003735   0.213868  -0.017 0.9861    
#tchiyr.std:hsz.std            -0.131678   0.255391  -0.516 0.6061    
#tchiyr.std:nsk.std             0.083365   0.121167   0.688 0.4914
# DIFFERENCES?
# very similar to nb model

# test the other two-way effects of time of day
o_c.tpm.tt.gaus.v2 <- glmmTMB(log(n.o_c.tpm+2) ~
                              tchiyr.std +
                              stthr.tri.a +
                              hsz.std +
                              nsk.std +
                              tchiyr.std:stthr.tri.a +
                              tchiyr.std:hsz.std +
                              tchiyr.std:nsk.std +
                              (1|aclew_child_id),
                            data=turn.transitions.tt)
#summary(o_c.tpm.tt.gaus.v2)
#stthr.tri.amidday             -0.035550   0.223546  -0.159              0.8736    
#stthr.tri.amorning            -0.127097   0.215122  -0.591              0.5546    
#tchiyr.std:stthr.tri.amidday   0.003734   0.213868   0.017              0.9861    
#tchiyr.std:stthr.tri.amorning  0.169912   0.227402   0.747              0.4550    
# DIFFERENCES?
# very similar to nb model

# Write model results out for input to shiny
o_c.tpm.models <- bind_rows(
  broom.mixed::tidy(o_c.tpm.rand.zinb) %>%
    mutate(model = "o_c.tpm_random_z-inb"),
  broom.mixed::tidy(o_c.tpm.rand.zinb.v2) %>%
    mutate(model = "o_c.tpm_random_z-inb.v2"),
  broom.mixed::tidy(o_c.tpm.rand.gaus) %>%
    mutate(model = "o_c.tpm_random_log_gaus"),
  broom.mixed::tidy(o_c.tpm.rand.gaus.v2) %>%
    mutate(model = "o_c.tpm_random_log_gaus.v2"),
  broom.mixed::tidy(o_c.tpm.tt.nb) %>%
    mutate(model = "o_c.tpm_turntaking_nb"),
  broom.mixed::tidy(o_c.tpm.tt.nb.v2) %>%
    mutate(model = "o_c.tpm_turntaking_nb.v2"),
  broom.mixed::tidy(o_c.tpm.tt.gaus) %>%
    mutate(model = "o_c.tpm_turntaking_log_gaus"),
  broom.mixed::tidy(o_c.tpm.tt.gaus.v2) %>%
    mutate(model = "o_c.tpm_turntaking_log_gaus.v2"))

## Sequence duration random sample ####
seqdur.random.distribution <- ggplot(turn.sequences.rand,
                            aes(round((seq.dur*60),0))) +
  geom_histogram() +
  ylab("# of clips") +
  xlab ("Sequence duration (sec)") +
  basic.theme
ggsave(paste0(shiny.input.img.path, "seqdur_random_distribution.png"),
       plot = seqdur.random.distribution,
       width = 8, height = 6, dpi = 300)
#sd(round((turn.sequences.rand$seq.dur*60),0))^2
#mean(round((turn.sequences.rand$seq.dur*60),0))
# mean is much smaller than variance
# non-zero values -- no zero inflation
turn.sequences.rand$uniq.segment <- paste0(turn.sequences.rand$aclew_child_id, "-",
                                          turn.sequences.rand$segment)
seqdur.sec.rand.nb <- glmmTMB(round((seq.dur*60),0) ~
                                tchiyr.std +
                                stthr.tri +
                                hsz.std +
                                nsk.std +
                                tchiyr.std:stthr.tri +
                                tchiyr.std:hsz.std +
                                tchiyr.std:nsk.std +
                                (1|uniq.segment) +
                                (1|aclew_child_id),
                              data=turn.sequences.rand,
                              family="nbinom1")
#seqdur.sec.rand.nb.res = simulateResiduals(seqdur.sec.rand.nb)
#plot(seqdur.sec.rand.nb.res, rank = T) # (manually saved)
#summary(seqdur.sec.rand.nb) # (no significant effects)
#Conditional model:
#                               Estimate Std. Error z value Pr(>|z|)    
#(Intercept)                    2.216674   0.142451  15.561 <0.0001 ***
#tchiyr.std                     0.105663   0.193521   0.546 0.585    
#stthr.trimorning               0.136008   0.174476   0.780 0.436    
#stthr.triafternoon             0.128729   0.160016   0.804 0.421    
#hsz.std                        0.010461   0.083797   0.125 0.901    
#nsk.std                        0.006317   0.048557   0.130 0.896    
#tchiyr.std:stthr.trimorning    0.196619   0.188012   1.046 0.296    
#tchiyr.std:stthr.triafternoon  0.045138   0.176905   0.255 0.799    
#tchiyr.std:hsz.std             0.142867   0.113612   1.257 0.209    
#tchiyr.std:nsk.std            -0.032344   0.063522  -0.509 0.611 
# no results to save for reporting

# test the other two-way effects of time of day
seqdur.sec.rand.nb.v2 <- glmmTMB(round((seq.dur*60),0) ~
                                   tchiyr.std +
                                   stthr.tri.a +
                                   hsz.std +
                                   nsk.std +
                                   tchiyr.std:stthr.tri.a +
                                   tchiyr.std:hsz.std +
                                   tchiyr.std:nsk.std +
                                   (1|uniq.segment) +
                                   (1|aclew_child_id),
                                 data=turn.sequences.rand,
                                 family="nbinom1")
#summary(seqdur.sec.rand.nb.v2) # (no significant effects)
#stthr.tri.amidday             -0.128729   0.160016  -0.804  0.421    
#stthr.tri.amorning             0.007280   0.116200   0.063  0.950    
#tchiyr.std:stthr.tri.amidday  -0.045138   0.176905  -0.255  0.799    
#tchiyr.std:stthr.tri.amorning  0.151481   0.129067   1.174  0.241    
# no results to save for reporting

## Sequence duration tt sample ####
seqdur.tt.distribution <- ggplot(turn.sequences.tt,
                            aes(round((seq.dur*60),0))) +
  geom_histogram() +
  ylab("# of clips") +
  xlab ("Sequence duration (sec)") +
  basic.theme
ggsave(paste0(shiny.input.img.path, "seqdur_turntaking_distribution.png"),
       plot = seqdur.tt.distribution,
       width = 8, height = 6, dpi = 300)
#sd(round((turn.sequences.tt$seq.dur*60),0))^2
#mean(round((turn.sequences.tt$seq.dur*60),0))
# mean is much smaller than variance
# non-zero values
turn.sequences.tt$uniq.segment <- paste0(turn.sequences.tt$aclew_child_id, "-",
                                          turn.sequences.tt$segment)
seqdur.sec.tt.nb <- glmmTMB(round((seq.dur*60),0) ~
                                tchiyr.std +
                                stthr.tri +
                                hsz.std +
                                nsk.std +
                                tchiyr.std:stthr.tri +
                                tchiyr.std:hsz.std +
                                tchiyr.std:nsk.std +
                                (1|uniq.segment) +
                                (1|aclew_child_id),
                              data=turn.sequences.tt,
                              family="nbinom1")
#seqdur.sec.tt.nb.res = simulateResiduals(seqdur.sec.tt.nb)
#plot(seqdur.sec.tt.nb.res, rank = T) # (manually saved)
#summary(seqdur.sec.tt.nb)
#Conditional model:
#                              Estimate Std. Error z value  Pr(>|z|)    
#(Intercept)                    2.24687    0.13586  16.539  < 0.0001 ***
#tchiyr.std                    -0.18273    0.12115  -1.508  0.13149    
#stthr.trimorning               0.05915    0.16101   0.367  0.71336    
#stthr.triafternoon             0.37875    0.14522   2.608  0.00911 ** 
#hsz.std                       -0.17315    0.09966  -1.737  0.08233 .  
#nsk.std                       -0.01070    0.05785  -0.185  0.85333    
#tchiyr.std:stthr.trimorning   -0.02062    0.16977  -0.121  0.90331    
#tchiyr.std:stthr.triafternoon  0.01937    0.14365   0.135  0.89271    
#tchiyr.std:hsz.std            -0.17535    0.12841  -1.366  0.17207    
#tchiyr.std:nsk.std             0.03003    0.07884   0.381  0.70327
# save for reporting
seqdur.sec.tt.nb.COEF.morn <-
  coef(summary(seqdur.sec.tt.nb))[[1]]["stthr.trimorning",]
seqdur.sec.tt.nb.COEF.aft <-
  coef(summary(seqdur.sec.tt.nb))[[1]]["stthr.triafternoon",]
seqdur.sec.tt.nb.COEF.hsz <-
  coef(summary(seqdur.sec.tt.nb))[[1]]["hsz.std",]

# test the other two-way effects of time of day
seqdur.sec.tt.nb.v2 <- glmmTMB(round((seq.dur*60),0) ~
                                 tchiyr.std +
                                 stthr.tri.a +
                                 hsz.std +
                                 nsk.std +
                                 tchiyr.std:stthr.tri.a +
                                 tchiyr.std:hsz.std +
                                 tchiyr.std:nsk.std +
                                 (1|uniq.segment) +
                                 (1|aclew_child_id),
                               data=turn.sequences.tt,
                               family="nbinom1")
#summary(seqdur.sec.tt.nb.v2)
#stthr.tri.amidday             -0.37875    0.14522  -2.608              0.00911 ** 
#stthr.tri.amorning            -0.31960    0.15110  -2.115              0.03442 *  
#tchiyr.std:stthr.tri.amidday  -0.01937    0.14365  -0.135              0.89271    
#tchiyr.std:stthr.tri.amorning -0.04000    0.16605  -0.241              0.80965    
# save for reporting
seqdur.sec.tt.nb.v2.COEF.morn <-
  coef(summary(seqdur.sec.tt.nb.v2))[[1]]["stthr.tri.amorning",]
seqdur.sec.tt.nb.v2.COEF.age.morn <-
  coef(summary(seqdur.sec.tt.nb.v2))[[1]]["tchiyr.std:stthr.tri.amorning",]

## Sequence duration random sample ####
seqdur.sec.rand.gaus <- glmmTMB(log(seq.dur) ~
                                tchiyr.std +
                                stthr.tri +
                                hsz.std +
                                nsk.std +
                                tchiyr.std:stthr.tri +
                                tchiyr.std:hsz.std +
                                tchiyr.std:nsk.std +
                                (1|uniq.segment) +
                                (1|aclew_child_id),
                              data=turn.sequences.rand)
#seqdur.sec.rand.gaus.res = simulateResiduals(seqdur.sec.rand.gaus)
#plot(seqdur.sec.rand.gaus.res, rank = T) # (manually saved)
#summary(seqdur.sec.rand.gaus) # (no significant effects)
#Conditional model:
#                               Estimate Std. Error z value Pr(>|z|)    
#(Intercept)                   -2.238598   0.179612 -12.464 <0.0001 ***
#tchiyr.std                     0.101763   0.247126   0.412 0.680    
#stthr.trimorning               0.079977   0.222088   0.360 0.719    
#stthr.triafternoon             0.132569   0.203187   0.652 0.514    
#hsz.std                        0.007682   0.107943   0.071 0.943    
#nsk.std                        0.010111   0.060352   0.168 0.867    
#tchiyr.std:stthr.trimorning    0.334872   0.241382   1.387 0.165    
#tchiyr.std:stthr.triafternoon  0.109382   0.223645   0.489 0.625    
#tchiyr.std:hsz.std             0.166068   0.144025   1.153 0.249    
#tchiyr.std:nsk.std            -0.059959   0.080566  -0.744 0.457
# DIFFERENCES?
# very similar to nb model

# test the other two-way effects of time of day
seqdur.sec.rand.gaus.v2 <- glmmTMB(log(seq.dur) ~
                                   tchiyr.std +
                                   stthr.tri.a +
                                   hsz.std +
                                   nsk.std +
                                   tchiyr.std:stthr.tri.a +
                                   tchiyr.std:hsz.std +
                                   tchiyr.std:nsk.std +
                                   (1|uniq.segment) +
                                   (1|aclew_child_id),
                                 data=turn.sequences.rand)
#summary(seqdur.sec.rand.gaus) # (no significant effects)
#stthr.tri.amidday             -0.132569   0.203187  -0.652  0.514    
#stthr.tri.amorning            -0.052592   0.146310  -0.359  0.719    
#tchiyr.std:stthr.tri.amidday  -0.109382   0.223645  -0.489  0.625    
#tchiyr.std:stthr.tri.amorning  0.225490   0.167873   1.343  0.179    
# DIFFERENCES?
# very similar to nb model


## Sequence duration tt sample ####
seqdur.sec.tt.gaus <- glmmTMB(log(seq.dur) ~
                                tchiyr.std +
                                stthr.tri +
                                hsz.std +
                                nsk.std +
                                tchiyr.std:stthr.tri +
                                tchiyr.std:hsz.std +
                                tchiyr.std:nsk.std +
                                (1|uniq.segment) +
                                (1|aclew_child_id),
                              data=turn.sequences.tt)
#seqdur.sec.tt.gaus.res = simulateResiduals(seqdur.sec.tt.gaus)
#plot(seqdur.sec.tt.gaus.res, rank = T) # (manually saved)
#summary(seqdur.sec.tt.gaus)
#Conditional model:
#                              Estimate Std. Error z value Pr(>|z|)    
#(Intercept)                   -2.32810    0.15565 -14.957 <0.0001 ***
#tchiyr.std                    -0.20332    0.14515  -1.401 0.1613    
#stthr.trimorning               0.07604    0.19411   0.392 0.6953    
#stthr.triafternoon             0.57536    0.18525   3.106 0.0019 ** 
#hsz.std                       -0.23223    0.11533  -2.014 0.0441 *  
#nsk.std                       -0.02348    0.07391  -0.318 0.7507    
#tchiyr.std:stthr.trimorning    0.01534    0.20032   0.077 0.9390    
#tchiyr.std:stthr.triafternoon -0.01424    0.18700  -0.076 0.9393    
#tchiyr.std:hsz.std            -0.19712    0.14911  -1.322 0.1862    
#tchiyr.std:nsk.std             0.04998    0.10153   0.492 0.6226
# DIFFERENCES?
# broadly similar to nb model

# test the other two-way effects of time of day
seqdur.sec.tt.gaus.v2 <- glmmTMB(log(seq.dur) ~
                                 tchiyr.std +
                                 stthr.tri.a +
                                 hsz.std +
                                 nsk.std +
                                 tchiyr.std:stthr.tri.a +
                                 tchiyr.std:hsz.std +
                                 tchiyr.std:nsk.std +
                                 (1|uniq.segment) +
                                 (1|aclew_child_id),
                               data=turn.sequences.tt)
#summary(seqdur.sec.tt.gaus.v2)
#stthr.tri.amidday             -0.57536    0.18525  -3.106              0.0019 ** 
#stthr.tri.amorning            -0.49932    0.19388  -2.575              0.0100 *  
#tchiyr.std:stthr.tri.amidday   0.01424    0.18700   0.076              0.9393    
#tchiyr.std:stthr.tri.amorning  0.02958    0.21021   0.141              0.8881    
# DIFFERENCES?
# very similar to nb model

# Write model results out for input to shiny
seqdur.sec.models <- bind_rows(
  broom.mixed::tidy(seqdur.sec.rand.nb) %>%
    mutate(model = "seqdur.sec_random_nb"),
  broom.mixed::tidy(seqdur.sec.rand.nb.v2) %>%
    mutate(model = "seqdur.sec_random_nb.v2"),
  broom.mixed::tidy(seqdur.sec.rand.gaus) %>%
    mutate(model = "seqdur.sec_random_log_gaus"),
  broom.mixed::tidy(seqdur.sec.rand.gaus.v2) %>%
    mutate(model = "seqdur.sec_random_log_gaus.v2"),
  broom.mixed::tidy(seqdur.sec.tt.nb) %>%
    mutate(model = "seqdur.sec_turntaking_nb"),
  broom.mixed::tidy(seqdur.sec.tt.nb.v2) %>%
    mutate(model = "seqdur.sec_turntaking_nb.v2"),
  broom.mixed::tidy(seqdur.sec.tt.gaus) %>%
    mutate(model = "seqdur.sec_turntaking_log_gaus"),
  broom.mixed::tidy(seqdur.sec.tt.gaus.v2) %>%
    mutate(model = "seqdur.sec_turntaking_log_gaus.v2"))

# get average overall turn transition rate
ttrate.all.avg.bychild <- turn.transitions %>%
  filter(grepl("tt", segment)) %>%
  replace_na(list(tm1.tier = 0, tp1.tier = 0)) %>%
  mutate(tm1.bin = ifelse(tm1.tier == 0, 0, 1),
         tp1.bin = ifelse(tp1.tier == 0, 0, 1)) %>%
  group_by(aclew_child_id, segment) %>%
  summarise(oct.raw = sum(tm1.bin), cot.raw = sum(tp1.bin)) %>%
  mutate(segdur = ifelse(grepl('ext', segment), 5, 1)) %>%
  group_by(aclew_child_id) %>%
  summarise (alltts = sum(oct.raw + cot.raw), sampledur = sum(segdur)) %>%
  mutate(peakttrate = alltts/sampledur)
ttrate.all.avg <- mean(ttrate.all.avg.bychild$peakttrate)

# set up a table of tt onset times within the randomly sampled segments
tts.all.rand <- turn.transitions %>%
  filter(grepl("rand", segment)) %>%
  select(aclew_child_id, segment, tm1.stop, tp1.start)
rand.seg.starts <- seg.info %>%
  filter(grepl("rand", CodeName)) %>%
  select(aclew_id, CodeName, clipoffset.sec)
rand.seg.secs <- tibble()
for (i in 1:nrow(rand.seg.starts)) {
  seg.secs <- tibble(
    aclew_child_id = rand.seg.starts$aclew_id[i],
    segment = rand.seg.starts$CodeName[i],
    segoffset.sec = c(0:299), # seconds in a (5-min) random clip
    segtime.sec = seq(rand.seg.starts$clipoffset.sec[i],
                      (rand.seg.starts$clipoffset.sec[i] + 299), 1),
    bin.tt = 0
  )
  tts.all.rand.chi <- tts.all.rand %>%
    filter(aclew_child_id == rand.seg.starts$aclew_id[i],
           segment == rand.seg.starts$CodeName[i])
  # add in tm1s
  tts.all.rand.chi.tm1 <- tts.all.rand.chi$tm1.stop[!is.na(tts.all.rand.chi$tm1.stop)]
  if (length(tts.all.rand.chi.tm1) > 0) {
    tm1.onsets.sec <- tts.all.rand.chi.tm1/1000
    # REMINDER: example use of findInterval:
    # findInterval(c(2,5,7,8), c(4,5,6,8,9,10)) => [1] 0 2 3 4
    tt.idx <- findInterval(tm1.onsets.sec, seg.secs$segtime.sec)
    seg.secs$bin.tt[tt.idx] <- seg.secs$bin.tt[tt.idx] + 1
  }
  # add in tp1s
  tts.all.rand.chi.tp1 <- tts.all.rand.chi$tp1.start[!is.na(tts.all.rand.chi$tp1.start)]
  if (length(tts.all.rand.chi.tp1) > 0) {
    tp1.onsets.sec <- tts.all.rand.chi.tp1/1000
    tt.idx <- findInterval(tp1.onsets.sec, seg.secs$segtime.sec)
    seg.secs$bin.tt[tt.idx] <- seg.secs$bin.tt[tt.idx] + 1
  }
  rand.seg.secs <- bind_rows(rand.seg.secs, seg.secs)
}

# get average tt rates in 1-min windows
rand.window.tts <- rand.seg.secs %>%
  filter(segoffset.sec >= 59) %>%
  select(-bin.tt) %>%
  mutate(ttr.min = 0)
for (chi in unique(rand.seg.secs$aclew_child_id)) {
  for (seg in unique(rand.seg.secs$segment)) {
    for (minstart in c(1:241)) {
      rand.window.tts$ttr.min[which(rand.window.tts$aclew_child_id == chi &
                              rand.window.tts$segment == seg &
                              rand.window.tts$segoffset.sec == minstart+59-1)] <-
        sum(subset(rand.seg.secs, aclew_child_id == chi &
                     segment == seg)$bin.tt[(minstart:(minstart+59))])
    }
  }
}
rand.window.tts <- rand.window.tts %>%
  left_join(select(ttrate.all.avg.bychild, c(aclew_child_id, peakttrate))) %>%
  mutate(GE.chipeak = ifelse(ttr.min >= peakttrate, 1, 0),
         GE.allpeak = ifelse(ttr.min >= ttrate.all.avg, 1, 0))
write_csv(rand.window.tts,
          paste0(metadata.path, "ttr.random.1minwindow.csv"))

# calculate # seconds in "peak" tt rates and duration of peaks
ttrGEpeak.segs <- rand.window.tts %>%
  filter(GE.allpeak == 1)
ttrGEpeak.segs$prevpeak <- c(-1, ttrGEpeak.segs$segoffset.sec[1:nrow(ttrGEpeak.segs)-1])
ttrGEpeak.segs$newstreak <- ifelse(ttrGEpeak.segs$prevpeak ==
                                        (ttrGEpeak.segs$segoffset.sec - 1), 0, 1)
ttrGEpeak.segs <- select(ttrGEpeak.segs, -prevpeak)
ttrGEpeaks <- filter(ttrGEpeak.segs, newstreak == 1)
ttrGEpeaks$end.sec <- 0
for (i in 1:nrow(ttrGEpeaks)) {
  chi <- ttrGEpeaks$aclew_child_id[i]
  seg <- ttrGEpeaks$segment[i]
  start <- ttrGEpeaks$segoffset.sec[i]
  start.idx <- which(ttrGEpeak.segs$aclew_child_id == chi &
                       ttrGEpeak.segs$segment == seg &
                       ttrGEpeak.segs$segoffset.sec == start)
  newstreaks <- which(ttrGEpeak.segs$aclew_child_id == chi &
                        ttrGEpeak.segs$segment == seg &
                        ttrGEpeak.segs$segoffset.sec > start &
                        ttrGEpeak.segs$newstreak == 1)
  if (length(newstreaks) > 0) {
    end.sec <- ttrGEpeak.segs$segoffset.sec[newstreaks[1] - 1]
  } else {
    end.sec <- ttrGEpeak.segs$segoffset.sec[max(
      which(ttrGEpeak.segs$aclew_child_id == chi &
              ttrGEpeak.segs$segment == seg &
              ttrGEpeak.segs$segoffset.sec >= start))]
  }
  ttrGEpeaks$end.sec[i] <- end.sec + 1
}
ttrGEpeaks <- ttrGEpeaks %>%
  mutate(start.sec = segoffset.sec - 59,
         peak.dur = end.sec - start.sec)
ttrGEpeaks.summary <- ttrGEpeaks %>%
  group_by(aclew_child_id) %>%
  summarise(npeaks = n(),
            pkdur.mean = mean(peak.dur),
            pkdur.sum = sum(peak.dur),
            pkdur.mph = ((pkdur.sum/60)/45) * 60) %>% # (peak mns/45 poss mns) * 60 mn/hr
  # add back zero estimates for children with no tt peaks in the random data
  full_join(select(ptcp.info, aclew_child_id)) %>%
  replace_na(list(npeaks = 0, pkdur.mean = 0, pkdur.sum = 0, pkdur.mph = 0))

ttrGEpeak.segs.chi <- rand.window.tts %>%
  filter(GE.chipeak == 1)

avg.random.peak.durs <- ttrGEpeaks %>%
  group_by(aclew_child_id) %>%
  summarise(meandur = mean(peak.dur), 
            meddur = median(peak.dur))

rand.window.tts <- rand.window.tts %>%
  mutate(endoffsetsec = segoffset.sec + 1,
         delta.allpeak = ttr.min - ttrate.all.avg) %>%
  left_join(ptcp.info)
write_csv(rand.window.tts, paste0(shiny.input.path, "random_window_tts.csv"))
ttr.random <- ggplot(rand.window.tts,
                          aes(x = as.factor(age_mo_round), y = ttr.min)) +
  geom_jitter(color = "gray80", size = 0.3) +
  geom_violin(color = "black", fill = "black") + 
  ylab("High-tt secs/min") + xlab("Age (months)")	+
  scale_y_continuous(limits=c(-1, round(max(rand.window.tts$ttr.min) + 5, -1)),
                     breaks=seq(0, round(max(rand.window.tts$ttr.min) + 5, -1), 10)) +
  coord_cartesian(ylim=c(0, round(max(rand.window.tts$ttr.min) + 5, -1))) +
  geom_hline(yintercept = ttrate.all.avg, color = "red") +
  geom_point(color = "red", size = 2.5, aes(y = peakttrate)) +
#  scale_color_manual(values = viridis(3)) +
#  scale_fill_manual(values = viridis(3)) +
  theme_apa() +
  theme(legend.position="none",
        axis.line = element_line(color="black", size = 0.4))

## Vocal maturity
# all vocalization types
chi.vm.lx.utts <- all.data %>%
  filter((tier == "vcm@CHI" | tier == "lex@CHI" | tier == "mwu@CHI") & !is.na(val)) %>%
  mutate(voc.rating = ifelse(val == "M", 4,
                             ifelse((val == "1" | val == "W"), 3,
                                    ifelse((val == "0" | val == "C"), 2,
                                           ifelse(val == "N", 1, 0))))) %>%
  filter(voc.rating > 0) %>%
  group_by(aclew_child_id, segment, sample, start) %>%
  summarise(max_voc.rtg = max(voc.rating))
all.voc.types.per.child <- tibble(
  aclew_child_id = rep(ptcp.info$aclew_child_id, 4),
  max_voc.rtg = c(rep(1, length(ptcp.info$aclew_child_id)),
                  rep(2, length(ptcp.info$aclew_child_id)),
                  rep(3, length(ptcp.info$aclew_child_id)),
                  rep(4, length(ptcp.info$aclew_child_id)))
)
chi.nvocs <- chi.vm.lx.utts %>%
  group_by(aclew_child_id) %>%
  summarise(n_vocs = n())
chi.vm.lx.voc.type.props <- chi.vm.lx.utts %>%
  group_by(aclew_child_id, max_voc.rtg) %>%
  summarise(n_voc.type = n()) %>%
  full_join(all.voc.types.per.child, by = c("aclew_child_id", "max_voc.rtg")) %>%
  replace_na(list(n_voc.type = 0)) %>%
  full_join(chi.nvocs, by = "aclew_child_id") %>%
  mutate(prop_voc.type = round(n_voc.type/n_vocs, 3)) %>%
  arrange(aclew_child_id, max_voc.rtg) %>%
  full_join(ptcp.info, by = "aclew_child_id")

# speech-like vs. non-speech-like only, only under 19mo
chi.vm.lx.utts.all <- all.data %>%
  filter((tier == "vcm@CHI" | tier == "lex@CHI" | tier == "mwu@CHI") &
           !is.na(val) & age_mo_round < 19) %>%
  mutate(voc.rating = ifelse(val == "M", 1,
                             ifelse((val == "1" | val == "W"), 1,
                                    ifelse((val == "0" | val == "C"), 1, 0)))) %>%
  group_by(aclew_child_id, segment, sample, start) %>%
  summarise(max_voc.rtg = max(voc.rating))
all.voc.types.per.child.all <- tibble(
  aclew_child_id = rep(ptcp.info$aclew_child_id, 2),
  max_voc.rtg = c(rep(0, length(ptcp.info$aclew_child_id)),
                  rep(1, length(ptcp.info$aclew_child_id)))
)
chi.nvocs.all <- chi.vm.lx.utts.all %>%
  group_by(aclew_child_id) %>%
  summarise(n_vocs = n())
chi.vm.lx.voc.type.props.bin <- chi.vm.lx.utts.all %>%
  group_by(aclew_child_id, max_voc.rtg) %>%
  summarise(n_voc.type = n()) %>%
  left_join(all.voc.types.per.child.all, by = c("aclew_child_id", "max_voc.rtg")) %>%
  replace_na(list(n_voc.type = 0)) %>%
  left_join(chi.nvocs.all, by = "aclew_child_id") %>%
  mutate(prop_voc.type = round(n_voc.type/n_vocs, 3)) %>%
  arrange(aclew_child_id, max_voc.rtg) %>%
  left_join(ptcp.info, by = "aclew_child_id")

chi.vm.lx.voc.type.props <- chi.vm.lx.voc.type.props %>%
  mutate(voc.type = factor(as.factor(max_voc.rtg),
                           labels = c("NCB", "CB", "SW", "MW")))
write_csv(chi.vm.lx.voc.type.props,
          paste0(shiny.input.path,
                 "all_vocmat-types_proportions.csv"))
voc.mat.by.age <- ggplot(
  data = chi.vm.lx.voc.type.props,
  aes(x = age_mo_round, y = prop_voc.type, group = as.factor(voc.type))) +
  geom_point(aes(color = as.factor(voc.type))) +
  geom_smooth(aes(color = as.factor(voc.type),
                  fill = as.factor(voc.type)), method = "loess") +
  ylab("Prop of linguistic vocs") + xlab("Child age (mo)") +
  labs(fill='Voc type') +
  labs(color='Voc type') +
  scale_color_manual(values = c("gray80", "gray54",
                                  "gray27", "black")) +
  scale_fill_manual(values = c("gray80", "gray54",
                                  "gray27", "black")) +
  scale_y_continuous(limits=c(-0.5,1.5),
                     breaks=seq(0,1,0.2)) +
  scale_x_continuous(limits=c(0,38),
                     breaks=seq(0,38,6)) +
  coord_cartesian(ylim=c(0,1), xlim=c(0,38)) +
  theme_apa() + theme(
    legend.position = c(0.9, 0.85),
    legend.background = element_rect(fill="transparent"),
    axis.line = element_line(color="black", size = 0.4))
voc.mat.by.age

all.models <- bind_rows(TCDS.models, ODS.models,
                        c_o.tpm.models, o_c.tpm.models,
                        seqdur.sec.models)
write_csv(all.models, paste0(shiny.input.path, "all_model_tables.csv"))
