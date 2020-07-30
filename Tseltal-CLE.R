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
         "Mother's age" = mat_age_rd, "Level of maternal education" = mat_ed,
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
  scale_color_manual(values = c(viridis(3)[1], viridis(3)[2], "blue"))
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

# Add in speaker sex
all.data$SpkrSex <- "Not known"
all.data$SpkrSex[grepl("FA|FC|FU", all.data$speaker)] <- "Female"
all.data$SpkrSex[grepl("MA|MC|MU", all.data$speaker)] <- "Male"
all.rand.segments.age.sx <- tibble(
  aclew_child_id = rep(unique(all.data$aclew_child_id),
                4*n.unique.rand.segs),
  segment = rep(unique(all.data$segment[grepl("random", all.data$segment)]),
                4*n.unique.recs),
  SpkrAge = c(rep("Adult", (n.unique.rand.segs * n.unique.recs)),
              rep("Child", (n.unique.rand.segs * n.unique.recs)),
              rep("Adult", (n.unique.rand.segs * n.unique.recs)),
              rep("Child", (n.unique.rand.segs * n.unique.recs))),
  SpkrSex = c(rep("Female", (n.unique.rand.segs * (2 * n.unique.recs))),
              rep("Male", (n.unique.rand.segs * (2 * n.unique.recs)))))

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
quantity.rand$stthr.tri.o <- factor(quantity.rand$stthr.tri,
                                  levels = c("morning", "afternoon", "midday"))

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
quantity.rand.sa$stthr.tri.o <- factor(quantity.rand.sa$stthr.tri,
                                       levels = c("morning", "afternoon", "midday"))


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
quantity.nonrand.tt$stthr.tri.o <- factor(quantity.nonrand.tt$stthr.tri,
                                          levels = c("morning", "afternoon", "midday"))

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
quantity.nonrand.tt.sa$stthr.tri.o <- factor(quantity.nonrand.tt.sa$stthr.tri,
                                             levels = c("morning", "afternoon", "midday"))


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
  ylab("TCDS (min/hr)") + xlab("")	+
  ggtitle("Random") +
  labs(alpha = "Age") +
  scale_y_continuous(limits=c(-10,40),
                     breaks=c(0,10,20, 30, 40)) +
  scale_color_manual(guide = FALSE, values = viridis(3)) +
  scale_fill_manual(guide = FALSE, values = viridis(3)) +
  scale_alpha_discrete(range = c(0.2, 0.8)) +
  coord_cartesian(ylim=c(0,40)) +
  theme_apa() +
  theme(legend.position="none",
        axis.line = element_line(color="black", size = 0.4),
        plot.title = element_text(hjust = 0.5))

tod.ods.rand <- ggplot(
  data = subset(quantity.rand_and_tt.all, Sample == "Random"), aes(
  y = round(ods_mph, 0), x = stthr.tri.centered,
  group = interaction(SplitAge, stthr.tri.centered),
  color = Sample,
  fill = Sample,
  alpha = SplitAge)) +
  geom_jitter() +
  geom_boxplot(color = "black") +
  ylab("ODS (min/hr)") + xlab("")	+
  ggtitle("Random") +
  labs(alpha = "Age") +
  scale_y_continuous(limits=c(-10,80),
                     breaks=seq(0,80,20)) +
  scale_color_manual(guide = FALSE, values = viridis(3)) +
  scale_fill_manual(guide = FALSE, values = viridis(3)) +
  scale_alpha_discrete(range = c(0.2, 0.8)) +
  coord_cartesian(ylim=c(0,80)) +
  theme_apa() +
  theme(legend.position="none",
        axis.line = element_line(color="black", size = 0.4),
        plot.title = element_text(hjust = 0.5))

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
  ggtitle("Turn taking") +
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
        axis.line = element_line(color="black", size = 0.4),
        plot.title = element_text(hjust = 0.5))

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
  ggtitle("Turn taking") +
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
        axis.line = element_line(color="black", size = 0.4),
        plot.title = element_text(hjust = 0.5))


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
  geom_histogram(binwidth = 2) +
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
                           tchiyr.std:nsk.std +
                           (1|aclew_child_id),
                         data=quantity.rand,
                         ziformula=~nsk.std,
                         family="nbinom1")
#tds.rand.zinb.res = simulateResiduals(tds.rand.zinb)
#plot(tds.rand.zinb.res, rank = T) # (manually saved)
#summary(tds.rand.zinb)
# Data: quantity.rand
# 
#      AIC      BIC   logLik deviance df.resid 
#    397.9    430.4   -186.0    371.9       77 
# 
# Random effects:
# 
# Conditional model:
#  Groups         Name        Variance Std.Dev.
#  aclew_child_id (Intercept) 0.09779  0.3127  
# Number of obs: 90, groups:  aclew_child_id, 10
# 
# Overdispersion parameter for nbinom1 family (): 4.39 
# 
# Conditional model:
#                               Estimate Std. Error z value Pr(>|z|)   
# (Intercept)                    0.91432    0.36101   2.533  0.01132 * 
# tchiyr.std                     0.59969    0.35796   1.675  0.09387 . 
# stthr.trimorning               0.82768    0.39530   2.094  0.03628 * 
# stthr.triafternoon             0.48838    0.37341   1.308  0.19091   
# hsz.std                        0.00822    0.22447   0.037  0.97079   
# nsk.std                       -0.12070    0.16013  -0.754  0.45098   
# tchiyr.std:stthr.trimorning   -0.28375    0.38955  -0.728  0.46636   
# tchiyr.std:stthr.triafternoon -0.85346    0.37781  -2.259  0.02389 * 
# tchiyr.std:nsk.std             0.57376    0.19449   2.950  0.00318 **
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Zero-inflation model:
#             Estimate Std. Error z value Pr(>|z|)
# (Intercept)   -57.43   15426.18  -0.004    0.997
# nsk.std       -55.68   15691.06  -0.004    0.997
# save for reporting
tds.rand.zinb.disp <- round(sigma(tds.rand.zinb), 2)
tds.rand.zinb.COEF.age <-
  coef(summary(tds.rand.zinb))[[1]]["tchiyr.std",] 
tds.rand.zinb.COEF.midd.morn <-
  coef(summary(tds.rand.zinb))[[1]]["stthr.trimorning",] 
tds.rand.zinb.COEF.midd.aft <-
  coef(summary(tds.rand.zinb))[[1]]["stthr.triafternoon",] 
tds.rand.zinb.COEF.age.midd.morn <-
  coef(summary(tds.rand.zinb))[[1]]["tchiyr.std:stthr.trimorning",]
tds.rand.zinb.COEF.age.midd.aft <-
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
                           tchiyr.std:nsk.std +
                           (1|aclew_child_id),
                         data=quantity.rand,
                         ziformula=~nsk.std,
                         family="nbinom1")
#summary(tds.rand.zinb.v2)
# stthr.tri.amorning             0.339309   0.268777   1.262  0.20680    
# tchiyr.std:stthr.tri.amorning  0.569710   0.300505   1.896  0.05798 .  
# save for reporting
tds.rand.zinb.v2.COEF.aft.morn <-
  coef(summary(tds.rand.zinb.v2))[[1]]["stthr.tri.amorning",] 
tds.rand.zinb.v2.COEF.age.aft.morn <-
  coef(summary(tds.rand.zinb.v2))[[1]]["tchiyr.std:stthr.tri.amorning",] 

## TCDS tt sample ####
tcds.tt.distribution <- ggplot(quantity.nonrand.tt,
                            aes(round(tds_mph,0))) +
  geom_histogram(binwidth = 2) +
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
                       tchiyr.std:nsk.std +
                       (1|aclew_child_id),
                     data=quantity.nonrand.tt,
                     family="nbinom1")
#tds.tt.nb.res = simulateResiduals(tds.tt.nb)
#plot(tds.tt.nb.res, rank = T) # (manually saved)
#summary(tds.tt.nb)
# Data: quantity.nonrand.tt
# 
#      AIC      BIC   logLik deviance df.resid 
#    432.2    455.0   -205.1    410.2       48 
# 
# Random effects:
# 
# Conditional model:
#  Groups         Name        Variance Std.Dev.
#  aclew_child_id (Intercept) 0.03781  0.1944  
# Number of obs: 59, groups:  aclew_child_id, 10
# 
# Overdispersion parameter for nbinom1 family (): 4.98 
# 
#                                Estimate Std. Error z value Pr(>|z|)   
# (Intercept)                    2.515650   0.222213  11.321 <0.001 ***
# tchiyr.std                     0.081027   0.212018   0.382  0.702    
# stthr.trimorning               0.138024   0.287574   0.480  0.631    
# stthr.triafternoon             0.062270   0.267289   0.233  0.816    
# hsz.std                        0.117579   0.136842   0.859  0.390    
# nsk.std                       -0.129227   0.105435  -1.226  0.220    
# tchiyr.std:stthr.trimorning   -0.133363   0.285654  -0.467  0.641    
# tchiyr.std:stthr.triafternoon  0.001369   0.236523   0.006  0.995    
# tchiyr.std:nsk.std             0.059666   0.131088   0.455  0.649  
# save for reporting
tds.tt.nb.disp <- round(sigma(tds.tt.nb), 2)

# test the other two-way effects of time of day
tds.tt.nb.v2 <- glmmTMB(round(tds_mph,0) ~
                       tchiyr.std +
                       stthr.tri.a +
                       hsz.std +
                       nsk.std +
                       tchiyr.std:stthr.tri.a +
                       tchiyr.std:nsk.std +
                       (1|aclew_child_id),
                     data=quantity.nonrand.tt,
                     family="nbinom1")
#summary(tds.tt.nb.v2)
# stthr.tri.amorning             0.07575    0.22426   0.338   0.736    
# tchiyr.std:stthr.tri.amorning -0.13474    0.26270  -0.513   0.608    
# no results to save for reporting


## TCDS random sample ####
tds.rand.gaus <- glmmTMB(log(tds_mph+1) ~
                           tchiyr.std +
                           stthr.tri +
                           hsz.std +
                           nsk.std +
                           tchiyr.std:stthr.tri +
                           tchiyr.std:nsk.std +
                           (1|aclew_child_id),
                         data=quantity.rand)
#tds.rand.gaus.res = simulateResiduals(tds.rand.gaus)
#plot(tds.rand.gaus.res, rank = T) # (manually saved)
#summary(tds.rand.gaus)
# Data: quantity.rand
# 
#      AIC      BIC   logLik deviance df.resid 
#    250.4    277.9   -114.2    228.4       79 
# 
# Random effects:
# 
# Conditional model:
#  Groups         Name        Variance Std.Dev.
#  aclew_child_id (Intercept) 0.04233  0.2057  
#  Residual                   0.70589  0.8402  
# Number of obs: 90, groups:  aclew_child_id, 10
# 
# Dispersion estimate for gaussian family (sigma^2): 0.706 
# 
# Conditional model:
#                               Estimate Std. Error z value  Pr(>|z|)    
# (Intercept)                     0.8188     0.1891   4.329  0.000015 ***
# tchiyr.std                      0.5404     0.2235   2.418  0.01560 *
# stthr.trimorning                0.5048     0.2505   2.015  0.04386 *  
# stthr.triafternoon              0.2897     0.2211   1.310  0.19012    
# hsz.std                        -0.1547     0.1565  -0.989  0.32269    
# nsk.std                         0.2290     0.1185   1.933  0.05328 .
# tchiyr.std:stthr.trimorning    -0.1737     0.2676  -0.649  0.51642    
# tchiyr.std:stthr.triafternoon  -0.6810     0.2388  -2.851  0.00436 ** 
# tchiyr.std:nsk.std              0.2308     0.1392   1.658  0.09726 .  
# DIFFERENCES between z-i nb and logged gaussian?
# broadly similar

# test the other two-way effects of time of day
tds.rand.gaus.v2 <- glmmTMB(log(tds_mph+1) ~
                           tchiyr.std +
                           stthr.tri.a +
                           hsz.std +
                           nsk.std +
                           tchiyr.std:stthr.tri.a +
                           tchiyr.std:nsk.std +
                           (1|aclew_child_id),
                         data=quantity.rand)
#summary(tds.rand.gaus.v2)
# stthr.tri.amorning              0.2151     0.2193   0.981  0.32662    
# tchiyr.std:stthr.tri.amorning   0.5073     0.2298   2.207  0.02729 *  
# DIFFERENCES between z-i nb and logged gaussian?
# broadly similar

## TCDS tt sample ####
tds.tt.gaus <- glmmTMB(log(tds_mph+1) ~
                       tchiyr.std +
                       stthr.tri +
                       hsz.std +
                       nsk.std +
                       tchiyr.std:stthr.tri +
                       tchiyr.std:nsk.std +
                       (1|aclew_child_id),
                     data=quantity.nonrand.tt)
#tds.tt.gaus.res = simulateResiduals(tds.tt.gaus)
#plot(tds.tt.gaus.res, rank = T) # (manually saved)
#summary(tds.tt.gaus)
# Data: quantity.nonrand.tt
# 
#      AIC      BIC   logLik deviance df.resid 
#    153.2    176.1    -65.6    131.2       48 
# 
# Random effects:
# 
# Conditional model:
#  Groups         Name        Variance Std.Dev.
#  aclew_child_id (Intercept) 0.04876  0.2208  
#  Residual                   0.50110  0.7079  
# Number of obs: 59, groups:  aclew_child_id, 10
# 
# Dispersion estimate for gaussian family (sigma^2): 0.501 
# 
# Conditional model:
#                               Estimate Std. Error z value  Pr(>|z|)    
# (Intercept)                    2.40519    0.25565   9.408 <0.001 ***
# tchiyr.std                     0.08664    0.23371   0.371  0.711    
# stthr.trimorning               0.12696    0.33565   0.378  0.705    
# stthr.triafternoon             0.05228    0.30508   0.171  0.864    
# hsz.std                        0.13424    0.15072   0.891  0.373    
# nsk.std                       -0.13719    0.11791  -1.164  0.245    
# tchiyr.std:stthr.trimorning   -0.16565    0.31907  -0.519  0.604    
# tchiyr.std:stthr.triafternoon  0.04031    0.27258   0.148  0.882    
# tchiyr.std:nsk.std             0.07303    0.14829   0.492  0.622    
# DIFFERENCES between z-i nb and logged gaussian?
# broadly similar


# test the other two-way effects of time of day
tds.tt.gaus.v2 <- glmmTMB(log(tds_mph+1) ~
                       tchiyr.std +
                       stthr.tri.a +
                       hsz.std +
                       nsk.std +
                       tchiyr.std:stthr.tri.a +
                       tchiyr.std:nsk.std +
                       (1|aclew_child_id),
                     data=quantity.nonrand.tt)
#summary(tds.tt.gaus.v2)
# stthr.tri.amorning             0.07467    0.25899   0.288   0.773    
# tchiyr.std:stthr.tri.amorning -0.20595    0.29278  -0.703   0.482    
# DIFFERENCES between z-i nb and logged gaussian?
# broadly similar

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

tds.per.seg.rand.age.sx <- all.data %>%
  filter(sample == "random" & speaker != "CHI" & SpkrAge != "Not known" &
           grepl("xds@", tier) & val == "T") %>%
  group_by(aclew_child_id, SpkrAge, SpkrSex, segment, segment_dur) %>%
  summarise(tds_min.age.sx = round(sum(dur)/60000,3)) %>%
  full_join(all.rand.segments.age.sx, by = c("aclew_child_id", "segment", "SpkrAge", "SpkrSex")) %>%
  replace_na(list(segment_dur = 5, tds_min.age.sx = 0)) %>%
  mutate(tds_mph.age.sx = (tds_min.age.sx/segment_dur)*60) %>%
  arrange(aclew_child_id, segment, SpkrAge, SpkrSex) %>%
  group_by(aclew_child_id, SpkrAge, SpkrSex) %>%
  summarize(m_tds_mph.age.sx = mean(tds_mph.age.sx))
tds.per.seg.rand.age.sx.FA <- filter(tds.per.seg.rand.age.sx,
                                     SpkrAge == "Adult" & SpkrSex == "Female") %>%
  rename(m_tds_mph.age.sx.FA = m_tds_mph.age.sx)
tds.per.seg.rand.age.sx.MA <- filter(tds.per.seg.rand.age.sx,
                                     SpkrAge == "Adult" & SpkrSex == "Male") %>%
  rename(m_tds_mph.age.sx.MA = m_tds_mph.age.sx)
adu.sx.ratio <- full_join(tds.per.seg.rand.age.sx.FA, tds.per.seg.rand.age.sx.MA,
                          by = "aclew_child_id") %>%
  mutate(ratioFAMA = m_tds_mph.age.sx.FA/m_tds_mph.age.sx.MA) %>%
  dplyr::select(aclew_child_id, ratioFAMA)


## ODS random sample ####
ods.random.distribution <- ggplot(quantity.rand,
                            aes(round(ods_mph,0))) +
  geom_histogram(binwidth = 2) +
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
                           tchiyr.std:nsk.std +
                           (1|aclew_child_id),
                         data=quantity.rand,
                         ziformula=~nsk.std,
                         family="nbinom1")
#ods.rand.zinb.res = simulateResiduals(ods.rand.zinb)
#plot(ods.rand.zinb.res, rank = T) # (manually saved)
#summary(ods.rand.zinb)
# Data: quantity.rand
# 
#      AIC      BIC   logLik deviance df.resid 
#    542.9    575.4   -258.4    516.9       77 
# 
# Random effects:
# 
# Conditional model:
#  Groups         Name        Variance                      Std.Dev.         
#  aclew_child_id (Intercept) 0.000000000000000000000003448 0.000000000001857
# Number of obs: 90, groups:  aclew_child_id, 10
# 
# Overdispersion parameter for nbinom1 family ():  6.5 
# 
# Conditional model:
#                               Estimate Std. Error z value   Pr(>|z|)    
# (Intercept)                    2.70946    0.16060  16.871 < 0.0001 ***
# tchiyr.std                    -0.38657    0.15935  -2.426   0.0153 *  
# stthr.trimorning               0.44961    0.18051   2.491   0.0127 *  
# stthr.triafternoon             0.32721    0.16402   1.995   0.0461 *  
# hsz.std                       -0.12098    0.07965  -1.519   0.1288    
# nsk.std                        0.68245    0.09360   7.291   0.0001 ***
# tchiyr.std:stthr.trimorning    0.26334    0.20038   1.314   0.1888    
# tchiyr.std:stthr.triafternoon  0.42014    0.17366   2.419   0.0155 *  
# tchiyr.std:nsk.std             0.14446    0.11210   1.289   0.1975    
# -
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Zero-inflation model:
#             Estimate Std. Error z value Pr(>|z|)
# (Intercept)   -51.51   13502.22  -0.004    0.997
# nsk.std       -55.02   13734.07  -0.004    0.997
# save for reporting
ods.rand.zinb.disp <- round(sigma(ods.rand.zinb), 2)
ods.rand.zinb.COEF.age <-
  coef(summary(ods.rand.zinb))[[1]]["tchiyr.std",] 
ods.rand.zinb.COEF.midd.morn <-
  coef(summary(ods.rand.zinb))[[1]]["stthr.trimorning",] 
ods.rand.zinb.COEF.midd.aft <-
  coef(summary(ods.rand.zinb))[[1]]["stthr.triafternoon",] 
ods.rand.zinb.COEF.nsk <-
  coef(summary(ods.rand.zinb))[[1]]["nsk.std",] 
ods.rand.zinb.COEF.age.midd.morn <-
  coef(summary(ods.rand.zinb))[[1]]["tchiyr.std:stthr.trimorning",]
ods.rand.zinb.COEF.age.midd.aft <-
  coef(summary(ods.rand.zinb))[[1]]["tchiyr.std:stthr.triafternoon",]

# test the other two-way effects of time of day
# !! singular convergence !!
# likely due to insufficient data for the afternoon vs. morning contrast
ods.rand.zinb.v2 <- glmmTMB(round(ods_mph,0) ~
                              tchiyr.std +
                              stthr.tri.a +
                              hsz.std +
                              nsk.std +
                              tchiyr.std:stthr.tri.a +
                              tchiyr.std:nsk.std +
                              (1|aclew_child_id),
                            data=quantity.rand,
                            ziformula=~nsk.std,
                            family="nbinom1")
#ods.rand.zinb.v2.res = simulateResiduals(ods.rand.zinb.v2)
#plot(ods.rand.zinb.v2.res, rank = T) # (manually saved)
# residuals look okay and the other estimates look similar
# to those from the converging model above, so I will go with these
# results for now (hopefully more data in the future will help here)
# save for reporting
#summary(ods.rand.zinb.v2)
# stthr.tri.amorning             0.12240    0.14739   0.830               0.4063    
# tchiyr.std:stthr.tri.amorning -0.15681    0.16012  -0.979               0.3274    
ods.rand.zinb.v2.COEF.aft.morn <-
  coef(summary(ods.rand.zinb.v2))[[1]]["stthr.tri.amorning",] 
ods.rand.zinb.v2.COEF.age.aft.morn <-
  coef(summary(ods.rand.zinb.v2))[[1]]["tchiyr.std:stthr.tri.amorning",] 

## ODS tt sample ####
ods.tt.distribution <- ggplot(quantity.nonrand.tt,
                            aes(round(ods_mph,0))) +
  geom_histogram(binwidth = 2) +
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
# !! non-convergence w/ standard stthr.tri predictor, so
# used the alternative stthr.tri.o to get both two-way
# time-of-day comparisons
ods.tt.zinb <- glmmTMB(round(ods_mph,0) ~
                         tchiyr.std +
                         stthr.tri.o +
                         hsz.std +
                         nsk.std +
                         tchiyr.std:stthr.tri.o +
                         tchiyr.std:nsk.std +
                         (1|aclew_child_id),
                       data=quantity.nonrand.tt,
                       ziformula=~nsk.std,
                       family="nbinom1")
#ods.tt.zinb.res = simulateResiduals(ods.tt.zinb)
#plot(ods.tt.zinb.res, rank = T) # (manually saved)
#summary(ods.tt.zinb)
# Data: quantity.nonrand.tt
# 
#      AIC      BIC   logLik deviance df.resid 
#    365.7    392.7   -169.8    339.7       46 
# 
# Random effects:
# 
# Conditional model:
#  Groups         Name        Variance                     Std.Dev.         
#  aclew_child_id (Intercept) 0.00000000000000000000004144 0.000000000006437
# Number of obs: 59, groups:  aclew_child_id, 10
# 
# Overdispersion parameter for nbinom1 family (): 4.22 
# 
# Conditional model:
#                                  Estimate Std. Error z value   Pr(>|z|)    
# (Intercept)                      2.642703   0.164931  16.023 < 0.000001 ***
# tchiyr.std                      -0.803159   0.234119  -3.431   0.000602 ***
# stthr.tri.oafternoon            -0.608341   0.251958  -2.414   0.015759 *  
# stthr.tri.omidday               -0.002503   0.263993  -0.009   0.992434    
# hsz.std                         -0.184692   0.087090  -2.121   0.033948 *  
# nsk.std                          0.627375   0.097367   6.443   0.000001 ***
# tchiyr.std:stthr.tri.oafternoon  0.475480   0.294358   1.615   0.106243    
# tchiyr.std:stthr.tri.omidday     0.536286   0.302352   1.774   0.076111 .  
# tchiyr.std:nsk.std              -0.012219   0.134777  -0.091   0.927763    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Zero-inflation model:
#             Estimate Std. Error z value Pr(>|z|)
# (Intercept)   -31.97   11304.01  -0.003    0.998
# nsk.std       -31.33   11122.86  -0.003    0.998
# save for reporting
ods.tt.zinb.disp <- round(sigma(ods.tt.zinb), 2)
ods.tt.zinb.COEF.age <-
  coef(summary(ods.tt.zinb))[[1]]["tchiyr.std",] 
ods.tt.zinb.COEF.morn.aft <-
  coef(summary(ods.tt.zinb))[[1]]["stthr.tri.oafternoon",] 
ods.tt.zinb.COEF.morn.midd <-
  coef(summary(ods.tt.zinb))[[1]]["stthr.tri.omidday",] 
ods.tt.zinb.COEF.hsz <-
  coef(summary(ods.tt.zinb))[[1]]["hsz.std",] 
ods.tt.zinb.COEF.nsk <-
  coef(summary(ods.tt.zinb))[[1]]["nsk.std",] 
ods.tt.zinb.COEF.age.morn.aft <-
  coef(summary(ods.tt.zinb))[[1]]["tchiyr.std:stthr.tri.oafternoon",] 
ods.tt.zinb.COEF.age.morn.midd <-
  coef(summary(ods.tt.zinb))[[1]]["tchiyr.std:stthr.tri.omidday",] 

# test the other two-way effects of time of day
ods.tt.zinb.v2 <- glmmTMB(round(ods_mph,0) ~
                            tchiyr.std +
                            stthr.tri.a +
                            hsz.std +
                            nsk.std +
                            tchiyr.std:stthr.tri.a +
                            tchiyr.std:nsk.std +
                            (1|aclew_child_id),
                          data=quantity.nonrand.tt,
                          ziformula=~nsk.std,
                          family="nbinom1")
#summary(ods.tt.zinb.v2)
# save for reporting
# stthr.tri.amidday              0.60584    0.29290   2.068   0.0386 *  
# tchiyr.std:stthr.tri.amidday   0.06081    0.31055   0.196   0.8448    
ods.tt.zinb.v2.COEF.aft.midd <-
  coef(summary(ods.tt.zinb.v2))[[1]]["stthr.tri.amidday",] 
ods.tt.zinb.v2.COEF.age.aft.midd <-
  coef(summary(ods.tt.zinb.v2))[[1]]["tchiyr.std:stthr.tri.amidday",] 


## ODS random sample ####
ods.rand.gaus <- glmmTMB(log(ods_mph+1) ~
                           tchiyr.std +
                           stthr.tri +
                           hsz.std +
                           nsk.std +
                           tchiyr.std:stthr.tri +
                           tchiyr.std:nsk.std +
                           (1|aclew_child_id),
                         data=quantity.rand)
#ods.rand.gaus.res = simulateResiduals(ods.rand.gaus)
#plot(ods.rand.gaus.res, rank = T) # (manually saved)
#summary(ods.rand.gaus)
# Data: quantity.rand
# 
#      AIC      BIC   logLik deviance df.resid 
#    220.6    248.1    -99.3    198.6       79 
# 
# Random effects:
# 
# Conditional model:
#  Groups         Name        Variance        Std.Dev.  
#  aclew_child_id (Intercept) 0.0000000008073 0.00002841
#  Residual                   0.5317946982101 0.72924255
# Number of obs: 90, groups:  aclew_child_id, 10
# 
# Dispersion estimate for gaussian family (sigma^2): 0.532 
# 
# Conditional model:
#                               Estimate Std. Error z value   Pr(>|z|)    
# (Intercept)                    2.04376    0.15287  13.369 < 0.000001 ***
# tchiyr.std                    -0.26450    0.17715  -1.493   0.135415    
# stthr.trimorning               0.23286    0.21438   1.086   0.277397    
# stthr.triafternoon             0.35269    0.18999   1.856   0.063407 .  
# hsz.std                       -0.37882    0.11236  -3.371   0.000748 ***
# nsk.std                        1.55757    0.09555  16.302 < 0.000001 ***
# tchiyr.std:stthr.trimorning    0.07095    0.22658   0.313   0.754176    
# tchiyr.std:stthr.triafternoon  0.42755    0.20507   2.085   0.037084 *  
# tchiyr.std:nsk.std             0.17956    0.11344   1.583   0.113435    
# DIFFERENCES between z-i nb and logged gaussian?
# some non-sig <-> sig (age, morn-vs-midd, hsz), but all going in the same direction

# test the other two-way effects of time of day
ods.rand.gaus.v2 <- glmmTMB(log(ods_mph+1) ~
                              tchiyr.std +
                              stthr.tri.a +
                              hsz.std +
                              nsk.std +
                              tchiyr.std:stthr.tri.a +
                              tchiyr.std:nsk.std +
                              (1|aclew_child_id),
                            data=quantity.rand)
#summary(ods.rand.gaus.v2)
# stthr.tri.amorning            -0.11983    0.18738  -0.640   0.522494    
# tchiyr.std:stthr.tri.amorning -0.35659    0.19600  -1.819   0.068861 .  
# save for reporting
# DIFFERENCES between z-i nb and logged gaussian?
# broadly similar

## ODS tt sample ####
ods.tt.gaus <- glmmTMB(log(ods_mph+1) ~
                         tchiyr.std +
                         stthr.tri +
                         hsz.std +
                         nsk.std +
                         tchiyr.std:stthr.tri +
                         tchiyr.std:nsk.std +
                         (1|aclew_child_id),
                       data=quantity.nonrand.tt)
#ods.tt.gaus.res = simulateResiduals(ods.tt.gaus)
#plot(ods.tt.gaus.res, rank = T) # (manually saved)
#summary(ods.tt.gaus)
# Random effects:
# 
# Conditional model:
#  Groups         Name        Variance        Std.Dev.  
#  aclew_child_id (Intercept) 0.0000000001595 0.00001263
#  Residual                   0.4803202867131 0.69305143
# Number of obs: 59, groups:  aclew_child_id, 10
# 
# Dispersion estimate for gaussian family (sigma^2): 0.48 
# 
# Conditional model:
#                               Estimate Std. Error z value   Pr(>|z|)    
# (Intercept)                    2.34259    0.21441  10.926 < 0.00001 ***
# tchiyr.std                    -0.14993    0.20131  -0.745   0.45642    
# stthr.trimorning              -0.24611    0.28518  -0.863   0.38813    
# stthr.triafternoon            -0.71950    0.25935  -2.774   0.00553 ** 
# hsz.std                       -0.26907    0.11809  -2.278   0.02270 *  
# nsk.std                        1.09149    0.10899  10.015 < 0.00001 ***
# tchiyr.std:stthr.trimorning   -0.24622    0.28606  -0.861   0.38940    
# tchiyr.std:stthr.triafternoon  0.05821    0.25601   0.227   0.82013    
# tchiyr.std:nsk.std            -0.08266    0.13667  -0.605   0.54531  
# DIFFERENCES between z-i nb and logged gaussian?
# some non-sig <-> sig (age), but all going in the same direction

# test the other two-way effects of time of day
ods.tt.gaus.v2 <- glmmTMB(log(ods_mph+1) ~
                            tchiyr.std +
                            stthr.tri.a +
                            hsz.std +
                            nsk.std +
                            tchiyr.std:stthr.tri.a +
                            tchiyr.std:nsk.std +
                            (1|aclew_child_id),
                          data=quantity.nonrand.tt)
#summary(ods.tt.gaus.v2)
# stthr.tri.amorning             0.47338    0.24385   1.941              0.05222 .  
# tchiyr.std:stthr.tri.amorning -0.30443    0.26895  -1.132              0.25767    
# save for reporting
# DIFFERENCES between z-i nb and logged gaussian?
# broadly similar

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


all.models <- bind_rows(TCDS.models, ODS.models)
write_csv(all.models, paste0(shiny.input.path, "all_model_tables.csv"))

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

#voc.mat.by.age
