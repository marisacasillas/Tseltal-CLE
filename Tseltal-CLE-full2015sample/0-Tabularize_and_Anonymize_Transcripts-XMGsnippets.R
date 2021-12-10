# Use this file to anonymize and tabularize transcript data
library(tidyverse)
library(phonfieldwork)

raw.data.path <- "transcripts/raw/eaf/" # text files exported from ELAN

files <- list.files(path = raw.data.path, pattern = "*.eaf")
all.utterances <- data.frame()

for (file in files) {
  file.contents <- eaf_to_df(paste0(raw.data.path, file))
  if (!(is.null(file.contents))) {
    spch.tbl <- file.contents %>%
      filter(tier_name %in% c("CHI", "OTR")) %>%
      select(tier_name, time_start, time_end) %>%
      mutate(
        dur.sec = time_end - time_start,
        recording = str_extract(file,
                                "M-F\\d{2}-C\\d{2}(-rec\\d)?"),
        clip.onset = str_extract_all(file, "\\d{6}")[[1]][1],
        tod.clip.onset = str_extract_all(file, "\\d{6}")[[1]][2]
      ) %>%
      rename("tier.name" = tier_name,
             "time.start" = time_start,
             "time.end" = time_end)
    if ("OTR" %in% spch.tbl$tier.name) {
      quien.tbl <- file.contents %>%
        filter(tier_name %in% "quien@OTR") %>%
        select(content, time_start, time_end) %>%
        mutate(
          recording = str_extract(file,
                                  "M-F\\d{2}-C\\d{2}(-rec\\d)?"),
          clip.onset = str_extract_all(file, "\\d{6}")[[1]][1],
        ) %>%
        rename("time.start" = time_start,
               "time.end" = time_end,
               "spkr.type" = content) %>%
        mutate(tier.name = "OTR")
      spch.tbl <- spch.tbl %>%
        left_join(quien.tbl,
                  by = c("tier.name", "time.start", "time.end",
                         "recording", "clip.onset")) %>%
        mutate(spkr.type = case_when(
          !is.na(spkr.type) ~ spkr.type,
          is.na(spkr.type) & tier.name == "OTR" ~ "missing",
          tier.name == "CHI" ~ "CHI",
          TRUE ~ "uh oh"
        ))
    } else {
      spch.tbl <- spch.tbl %>%
        mutate(spkr.type = case_when(
          tier.name == "CHI" ~ "CHI",
          TRUE ~ "uh oh"
          ))
    }
    all.utterances <- bind_rows(all.utterances, spch.tbl)
  }
}

all.utterances <- all.utterances %>%
  mutate(spkr.type.corr = case_when(
    tier.name == "CHI" ~ "CHI",
    grepl("(hombre)|(señor)", spkr.type) ~ "hombre",
    grepl("([mM]uj+er)|(muchacha)", spkr.type) ~ "mujer",
    grepl("niña", spkr.type) ~ "niña",
    grepl("n[iu]ño", spkr.type) ~ "niño",
    TRUE ~ "other/unknown"
  ))

write_csv(all.utterances, "all.CHI.OTR.utterances.csv")

# Write out list of total annotated clips
all.annotated.clips <- tibble(
  recording = str_extract(files,
                          "M-F\\d{2}-C\\d{2}(-rec\\d)?"),
  filename = files,
  clip.onset = unlist(map(str_extract_all(files, "\\d{6}"), 1)),
  segment.num = 0
)
for (rec in unique(all.annotated.clips$recording)) {
  rec.clips <- which(all.annotated.clips$recording == rec)
  n <- 1
  for (idx in rec.clips) {
    all.annotated.clips$segment.num[idx] <- n
    n <- n + 1
  }
}

write_csv(all.annotated.clips, "all.annotated.clips.csv")
