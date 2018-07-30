# Use this file to anonymize and check transcript data
library(tidyverse)

raw.data.path <- "transcripts/raw/" # text files exported from ELAN
anon.data.path <- "transcripts/anon/" # text files to use as input to analysis

files <- list.files(path=raw.data.path,pattern="*.txt")
all.data <- data.frame()
for (i in 1:length(files)) {
  print(files[i])
  newfile <- read.table(paste0(raw.data.path, files[i]), quote = "",
                        sep="\t", header=FALSE, stringsAsFactors = FALSE)
  names(newfile)[1:6] <- c("tier", "speaker", "start", "stop", "dur", "val")
  transcription.tiers <- which(newfile$tier == newfile$speaker)
  newfile$val[transcription.tiers] <- "0."
  uniq.tiers <- unique(newfile$tier)
  if ("code" %in% uniq.tiers) {
    newfile <- newfile %>% filter(!(str_detect(tier, "AudioOnly|Notes|context|code_num|on_off")))
    write_csv(newfile, paste0(anon.data.path, files[i]))
  } else {
    print("Tiers missing!")
  }
}