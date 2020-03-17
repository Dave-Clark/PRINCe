### check dna extract labels for any duplicates
library(data.table)

dna <- fread("PRINCE_dna_extracts_raw.csv")

# coerce labels to upper case for consistency
dna[, sample := toupper(sample)]

# extract site name, core number, and depth
dna[, ":="(site = sapply(strsplit(sample, "-"), "[[", 1),
  core = as.integer(sapply(strsplit(dna$sample, "-"), function(x)
    tryCatch(x[[2]],
      warning = function(w) return(NA),
      error = function(e) return(NA)))),
  depth = sapply(strsplit(dna$sample, "-"), function(x)
    tryCatch(x[[3]],
      warning = function(w) return(NA),
      error = function(e) return(NA))),
  season = sapply(strsplit(dna$sample, "-"), function(x)
    tryCatch(x[[4]],
      warning = function(w) return(NA),
      error = function(e) return(NA)))
)]

# reformat season label
dna[, labelType := ifelse(grepl("PL", season), "paper", "original")]
dna[, season := ifelse(grepl("S", season), "Summer", "Winter")]

# need to check missing samples
# C7 summer samples are actually C8.
# S5 didn't get done in winter, so can be dropped.
# C5-5-W needs to be replaced by C5-7-W due to coring through wood

# update 24/01/19
# relabelled sample, recorded in lab book
# relabelled google sheet
# looks good now
sum(dna[use == "yes", table(site, season)])
# 147 prince samples including negs
# plus 3 chalk samples for prelim data (qPCR only)
# = 150 samples in total

# subset out useful samples
useDna <- dna[use == "yes"]

# create vector of well positions
wellPositions <- paste(rep(letters[1:8], each = 12), rep(1:12, times = 8),
  sep = "-")

# concatenate so that enough wells are available for all samples
wellPositions <- rep(wellPositions, times = 2)[1:nrow(useDna)]

# assign each sample a random well position
set.seed(001)
useDna <- useDna[sample(nrow(useDna)), ]

useDna[, well := wellPositions]

# assign plate to duplicated well positions
useDna[, plate := 1:.N, by = well]

fwrite(useDna, "randomised_plate_positions.csv", row.names = F)
