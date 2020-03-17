# randomise order of DNA extractions for reactor exp
dat <- read.csv("sampleCodes.csv")

# randomly shuffle rows
dat <- dat[sample(nrow(dat)), ]

# add tube number/index column
dat$tube <- 1:24

write.csv(dat, "randomisedDnaExtracts.csv")
