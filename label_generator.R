siteNames <- c(paste0("C", 1:10), paste0("S", 1:5))

siteNames <- rep(siteNames, each = 6)

replicate <- rep(1:3, each = 2, times = 15)

depth <- rep(c("5cm", "10cm"), times = 45)

labels <- paste(siteNames, replicate, depth, sep = "-")

labels <- paste0(rep(labels, each = 2), "_",
  rep(c("DNA+Nut", "RNA"), times = 45))

write.table(labels, "labels.txt", sep = "\n", row.names = F, quote = F)
