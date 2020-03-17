library(data.table)
library(stringr)
library(ggplot2)
library(gtools)

dat <- fread("PRINCE_porewater_total_P_raw.csv")

dat[, sample := toupper(sample)]

# extract data from sample label
dat[sample != "STD", ":="(
  site = sapply(strsplit(sample, "-"), "[[", 1),
  probe = sapply(str_extract_all(sample, "\\d+"), "[", 2),
  treatment = ifelse(grepl("N-", sample), "Nitrite", "Ammonium"),
  season = ifelse(grepl("-S", sample), "Summer", "Winter")
)]

# visualise all standard curves
ggplot(dat[sample == "STD",], aes(
    x = conc, y = absorbance, col = as.factor(batch))) +
  geom_point() +
  facet_wrap(~ batch)

# drop smallest standard as this screws up linearity and standard 7
stds <- dat[sample == "STD" & conc != 0.06125]

ggplot(dat[sample == "STD",], aes(
    x = absorbance, y = conc, col = as.factor(batch))) +
  geom_point()
stds[, batch := as.factor(batch)]
stdMod <- lm(conc ~ absorbance * batch, data = stds)

dat$batch <- as.factor(dat$batch)

# predict P from std curve
dat[sample != "STD", conc := predict(stdMod, dat[sample != "STD"])]

# exclude standards
totalP <- dat[!is.na(site) & site != "STS"]

# multiply by dilution factor
totalP[, conc := conc * dilution]

# remove outliers
totalP[conc <= 0, conc := NA]

fwrite(totalP, "PRINCe_total_P.csv")

# read in location data to merge site names
coords <- fread("coords.csv")

# merge with P data
totalP[coords, siteName := i.site, on = c(site = "siteCode")]

# make site labels
totalP[, label := paste0(site, " - ", siteName)]

# reorder sites
totalP[, label := factor(label,
  levels = mixedsort(unique(label)))]

phosPlot <- ggplot(totalP, aes(x = label, y = conc, col = season)) +
  geom_boxplot() +
  labs(x = "Site", y = expression(Total~phosphorous~(mg~L^1)),
    col = "Season") +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    panel.grid = element_blank())

ggsave("../graphics/total_phosphorous.pdf", phosPlot, width = 6, height = 5,
  device = "pdf")
