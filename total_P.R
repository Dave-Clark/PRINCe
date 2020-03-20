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

totalP[, probe := as.integer(probe)]
totalP[, depth := ifelse(probe <= 8, "5cm", "10cm")]

totalP <- totalP[depth == "5cm"]

totalP[treatment == "Ammonium", microbialCore := ifelse(
  probe %in% c(1, 2, 3), 1, ifelse(
    probe %in% c(4, 5, 6), 2, 3))]

totalP[treatment == "Nitrite", microbialCore := ifelse(
  probe %in% c(1, 2, 3), 4, ifelse(
    probe %in% c(4, 5, 6), 5, 6))]

# calculate mean P by site, season, and microb core
microP <- totalP[, mean(conc), by = c(
  "site", "season", "treatment", "microbialCore")]

names(microP)[5] <- "total_P"

p_5cm <- ggplot(microP, aes(x = site, y = total_P, col = season)) + 
  geom_boxplot(width = 0.6) +
  labs(x = "Site", y = expression(Total~P~(mg~L^-1))) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    panel.grid = element_blank())

ggsave("../graphics/p_5cm.pdf", p_5cm, height = 5, width = 7, device = "pdf")

qpcr <- fread("qPCR_raw/mergedqPCRData.csv")

qpcr[, label := paste(site, core, season, sep = "_")]

microP[, label := paste(site, microbialCore, season, sep = "_")]

usefulCols <- c("AOA_amoA", "AOB_amoA", "Bac_16S", "Anammox_hzo", "Anammox_hzsA", "nirK_C1", "nirK_C2", "nirK_C3", "nirS_C1", "nirS_C2", "nosZ_C1", "Nitrospira_nxrB", "Thaum_ureC", "label")

allDat <- merge(qpcr[, .SD, .SDcols = usefulCols], microP, by = "label")

meltDat <- melt(allDat[, !c("label", "treatment", "microbialCore")],
  id.vars = c("site", "season", "total_P"))

qpcrP <- ggplot(meltDat,
    aes(x = total_P, y = value, col = site, shape = season)) +
  geom_point(size = 2, alpha = 0.6) +
  facet_wrap(~ variable, scales = "free_y") +
  labs(x = expression(Total~P~(mg~L^-1)),
    y = expression(Gene~copies~(g^-1~dry~sediment))) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    strip.text.x = element_text(size = 14),
    panel.grid = element_blank())

ggsave("../graphics/qpcr_panel.pdf", qpcrP, height = 10, width = 13, device = "pdf")

glmer.nb()
