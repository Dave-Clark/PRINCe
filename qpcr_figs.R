library(data.table)
library(ggplot2)
library(patchwork)

dat <- fread("mergedqPCRData.csv")
dat <- dat[!site %in% c("NEG1", "NEG2", "NEG3")]

# amoA ratio
aoaRatio <- ggplot(dat, aes(x = site, col = season, y = AOA_amoA/AOB_amoA)) +
  geom_hline(yintercept = 1, linetype = 2, col = "grey") +
  geom_boxplot() +
  labs(x = "Site", y = expression(AOA~amoA:AOB~amoA~ratio)) +
  theme_bw() +
  facet_wrap(~site, scales = "free") +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_blank(),
    strip.text.x = element_text(size = 14),
    panel.grid = element_blank())

ggsave("../../graphics/amoA_ratio.pdf", aoaRatio, height = 8, width = 12,
  device = "pdf")

ammoniaOx <- melt(
  dat[, .SD, .SDcols = c("site", "season", "AOA_amoA", "AOB_amoA")],
  id.vars = c("site", "season"))

amoaPlot <- ggplot(ammoniaOx, aes(x = season, y = value, fill = variable)) +
  geom_boxplot() +
  facet_wrap(~site, scales = "free_y") +
  labs(x = "Season", y = expression(Gene~copies~(g^-1))) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = NULL,
    strip.text.x = element_text(size = 14),
    panel.grid = element_blank())

ggsave("../../graphics/ammonia_oxidisers.pdf", amoaPlot, height = 8,
  width = 12, device = "pdf")

# look at correlation between anammox genes
ggplot(dat,
  aes(x = Anammox_hzo, y = Anammox_hzsA, col = site, shape = season)) +
  geom_point()

anammox <- melt(
  dat[, .SD, .SDcols = c("site", "season", "Anammox_hzo", "Anammox_hzsA")],
  id.vars = c("site", "season"))

anammoxPlot <- ggplot(anammox, aes(x = site, y = value, col = season)) +
  geom_boxplot() +
  facet_wrap(~variable, scales = "free_y") +
  labs(x = "Site", y = expression(Anammox~functional~gene~abundance~(g^-1))) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_blank(),
    strip.text.x = element_text(size = 14),
    panel.grid = element_blank())

ggsave("../../graphics/anammox_abundance.pdf", anammoxPlot, height = 5, width = 12, device = "pdf")


# nirK dynamics
nirK <- melt(
  dat[, .SD, .SDcols = c("site", "season", "nirK_C1", "nirK_C2", "nirK_C3")],
  id.vars = c("site", "season"))

nirS <- melt(
    dat[, .SD, .SDcols = c("site", "season", "nirS_C1", "nirS_C2")],
    id.vars = c("site", "season"))

nirKPlot <- ggplot(nirK, aes(x = season, fill = variable, y = value)) +
  geom_boxplot() +
  facet_wrap(~site) +
  labs(x = "Site", y = expression(nirK~gene~abundance~(g^-1))) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_blank(),
    strip.text.x = element_text(size = 14),
    panel.grid = element_blank())

ggsave("../../graphics/nirK_abundance.pdf", nirKPlot, height = 8, width = 10)

nirSPlot <- ggplot(nirS, aes(x = season, fill = variable, y = value)) +
  geom_boxplot() +
  facet_wrap(~site) +
  labs(x = "Site", y = expression(nirK~gene~abundance~(g^-1))) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_blank(),
    strip.text.x = element_text(size = 14),
    panel.grid = element_blank())

ggsave("../../graphics/nirS_abundance.pdf", nirSPlot, height = 8, width = 10)

nxrBPlot <- ggplot(dat, aes(x = site, col = season, y = Nitrospira_nxrB)) +
  geom_boxplot() +
  labs(x = "Site", y = expression(Nitrospira~nxrB~gene~abundance~(g^-1))) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_blank(),
    strip.text.x = element_text(size = 14),
    panel.grid = element_blank())

ggsave("../../nxrB_abundance.pdf", nxrBPlot, height = 5, width = 8,
  device = "pdf")

ureCPlot <- ggplot(dat, aes(x = site, fill = season, y = Thaum_ureC)) +
    geom_boxplot() +
    scale_y_log10() +
    labs(x = "Site",
      y = expression(atop(Thaumarchaeal~italic(ureC),  gene~copies~(g^-1)))) +
    theme_bw() +
    theme(axis.text = element_text(size = 18),
      axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.title = element_text(size = 20),
      legend.text = element_text(size = 16),
      legend.position = c(0.8, 0.9),
      legend.title = element_blank(),
      panel.grid = element_blank(),
      aspect.ratio = 1,
      plot.title = element_text(size = 20))

 ureCProp <- ggplot(dat,
      aes(x = site, fill = season, y = Thaum_ureC/AOA_amoA)) +
    geom_boxplot() +
    labs(x = "Site", y = expression(Archaeal~italic(ureC):italic(amoA))) +
    theme_bw() +
    theme(axis.text = element_text(size = 18),
      axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.title = element_text(size = 20),
      legend.position = "none",
      panel.grid = element_blank(),
      aspect.ratio = 1,
      plot.title = element_text(size = 20))

ureCPanel <- ureCPlot + ureCProp + plot_layout(ncol = 2)

ggsave("/media/dave/storage/NERC_fellowship/ureC_plot.pdf", ureCPanel,
  height = 4, width = 9, device = "pdf")


ggsave("/media/storage/Dropbox/cv/QMUL_lectureship_19_06/Prelim_ureC.pdf",
  ureCPlot, height = 5, width = 5.5, device = "pdf")

ureCRatio <- ggplot(dat, aes(x = site, col = season, y = AOA_amoA/Thaum_ureC)) +
  geom_boxplot() +
  labs(x = "Site", y = expression(AOA~amoA:Thaumarchaeal~ureC~ratio)) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_blank(),
    strip.text.x = element_text(size = 14),
    panel.grid = element_blank())

    ggplot(dat, aes(x = site, col = season, y = Thaum_ureC/AOA_amoA)) +
      geom_boxplot() +
      labs(x = "Site", y = expression(AOA~amoA:Thaumarchaeal~ureC~ratio)) +
      theme_bw() +
      theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        strip.text.x = element_text(size = 14),
        panel.grid = element_blank())


bac16S <- ggplot(dat, aes(x = site, col = season, y = Bac_16S)) +
    geom_boxplot() +
    labs(x = "Site", y = "Bacterial 16S rRNA copies (g-1)") +
    theme_bw() +
    theme(axis.text = element_text(size = 16),
      axis.title = element_text(size = 18),
      legend.text = element_text(size = 14),
      legend.title = element_blank(),
      strip.text.x = element_text(size = 14),
      panel.grid = element_blank())

  ggsave("../../graphics/Bac_16S_abundance.pdf", bac16S, height = 5, width = 8,
    device = "pdf")
