library(data.table)
library(ggplot2)
library(patchwork)

dat <- fread("reactor_qPCR.csv")

genes <- names(dat)[10:19]

# normalise counts to 1g wet weight
dat[, (genes) := lapply(.SD, function(x)
  (1/wet) * 100 * x), .SDcols = genes]

dat[, experiment := paste0("Experiment ", experiment)]

arch <- ggplot(na.omit(dat), aes(x = days, y = arch_16S, col = treatment)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_y_continuous(labels = function(x) format(x, scientific=T)) +
  facet_wrap(~ experiment) +
  labs(x = "Days after P addition",
    y = expression(Gene~copies~g^-1~sediment), col = NULL) +
  ggtitle("Archaeal 16S rRNA") +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text.x = element_text(size = 14),
    panel.grid = element_blank())

bac <- ggplot(na.omit(dat), aes(x = days, y = bac_16S, col = treatment)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_y_continuous(labels = function(x) format(x, scientific=T)) +
  facet_wrap(~ experiment) +
  labs(x = "Days after P addition",
    y = expression(Gene~copies~g^-1~sediment), col = NULL) +
  ggtitle("Bacterial 16S rRNA") +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text.x = element_text(size = 14),
    panel.grid = element_blank())

phylo <- arch / bac

ggsave("../../graphics/PRINCe_reactors_16S_qPCR.pdf", phylo, height = 8,
  width = 9, device = "pdf")

AOA <- ggplot(na.omit(dat), aes(x = days, y = AOA, col = treatment)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_y_continuous(labels = function(x) format(x, scientific=T)) +
  facet_wrap(~ experiment) +
  labs(x = "Days after P addition",
    y = expression(Gene~copies~g^-1~sediment), col = NULL) +
  ggtitle(expression(Archaeal~italic(amoA))) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text.x = element_text(size = 14),
    panel.grid = element_blank())

AOB <- ggplot(na.omit(dat), aes(x = days, y = AOB, col = treatment)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_y_continuous(labels = function(x) format(x, scientific=T)) +
  facet_wrap(~ experiment) +
  labs(x = "Days after P addition",
    y = expression(Gene~copies~g^-1~sediment), col = NULL) +
  ggtitle(expression(Bacterial~italic(amoA))) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text.x = element_text(size = 14),
    panel.grid = element_blank())

comammox <- ggplot(na.omit(dat),
    aes(x = days, y = comammox_amoA, col = treatment)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_y_continuous(labels = function(x) format(x, scientific=T)) +
  facet_wrap(~ experiment) +
  labs(x = "Days after P addition",
    y = expression(Gene~copies~g^-1~sediment), col = NULL) +
  ggtitle(expression(Comammox~bacterial~italic(amoA))) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text.x = element_text(size = 14),
    panel.grid = element_blank())

amoaPlot <- AOA / AOB / comammox

ggsave("../../graphics/PRINCe_reactors_amoA_qPCR.pdf", amoaPlot, height = 11,
  width = 9, device = "pdf")

hzo <- ggplot(na.omit(dat), aes(x = days, y = hzo, col = treatment)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_y_continuous(labels = function(x) format(x, scientific=T)) +
  facet_wrap(~ experiment) +
  labs(x = "Days after P addition",
    y = expression(Gene~copies~g^-1~sediment), col = NULL) +
  ggtitle(expression(Anammox~bacterial~italic(hzo))) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text.x = element_text(size = 14),
    panel.grid = element_blank())

nxrb <- ggplot(na.omit(dat), aes(x = days, y = nxrB, col = treatment)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_y_continuous(labels = function(x) format(x, scientific=T)) +
  facet_wrap(~ experiment) +
  labs(x = "Days after P addition",
    y = expression(Gene~copies~g^-1~sediment), col = NULL) +
  ggtitle(expression(italic(Nitrospira~nxrB))) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text.x = element_text(size = 14),
    panel.grid = element_blank())

nitrif <- hzo / nxrb

ggsave("../../graphics/PRINCe_reactors_nxrb_hzo_qPCR.pdf", nitrif, height = 8,
  width = 9, device = "pdf")

nirk <- ggplot(na.omit(dat), aes(x = days, y = nirK_C3, col = treatment)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_y_continuous(labels = function(x) format(x, scientific=T)) +
  facet_wrap(~ experiment) +
  labs(x = "Days after P addition",
    y = expression(Gene~copies~g^-1~sediment), col = NULL) +
  ggtitle(expression(italic(nirK)~Clade~III)) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text.x = element_text(size = 14),
    panel.grid = element_blank())

nirs <- ggplot(na.omit(dat), aes(x = days, y = nirS_C1, col = treatment)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_y_continuous(labels = function(x) format(x, scientific=T)) +
  facet_wrap(~ experiment) +
  labs(x = "Days after P addition",
    y = expression(Gene~copies~g^-1~sediment), col = NULL) +
  ggtitle(expression(italic(nirS)~Clade~I)) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text.x = element_text(size = 14),
    panel.grid = element_blank())

nosz <- ggplot(na.omit(dat),
    aes(x = days, y = nosZ_C1, col = treatment)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_y_continuous(labels = function(x) format(x, scientific=T)) +
  facet_wrap(~ experiment) +
  labs(x = "Days after P addition",
    y = expression(Gene~copies~g^-1~sediment), col = NULL) +
  ggtitle(expression(italic(nosZ)~Clade~I)) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text.x = element_text(size = 14),
    panel.grid = element_blank())

denit <- nirs / nirk / nosz

ggsave("../../graphics/PRINCe_reactors_denitrifiers_qPCR.pdf", denit,
  height = 11, width = 9, device = "pdf")
