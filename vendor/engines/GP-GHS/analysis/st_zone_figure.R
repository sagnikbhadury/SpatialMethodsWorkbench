## Figure of the spatial transcriptomics slide and its six pathologist zones.
## Left: the H&E image stored with the data. Right: the same image with every
## analysed spot coloured by the zone that the pathologist assigned to it.
## Nothing is estimated here; the zones are read from the data file.
##
## Usage (R 4.4.1): Rscript st_zone_figure.R   (after st_prepare.R)

suppressPackageStartupMessages({ library(hdf5r); library(ggplot2); library(gridExtra); library(grid) })

DATA_DIR <- Sys.getenv("ST_DATA_DIR", file.path(
  "C:/Users/bhadury/University of Michigan Dropbox/Sagnik Bhadury",
  "Data for Analysis/CRC_CMS_Visium_cellxgene"))
OUT_DIR <- Sys.getenv("ST_OUT_DIR", "st_results")
prep <- readRDS(file.path(OUT_DIR, "st_prepared.rds")); obs <- prep$obs; zt <- prep$zone_tab
zone_cols <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#4a3aa7")[seq_len(nrow(zt))]

f <- H5File$new(file.path(DATA_DIR, paste0(prep$slide, ".h5ad")), "r")
lib <- setdiff(names(f[["uns/spatial"]]), "is_single")[1]
g <- f[["uns/spatial"]][[lib]]
res <- if ("hires" %in% names(g[["images"]])) "hires" else names(g[["images"]])[1]
img <- g[["images"]][[res]]$read()
sf_name <- grep(res, names(g[["scalefactors"]]), value = TRUE)[1]
sf <- if (is.na(sf_name)) 1 else g[["scalefactors"]][[sf_name]]$read()
cat("image:", res, "| array dims as read:", dim(img), "| scale factor:", sf_name, sf, "\n")
f$close_all()
## hdf5r returns the array with reversed dimensions (channel, column, row)
if (dim(img)[1] %in% 3:4) img <- aperm(img, c(3, 2, 1))
if (max(img) > 1) img <- img / 255
img <- img[, , 1:3]
H <- dim(img)[1]; W <- dim(img)[2]
obs$px <- obs$x * sf; obs$py <- obs$y * sf
pad <- 0.04 * max(diff(range(obs$px)), diff(range(obs$py)))
xl <- c(max(1, min(obs$px) - pad), min(W, max(obs$px) + pad))
yl <- c(max(1, min(obs$py) - pad), min(H, max(obs$py) + pad))
ras <- rasterGrob(as.raster(img), interpolate = TRUE)
base <- theme_void(base_size = 10) +
  theme(plot.title = element_text(size = 10.5, face = "bold"), plot.subtitle = element_text(size = 9),
        legend.position = "bottom", plot.margin = margin(4, 6, 4, 6))
he <- function() annotation_custom(ras, xmin = 0.5, xmax = W + 0.5, ymin = -(H + 0.5), ymax = -0.5)
p1 <- ggplot(obs, aes(px, -py)) + he() + geom_blank() +
  coord_equal(xlim = xl, ylim = -rev(yl), expand = FALSE) +
  labs(title = "(a) Hematoxylin and eosin image",
       subtitle = "One Visium slide of a rectal adenocarcinoma") + base
obs$zone_lab <- factor(sprintf("%s (%d spots)", zt$name[obs$zone], zt$n_spots[obs$zone]),
                       sprintf("%s (%d spots)", zt$name, zt$n_spots))
p2 <- ggplot(obs, aes(px, -py, color = zone_lab)) + he() +
  geom_point(size = 1.05, alpha = 0.95) +
  scale_color_manual(values = setNames(zone_cols, levels(obs$zone_lab)), name = NULL) +
  coord_equal(xlim = xl, ylim = -rev(yl), expand = FALSE) +
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 2.6))) +
  labs(title = "(b) Zones annotated by the pathologist",
       subtitle = sprintf("%s spots under tissue, each assigned to one of %d zones",
                          format(nrow(obs), big.mark = ","), nrow(zt))) + base
## the left panel gets a bottom margin as tall as the legend so that the two images align
ggsave(file.path("manuscript", "Image_s", "st_zones.png"),
       arrangeGrob(p1 + theme(plot.margin = margin(4, 6, 58, 6)), p2, nrow = 1),
       width = 11, height = 5.6, dpi = 200)
cat("figure written: manuscript/Image_s/st_zones.png | image", W, "x", H, "\n")
