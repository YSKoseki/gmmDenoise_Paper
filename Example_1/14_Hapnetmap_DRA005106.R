# 14_Hapnetmap_DRA005106.R
# Last updated on 2023.5.26 by YK
# An R script to create haplotype networks and haplotype frequency maps
# R 4.1.2
# 
# Packages required
library(treeio); packageVersion("treeio") # ‘1.18.1’
library(speedyseq); packageVersion("speedyseq") # ‘0.5.3.9018’
library(DECIPHER); packageVersion("DECIPHER") # ‘2.22.0’
library(pegas); packageVersion("pegas") # ‘1.1’
library(wesanderson); packageVersion("wesanderson") # ‘0.3.6’
library(svglite); packageVersion("svglite") # 2.1.0
library(sf); packageVersion("sf") # ‘1.0.6’
library(rmapshaper); packageVersion("rmapshaper") # ‘0.4.6’
library(ggspatial); packageVersion("ggspatial") # ‘1.1.7’
library(cowplot); packageVersion("cowplot") # ‘1.1.1’
library(tidyverse); packageVersion("tidyverse") # ‘1.3.1’

# Paths
## List of phyloseq objects
path_node <- "./13-DNfocal_DRA005106/01-Saved_object/focalnode.obj"
path_phylo <- "./13-DNfocal_DRA005106/01-Saved_object/phylo_pre4.obj"
path_shape <- "./14-shapefiles"

## Output directory
path_output <- "./14-Hapnetmap_DRA005106.R"
dir.create(path_output, recursive=TRUE)

# Load the phyloseq objects
focalnode <- readRDS(path_node)
phylo_pre4 <- readRDS(path_phylo)

# Haplotype networks
## Create objects for plotting
get.hapnet <- function(x, y) {
  label <- x %>% phy_tree() %>%
    treeio::tree_subset(node = y, levels_back = 0) %>%
    `[[`("tip.label")
  phyfocal <- x %>% prune_taxa(label, .) %>% # The selected species only
    filter_tax_table(denoised == "no") %>% # Retained ASVs only
    prune_samples(sample_sums(.) > 0, .) # Non-zero samples only
  taxa_names(phyfocal) <- phyfocal %>% tax_table() %>%
    `[`(, "glob") %>% as.vector()
  nsamp <- phyfocal %>% otu_table() %>%
    apply(2, function(x) {which(x > 0) %>% length()}) 
  reads <- phyfocal %>% otu_table() %>% colSums()
  dnabin <- phyfocal %>% refseq() %>%
    DECIPHER::AlignSeqs(processors = NULL, verbose = FALSE) %>%
    DECIPHER::StaggerAlignment(processors = NULL, verbose = FALSE) %>%
    rep(times = reads) %>% as.DNAbin()
  hap <- pegas::haplotype(dnabin, strict = TRUE,
                          labels = names(reads))
  net <- pegas::haploNet(hap)
  # See if I haven't made any mistakes about labeling
  read_net <- attr(net, "freq") %>% as.numeric()
  names(read_net) <- attr(net, "labels")
  if (identical(reads, read_net)) {
    cat("The ASV read counts are arranged consistent with the ASV labels.\n")
  } else {
    stop("The ASV read counts are NOT arranged consistent with the ASV labels!\n")
  }
  if (identical(names(nsamp), names(read_net))) {
    cat("The ASV-detected sample sizes are arranged consistent with the ASV labels.\n")
  } else {
    stop("The ASV-detected sample sizes are NOT arranged consistent with the ASV labels!\n")
  }
  z <- list(net = net, reads = reads, nsamp = nsamp)
  return(z)
}

hapnet <- get.hapnet(phylo_pre4[["dada"]], focalnode[["dada"]])

## Test plot
replot.hapnet <- function(hapnet, size = "reads", labels = TRUE,
                          scale.ratio = 100, cex = 1, lwd = 1,
                          col = NA, bg = NA, alpha = 1) {
  net <- hapnet[["net"]]
  if (size == "reads") {N <- hapnet[["reads"]]}
  if (size == "sqrtreads") {N <- sqrt(hapnet[["reads"]])}
  if (size == "logreads") {N <- log10(hapnet[["reads"]])}
  if (size == "nsamp") {N <- hapnet[["nsamp"]]}
  xy <- hapnet[["xy"]]
  nhap <- length(N)
  if (is.na(bg)) {
    if (nhap <= 5) {
      par_bg <- wes_palette("Darjeeling1", n = nhap) %>%
        alpha(alpha)
    }
    else {
      if (nhap <= 10) {
        par_bg <- c(wes_palette("Darjeeling1"),
                    wes_palette("Darjeeling2", n = nhap - 5)) %>%
          alpha(alpha)
      }
      else {
        par_bg <- c(wes_palette("Darjeeling1"),
                    wes_palette("Darjeeling2"),
                    wes_palette("Cavalcanti1", n = nhap - 10)) %>%
          alpha(alpha)
      }
    }
  }
  else {
    par_bg <- bg
  }
  if (is.na(col)) {
    par_col <- "black"
  }
  else {
    par_col <- col
  }
  plot(net, xy = xy, labels = labels, scale.ratio = scale.ratio, cex = cex,
       lwd = lwd, size = N, threshold = 0, show.mutation = 1, legend = FALSE,
       col = par_col, bg = par_bg)
  xy <- replot()
  return(xy)
}
tmp <- replot.hapnet(hapnet,
                     size = "sqrtreads",
                     scale.ratio = 3, cex = 1.5, alpha = .7, lwd = 3)
##Notrun hapnet[["dada"]][["xy"]] <- tmp
rm(tmp)
## Final plot
svglite(paste0(path_output, "/02-Hapnet.svg"),
        width = 5, height = 3)
replot.hapnet(hapnet,
              size = "sqrtreads",
              scale.ratio = 10, cex = .8, alpha = .7, lwd = 1)
dev.off()


# Geographical maps
## Lake Biwa obj
lake_biwa <- read_sf(paste0(path_shape, "/W09-05_GML/W09-05-g_Lake.shp")) %>%
  filter(W09_001 == "琵琶湖")
st_crs(lake_biwa) <- 4612 # JGD2000, EPSG:4612
## Land obj
land <- read_sf(paste0(path_shape, "/N03-140401_GML/N03-14_140401.shp")) %>%
  aggregate(list(substr(.$N03_007, 1, 2)), head, n = 1) %>%
  select(ID = Group.1, name = N03_001, geometry)
land.01 <- land %>% ms_simplify(keep = .01)
land.001 <- land %>% ms_simplify(keep = .001)
land.01focal <- land.01 %>%
  filter(name %in% c("福井県", "岐阜県", "三重県", "滋賀県", "京都府", "大阪府", "兵庫県")) %>%
  st_union()
land.001japan <- land.001 %>% aggregate(list(rep(1, 47)), head, n = 1)
land.001japan[1, "name"] <- "日本"
## Mapping
range_lon <- c(135.45, 136.75)
range_lat <- c(34.85, 35.95)
range_lon2 <- c(126, 149)
range_lat2 <- c(26, 46)
annot_lake <- c(136.16, 35.33)
map_main <- ggplot() +
  geom_sf(data = land.01focal, size = .4, color = "black") +  
  geom_sf(data = lake_biwa, fill = "white", size = .4, color = "black") +
  coord_sf(xlim = range_lon, ylim = range_lat, crs = st_crs(4612)) +
  annotate("text", x = annot_lake[1], y = annot_lake[2],
           label = "Lake\nBiwa", size = 6) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        axis.text = element_text(color = "black", size = 13),
        axis.title = element_blank()) +
  annotation_north_arrow(location = "bl", which_north = "true",
                         height = unit(2, "cm"), width = unit(2, "cm"),
                         pad_x = unit(0.5, "cm"), pad_y = unit(1, "cm"),
                         style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "bl", width_hint = 0.25) +
  panel_border(color = "black")
map_inset <- ggplot() +
  geom_sf(data = land.001japan, size = .3, color = "black") +
  geom_sf(data = lake_biwa, fill = "white", size = .03, color = "black") +
  coord_sf(xlim = range_lon2, ylim = range_lat2, crs = st_crs(4612)) +
  geom_rect(aes(xmin=range_lon[1], xmax=range_lon[2],
                ymin=range_lat[1], ymax=range_lat[2]),
            fill = NA, color = "red", size = .6) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0, 0, -.1, -.1), "null")) +
  panel_border(color = "black")


# Haplotype frequency data
## Create data table (tibble)
get.hapfreq <- function(x, y) {
  label <- x %>% phy_tree() %>%
    treeio::tree_subset(node = y, levels_back = 0) %>%
    `[[`("tip.label")
  phyfocal <- x %>%
    prune_taxa(label, .) %>% # The selected species only
    filter_tax_table(denoised == "no") %>% # Retained ASVs only
    prune_samples(sample_sums(.) > 0, .) # Non-zero samples only
  taxa_names(phyfocal) <- phyfocal %>% tax_table() %>%
    `[`(, "glob") %>% as.vector()
  otutab <- phyfocal %>% otu_table() %>% as.data.frame()
  otutab <- otutab %>%
    rowSums() %>% cbind(otutab, tot = .) %>%
    rownames_to_column(var = "sample")
  sampdat <- phyfocal %>% sample_data() %>% as.data.frame()
  class(sampdat) <- "data.frame" # Transform from phyloseq sample_data to data.frame
  sampdat <- sampdat %>%
    rownames_to_column(var = "sample")
  df <- right_join(sampdat, otutab, by = "sample")
  return(df)
}
hap_tib <- get.hapfreq(phylo_pre4[["dada"]], focalnode[["dada"]]) %>%
  group_by(lat, lon) %>%
  summarize(across(c(sample, river), first),
            across(c(starts_with("ASV_"), tot), sum))
hapid <- hap_tib %>% colnames() %>%
  str_detect(pattern = "^ASV") %>%
  subset(colnames(hap_tib), .)
## Create pie chart by sampling location
hap_tib2 <- hap_tib %>%
  tidyr::pivot_longer(cols = all_of(hapid),
                      names_to = "haplo", values_to = "reads") %>%
  tidyr::nest(data = c(haplo, reads)) %>%
  # Make a pie chart from each row and convert to grob
  mutate(pie.grob = map(data,
                        function(d) {
                          ggplotGrob(
                            ggplot(d, aes(x = 1, y = reads, fill = haplo)) +
                              geom_col(
                                color = "black", lwd = .4, alpha = .7,
                                show.legend = FALSE
                              ) +
                              scale_fill_manual(
                                values = wes_palette(
                                  "Darjeeling1",
                                  length(hapid)
                                )
                              ) +
                              coord_polar(theta = "y") +
                              theme_void()
                          )
                        }
  )) %>%
  # convert each grob to an annotation_custom layer. Also adjust the radius
  #   value to a reasonable size.
  rowwise() %>%
  mutate(radius = .07) %>%
  mutate(subgrob = list(
    annotation_custom(grob = pie.grob,
                      xmin = lon - radius, xmax = lon + radius,
                      ymin = lat - radius, ymax = lat + radius)
  ))


# Draw pie charts on maps
hapmap <- map_main + hap_tib2$subgrob
(map_final <- hapmap %>%
    ggdraw() +
    draw_plot(map_inset, x = 0.123, y = 0.67, width = 0.3, height = 0.3))
save_plot(paste0(path_output, "/03-Hapmap.svg"), map_final,
          base_height = 7, base_asp = 1 / 1)


# Save data
## Save R objects
path_saved_object <- paste0(path_output, "/01-Saved_object")
dir.create(path_saved_object, recursive=TRUE)
saveRDS(land, paste0(path_saved_object, "/land", ".obj")); rm(land)

# Save workspace and session info
save.image(paste0(path_output, "/hapnetmap.RData"))
writeLines(capture.output(sessionInfo()), paste0(path_output, "/hapnetmap.info"))
