# 6_sp_correlation.R: this script gets species correlations.

rm(list = ls())

library(tidyverse)
library(spOccupancy)
library(ggcorrplot)
library(ggdist)
library(posterior)
library(here)

# Species correlations

load(here('data', 'data_bundles.rdata'))

sp.codes <- c("Cattle","Brocket deer", 
              "Collared peccary", 
              "Lowland tapir","Poaching",
              "White-lipped peccary", "White-tailed deer")


load(here("results", "spOccupancy_fits", "julian_dry_rainy_spatial", 
          "dry_rainy_ct_web-4-chain-2023-07-25.R"))

dry_rainy_web_ct <- out_web; rm(out_web)

load(here("results", "spOccupancy_fits", "julian_dry_rainy_spatial", 
          "dry_rainy_au_web-4-chain-2023-07-24.R"))
dry_rainy_web_au <- out_web; rm(out_web)

source(here("functions", "get_correlations_spOcc.R"))

N <- 7 # number of species
n.factors <- 2


mod_obj <- list(dry_rainy_web_ct, dry_rainy_web_au)


mod_list <- list(dry_rainy_web_ct, dry_rainy_web_au) 

corr_array_list <- lapply(mod_list, function(x) get.corr.spOcc(x, N, n.factors)$correlation_array) #first argument is model object from spOccupancy, second arg is number of species, third arg is number of latent variables
# 

names(corr_array_list) <- c('CTs', 'ARUs')

corr_ppd_pairwise <- function(species.1 = 'Poaching', species.2 = "Lowland tapir") {
  
  sp1.idx <- which(sp.codes == species.1)
  sp2.idx <- which(sp.codes == species.2)

  out <- data.frame(do.call(cbind, lapply(corr_array_list, function(x) x[sp1.idx, sp2.idx, ])), Species = paste0(species.1, ' - ', species.2)) %>% gather(Model, Corr, -Species)
  
  return(out)
}


corr_pairwise <- list(
  corr_ppd_pairwise('Cattle', 'Collared peccary'),
  corr_ppd_pairwise('Poaching', 'Collared peccary'),
  corr_ppd_pairwise('Cattle', 'Lowland tapir'),
  corr_ppd_pairwise('Poaching', 'Lowland tapir'),
  corr_ppd_pairwise('Cattle', 'Brocket deer'),
  corr_ppd_pairwise('Poaching', 'Brocket deer'),
  corr_ppd_pairwise('Cattle', 'White-lipped peccary'),
  corr_ppd_pairwise('Poaching', 'White-lipped peccary'),
  corr_ppd_pairwise('Cattle', 'White-tailed deer'),
  corr_ppd_pairwise('Poaching', 'White-tailed deer'),
  corr_ppd_pairwise('Cattle', 'Poaching')
)


# Correlation density plots -----------------------------------------------------------


corr_dist_plots_density <- corr_pairwise %>% 
  map(~.x %>% 
        mutate(Corr = signif(x = Corr, digits = 2)) %>% 
  ggplot(aes(x = Corr)) +
  geom_density(aes(fill=factor(Model)), 
               linewidth = 0.2, 
               alpha = 0.4) + 
  labs(x = 'Correlation', y = 'Density', 
       title = .x$Species) +
  scale_fill_viridis_d(name = 'Method', begin = 0, end = 0.9, labels = c('Audio', 'CT')) +
    theme(axis.text.y = element_text(size=13),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size=13),
          axis.title.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.position = "none",
          title = element_text(size = 13))

)


# Correlation point estimates plot ---------------------------------------------

corr_dist_plots_pointint <- corr_pairwise %>% 
  map(~.x %>% 
        ggplot(aes(Species, Corr, 
                   fill = Model, shape = Model)) +
        coord_flip(ylim = c(-1, 1)) +
        stat_pointinterval(position = "dodge",
                           point_interval = "median_qi",
                           .width = c(0.5, 0.89),
                           point_size = 3.5)  +
        scale_shape_manual(values = c(8,16, 24, 25)) +
        geom_hline(aes(yintercept = 0), linetype = 2)  +
        labs(y = .x$Species) +
        theme(axis.text.y = element_blank(),
              axis.title.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.x = element_text(size=13),
              axis.title.x = element_text(size=15),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black",),
              legend.position = "none",
              title = element_text(size = 13))
  )


# Correlation halfeye plot ----------------------------------------------------


corr_dist_plots_halfeye <- corr_pairwise %>% 
  map(~.x %>% 
        ggplot(aes(Species, Corr, 
                   fill = Model),
                   shape = Model) +
        coord_flip(ylim = c(-1, 1)) +
        stat_halfeye(position = "dodgejust",
                     aes(fill = after_stat(cut_cdf_qi(cdf, 
                                                .width = c(0.5, 0.89, 1))),
                         shape = Model),
                     point_interval = "median_qi",
                     .width = c(0.5, 0.89),
                     point_size = 3.5) +
        scale_fill_manual(values = c("lightseagreen", "paleturquoise", "white")) +
        scale_shape_manual(values = c(8,16)) +
        geom_hline(aes(yintercept = 0), linetype = 2)  +
        labs(title = .x$Species) +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=13),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black",),
        legend.position = "none",
        title = element_text(size = 13))


)

ci_corr <- corr_pairwise %>% 
  map(~.x %>%
        group_by(Species, Model) %>% median_qi(Corr, .width = 0.89)
  ) %>% 
  bind_rows()


# Save data ---------------------------------------------------------------


saveRDS(corr_dist_plots_pointint, "figures/corr_dist_plots_point.RDS")
saveRDS(corr_dist_plots_density, "figures/corr_dist_plots_density.RDS")
saveRDS(corr_dist_plots_halfeye, "figures/corr_dist_plots_halfeye.RDS")
saveRDS(ci_corr, "results/estimates/ci_corr.RDS")





