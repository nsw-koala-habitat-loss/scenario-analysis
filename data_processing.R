# get required libraries
library(foreign)
library(tidyverse)
library(terra)
library(tidyterra)
library(sf)
library(exactextractr)
library(qs)
library(INLA)
library(wesanderson)

# read functions
source("functions.R")

# read in woody vegetation extent data for 2011
woody <- rast(paste(getwd(), "/input/woody_cover/woody_nsw_2011.tif", sep = "")) %>% round()
woody_mask <- woody %>% classify(cbind(c(0, 1), c(NA, 1)))

# read in koala habitat data and reclassify so 1 = koala habitat, 0 = non-habitat then aggregate, reproject, and snap to woody cover layer
# then aggregate, reproject, and snap to woody cover layer
khab <- rast(paste(getwd(), "/input/koala_habitat/KoalaHabitatSuitabilityModelClasses.tif", sep = "")) %>% aggregate(5, fun = "modal") %>% round() %>% classify(cbind(c(1, 2, 3, 4, 5, 6), c(0, 0, 0, 1, 1, 1))) %>% project(woody)

# read in properties
props1 <- vect("input/risk_maps/Pred_Ag.shp") %>% as_sf()
props <- vect("input/risk_maps/Khab_risk_Ag.shp") %>% as_sf() %>% cbind(props1$KMR)
names(props) <- c("Woody", "Risk", "Khab_P", "KhabRisk", "KMR", "geometry")
rm(props1)
gc()

# read in nvr data, reclassify for treatment (category 1 and category 2 regulated) and control (category 2 - vulnerable and sensitive) then aggregate, reproject, and snap to woody cover layer
nvr_incl <- rast(paste(getwd(), "/input/nvr_mapping_v4/naluma_nsw_2017_abel0_c20221212_u9.tif", sep = "")) %>% aggregate(5, fun = "modal") %>% round() %>% classify(cbind(c(NA, 1), c(1, NA))) %>% project(woody)
gc()
# get treatment layer (1 = treatment, 0 = control)
nvr_treat <- rast(paste(getwd(), "/input/nvr_mapping_v4/naluma_nsw_2017_abkl0_c20221212_u9.tif", sep = "")) %>% aggregate(5, fun = "modal") %>% round() %>% classify(cbind(c(NA, 3, 4, 5, 6), c(1, 0, 0, NA, 0))) %>% project(woody)
gc()
# get control layer (1 = control, 0 = treatment)
nvr_contr <- rast(paste(getwd(), "/input/nvr_mapping_v4/naluma_nsw_2017_abkl0_c20221212_u9.tif", sep = "")) %>% aggregate(5, fun = "modal") %>% round() %>% classify(cbind(c(NA, 3, 4, 5, 6), c(0, 1, 1, NA, 1))) %>% project(woody)
gc()

# generate include layer
# multiply by the woody layer as we only care about raster cells that were woody in 2011
include <- nvr_incl * woody_mask
# generate treatment and control layers
treatment <- include * nvr_treat
treatment <- treatment %>% round()
control <- include * nvr_contr
control <- control %>% round()

# free up memory
rm(nvr_incl, nvr_treat, nvr_contr, include)
gc()

# get zonal statistics for baseline woody vegetation extent and koala habitat
woody_treat <- woody * treatment
woody_contr <- woody * control
khab_treat <- khab * woody_treat
khab_contr <- khab * woody_contr
# zonal statistics for koala habitat only
zonal_koala_treat <- as_tibble(props)
zonal_koala_contr <- as_tibble(props)
zonal_koala_treat <- bind_cols(zonal_koala_treat, khab_treat %>% exact_extract(props, "sum") %>% ceiling() %>% as_tibble())
zonal_koala_contr <- bind_cols(zonal_koala_contr, khab_contr %>% exact_extract(props, "sum") %>% ceiling() %>% as_tibble()) 

rm(woody_treat, woody_contr, khab_contr)
gc()

# get the new nvr layer that splits self-assessment and exempt areas
nvr_new_cat1cat2 <- rast(paste(getwd(), "/input//nvr_mapping_v11/naluma_nsw_2017_abpl0_c20240321_u11.tif", sep = "")) %>% aggregate(5, fun = "modal") %>% round() %>% classify(cbind(c(1, 2, 3, 4, 5, 6), c(1, 1, 0, 0, 0, 0))) %>% project(woody)
nvr_new_cat2 <- rast(paste(getwd(), "/input//nvr_mapping_v11/naluma_nsw_2017_abpl0_c20240321_u11.tif", sep = "")) %>% aggregate(5, fun = "modal") %>% round() %>% classify(cbind(c(1, 2, 3, 4, 5, 6), c(0, 1, 0, 0, 0, 0))) %>% project(woody)

# get the zonal statistics for treated koala habitat for the new nvr categories
khab_cat1cat2_treat <- khab_treat * nvr_new_cat1cat2
khab_cat2_treat <- khab_treat * nvr_new_cat2
zonal_koala_cat1cat2_treat <- as_tibble(props)
zonal_koala_cat2_treat <- as_tibble(props)
zonal_koala_cat1cat2_treat <- bind_cols(zonal_koala_cat1cat2_treat, khab_cat1cat2_treat %>% exact_extract(props, "sum") %>% ceiling() %>% as_tibble())
zonal_koala_cat2_treat <- bind_cols(zonal_koala_cat2_treat, khab_cat2_treat %>% exact_extract(props, "sum") %>% ceiling() %>% as_tibble())
khab_combined <- zonal_koala_treat %>% rename(treat = value) %>% bind_cols(zonal_koala_contr %>% rename(contr = value) %>% select(contr)) %>% bind_cols(zonal_koala_cat1cat2_treat %>% rename(treatc1c2 = value) %>% select(treatc1c2)) %>% bind_cols(zonal_koala_cat2_treat %>% rename(treatc2 = value) %>% select(treatc2))
khab_combined_matched <- khab_combined %>% filter((treat > 0 & contr == 0) | (treat == 0 & contr > 0))
khab_combined_mixed <- khab_combined %>% filter(treat > 0 & contr > 0)

# free up memory
rm(khab_cat1cat2_treat, khab_cat2_treat, nvr_new_cat1cat2, nvr_new_cat2)
gc()

# select properties that have cat 1 and cat 2 land
khab_combined_matched_cat1cat2 <- khab_combined_matched %>% filter(treatc1c2 > 0) %>% select(KMR, KhabRisk, treatc1c2)
khab_combined_mixed_cat1cat2 <- khab_combined_mixed %>% filter(treatc1c2 > 0) %>% select(KMR, KhabRisk, treatc1c2)

# select properties that have cat 2 land
khab_combined_matched_cat2 <- khab_combined_matched %>% filter(treatc2 > 0) %>% select(KMR, KhabRisk, treatc2)
khab_combined_mixed_cat2 <- khab_combined_mixed %>% filter(treatc2 > 0) %>% select(KMR, KhabRisk, treatc2)

# run scenarios

# load models if needed
Model_Matched_Koala <- readRDS("input/models/model_matched_koala.rds")
Model_Mixed_Koala <- readRDS("input/models/model_mixed_koala.rds")

# get posterior samples
# matched
Post_Matched_Koala <- inla.posterior.sample(10000, Model_Matched_Koala, selection = list("ba1:ci1" = 1, "time:ba1:ci1" = 1))
Post_Matched_Koala_BACI <- lapply(Post_Matched_Koala, function(x) {x$latent["ba1:ci1:1",1]}) %>% unlist() %>% as.matrix()
Post_Matched_Koala_BACITIME <- lapply(Post_Matched_Koala, function(x) {x$latent["time:ba1:ci1:1",1]}) %>% unlist() %>% as.matrix()
# mixed
Post_Mixed_Koala <- inla.posterior.sample(10000, Model_Mixed_Koala, selection = list("ba1:ci1" = 1, "time:ba1:ci1" = 1))
Post_Mixed_Koala_BACI <- lapply(Post_Mixed_Koala, function(x) {x$latent["ba1:ci1:1",1]}) %>% unlist() %>% as.matrix()
Post_Mixed_Koala_BACITIME <- lapply(Post_Mixed_Koala, function(x) {x$latent["time:ba1:ci1:1",1]}) %>% unlist() %>% as.matrix()

# matched cat 1 and cat 2
# simulate scanarios for each KMR
Matched_Scenario_Cat1Cat2 <- list()
index <- 1
# matched cat 1 and cat 2
for (i in unique(khab_combined_matched_cat1cat2$KMR)) {
    # get the risk, amount, and params
    Risk <- khab_combined_matched_cat1cat2 %>% filter(KMR == i) %>% pull(KhabRisk)
    Amount <- khab_combined_matched_cat1cat2 %>% filter(KMR == i) %>% pull(treatc1c2)
    Params <- cbind(Post_Matched_Koala_BACI, Post_Matched_Koala_BACITIME) %>% as.matrix()
    # simulate the scenario
    Matched_Scenario_Cat1Cat2[[index]] <- get_scenario(Risk, Amount, Params, 6)
    # increment index
    index <- index + 1
}
# give the list names
names(Matched_Scenario_Cat1Cat2) <- unique(khab_combined_matched_cat1cat2$KMR)
# rename the AREA column to V1 and add KMR field
index <- 1
for (i in unique(khab_combined_matched_cat1cat2$KMR)) {
    Matched_Scenario_Cat1Cat2[[index]] <- Matched_Scenario_Cat1Cat2[[index]] %>% mutate(KMR = i) %>% rename(AREA = V1)
    # increment index
    index <- index + 1
}
# bind all rows together
Matched_Scenario_Cat1Cat2_Combined <- do.call(rbind, Matched_Scenario_Cat1Cat2)
# save output
saveRDS(Matched_Scenario_Cat1Cat2_Combined, "output/matched_scenario_cat1cat2.rds")

# mixed cat 1 and cat 2
# simulate scanarios for each KMR
Mixed_Scenario_Cat1Cat2 <- list()
index <- 1
# mixed cat 1 and cat 2
for (i in unique(khab_combined_mixed_cat1cat2$KMR)) {
    # get the risk, amount, and params
    Risk <- khab_combined_mixed_cat1cat2 %>% filter(KMR == i) %>% pull(KhabRisk)
    Amount <- khab_combined_mixed_cat1cat2 %>% filter(KMR == i) %>% pull(treatc1c2)
    Params <- cbind(Post_Mixed_Koala_BACI, Post_Mixed_Koala_BACITIME) %>% as.matrix()
    # simulate the scenario
    Mixed_Scenario_Cat1Cat2[[index]] <- get_scenario(Risk, Amount, Params, 6)
    # increment index
    index <- index + 1
}
# give the list names
names(Mixed_Scenario_Cat1Cat2) <- unique(khab_combined_mixed_cat1cat2$KMR)
# rename the AREA column to V1 and add KMR field
index <- 1
for (i in unique(khab_combined_mixed_cat1cat2$KMR)) {
    Mixed_Scenario_Cat1Cat2[[index]] <- Mixed_Scenario_Cat1Cat2[[index]] %>% mutate(KMR = i) %>% rename(AREA = V1)
    # increment index
    index <- index + 1
}
# bind all rows together
Mixed_Scenario_Cat1Cat2_Combined <- do.call(rbind, Mixed_Scenario_Cat1Cat2)
# save output
saveRDS(Mixed_Scenario_Cat1Cat2_Combined, "output/mixed_scenario_cat1cat2.rds")

# matched cat 2
# simulate scanarios for each KMR
Matched_Scenario_Cat2 <- list()
index <- 1
# matched cat 2
for (i in unique(khab_combined_matched_cat2$KMR)) {
    # get the risk, amount, and params
    Risk <- khab_combined_matched_cat2 %>% filter(KMR == i) %>% pull(KhabRisk)
    Amount <- khab_combined_matched_cat2 %>% filter(KMR == i) %>% pull(treatc2)
    Params <- cbind(Post_Matched_Koala_BACI, Post_Matched_Koala_BACITIME) %>% as.matrix()
    # simulate the scenario
    Matched_Scenario_Cat2[[index]] <- get_scenario(Risk, Amount, Params, 6)
    # increment index
    index <- index + 1
}
# give the list names
names(Matched_Scenario_Cat2) <- unique(khab_combined_matched_cat2$KMR)
# rename the AREA column to V1 and add KMR field
index <- 1
for (i in unique(khab_combined_matched_cat2$KMR)) {
    Matched_Scenario_Cat2[[index]] <- Matched_Scenario_Cat2[[index]] %>% mutate(KMR = i) %>% rename(AREA = V1)
    # increment index
    index <- index + 1
}
# bind all rows together
Matched_Scenario_Cat2_Combined <- do.call(rbind, Matched_Scenario_Cat2)
# save output
saveRDS(Matched_Scenario_Cat2_Combined, "output/matched_scenario_cat2.rds")

# mixed cat 2
# simulate scanarios for each KMR
Mixed_Scenario_Cat2 <- list()
index <- 1
# mixed cat 2
for (i in unique(khab_combined_mixed_cat2$KMR)) {
    # get the risk, amount, and params
    Risk <- khab_combined_mixed_cat2 %>% filter(KMR == i) %>% pull(KhabRisk)
    Amount <- khab_combined_mixed_cat2 %>% filter(KMR == i) %>% pull(treatc2)
    Params <- cbind(Post_Mixed_Koala_BACI, Post_Mixed_Koala_BACITIME) %>% as.matrix()
    # simulate the scenario
    Mixed_Scenario_Cat2[[index]] <- get_scenario(Risk, Amount, Params, 6)
    # increment index
    index <- index + 1
}
# give the list names
names(Mixed_Scenario_Cat2) <- unique(khab_combined_mixed_cat2$KMR)
# rename the AREA column to V1 and add KMR field
index <- 1
for (i in unique(khab_combined_mixed_cat2$KMR)) {
    Mixed_Scenario_Cat2[[index]] <- Mixed_Scenario_Cat2[[index]] %>% mutate(KMR = i) %>% rename(AREA = V1)
    # increment index
    index <- index + 1
}
# bind all rows together
Mixed_Scenario_Cat2_Combined <- do.call(rbind, Mixed_Scenario_Cat2)
# save output
saveRDS(Mixed_Scenario_Cat2_Combined, "output/mixed_scenario_cat2.rds")

# get the full matched plus mixed scanarios
Scenario_Cat1Cat2 <- Matched_Scenario_Cat1Cat2_Combined %>% as_tibble()
Scenario_Cat1Cat2$AREA <- Matched_Scenario_Cat1Cat2_Combined$AREA + Mixed_Scenario_Cat1Cat2_Combined$AREA
Scenario_Cat2 <- Matched_Scenario_Cat2_Combined %>% as_tibble()
Scenario_Cat2$AREA <- Matched_Scenario_Cat2_Combined$AREA + Mixed_Scenario_Cat2_Combined$AREA

# calculate some summary measures for each KMR
SummaryCat1Cat2 <- Scenario_Cat1Cat2 %>% group_by(KMR) %>% summarise(Median = median(AREA), Mean = mean(AREA), Lower = quantile(AREA, 0.025), Upper = quantile(AREA, 0.975), P = sum(ifelse(AREA >= 0, 1, 0)) / n())
SummaryCat2 <- Scenario_Cat2 %>% group_by(KMR) %>% summarise(Median = median(AREA), Mean = mean(AREA), Lower = quantile(AREA, 0.025), Upper = quantile(AREA, 0.975), , P = sum(ifelse(AREA >= 0, 1, 0)) / n())

# now make some plots

# Cat 1 and Cat 2 scenario
Plot <- ggplot(Scenario_Cat1Cat2, aes(x = KMR, y = AREA, fill = KMR)) + geom_violin(color = NA) + theme_minimal() + scale_y_continuous(limits = c(-500, 5000), breaks = seq(-500, 5000, by = 500)) + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x=element_blank()) + geom_hline(yintercept = 0) + labs(y = "Reduction in clearing (ha / year)") + theme(legend.title = element_blank(), legend.position="bottom", axis.title.y = element_text(size = 18)) 

ggsave(Plot, file = "output/figures/scenario_cat1cat2.jpg", width = 30, height = 20, units = "cm", dpi = 300)

# Cat 2 scenario
Plot <- ggplot(Scenario_Cat2, aes(x = KMR, y = AREA, fill = KMR)) + geom_violin(color = NA) + theme_minimal() + scale_y_continuous(limits = c(-500, 5000), breaks = seq(-500, 5000, by = 500)) + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x=element_blank()) + geom_hline(yintercept = 0) + labs(y = "Reduction in clearing (ha / year)") + theme(legend.title = element_blank(), legend.position="bottom", axis.title.y = element_text(size = 18))

ggsave(Plot, file = "output/figures/scenario_cat2.jpg", width = 30, height = 20, units = "cm", dpi = 300)

