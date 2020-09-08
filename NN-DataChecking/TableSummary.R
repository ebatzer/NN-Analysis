library(plyr)

full_cover <- read.csv("Data/full-cover-23-april-2020.csv")
taxonomy <- read.csv("Data/site-taxonomy-2019-10-07.csv")

source(file = "NN-DataChecking/taxonomic_adjustments.R")

full_cover_adjusted <- Taxonomic.Adjustments(full_cover)

dim(full_cover_adjusted)

detach("package:plyr", unload=TRUE)

library(tidyverse)

changed_species <- full_cover_adjusted %>% select(site_code, Taxon, old_taxon) %>%
  filter(! Taxon == old_taxon) %>%
  rename("new_taxon" = Taxon) %>%
  distinct()

taxon_tally <- full_cover %>% 
  filter(live == 1 & 
         !is.na(max_cover)) %>%
  group_by(site_code, Taxon) %>%
  summarise(mean_cover = mean(max_cover),
            obs_freq = n())

# Generate table of proposed updates
proposed_updates <- taxon_tally %>% 
  left_join(changed_species, by = c("site_code" = "site_code", "Taxon" = "old_taxon")) %>%
  rowwise() %>% 
  mutate(new_taxon = if_else(is.na(new_taxon), Taxon, new_taxon)) %>%
  select(site_code, Taxon, new_taxon) %>%
  filter(paste(site_code, new_taxon) %in% paste(changed_species$site_code, changed_species$new_taxon)) %>%
  rename("old_taxon" = "Taxon")
  
# write.csv(proposed_updates, file = "NN-DataChecking/data/proposed_adjustments.csv")

# Generate table of frequencies, average cover
site_block_features <- full_cover %>% group_by(site_code) %>%
  summarise(n_years = max(year_trt + 1),
            first_year = min(year),
            last_year = max(year),
            n_blocks = length(unique(block)),
            n_plots = length(unique(plot)))

spec_summary <- full_cover %>% 
  filter(live == 1 & 
           !is.na(max_cover)) %>%
  group_by(site_code, Taxon, year) %>%
  summarise(mean_cover = mean(max_cover),
            obs_freq = n()) %>%
  left_join(site_block_features, by = "site_code") %>%
  mutate(mean_cover = mean_cover * (obs_freq / n_plots)) %>%
  select(site_code, Taxon, year, mean_cover, obs_freq) %>%
  ungroup() %>%
  mutate(Taxon = tolower(Taxon))

spec_summary_bytrt <- full_cover %>% 
  filter(live == 1 & 
           !is.na(max_cover)) %>%
  mutate(trt = replace_na(trt, "Control")) %>%
  group_by(site_code, Taxon, year, trt) %>%
  summarise(mean_cover = mean(max_cover),
            obs_freq = n()) %>%
  left_join(site_block_features, by = "site_code") %>%
  mutate(mean_cover = mean_cover * (obs_freq / n_plots)) %>%
  select(site_code, Taxon, year, mean_cover, obs_freq, trt) %>%
  ungroup() %>%
  mutate(Taxon = tolower(Taxon))

spec_summary <- bind_rows(spec_summary, spec_summary_bytrt) %>%
  mutate(trt = replace_na(as.character(trt), "ALL"))

taxonomy_cleaned <- taxonomy %>% 
  mutate(local_name = tolower(gsub("sp\\. .*", "sp\\.", local_name)),
         standard_taxon = tolower(standard_taxon)) %>%
  select(site_code, Family, local_name, standard_taxon) %>%
  distinct()

spec_summary <- spec_summary %>%
  left_join(taxonomy_cleaned,
            by = c("Taxon" = "standard_taxon", "site_code"))
                    
# Writing new dataframes for app
# write.csv(site_block_features, "NN-DataChecking/data/sb_features.csv")
write.csv(spec_summary, "NN-DataChecking/data/spec_summary.csv")
