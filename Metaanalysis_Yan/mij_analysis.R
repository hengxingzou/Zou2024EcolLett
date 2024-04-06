library(tidyverse)

# setwd("C:/Users/xy4279/Box/UT/4th year/priority effect collab HXZou/meta-analysis-mij")
# read in meta-analysis dataset 
live_ref_data = read_csv("Metaanalysis_Yan/data_live_ref.csv")
sterile_ref_data = read_csv("Metaanalysis_Yan/data_sterile_ref.csv")


### for dataset that use live reference soil 
# (look at effect of PSF, the part of microbiome that dynamically 
# responds to changes in the plant community)
live_ref = live_ref_data %>% 
  mutate(mAA = log(`AinA mean`)-log(`AinR mean`),
         mAA_var = (`AinA se`/`AinA mean`)^2 + (`AinR se`/`AinR mean`)^2, # for meta-analysis method, no need to worry for now
         mAB = log(`AinB mean`)-log(`AinR mean`),
         mAB_var = (`AinB se`/`AinB mean`)^2 + (`AinR se`/`AinR mean`)^2,
         mBA = log(`BinA mean`)-log(`BinR mean`),
         mBA_var = (`BinA se`/`BinA mean`)^2 + (`BinR se`/`BinR mean`)^2,
         mBB = log(`BinB mean`)-log(`BinR mean`),
         mBB_var = (`BinB se`/`BinB mean`)^2 + (`BinR se`/`BinR mean`)^2) %>% view

# we first combine all conspecific effects then remove duplicates
live_sp_A = live_ref %>% select(Experiment, `Species A`, mAA) %>%  # can include more columns for additional stats
  rename(species = `Species A`,
         mii = mAA)
live_sp_B = live_ref %>% select(Experiment, `Species B`, mBB) %>% 
  rename(species = `Species B`,
         mii = mBB)

mii_live = rbind(live_sp_A, live_sp_B) %>% 
  distinct(Experiment, species, .keep_all = T) %>% 
  pull(mii)

    
hist(mii_live)
mean(mii_live) 
# a negative value means a negative effect of the cultivated microbiome (PSF)
# on plant growth compared to live soil with "naive" microbiome
sd(mii_live)
# this is the (sample) distribution of reported mean mii of 56 species across 18 experiments 
# not the estimated overall mean of mii across species/experiments

# extract heterospecific microbial effect
mij_live = live_ref %>% 
  select(mAB,mBA) %>% unlist
hist(mij_live)
mean(mij_live)
# a less negative value compared to mii means a less negative heterospecific effect 
# of PSF than conspecific effect
sd(mij_live)



### for dataset that use sterile reference soil 
# (look at effect of the whole microbial community)
sterile_ref = sterile_ref_data %>% 
  mutate(mAA = log(`AinA mean`)-log(`AinR mean`),
         mAA_var = (`AinA se`/`AinA mean`)^2 + (`AinR se`/`AinR mean`)^2,
         mAB = log(`AinB mean`)-log(`AinR mean`),
         mAB_var = (`AinB se`/`AinB mean`)^2 + (`AinR se`/`AinR mean`)^2,
         mBA = log(`BinA mean`)-log(`BinR mean`),
         mBA_var = (`BinA se`/`BinA mean`)^2 + (`BinR se`/`BinR mean`)^2,
         mBB = log(`BinB mean`)-log(`BinR mean`),
         mBB_var = (`BinB se`/`BinB mean`)^2 + (`BinR se`/`BinR mean`)^2) %>% 
  # removing outliers (reported 0 biomass for non-germination, see analysis code for Yan et al. 2021)
  filter(!(species_pair %in% c("colubrina_spinosa_iriartea_deltoidea", 
                               "rumex_occidentalis_rumex_salicifolius"))) %>% view


sterile_sp_A = sterile_ref %>% select(Experiment, `Species A`, mAA) %>%  # can include more columns for additional stats
  rename(species = `Species A`,
         mii = mAA)
sterile_sp_B = sterile_ref %>% select(Experiment, `Species B`, mBB) %>% 
  rename(species = `Species B`,
         mii = mBB)

mii_sterile = rbind(sterile_sp_A, sterile_sp_B) %>% 
  distinct(Experiment, species, .keep_all = T) %>% 
  pull(mii)
hist(mii_sterile)
mean(mii_sterile) 
# a positive value means the entire soil microbiome of one plant species 
# has an overal positive effect on conspecific, compared to no microbes
sd(mii_sterile)


mij_sterile = sterile_ref %>% 
  select(mAB,mBA) %>% unlist
hist(mij_sterile)
mean(mij_sterile)
# a more positive value than mij means the entire soil microbiome of one plant species 
# has an overal more positive effect on heterospecific than conspecifics 
sd(mij_sterile)


summary = data.frame(live_ref = c(paste(round(mean(mii_live), 2), round(sd(mii_live), 2), sep=","),
                          paste(round(mean(mij_live), 2), round(sd(mij_live), 2), sep=",")),
             sterile_ref = c(paste(round(mean(mii_sterile), 2), round(sd(mii_sterile), 2), sep=","),
                             paste(round(mean(mij_sterile), 2), round(sd(mij_sterile), 2), sep=",")))
rownames(summary) = c("conspecific", "heterospecific")
summary

# either way, conspecific interactions through PSF (soil microbes) tend to be 
# more negative (less positive) than heterospecific interactions


