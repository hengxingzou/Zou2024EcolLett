source("Functions_PSF.R")

data_dir = "Data/"

# Read data from Yan et al. 2022 PNAS
# For more information on data collection and calculation, see https://www.pnas.org/doi/abs/10.1073/pnas.2122088119
# Only using data from experiments with sterile reference soil

sterile_ref_data = read_csv(paste0(data_dir, "/Metaanalysis_Yan/data_sterile_ref.csv"))

# Dataset that uses sterile reference soil

sterile_ref = sterile_ref_data %>% 
  mutate(mAA = log(`AinA mean`)-log(`AinR mean`),
         mAA_var = (`AinA se`/`AinA mean`)^2 + (`AinR se`/`AinR mean`)^2,
         mAB = log(`AinB mean`)-log(`AinR mean`),
         mAB_var = (`AinB se`/`AinB mean`)^2 + (`AinR se`/`AinR mean`)^2,
         mBA = log(`BinA mean`)-log(`BinR mean`),
         mBA_var = (`BinA se`/`BinA mean`)^2 + (`BinR se`/`BinR mean`)^2,
         mBB = log(`BinB mean`)-log(`BinR mean`),
         mBB_var = (`BinB se`/`BinB mean`)^2 + (`BinR se`/`BinR mean`)^2) %>% 
  # removing outliers (reported 0 biomass for non-germination, see analysis code for Yan et al. 2022)
  filter(!(species_pair %in% c("colubrina_spinosa_iriartea_deltoidea", 
                               "rumex_occidentalis_rumex_salicifolius"))) %>% 
  distinct()

# Select those with negative m_xx and m_yx

m_sterile_neg = sterile_ref %>% 
  select(Experiment, `Species A`, `Species B`, mAA, mBB, mAB, mBA) %>% 
  filter(if_all(starts_with("m"), ~ . < 0)) %>% 
  mutate(Reference = "Sterile")

# Extract m_xx and m_yx values, then calculate the ratio of m_xx/m_yx

m_sterile_neg %<>% mutate(Self_Other_Ratio_A = mAA/mBA, 
                          Self_Other_ratio_B = mBB/mAB)

self_other_ratios = c(m_sterile_neg$Self_Other_Ratio_A, m_sterile_neg$Self_Other_ratio_B)

# Examine the distribution (log normal)

hist(log(self_other_ratios), breaks = 50)

# Fit a log normal distribution

fit_params = MASS::fitdistr(self_other_ratios, "lognormal")

fit_params

# meanlog: 0.109, sdlog: 1.168

