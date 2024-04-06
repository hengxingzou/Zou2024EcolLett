source("Functions_PSF.R")

data_dir = "Data/m_ii_small/"


########## Shared Parameters Among All Simulations ##########


nspp = 3

# Note that to keep the dimension of data structures constant, two-plant communities are
# just three-plant communities with initial density of plant k, its associated
# microbiome, and any parameters related to k set to 0; hence nspp is set to 3
# for all simulations

# Initial Conditions

plants_ini = rep(0, nspp)
seeds_ini_single = 1
microbes_ini_single = 0.1

# Density gradient for response surface

seeds_grad = seq(1, 6, 1)

# Model parameters

para = tibble(d = rep(0.1, nspp), # plant intraspecific competition
              mu = rep(0.1, nspp), # microbe intraspecific competition
              surv_y = rep(0.3, nspp), # end-of-year microbe survival
              lambda = rep(3, nspp), # plant fecundity
              K = rep(20, nspp) # plant carrying capacity; not used in the model
)

# Time

tau = 12
lc_len = 6
max_arriv = tau-lc_len-1

# Cultivation and feedback rates

m = -0.3
v = 0.3

# All equal feedback rates

# m_matrix = matrix(rep(m, nspp*nspp), nrow = nspp, byrow = T)

# Microbiome specializes on its host (m_ij_small)

# m_matrix = matrix(rep(-0.1, nspp*nspp), nrow = nspp, byrow = T)
# diag(m_matrix) = rep(m, nspp)

# Microbiome specializes on others (m_ii_small)

m_matrix = matrix(rep(m, nspp*nspp), nrow = nspp, byrow = T)
diag(m_matrix) = -0.1


# For "Standard", use percent_ic = 0.3 for low and percent_ic = 0.6 for high overlap


########## Two Plants, Shared Parameters ##########


# Initial conditions

seeds_ini = c(seeds_ini_single, seeds_ini_single, 0)
microbes_ini = c(microbes_ini_single, microbes_ini_single, 0)

# Density gradient for response surface

initpop_combns = expand_grid(seeds_grad, seeds_grad, 0)

arriv_combns_ij = expand_grid(seq(0, max_arriv, 1), seq(0, max_arriv, 1), 0)


########## Two Plants, No Overlap ##########


# To generate v_matrix that is compatible to dimensions of other data structures
# in two-plant models, first generate 2*2 matrix, then append 0 to the third row and column

percent_ic = 0
v_matrix = diag(nspp-1)*rep(v, nspp-1)*(1-percent_ic)
v_matrix[lower.tri(v_matrix)] = percent_ic*v/(nspp-2)
v_matrix[upper.tri(v_matrix)] = percent_ic*v/(nspp-2)
v_matrix = rbind(v_matrix, rep(0, nspp-1)) %>% 
  cbind(rep(0, nspp))

# long time series

t = 50*tau

ij_nc_ts = foreach(row = 1:nrow(arriv_combns_ij), .combine = rbind) %dopar% {
  timeseq = arriv_combns_ij %>% slice(row) %>% as.numeric()
  model(plants_ini, seeds_ini, microbes_ini, para, m_matrix, v_matrix, timeseq, tau, t, lc_len) %>% 
    pivot_longer(cols = 1:3, values_to = "Population", names_to = "Group") %>% 
    mutate(pi = timeseq[1], pj = timeseq[2], pk = timeseq[3])
}

# response surface

t = 1*tau # only one year

ij_nc_rs = tibble()

for (row in 1:nrow(arriv_combns_ij)) {
  timeseq = arriv_combns_ij %>% slice(row) %>% as.numeric()
  
  dta = foreach(inits = 1:nrow(initpop_combns), .combine = rbind) %dopar% {
    seeds_ini = initpop_combns %>% slice(inits) %>% as.numeric()
    model(plants_ini, seeds_ini, microbes_ini, para, m_matrix, v_matrix, timeseq, tau, t, lc_len) %>% 
      pivot_longer(cols = 1:3, values_to = "Population", names_to = "Group") %>% 
      mutate(pi = timeseq[1], pj = timeseq[2], pk = timeseq[3], 
             init1 = seeds_ini[1], init2 = seeds_ini[2], init3 = seeds_ini[3])
  }
  
  ij_nc_rs = rbind(ij_nc_rs, dta)
}

# write data

write_csv(ij_nc_ts, paste0(data_dir, "ij_nc_ts.csv"))
write_csv(ij_nc_rs, paste0(data_dir, "ij_nc_rs.csv"))


########## Two Plants, Low Overlap ##########


# To generate v_matrix that is compatible to dimensions of other data structures
# in two-plant models, first generate 2*2 matrix, then append 0 to the third row and column

percent_ic = 0.25
v_matrix = diag(nspp-1)*rep(v, nspp-1)*(1-percent_ic)
v_matrix[lower.tri(v_matrix)] = percent_ic*v/(nspp-2)
v_matrix[upper.tri(v_matrix)] = percent_ic*v/(nspp-2)
v_matrix = rbind(v_matrix, rep(0, nspp-1)) %>% 
  cbind(rep(0, nspp))

# long time series

t = 50*tau

ij_al_ts = foreach(row = 1:nrow(arriv_combns_ij), .combine = rbind) %dopar% {
  timeseq = arriv_combns_ij %>% slice(row) %>% as.numeric()
  model(plants_ini, seeds_ini, microbes_ini, para, m_matrix, v_matrix, timeseq, tau, t, lc_len) %>% 
    pivot_longer(cols = 1:3, values_to = "Population", names_to = "Group") %>% 
    mutate(pi = timeseq[1], pj = timeseq[2], pk = timeseq[3])
}

# response surface

t = 1*tau # only one year

ij_al_rs = tibble()

for (row in 1:nrow(arriv_combns_ij)) {
  timeseq = arriv_combns_ij %>% slice(row) %>% as.numeric()
  
  dta = foreach(inits = 1:nrow(initpop_combns), .combine = rbind) %dopar% {
    seeds_ini = initpop_combns %>% slice(inits) %>% as.numeric()
    model(plants_ini, seeds_ini, microbes_ini, para, m_matrix, v_matrix, timeseq, tau, t, lc_len) %>% 
      pivot_longer(cols = 1:3, values_to = "Population", names_to = "Group") %>% 
      mutate(pi = timeseq[1], pj = timeseq[2], pk = timeseq[3], 
             init1 = seeds_ini[1], init2 = seeds_ini[2], init3 = seeds_ini[3])
  }
  
  ij_al_rs = rbind(ij_al_rs, dta)
}

# write data

write_csv(ij_al_ts, paste0(data_dir, "ij_al_ts.csv"))
write_csv(ij_al_rs, paste0(data_dir, "ij_al_rs.csv"))


########## Two Plants, High Overlap ##########


# To generate v_matrix that is compatible to dimensions of other data structures
# in two-plant models, first generate 2*2 matrix, then append 0 to the third row and column

percent_ic = 0.5
v_matrix = diag(nspp-1)*rep(v, nspp-1)*(1-percent_ic)
v_matrix[lower.tri(v_matrix)] = percent_ic*v/(nspp-2)
v_matrix[upper.tri(v_matrix)] = percent_ic*v/(nspp-2)
v_matrix = rbind(v_matrix, rep(0, nspp-1)) %>% 
  cbind(rep(0, nspp))

# long time series

t = 50*tau

ij_co_ts = foreach(row = 1:nrow(arriv_combns_ij), .combine = rbind) %dopar% {
  timeseq = arriv_combns_ij %>% slice(row) %>% as.numeric()
  model(plants_ini, seeds_ini, microbes_ini, para, m_matrix, v_matrix, timeseq, tau, t, lc_len) %>% 
    pivot_longer(cols = 1:3, values_to = "Population", names_to = "Group") %>% 
    mutate(pi = timeseq[1], pj = timeseq[2], pk = timeseq[3])
}

# response surface

t = 1*tau # only one year

ij_co_rs = tibble()

for (row in 1:nrow(arriv_combns_ij)) {
  timeseq = arriv_combns_ij %>% slice(row) %>% as.numeric()
  
  dta = foreach(inits = 1:nrow(initpop_combns), .combine = rbind) %dopar% {
    seeds_ini = initpop_combns %>% slice(inits) %>% as.numeric()
    model(plants_ini, seeds_ini, microbes_ini, para, m_matrix, v_matrix, timeseq, tau, t, lc_len) %>% 
      pivot_longer(cols = 1:3, values_to = "Population", names_to = "Group") %>% 
      mutate(pi = timeseq[1], pj = timeseq[2], pk = timeseq[3], 
             init1 = seeds_ini[1], init2 = seeds_ini[2], init3 = seeds_ini[3])
  }
  
  ij_co_rs = rbind(ij_co_rs, dta)
}

# write data

write_csv(ij_co_ts, paste0(data_dir, "ij_co_ts.csv"))
write_csv(ij_co_rs, paste0(data_dir, "ij_co_rs.csv"))


########## Three Plants, Shared Parameters ##########


# Initial conditions

seeds_ini = rep(seeds_ini_single, nspp)
microbes_ini = rep(microbes_ini_single, nspp)

# Density gradient for response surface

initpop_combns = expand_grid(seeds_grad, seeds_grad, seeds_grad)

# Full combinations of arrival times

# arriv_combns = expand_grid(seq(0, max_arriv, 1), seq(0, max_arriv, 1), seq(0, max_arriv, 1))

# Three scenarios only (before, after, between, equal intervals; faster)

arriv_combns = rbind(tibble(i = 0, j = 0, k = seq(0, max_arriv, 1)),
                     tibble(i = max_arriv, j = max_arriv, k = seq(0, max_arriv, 1)),
                     tibble(i = 0, j = max_arriv, k = seq(0, max_arriv, 1))) %>%
  distinct()

# Longer year with fixed interval between i and j

# arriv_combns = rbind(tibble(i = 5, j = 11, k = seq(0, 11)))


########## Three Species, No Overlap ##########


percent_ic = 0
v_matrix = diag(nspp)*rep(v, nspp)*(1-percent_ic)
v_matrix[lower.tri(v_matrix)] = percent_ic*v/(nspp-1)
v_matrix[upper.tri(v_matrix)] = percent_ic*v/(nspp-1)

# long time series

t = 50*tau

ijk_nc_ts = foreach(row = 1:nrow(arriv_combns), .combine = rbind) %dopar% {
  timeseq = arriv_combns %>% slice(row) %>% as.numeric()
  model(plants_ini, seeds_ini, microbes_ini, para, m_matrix, v_matrix, timeseq, tau, t, lc_len) %>%
    pivot_longer(cols = 1:3, values_to = "Population", names_to = "Group") %>%
    mutate(pi = timeseq[1], pj = timeseq[2], pk = timeseq[3])
}

# response surface

t = 1*tau # only one year

ijk_nc_rs = tibble()

for (row in 1:nrow(arriv_combns)) {
  timeseq = arriv_combns %>% slice(row) %>% as.numeric()
  
  dta = foreach(inits = 1:nrow(initpop_combns), .combine = rbind) %dopar% {
    seeds_ini = initpop_combns %>% slice(inits) %>% as.numeric()
    model(plants_ini, seeds_ini, microbes_ini, para, m_matrix, v_matrix, timeseq, tau, t, lc_len) %>% 
      pivot_longer(cols = 1:3, values_to = "Population", names_to = "Group") %>% 
      mutate(pi = timeseq[1], pj = timeseq[2], pk = timeseq[3], 
             init1 = seeds_ini[1], init2 = seeds_ini[2], init3 = seeds_ini[3])
  }
  
  ijk_nc_rs = rbind(ijk_nc_rs, dta)
}

# write data

write_csv(ijk_nc_ts, paste0(data_dir, "ijk_nc_ts.csv"))
write_csv(ijk_nc_rs, paste0(data_dir, "ijk_nc_rs.csv"))


########## Three Species, Low Overlap ##########


percent_ic = 0.25
v_matrix = diag(nspp)*rep(v, nspp)*(1-percent_ic)
v_matrix[lower.tri(v_matrix)] = percent_ic*v/(nspp-1)
v_matrix[upper.tri(v_matrix)] = percent_ic*v/(nspp-1)

# long time series

t = 50*tau

ijk_al_ts = foreach(row = 1:nrow(arriv_combns), .combine = rbind) %dopar% {
  timeseq = arriv_combns %>% slice(row) %>% as.numeric()
  model(plants_ini, seeds_ini, microbes_ini, para, m_matrix, v_matrix, timeseq, tau, t, lc_len) %>%
    pivot_longer(cols = 1:3, values_to = "Population", names_to = "Group") %>%
    mutate(pi = timeseq[1], pj = timeseq[2], pk = timeseq[3])
}

# response surface

t = 1*tau # only one year

ijk_al_rs = tibble()

for (row in 1:nrow(arriv_combns)) {
  timeseq = arriv_combns %>% slice(row) %>% as.numeric()
  
  dta = foreach(inits = 1:nrow(initpop_combns), .combine = rbind) %dopar% {
    seeds_ini = initpop_combns %>% slice(inits) %>% as.numeric()
    model(plants_ini, seeds_ini, microbes_ini, para, m_matrix, v_matrix, timeseq, tau, t, lc_len) %>% 
      pivot_longer(cols = 1:3, values_to = "Population", names_to = "Group") %>% 
      mutate(pi = timeseq[1], pj = timeseq[2], pk = timeseq[3], 
             init1 = seeds_ini[1], init2 = seeds_ini[2], init3 = seeds_ini[3])
  }
  
  ijk_al_rs = rbind(ijk_al_rs, dta)
}

# write data

write_csv(ijk_al_ts, paste0(data_dir, "ijk_al_ts.csv"))
write_csv(ijk_al_rs, paste0(data_dir, "ijk_al_rs.csv"))


########## Three Species, High Overlap ##########


percent_ic = 0.5
v_matrix = diag(nspp)*rep(v, nspp)*(1-percent_ic)
v_matrix[lower.tri(v_matrix)] = percent_ic*v/(nspp-1)
v_matrix[upper.tri(v_matrix)] = percent_ic*v/(nspp-1)

# long time series

t = 50*tau

ijk_co_ts = foreach(row = 1:nrow(arriv_combns), .combine = rbind) %dopar% {
  timeseq = arriv_combns %>% slice(row) %>% as.numeric()
  model(plants_ini, seeds_ini, microbes_ini, para, m_matrix, v_matrix, timeseq, tau, t, lc_len) %>%
    pivot_longer(cols = 1:3, values_to = "Population", names_to = "Group") %>%
    mutate(pi = timeseq[1], pj = timeseq[2], pk = timeseq[3])
}

# response surface

t = 1*tau # only one year

ijk_co_rs = tibble()

for (row in 1:nrow(arriv_combns)) {
  timeseq = arriv_combns %>% slice(row) %>% as.numeric()
  
  dta = foreach(inits = 1:nrow(initpop_combns), .combine = rbind) %dopar% {
    seeds_ini = initpop_combns %>% slice(inits) %>% as.numeric()
    model(plants_ini, seeds_ini, microbes_ini, para, m_matrix, v_matrix, timeseq, tau, t, lc_len) %>% 
      pivot_longer(cols = 1:3, values_to = "Population", names_to = "Group") %>% 
      mutate(pi = timeseq[1], pj = timeseq[2], pk = timeseq[3], 
             init1 = seeds_ini[1], init2 = seeds_ini[2], init3 = seeds_ini[3])
  }
  
  ijk_co_rs = rbind(ijk_co_rs, dta)
}

# write data

write_csv(ijk_co_ts, paste0(data_dir, "ijk_co_ts.csv"))
write_csv(ijk_co_rs, paste0(data_dir, "ijk_co_rs.csv"))
