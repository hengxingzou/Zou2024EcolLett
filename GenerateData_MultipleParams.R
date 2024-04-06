source("Functions_PSF.R")
source("ModelFitting.R")

data_dir = "Data/MultipleParams"

rep = 5 # actual repetitions on the cluster: 100


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

# Model parameters, fixed

para = tibble(d = rep(0.1, nspp), # plant intraspecific competition
              mu = rep(0.1, nspp), # microbe intraspecific competition
              surv_y = rep(0.3, nspp), # end-of-year microbe survival
              lambda = rep(3, nspp) # plant fecundity
)

# Time

tau = 12
lc_len = 6
max_arriv = tau-lc_len-1
t = 1*tau # only one year for response surface

# Initial conditions

seeds_ini = rep(seeds_ini_single, nspp)
microbes_ini = rep(microbes_ini_single, nspp)

# Density gradient for response surface

initpop_combns = expand_grid(seeds_grad, seeds_grad, seeds_grad)

# Arrival time combinations

# arriv_combns = expand_grid(seq(0, max_arriv, 1), seq(0, max_arriv, 1), seq(0, max_arriv, 1))

arriv_combns = rbind(tibble(i = 0, j = 0, k = seq(0, max_arriv, 1)),
                     tibble(i = max_arriv, j = max_arriv, k = seq(0, max_arriv, 1)),
                     tibble(i = 0, j = max_arriv, k = seq(0, max_arriv, 1))) %>%
  distinct()

# Model fitting

inits = c(lambda = 3, alpha_intra = 0.2, alpha_inter1 = 0.2, alpha_inter2 = 0.2, 
          beta_112 = 0.2, beta_113 = 0.2, beta_123 = 0.2)
lowers = c(lambda = 1, alpha_intra = -3, alpha_inter1 = -3, alpha_inter2 = -3, 
           beta_112 = -3, beta_113 = -3, beta_123 = -3)
uppers = c(lambda = 5, alpha_intra = 3, alpha_inter1 = 3, alpha_inter2 = 3, 
           beta_112 = 3, beta_113 = 3, beta_123 = 3)

nls_ctrl = nls.control(maxiter = 1000, tol = 1e-8, minFactor = (1/10)*(1/1024),
                       printEval = FALSE, warnOnly = FALSE)


########## No Overlap ##########


alphabetas_nc_reps = foreach (r = 1:rep, .combine = rbind) %dopar% {
  
  # Cultivation and feedback rates, randomly generated from uniform distributions
  
  m = runif(1, -0.5, -0.1)
  v = runif(1, 0.1, 0.3)

  m_matrix = matrix(rep(m, nspp*nspp), nrow = nspp, byrow = T)
  
  # Randomly generated feedback rates, with the ratio of m_xx/m_yx informed by Yan et al. 2022
  # See FeedbackRatios.R for details
  # For simplification, draw only one ratio for all three species
  # In this scenario v is not randomly drawn to highlight the effect of m_xx/m_yx
  
  # m_ratio = rlnorm(1, meanlog = 0.142, sdlog = 1.087)
  # m_total = -0.5
  # v = 0.3
  
  # Percentage of interspecific feedback (m_yx in total m from plant x)
  # percent_im = 1 / (1+m_ratio)
  # m_matrix = diag(nspp)*rep(m_total, nspp)*(1-percent_im)
  # m_matrix[lower.tri(m_matrix)] = percent_im*m_total/(nspp-1)
  # m_matrix[upper.tri(m_matrix)] = percent_im*m_total/(nspp-1)
  
  # Construct v matrix
  
  percent_ic = 0
  v_matrix = diag(nspp)*rep(v, nspp)*(1-percent_ic)
  v_matrix[lower.tri(v_matrix)] = percent_ic*v/(nspp-1)
  v_matrix[upper.tri(v_matrix)] = percent_ic*v/(nspp-1)
  
  # Data generation
  
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
  
  # Model fitting
  
  alphabetas_nc = fit_all_hoi(arriv_combns, ijk_nc_rs) %>%
    mutate(pij = pj-pi, pjk = pk-pj, pik = pk-pi) %>% 
    mutate(rep = r)
  
}


########## Low Overlap ##########


alphabetas_al_reps = foreach (r = 1:rep, .combine = rbind) %dopar% {
  
  # Cultivation and feedback rates, randomly generated
  
  m = runif(1, -0.4, -0.2)
  v = runif(1, 0.2, 0.4)
  
  m_matrix = matrix(rep(m, nspp*nspp), nrow = nspp, byrow = T)
  
  # Randomly generated feedback rates, with the ratio of m_xx/m_yx informed by Yan et al. 2022
  # See FeedbackRatios.R for details
  # For simplification, draw only one ratio for all three species
  # In this scenario v is not randomly drawn to highlight the effect of m_xx/m_yx
  
  # m_ratio = rlnorm(1, meanlog = 0.142, sdlog = 1.087)
  # m_total = -0.5
  # v = 0.3
  
  # Percentage of interspecific feedback (m_yx in total m from plant x)
  # percent_im = 1 / (1+m_ratio)
  # m_matrix = diag(nspp)*rep(m_total, nspp)*(1-percent_im)
  # m_matrix[lower.tri(m_matrix)] = percent_im*m_total/(nspp-1)
  # m_matrix[upper.tri(m_matrix)] = percent_im*m_total/(nspp-1)
  
  # Construct v matrix
  percent_ic = 0.3
  v_matrix = diag(nspp)*rep(v, nspp)*(1-percent_ic)
  v_matrix[lower.tri(v_matrix)] = percent_ic*v/(nspp-1)
  v_matrix[upper.tri(v_matrix)] = percent_ic*v/(nspp-1)
  
  # Data generation
  
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
  
  # Model fitting
  
  alphabetas_al = fit_all_hoi(arriv_combns, ijk_al_rs) %>%
    mutate(pij = pj-pi, pjk = pk-pj, pik = pk-pi) %>% 
    mutate(rep = r)
  
}


########## High Overlap ##########


alphabetas_co_reps = foreach (r = 1:rep, .combine = rbind) %dopar% {
  
  # Cultivation and feedback rates, randomly generated
  
  m = runif(1, -0.4, -0.2)
  v = runif(1, 0.2, 0.4)
  
  m_matrix = matrix(rep(m, nspp*nspp), nrow = nspp, byrow = T)
  
  # Randomly generated feedback rates, with the ratio of m_xx/m_yx informed by Yan et al. 2022
  # See FeedbackRatios.R for details
  # For simplification, draw only one ratio for all three species
  # In this scenario v is not randomly drawn to highlight the effect of m_xx/m_yx
  
  # m_ratio = rlnorm(1, meanlog = 0.142, sdlog = 1.087)
  # m_total = -0.5
  # v = 0.3
  
  # Percentage of interspecific feedback (m_yx in total m from plant x)
  # percent_im = 1 / (1+m_ratio)
  # m_matrix = diag(nspp)*rep(m_total, nspp)*(1-percent_im)
  # m_matrix[lower.tri(m_matrix)] = percent_im*m_total/(nspp-1)
  # m_matrix[upper.tri(m_matrix)] = percent_im*m_total/(nspp-1)
  
  # Construct v matrix
  percent_ic = 0.6
  v_matrix = diag(nspp)*rep(v, nspp)*(1-percent_ic)
  v_matrix[lower.tri(v_matrix)] = percent_ic*v/(nspp-1)
  v_matrix[upper.tri(v_matrix)] = percent_ic*v/(nspp-1)
  
  # Data generation
  
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
  
  # Model fitting
  
  alphabetas_co = fit_all_hoi(arriv_combns, ijk_co_rs) %>%
    mutate(pij = pj-pi, pjk = pk-pj, pik = pk-pi) %>% 
    mutate(rep = r)
  
}


########## Write Data ##########


write_csv(alphabetas_al_reps, paste0(data_dir, "alphabetas_al_reps.csv"))
write_csv(alphabetas_nc_reps, paste0(data_dir, "alphabetas_nc_reps.csv"))
write_csv(alphabetas_co_reps, paste0(data_dir, "alphabetas_co_reps.csv"))