library(tidyverse)
library(magrittr)
library(minpack.lm)
library(patchwork)
library(doParallel)
library(foreach)
library(Metrics)
registerDoParallel(cores = 4)

# Color scheme inspired by Samurai Champloo, Ep2 "Redeye Reprisal", Dir. Watanabe, Shin'ichiro

color_scheme_3 = c("#A6453C", "#444366", "#CF8095")
color_scheme_5 = c("#A6453C", "#444366", "#CF8095", "#3A6870", "#8A8C45")

########## Supporting Functions ##########

# Get arrival times

get_timeseq = function(nspp, bound = NULL, var = 0, even = T) {
  if (even) {timeseq = seq(0, bound, length = nspp)}
  else {timeseq = sort(runif(nspp, 0, bound))}
  timeseq = sapply(timeseq, function(x) rnorm(1, x, sd = var)) # randomization, normal dist; runif(1, x-var, x+var)
  return(timeseq + abs(min(timeseq))) # normalize minimum to 0
}

# # Delta t matrix
# 
# init_deltap_matrix = function(mean_times) {
#   mean_times = (mean_times - min(mean_times))
#   deltat_matrix = matrix(0, nrow = length(mean_times), ncol = length(mean_times))
#   for (i in 1:length(mean_times)) {
#     for (j in 1:length(mean_times)) {
#       deltat_matrix[i, j] = mean_times[j]-mean_times[i]
#     }
#   }
#   return(deltat_matrix)
# }

# Get emergence vector (a) within a season

get_emergence = function(timeseq, tau, lc_len) { # lc_len: length of plant life cycle
  a_matrix = matrix(0, nrow = length(timeseq), ncol = tau)
  for (spp in 1:length(timeseq)) {
    a = rep(0, tau)
    present = 1+abs(timeseq[spp])
    a[present] = 1
    a_matrix[spp, ] = a
  }
  return(a_matrix)
}

# Get reproduction vector (c) within a season

get_reproduction = function(timeseq, tau, lc_len) {
  new_timeseq = timeseq+lc_len
  return(get_emergence(new_timeseq, tau, lc_len))
}


########## Simulation ##########

model = function(plants, seeds, microbes, para, m_matrix, v_matrix, timeseq, tau, t, lc_len, return_all = T) {
  
  if (t %% tau != 0) {
    t = (t %/% tau)*tau
    warning("length of simulation is not multiples of time steps within a year,
             \n coercing into nearest multiples")
  }
  
  # Initial density of microbes at the beginning of each season
  # microbes_ini = microbes
  
  b = rep(ifelse(1:tau %% tau == 0, 0, 1), times = t %/% tau)
  # all_pop = tibble(Plants = plants, Seeds = seeds, Microbes = microbes) %>% 
  #             mutate(Time = 0, Species = 1:length(plants)) # if filling population at time 0
  all_pop = tibble()
  
  for (p in 1:t) {
    
    counter = p %% tau # counting time steps within a season
    if (counter == 1) {
      a_matrix = get_emergence(timeseq, tau, lc_len)
      c_matrix = get_reproduction(timeseq, tau, lc_len)
    }
    if (counter == 0) {counter = tau}
    
    # Plants
    new_plants = 
      a_matrix[, counter]*seeds + # emergence from seed bank
      (1-c_matrix[, counter])*plants*(1 - para$d*plants + as.numeric(m_matrix %*% microbes)) # process within life cycle
    
    # Seeds
    new_seeds = c_matrix[, counter]*para$lambda*plants + #*(1-plants/para$K) + # reproduction
      (1-c_matrix[, counter])*(1-a_matrix[, counter])*seeds # maintenance of seed bank
    
    # Microbes
    # new_microbes = b[p]*microbes*(1 + as.numeric(v_matrix %*% plants) - para$mu * microbes) +
    #   (1-b[p])*microbes_ini # reset to initial density
    # new_microbes = microbes*(1 + as.numeric(v_matrix %*% plants) - para$mu * microbes) # density carries over
    new_microbes = b[p]*microbes*(1 + as.numeric(v_matrix %*% plants) - para$mu * microbes) +
      (1-b[p])*microbes*para$surv_y # density carries over with a mortality

    all_pop = rbind(all_pop, 
                    tibble(Plants = new_plants, Seeds = new_seeds, Microbes = new_microbes) %>% 
                      mutate(Time = p, Species = 1:length(plants)))
    
    plants = new_plants
    seeds = new_seeds
    microbes = new_microbes
    
  }
  
  if (return_all) {return(all_pop)}
  
  return(c(plants, seeds, microbes))
}
