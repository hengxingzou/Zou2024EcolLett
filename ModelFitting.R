########## Phenomenological Models ##########


# No HOI

ri_model_1 = function(N_self, para) {
  with(para, {
    t1 = alpha_intra*N_self
    lambda*exp(-t1)
  })
}

ri_model_2 = function(N_self, N_other, para) {
  with(para, {
    t1 = alpha_intra*N_self
    t2 = alpha_inter*N_other
    lambda*exp(-t1-t2)
  })
}

ri_model_3 = function(N_self, N_other1, N_other2, para) {
  with(para, {
    t1 = alpha_intra*N_self
    t2 = alpha_inter1*N_other1
    t3 = alpha_inter2*N_other2
    lambda*exp(-t1-t2-t3)
  })
}

# HOI

ri_model_3_hoi = function(N_self, N_other1, N_other2, para) {
  with(para, {
    
    t1 = alpha_intra*N_self
    t2 = alpha_inter1*N_other1
    t3 = alpha_inter2*N_other2
    
    b1 = beta_112*N_self*N_other1
    b2 = beta_113*N_self*N_other2
    b3 = beta_123*N_other1*N_other2
    
    lambda*exp(-t1-t2-t3-b1-b2-b3)
  })
}


########## Processing Data ##########


transform_dt = function(arriv_output, timeseq, ts = F) {
  # ts: T if time series
  
  if (ts) {
    
    init = arriv_output %>% 
      filter(Time == 1, Group == "Plants") # initial population
    transformed = bind_rows(init, 
                            arriv_output %>% filter(Time %% tau == 0, Group == "Seeds")) %>% 
      mutate(Season = Time %% tau) %>% 
      filter(pi == timeseq[1], pj == timeseq[2], pk == timeseq[3]) %>%
      pivot_wider(names_from = Species, values_from = Population) %>% 
      select(-(1:5))
    
    out_tp1 = transformed[-1, -1]
    out_t = transformed[-nrow(transformed), -1]
    transformed = cbind(out_tp1, out_t)
    colnames(transformed) = c("N1tp1", "N2tp1", "N3tp1", "N1t", "N2t", "N3t")
    
  } else {
    
    transformed = arriv_output %>% 
      filter(Time %% tau == 0, Group == "Seeds") %>% 
      mutate(Season = Time/tau) %>% 
      filter(pi == timeseq[1], pj == timeseq[2], pk == timeseq[3]) %>%
      pivot_wider(names_from = Species, values_from = Population) %>% 
      select(-(1:5), -Season) %>% 
      rename(N1tp1 = `1`, N2tp1 = `2`, N3tp1 = `3`, 
             N1t = init1, N2t = init2, N3t = init3)
  }
  
  
  return(transformed)
}


########## Fitting Pairwise Models ##########


# Model fitting, one species

fit_one_spp = function(arriv_combns, arriv_output, rep = 1) {
  
  fitted_a = tibble()
  
  for (i in 1:rep) {
    
    for (row in 1:nrow(arriv_combns)) {
      
      # Select one arrival time sequence
      timeseq = arriv_combns %>% slice(row) %>% as.numeric()
      
      # Preparing data
      sample_dta = transform_dt(arriv_output, timeseq)
      
      # Model fitting
      m1 = nlsLM(
        (N1tp1) ~ (N1t*ri_model_1(N1t, 
                                  para = list(lambda = lambda, 
                                              alpha_intra = alpha_intra))),
        data = sample_dta, 
        start = inits,
        lower = lowers,
        upper = uppers,
        algorithm = "port",
        control = nls_ctrl
      )
      
      # Extract alphas
      a_each = tibble(
        alphaii = summary(m1)$coefficients[2, 1],
        pi = timeseq[1], pj = timeseq[2], pk = timeseq[3], 
        rep = i
      )
      
      # Add to all results
      fitted_a = rbind(fitted_a, a_each)
      
    }
  }
  
  return(fitted_a)
  
}

# Model fitting, two species

fit_all_2 = function(arriv_combns, arriv_output, spp_pair, ts = F, rep = 1) {
  
  fitted_a = tibble()
  
  for (i in 1:rep) {
    
    for (row in 1:nrow(arriv_combns)) {
      
      # Select one arrival time sequence
      timeseq = arriv_combns %>% slice(row) %>% as.numeric()
      # print(timeseq)

      # Preparing data
      sample_dta = transform_dt(arriv_output, timeseq, ts)
      
      # Model fitting
      if (! 3 %in% spp_pair) { # species 1 and 2
        
        m1 = nlsLM(
          (N1tp1) ~ (N1t*ri_model_2(N1t, N2t,
                                    para = list(lambda = lambda, alpha_intra = alpha_intra, alpha_inter = alpha_inter))),
          data = sample_dta, 
          start = inits,
          lower = lowers,
          upper = uppers,
          algorithm = "port",
          control = nls_ctrl
        )
        
        m2 = nlsLM(
          (N2tp1) ~ (N2t*ri_model_2(N2t, N1t,
                                    para = list(lambda = lambda, alpha_intra = alpha_intra, alpha_inter = alpha_inter))),
          data = sample_dta, 
          start = inits,
          lower = lowers,
          upper = uppers,
          algorithm = "port",
          control = nls_ctrl
        )
        
        k = 3
      }
      
      if (! 1 %in% spp_pair) { # species 2 and 3
        
        m1 = nlsLM(
          (N2tp1) ~ (N2t*ri_model_2(N2t, N3t,
                                    para = list(lambda = lambda, alpha_intra = alpha_intra, alpha_inter = alpha_inter))),
          data = sample_dta, 
          start = inits,
          lower = lowers,
          upper = uppers,
          algorithm = "port",
          control = nls_ctrl
        )
        
        m2 = nlsLM(
          (N3tp1) ~ (N3t*ri_model_2(N3t, N2t,
                                    para = list(lambda = lambda, alpha_intra = alpha_intra, alpha_inter = alpha_inter))),
          data = sample_dta, 
          start = inits,
          lower = lowers,
          upper = uppers,
          algorithm = "port",
          control = nls_ctrl
        )
        
        k = 1
      }
      
      if (! 2 %in% spp_pair) { # species 1 and 3
        
        m1 = nlsLM(
          (N1tp1) ~ (N1t*ri_model_2(N1t, N3t,
                                    para = list(lambda = lambda, alpha_intra = alpha_intra, alpha_inter = alpha_inter))),
          data = sample_dta, 
          start = inits,
          lower = lowers,
          upper = uppers,
          algorithm = "port",
          control = nls_ctrl
        )
        
        m2 = nlsLM(
          (N3tp1) ~ (N3t*ri_model_2(N3t, N1t,
                                    para = list(lambda = lambda, alpha_intra = alpha_intra, alpha_inter = alpha_inter))),
          data = sample_dta, 
          start = inits,
          lower = lowers,
          upper = uppers,
          algorithm = "port",
          control = nls_ctrl
        )
        
        k = 2
      }
      
      # Extract alphas
      params1 = summary(m1)
      params2 = summary(m2)
      
      a_each = tibble(
        alphaii = params1$coefficients[2, 1], seii = params1$coefficients[2, 2],
        alphaij = params1$coefficients[3, 1], seij = params1$coefficients[3, 2],
        
        alphaji = params2$coefficients[3, 1], seji = params2$coefficients[3, 2],
        alphajj = params2$coefficients[2, 1], sejj = params2$coefficients[2, 2],
        
        pi = timeseq[spp_pair[1]], pj = timeseq[spp_pair[2]], pk = timeseq[k], 
        rep = i
      )
      
      # Add to all results
      fitted_a = rbind(fitted_a, a_each)
      
    }
  }
  
  return(fitted_a)
}


# Model fitting, three species

fit_all = function(arriv_combns, arriv_output, ts = F, rep = 1) {
  
  fitted_a = tibble()
  
  for (i in 1:rep) {
    
    fitted_each = foreach(row = 1:nrow(arriv_combns), .combine = rbind) %dopar% {
      
    # for (row in 1:nrow(arriv_combns)) {
      
      # Select one arrival time sequence
      timeseq = arriv_combns %>% slice(row) %>% as.numeric()

      # Preparing data
      sample_dta = transform_dt(arriv_output, timeseq, ts)
      
      # Model fitting: species 1, 2, 3
      m1 = nlsLM(
        (N1tp1) ~ (N1t*ri_model_3(N1t, N2t, N3t, 
                                  para = list(lambda = lambda, 
                                              alpha_intra = alpha_intra, alpha_inter1 = alpha_inter1, alpha_inter2 = alpha_inter2))),
        data = sample_dta, 
        start = inits,
        lower = lowers,
        upper = uppers,
        algorithm = "port",
        control = nls_ctrl
      )
      
      m2 = nlsLM(
        (N2tp1) ~ (N2t*ri_model_3(N2t, N1t, N3t, 
                                  para = list(lambda = lambda, 
                                              alpha_intra = alpha_intra, alpha_inter1 = alpha_inter1, alpha_inter2 = alpha_inter2))),
        data = sample_dta, 
        start = inits,
        lower = lowers,
        upper = uppers,
        algorithm = "port",
        control = nls_ctrl
      )
      
      m3 = nlsLM(
        (N3tp1) ~ (N3t*ri_model_3(N3t, N1t, N2t, 
                                  para = list(lambda = lambda, 
                                              alpha_intra = alpha_intra, alpha_inter1 = alpha_inter1, alpha_inter2 = alpha_inter2))),
        data = sample_dta, 
        start = inits,
        lower = lowers,
        upper = uppers,
        algorithm = "port",
        control = nls_ctrl
      )
      
      # Extract alphas
      params1 = summary(m1)
      params2 = summary(m2)
      params3 = summary(m3)
      
      tibble(
        alphaii = params1$coefficients[2, 1], seii = params1$coefficients[2, 2],
        alphaij = params1$coefficients[3, 1], seij = params1$coefficients[3, 2],
        alphaik = params1$coefficients[4, 1], seik = params1$coefficients[4, 2],
        
        alphaji = params2$coefficients[3, 1], seji = params2$coefficients[3, 2],
        alphajj = params2$coefficients[2, 1], sejj = params2$coefficients[2, 2],
        alphajk = params2$coefficients[4, 1], sejk = params2$coefficients[4, 2],
        
        alphaki = params3$coefficients[3, 1], seki = params3$coefficients[3, 2],
        alphakj = params3$coefficients[4, 1], sekj = params3$coefficients[4, 2],
        alphakk = params3$coefficients[2, 1], sekk = params3$coefficients[2, 2],
        
        pi = timeseq[1], pj = timeseq[2], pk = timeseq[3], 
        rep = i
      )
      
      # Add to all results
      # fitted_a = rbind(fitted_a, a_each)
      
    }
    
    fitted_a = rbind(fitted_a, fitted_each)
  }
  
  return(fitted_a)
}


########## Fitting HOI Models ##########


fit_all_hoi = function(arriv_combns, arriv_output, ts = F, rep = 1) {
  
  fitted_coeffs = tibble()
  
  for (i in 1:rep) {
    
    for (row in 1:nrow(arriv_combns)) {
      
      # Select one arrival time sequence
      timeseq = arriv_combns %>% slice(row) %>% as.numeric()
      
      # Preparing data
      sample_dta = transform_dt(arriv_output, timeseq, ts)

      # Model fitting: species 1, 2, 3
      m1 = nlsLM(
        (N1tp1) ~ (N1t*ri_model_3_hoi(N1t, N2t, N3t, 
                                      para = list(lambda = lambda, 
                                                  alpha_intra = alpha_intra, alpha_inter1 = alpha_inter1, alpha_inter2 = alpha_inter2, 
                                                  beta_112 = beta_112, beta_113 = beta_113, beta_123 = beta_123))),
        data = sample_dta, 
        start = inits,
        lower = lowers,
        upper = uppers,
        algorithm = "port",
        control = nls_ctrl
      )
      
      m2 = nlsLM(
        (N2tp1) ~ (N2t*ri_model_3_hoi(N2t, N1t, N3t, 
                                      para = list(lambda = lambda, 
                                                  alpha_intra = alpha_intra, alpha_inter1 = alpha_inter1, alpha_inter2 = alpha_inter2, 
                                                  beta_112 = beta_112, beta_113 = beta_113, beta_123 = beta_123))),
        data = sample_dta, 
        start = inits,
        lower = lowers,
        upper = uppers,
        algorithm = "port",
        control = nls_ctrl
      )
      
      m3 = nlsLM(
        (N3tp1) ~ (N3t*ri_model_3_hoi(N3t, N1t, N2t, 
                                      para = list(lambda = lambda, 
                                                  alpha_intra = alpha_intra, alpha_inter1 = alpha_inter1, alpha_inter2 = alpha_inter2, 
                                                  beta_112 = beta_112, beta_113 = beta_113, beta_123 = beta_123))),
        data = sample_dta, 
        start = inits,
        lower = lowers,
        upper = uppers,
        algorithm = "port",
        control = nls_ctrl
      )
      
      # Extract alphas
      params1 = summary(m1)
      params2 = summary(m2)
      params3 = summary(m3)
      
      coeff_each = tibble(
        
        # Pairwise alphas
        
        alphaii = params1$coefficients[2, 1], seii = params1$coefficients[2, 2],
        alphaij = params1$coefficients[3, 1], seij = params1$coefficients[3, 2],
        alphaik = params1$coefficients[4, 1], seik = params1$coefficients[4, 2],
        
        alphaji = params2$coefficients[3, 1], seji = params2$coefficients[3, 2],
        alphajj = params2$coefficients[2, 1], sejj = params2$coefficients[2, 2],
        alphajk = params2$coefficients[4, 1], sejk = params2$coefficients[4, 2],
        
        alphaki = params3$coefficients[3, 1], seki = params3$coefficients[3, 2],
        alphakj = params3$coefficients[4, 1], sekj = params3$coefficients[4, 2],
        alphakk = params3$coefficients[2, 1], sekk = params3$coefficients[2, 2],
        
        # HOI betas
        
        betaiij = params1$coefficients[5, 1], seiij = params1$coefficients[5, 2],
        betaiik = params1$coefficients[6, 1], seiik = params1$coefficients[6, 2],
        betaijk = params1$coefficients[7, 1], seijk = params1$coefficients[7, 2],
        
        betajji = params2$coefficients[5, 1], sejji = params1$coefficients[5, 2],
        betajjk = params2$coefficients[6, 1], sejjk = params1$coefficients[6, 2],
        betajik = params2$coefficients[7, 1], sejik = params1$coefficients[7, 2],
        
        betakki = params3$coefficients[5, 1], sekki = params1$coefficients[5, 2],
        betakkj = params3$coefficients[6, 1], sekkj = params1$coefficients[6, 2],
        betakij = params3$coefficients[7, 1], sekij = params1$coefficients[7, 2],
        
        pi = timeseq[1], pj = timeseq[2], pk = timeseq[3], 
        rep = i
      )
      
      # Add to all results
      fitted_coeffs = rbind(fitted_coeffs, coeff_each)
      
    }
  }
  
  return(fitted_coeffs)
}


########## Compare Pairwise and HOI Models ##########


compare_fit = function(sample_dta, spp) {
  
  # Initials for pairwise models
  inits = c(lambda = 3, alpha_intra = 0.2, alpha_inter1 = 0.2, alpha_inter2 = 0.2)
  lowers = c(lambda = 1, alpha_intra = -1, alpha_inter1 = -1, alpha_inter2 = -1)
  uppers = c(lambda = 5, alpha_intra = 1, alpha_inter1 = 1, alpha_inter2 = 1)
  
  if (spp == 1) {
    
    m_pairwise = nlsLM(
      (N1tp1) ~ (N1t*ri_model_3(N1t, N2t, N3t, 
                                para = list(lambda = lambda, 
                                            alpha_intra = alpha_intra, alpha_inter1 = alpha_inter1, alpha_inter2 = alpha_inter2))),
      data = sample_dta, 
      start = inits,
      lower = lowers,
      upper = uppers,
      algorithm = "port",
      control = nls_ctrl
    )
  }
  
  if (spp == 2) {
    
    m_pairwise = nlsLM(
      (N2tp1) ~ (N2t*ri_model_3(N2t, N1t, N3t, 
                                para = list(lambda = lambda, 
                                            alpha_intra = alpha_intra, alpha_inter1 = alpha_inter1, alpha_inter2 = alpha_inter2))),
      data = sample_dta, 
      start = inits,
      lower = lowers,
      upper = uppers,
      algorithm = "port",
      control = nls_ctrl
    )
  }
  
  # Initials for HOI models
  inits = c(lambda = 3, alpha_intra = 0.2, alpha_inter1 = 0.2, alpha_inter2 = 0.2, 
            beta_112 = 0.2, beta_113 = 0.2, beta_123 = 0.2)
  lowers = c(lambda = 1, alpha_intra = -1, alpha_inter1 = -1, alpha_inter2 = -1, 
             beta_112 = -1, beta_113 = -1, beta_123 = -1)
  uppers = c(lambda = 5, alpha_intra = 1, alpha_inter1 = 1, alpha_inter2 = 1, 
             beta_112 = 1, beta_113 = 1, beta_123 = 1)
  
  if (spp == 1) {
    
    m_hoi = nlsLM(
      (N1tp1) ~ (N1t*ri_model_3_hoi(N1t, N2t, N3t, 
                                    para = list(lambda = lambda, 
                                                alpha_intra = alpha_intra, alpha_inter1 = alpha_inter1, alpha_inter2 = alpha_inter2, 
                                                beta_112 = beta_112, beta_113 = beta_113, beta_123 = beta_123))),
      data = sample_dta, 
      start = inits,
      lower = lowers,
      upper = uppers,
      algorithm = "port",
      control = nls_ctrl
    )
  } 
  
  if (spp == 2) {
    
    m_hoi = nlsLM(
      (N2tp1) ~ (N2t*ri_model_3_hoi(N2t, N1t, N3t, 
                                    para = list(lambda = lambda, 
                                                alpha_intra = alpha_intra, alpha_inter1 = alpha_inter1, alpha_inter2 = alpha_inter2, 
                                                beta_112 = beta_112, beta_113 = beta_113, beta_123 = beta_123))),
      data = sample_dta, 
      start = inits,
      lower = lowers,
      upper = uppers,
      algorithm = "port",
      control = nls_ctrl
    )
  }
  
  return(list(m_pairwise, m_hoi, c(AIC(m_pairwise), AIC(m_hoi))))
}

calculate_rmse = function(sample_dta) {
  
  out = tibble()
  
  for (spp in 1:2) {
    
    models = compare_fit(sample_dta, spp)
    rmse_both = c(rmse(pull(sample_dta, spp+3), predict(models[[1]])), 
                  rmse(pull(sample_dta, spp+3), predict(models[[2]])))
    AIC_both = models[[3]]
    out = rbind(out, tibble(Species = spp, RMSE = rmse_both, AIC = AIC_both, Model = c("Pairwise", "HOI"), 
                            type = unique(sample_dta$type), case = unique(sample_dta$case)))
    
  }
  
  return(out)
}

