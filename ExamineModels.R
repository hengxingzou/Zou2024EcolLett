source("Functions_PSF.R")
source("ModelFitting.R")


########## All Pairwise Phenological Models ##########


bh_model_2 = function(N_self, N_other, para) {
  with(para, {
    t1 = alpha_intra*N_self
    t2 = alpha_inter*N_other
    lambda/(1+t1+t2)
  })
}

lv_model_2 = function(N_self, N_other, para) {
  with(para, {
    t1 = alpha_intra*N_self
    t2 = alpha_inter*N_other
    lambda*(1-t1-t2)
  })
}

ri_model_2 = function(N_self, N_other, para) {
  with(para, {
    t1 = alpha_intra*N_self
    t2 = alpha_inter*N_other
    lambda*exp(-t1-t2)
  })
}

lw_model_2 = function(N_self, N_other, para) {
  with(para, {
    t1 = N_self^alpha_intra
    t2 = N_other^alpha_inter
    lambda/(1+t1+t2)
  })
}


########## Read Data ##########


data_dir = "Data/Standard/"

# Response surface, for model fitting, use low overlap as an example

ij_data = read_csv(paste0(data_dir, "ij_al_rs.csv"))

# Arrival times of i and j

tau = 12
lc_len = 6

max_arriv = tau-lc_len-1

arriv_combns_ij = expand_grid(seq(0, max_arriv, 1), seq(0, max_arriv, 1), 0)

# Initial parameters for Lotka-Volterra, Beverton-Holt, Ricker, Law-Watkinson

inits = c(lambda = 3, alpha_intra = 0.2, alpha_inter = 0.2)
lowers = c(lambda = 1, alpha_intra = -3, alpha_inter = -3)
uppers = c(lambda = 5, alpha_intra = 3, alpha_inter = 3)

nls_ctrl = nls.control(maxiter = 1000, tol = 1e-8, minFactor = (1/10)*(1/1024),
                       printEval = FALSE, warnOnly = FALSE)

# Looping through the range of pj

output = tibble()

for (row in 1:nrow(arriv_combns_ij)) {

  # Select one arrival time sequence
  timeseq = arriv_combns_ij %>% slice(row) %>% as.numeric()
  
  # Preparing data
  sample_dta = transform_dt(ij_data, timeseq)
  
  # Lotka-Volterra
  
  lv_1 = nlsLM(
    (N1tp1) ~ (N1t*lv_model_2(N1t, N2t,
                              para = list(lambda = lambda, alpha_intra = alpha_intra, alpha_inter = alpha_inter))),
    data = sample_dta, 
    start = inits,
    lower = lowers,
    upper = uppers,
    algorithm = "port",
    control = nls_ctrl
  )
  
  # Beverton-Holt
  
  bh_1 = nlsLM(
    (N1tp1) ~ (N1t*bh_model_2(N1t, N2t,
                              para = list(lambda = lambda, alpha_intra = alpha_intra, alpha_inter = alpha_inter))),
    data = sample_dta, 
    start = inits,
    lower = lowers,
    upper = uppers,
    algorithm = "port",
    control = nls_ctrl
  )
  
  # Ricker
  
  ri_1 = nlsLM(
    (N1tp1) ~ (N1t*ri_model_2(N1t, N2t,
                              para = list(lambda = lambda, alpha_intra = alpha_intra, alpha_inter = alpha_inter))),
    data = sample_dta, 
    start = inits,
    lower = lowers,
    upper = uppers,
    algorithm = "port",
    control = nls_ctrl
  )
  
  # Law-Watkinson
  
  lw_1 = nlsLM(
    (N1tp1) ~ (N1t*lw_model_2(N1t, N2t,
                              para = list(lambda = lambda, alpha_intra = alpha_intra, alpha_inter = alpha_inter))),
    data = sample_dta, 
    start = inits,
    lower = lowers,
    upper = uppers,
    algorithm = "port",
    control = nls_ctrl
  )
  
  output = rbind(output, 
                 tibble(model = c("lv", "bh", "ri", "lw"), 
                        AIC = AIC(lv_1, bh_1, ri_1, lw_1)$AIC, 
                        pi = timeseq[1], pj = timeseq[2]))
  
}

# Visualize AIC values

AIC_comparison = 
  output %>% 
  filter(pi == 0 | pj == 0) %>% 
  ggplot(aes(x = (pj-pi), y = AIC, color = model)) + 
  geom_point(size = 3) + 
  scale_color_manual(name = "Model type",
                     values = color_scheme_5, 
                     labels = c("Beverton-Holt", "Lotka-Volterra", "Law-Watkinson", "Ricker")) + 
  scale_x_continuous(name = expression(paste(Delta, p[ij])), 
                     breaks = seq(-5, 5, 1)) +
  ylab("AIC") +
  theme_bw() + 
  theme(axis.text.y = element_text(size = 15), axis.text.x = element_text(size = 15), 
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), 
        legend.text = element_text(size = 17), legend.title = element_text(size = 20), 
        panel.grid.minor = element_blank(), 
        plot.title = element_text(size = 20),
        strip.text.x = element_blank(), strip.background = element_blank())

AIC_comparison

ggsave("Plots/AICComparison.png", plot = AIC_comparison,
       device = "png", width = 2400, height = 2400, units = "px")
