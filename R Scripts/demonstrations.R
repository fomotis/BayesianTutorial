library(rstan)
library(coda)
library(bayesplot)
library(bridgesampling)
library(tidyverse)
library(ggpubr)
library(nlme)
library(deSolve)
library(rbokeh)
library(plotly)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores(logical = FALSE))

################ Data ##########
combined_data <- read.csv("Data/combined_data.csv", header = T, na.strings = ".", 
                          stringsAsFactors = FALSE ) %>% 
  mutate(Temperature = factor(Temperature))
b5_mono_data <- combined_data %>% 
  dplyr::filter(Experiment == "Mono" & 
                  Strain == "B5" & Time > 1)
b5_mono_data <- b5_mono_data[order(b5_mono_data$Sample.ID2, b5_mono_data$Time), ]
#turn temperature to factor
b5_mono_data$Temperature <- factor(b5_mono_data$Temperature)

#creating numeric group for random effects
b5_mono_data <- b5_mono_data %>% dplyr::mutate(Groups2 = rep(1:length(unique(b5_mono_data$Sample.ID2)),
                                                             each = 8))

#indicator for treatment = 20 or 22
b5_mono_data <- b5_mono_data %>% 
  dplyr::mutate(T_2022 = ifelse(
    Treatment %in% c(20, 22), 1, 0)
  )
#remove NAs from the data
b5_mono_data_nona <- b5_mono_data[!is.na(b5_mono_data$m_fschlin), ]
#obtain the treatment design matrix
trt_mat <- as.matrix(b5_mono_data_nona[, c("T_18", "T_20", "T_22")])

################# function to fit model using stan
#mod is the stan model file, mod_data is the data needed to fit the model
stan_mod_func <- function(mod, mod_data, init = NULL, ...) {
  
  if(is.null(init)) {
    
    ret_mod <- stanc(mod) # formulate model from stan file
    sm <- stan_model(stanc_ret = ret_mod, verbose = FALSE)
    sm_fit <- sampling(sm, data = mod_data, iter = 3000, thin = 1, 
                       control = list(max_treedepth = 13, adapt_delta = 0.98), 
                       chains = 2)
    
  } else {
    
    ret_mod <- stanc(mod) # formulate model from stan file
    sm <- stan_model(stanc_ret = ret_mod, verbose = FALSE)
    sm_fit <- sampling(sm, data = mod_data, iter = 3000, thin = 1, 
                       control = list(max_treedepth = 13, adapt_delta = 0.98), init = init, 
                       chains = 2)
    
  }
  
  return(sm_fit)
}

### function to return the posterior mean of a parameter
par_postsums <- function(stan_mod, par = "yhat") {
  
  as.data.frame(summary(stan_mod, pars = par)$summary)
  
}

par_postsums2 <- function(stan_mod, pars = c("sigma","alpha", "beta", "lsigma")) {
  
  y <- par_postsums(stan_mod = stan_mod, par = pars)
  z1 <- as.data.frame(y[, c("mean", "2.5%", "97.5%")])
  names(z1) <- c("mean", "LCI", "UCI")
  cbind(z1, Parameters = row.names(z1))
  
}

########### function to extract residuals (should be a matrix) and plot them for bivariate normal model
error_plot <- function(stan_mod_obj, resid_name = "epsilon", 
                       number_rows, group = NULL, 
                       colnames = c("error_Trait", "error_abundance")) {
  
  errors <- matrix(par_postsums(stan_mod_obj, resid_name)[, "mean"], 
                   nrow = number_rows, byrow = TRUE)
  errors <- as.data.frame(errors)
  names(errors) <- colnames
  
  if(!is.null(group)) {
    
    errors$Group <- group
    
    p <- errors %>% ggplot(aes(x = get(colnames[2], errors), y = get(colnames[1], errors), 
                               group = Group, color = Group)) + 
      geom_point(size = 5) + 
      theme_bw() + geom_smooth(aes(x = get(colnames[2], errors), y = get(colnames[1], errors)), 
                               method = "lm", inherit.aes = FALSE, se = FALSE,
                               color = "black", data = errors) + 
      labs(x = expression(epsilon[abundance]), y = expression(epsilon[trait])) #+
    #geom_smooth(aes(x = get(colnames[2], errors), y = get(colnames[1], errors), 
    #                group = Group, color = Group), 
    #            method = "lm", inherit.aes = FALSE, se = FALSE, 
    #            data = errors)
    
  } else {
    
    p <- errors %>% ggplot(aes(x = get(colnames[2], errors), y = get(colnames[1], errors))) + 
      geom_point(size = 5) + 
      theme_bw() + geom_smooth(aes(x = get(colnames[2], errors), y = get(colnames[1], errors)), 
                               method = "lm", inherit.aes = FALSE, se = FALSE, color = "black", 
                               data = errors) + 
      labs(x = expression(epsilon[abundance]), y = expression(epsilon[trait]))
    
  }
  
  return(p)
}

### credible interval plots for parameter estimates
ciplot <- function(pdata, x, y, x_lab, y_lab, x_scale = NULL) {
  
  pd <- position_dodge(0.1)
  if(!is.null(x_scale)) {
    
    pdata %>% 
      ggplot(aes(x = get(x, pdata), y = get(y, pdata))) + 
      geom_errorbar(aes(ymin = LCI, ymax = UCI), 
                    color="black", width = 0.1, 
                    position = pd, data = pdata, 
                    inherit.aes = TRUE) +
      geom_point(size = 4, position = pd, shape = 21, 
                 fill = "red3", color = "black") + 
      theme_bw() + 
      labs(x = x_lab, y = y_lab) + 
      scale_x_discrete(labels = x_scale)
    
  } else {
    
    pdata %>% 
      ggplot(aes(x = get(x, pdata), y = get(y, pdata))) + 
      geom_point(size = 4, position = pd, shape = 21, 
                 fill = "red3", color = "black") + 
      geom_errorbar(aes(ymin = LCI, ymax = UCI), 
                    color="black", width = 0.1, 
                    position = pd, data = pdata, 
                    inherit.aes = TRUE) +
      theme_bw() + 
      labs(x = x_lab, y = y_lab)
  }
  
}

###### function to plot predicted values

pred_plot <- function(orig_data, predicted_values, x, y, x_lab, y_lab) {
  
  pdata <- cbind(orig_data, Pred = predicted_values)
  pdata2 <- pdata %>% filter(Replicate == 2)
  
  ggplot(data = pdata, aes(x = get(x, pdata), y = get(y, pdata), group = interaction(Temperature, Replicate), 
                           color = Temperature)) + 
    geom_point(size = 4) + 
    geom_line(aes(x = get(x, pdata2), y = get("Pred", pdata2), group = interaction(Temperature, Replicate), 
                  color = Temperature), size = 1.5, 
              linetype= "dashed", data = pdata2, inherit.aes = FALSE) +
    theme_bw() + color_palette("Dark2") +
    scale_x_continuous(breaks = seq(0, 30, by = 4)) +
    #scale_y_continuous(breaks = seq(4, 8, by = 0.5)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.9, vjust = 0.5), legend.position = "top") +
    labs(x = x_lab, y = y_lab)
}

#### abundance plot

abd_plot <- combined_data %>% filter(Strain == "B5" & Experiment == "Mono" & Time > 1) %>%
  ggplot(aes(x = Time, y = ldensity, group = interaction(Temperature, Replicate), color = Temperature)) +
  geom_point(size = 4) +
  geom_line(size = 1.5) +
  theme_minimal() +
  color_palette("Dark2") +
  scale_x_continuous(breaks = seq(0, 30, by = 4)) +
  scale_y_continuous(breaks = seq(0, 14, by = 1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.01, vjust = 0.5),
        legend.position = "right", 
        axis.line = element_line(size = 0.5)) +
  labs(x = "Time", y = "log(Abundance)")
abds <- ggplotly(abd_plot)
withr::with_dir("html", htmlwidgets::saveWidget(abds, file = "abd.html"))


########## Frequentist analysis for the Verhuls Model

SSLVE <- function(time, N0, r, K) {
  
  numerator <- K * N0
  denominator <- N0 + ((K - N0) * exp(-r * time))
  mu <- numerator / denominator
  return(mu)
  
}

ver_model <- gnls(ldensity ~ SSLVE(time = Time, N0, r, K), 
                data = b5_mono_data, 
                params = list(N0 ~ 1, r ~ Temperature, K ~ Temperature),
                start = c(4.2, 0.22, 0.23, 0.36, 9, 0, 0),
                correlation = NULL,
                na.action = na.omit
)

ver_model_coefs <- coef(ver_model)
ver_model_vacov <- vcov(ver_model)
ver_model_sum <- summary(ver_model)

names(ver_model_coefs) <- c("N0", "rInt", "r20", "r22", "KInt", "K20", "K22")
row.names(ver_model_vacov) <- c("N0", "rInt", "r20", "r22", "KInt", "K20", "K22")
colnames(ver_model_vacov) <- c("N0", "rInt", "r20", "r22", "KInt", "K20", "K22")
interestr <- c("rInt", paste("(rInt", c("r20)","r22)"), sep = "+"))
interestK <- c("KInt", paste("(KInt", c("K20)","K22)"), sep = "+"))
interestA <- paste(interestr,interestK, sep = "/")
interest <- c(interestr, interestK, interestA)
estimates <- t(sapply(interest, function(.x) {
  
  car::deltaMethod(ver_model_coefs, g. = .x, 
                   vcov. = ver_model_vacov, 
                   level = 0.95)
  
}))
colnames(estimates) <- c("Estimate", "SE", "LCI", "UCI")
  estimates <- as.data.frame(sapply(as.data.frame(estimates), as.numeric)) %>% 
  mutate(Parameters = c("r18", "r20", "r22", "K18", "K20", "K22", 
                        "A18", "A20", "A22"))

estimates$Type <- rep(c("growth rate", "carrying capacity", "intraspecific effect"), each = 3)

pr <- ciplot(estimates %>% filter(Type == "growth rate"), x = "Parameters", y = "Estimate",
             x_lab = "Parameters", y_lab = "Estimates", x_scale = c(expression(r[18]), expression(r[20]), 
                                                                    expression(r[22]))) +
  facet_wrap(.~Type, scales = "free_y")
pk <- ciplot(estimates %>% filter(Type == "carrying capacity"), x = "Parameters", y = "Estimate",
             x_lab = "Parameters", y_lab = "", x_scale = c(expression(K[18]), expression(K[20]), 
                                                                    expression(K[22]))) +
  facet_wrap(.~Type, scales = "free_y")
pa <- ciplot(estimates %>% filter(Type == "intraspecific effect"), x = "Parameters", y = "Estimate",
             x_lab = "Parameters", y_lab = "", x_scale = c(expression(A[18]), expression(A[20]), 
                                                                    expression(A[22]))) +
  facet_wrap(.~Type, scales = "free_y")

p1 <- ggarrange(plotlist = list(pr, pk, pa), ncol = 3, nrow = 1,
                align = "hv", common.legend = TRUE)
#ggplotly(p1)

###rbokeh
p1_list = vector("list", 3)
j<-1
for(i in c("growth rate", "carrying capacity", "intraspecific effect")) {
  
  p1_list[[j]] <- figure(xlim = estimates$Parameters[estimates$Type == i], 
         width = 220, height = 350, tools = c("pan", "wheel_zoom", "box_zoom", "box_select", "reset")) %>%
    ly_segments(Parameters, LCI, Parameters, UCI, data = estimates[estimates$Type == i,], 
                color = "black", width = 2) %>%
    ly_points(Parameters, Estimate, glyph = 16, data = estimates[estimates$Type == i,], 
              color = "red", size = 20, hover = c(Parameters, Estimate)) %>%
    y_axis(label = "Estimates")
  j <<- j+1
}

p1s <- grid_plot(p1_list, ncol = 3, same_axes = F, link_data = F)

withr::with_dir("html", htmlwidgets::saveWidget(p1s, file = "est_freq.html"))

##################### Bayesian Analysis
set.seed(1992)

#covraite matrix for the trait mean function
X_obs <- b5_mono_data_nona %>% 
  model.matrix(~I(Time^0.5):T_18 + I(Time^3):T_2022 + 
                 I((Time^3) * log(Time)):T_2022 +
                 T_18 + T_20 + T_22, data = . )

verhulst_data <- list(
  N_obs = nrow(X_obs),
  y_obs = b5_mono_data_nona[, "ldensity"],
  trt = b5_mono_data_nona[, c("T_18", "T_20", "T_22")],
  time = b5_mono_data_nona$Time,
  n_trt = ncol(trt_mat) 
)
stan_verhulst <- stan_mod_func(mod = "stan_models/verhulst.stan", mod_data = verhulst_data, 
                               init = list( 
                                 list(r = c(0.22, 0.23, 0.35)/2, K = c(9.72, 11.50, 10.99)/2), 
                                 list(r = c(0.22, 0.23, 0.35)*1.5, K = c(9.72, 11.50, 10.99)*1.5)
                               ))

ver_pdata <- par_postsums2(stan_verhulst, pars = c("r", "K", "A"))
ver_pdata$Parameters <- c(paste0("r", c(18, 20, 22)),
                          paste0("K", c(18, 20, 22)),
                          paste0("A", c(18, 20, 22)))

ver_pdata$Type <- rep(c("growth rate", "carrying capacity", "intraspecific effect"), each = 3)

### plotting Bayesian 
p2_list = vector("list", 3)
j<-1
for(i in c("growth rate", "carrying capacity", "intraspecific effect")) {
  
  p2_list[[j]] <- figure(xlim = ver_pdata$Parameters[ver_pdata$Type == i], 
                         width = 220, height = 350, tools = c("pan", "wheel_zoom", "box_zoom", "box_select", "reset")) %>%
    ly_segments(Parameters, LCI, Parameters, UCI, data = ver_pdata[ver_pdata$Type == i,], 
                color = "black", width = 2) %>%
    ly_points(Parameters, mean, glyph = 16, data = ver_pdata[ver_pdata$Type == i,], 
              color = "red", size = 20, hover = c(Parameters, mean)) %>%
    y_axis(label = "Estimates")
  j <<- j+1
}

p2s <- grid_plot(p2_list, ncol = 3, same_axes = F, link_data = F)
withr::with_dir("html", htmlwidgets::saveWidget(p2s, file = "est_bay.html"))

################ plotting both results together ##########
bay_fre_est <- rbind(estimates %>% dplyr::select(-SE), 
                     ver_pdata %>% dplyr::rename(Estimate = mean)
                     ) %>% 
              mutate(Method = rep(c("Frequentist", "Bayesian"), each = 9)
              )

pd <- position_dodge(1.0)
p3_list <- vector("list", 3)
j <- 1
for(i in c("growth rate", "carrying capacity", "intraspecific effect")) {
  
  if(j == 3) {
    
    p3_list[[j]] <- bay_fre_est %>% 
      filter(Type == i) %>%
      ggplot(aes(x = Parameters, y = Estimate, group = Method, colour = Method)) + 
      geom_errorbar(aes(ymin = LCI, ymax = UCI), colour = "black", 
                    width=.1, position = pd, size = 1.5) + 
      geom_point(position = pd, size = 9) + 
      theme_minimal() + 
      theme(axis.line = element_line(size = 0.5), 
            legend.position = "none", 
            legend.direction = "vertical")
    
  } else {
    
    p3_list[[j]] <- bay_fre_est %>% 
      filter(Type == i) %>%
      ggplot(aes(x = Parameters, y = Estimate, group = Method, colour = Method)) + 
      geom_errorbar(aes(ymin = LCI, ymax = UCI), colour = "black", 
                    width=.1, position = pd, size = 1.5) + 
      geom_point(position = pd, size = 9) + 
      theme_minimal() + 
      theme(legend.position = "none", 
            axis.line = element_line(size = 0.5)
              )
  }
  
  j <<- j+1
}
p3s_gga <- ggarrange(plotlist = p3_list, nrow = 1, ncol = 3)
p3s <- subplot(p3_list, nrows = 1)
withr::with_dir("html", htmlwidgets::saveWidget(p3s, file = "est_bayfreq.html"))


################# Lotka-Voltera for 2 species example

biculture_data <- combined_data %>% 
  dplyr::filter(Experiment == "BI") %>% 
  mutate(Temperature = factor(Temperature), 
         Temperature2 = case_when( 
           Temperature == 18 ~ 1,
           Temperature == 20 ~ 2,
           Temperature == 22 ~ 3
           
         ),
         Species = case_when(
           Strain == "B4" ~ 1,
           Strain == "B5" ~ 2
         )
         
  )

biculture_data_nona <- biculture_data %>% filter(Time > 1, !is.na(ldensity))

biculture_data_wide <- biculture_data_nona %>% pivot_wider(names_from = Strain, 
                                                      id_cols = c("Date2", "Sample.ID2", "Date", "Time",
                                                                  "Treatment", "Replicate", "Temperature",
                                                                  "Temperature2", "T_18", "T_20", "T_22", "Groups"
                                                      ),
                                                      values_from = c(ldensity, m_fschlin, m_yelbhlin, 
                                                                      m_redbhlin, v_fschlin, v_yelbhlin, 
                                                                      Number)
)

biculture_data %>% dplyr::filter( Time > 1) %>% 
  ggplot(aes(x = Time, y = ldensity, group = interaction(Strain, Temperature, Replicate), 
             color = Temperature)) + 
  geom_point(size = 4) + 
  geom_line(aes(linetype = Strain), size = 1.5) +
  theme_minimal() + 
  color_palette("Dark2") +
  scale_x_continuous(breaks = seq(0, 30, by = 4)) +
  scale_y_continuous(breaks = seq(0, 14, by = 1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.9, vjust = 0.5), 
        legend.position = "top") +
  labs(x = "Time", y = "Log(Abundance)")

### frequentist analysis

# 2 species LVE
LVE <- function(time, y, parms, ...) {
  
  dy <- rep(0, length(y))
  dy[1] <-  y[1] * parms$r1 * (1 - (parms$alpha11*y[1] + parms$alpha12*y[2]) )
  dy[2] <-  y[2] * parms$r2 * (1 - (parms$alpha21*y[1] + parms$alpha22*y[2]) )
  
  return(list(dy))
}

ts <- sort(unique(c(seq(min(biculture_data$Time), 50, by = 1), 
                    biculture_data$Time)))
#a test
out_solu <- ode(y = c(7, 6), 
    times = ts, 
    func = LVE,
    parms = list(r1 = 0.2, r2 = 0.3, alpha11 = 0.06, alpha21 = 0.02, alpha12 = 0.02, alpha22 = 0.06)
     )
lapply(c(18, 20, 22), function(i) {
  
  ddata <- biculture_data_wide %>% filter(Temperature == 18)
  # function to obtain sum of squares residual
  SSR <- function(parms) {
    
    # mapping the parameters to their respective positions
    r1 <- parms[1]
    r2 <- parms[2]
    alpha11 <- parms[3]
    alpha21 <- parms[4]
    alpha12 <- parms[5]
    alpha22 <- parms[6]
    sigma_b4 <- parms[7]
    sigma_b5 <- parms[8]
    rho <- 0 #parms[9]
    N0_1 <- 3.6 #parms[10]
    N0_2 <- 3.5 #parms[11]
    
    #a time vector
    time_seq <- ts
    # initial values
    y_inits <- c(N0_1, N0_2)
    
    #solve the ode
    out_solu <- ode(y = y_inits, times = time_seq, func = LVE, 
                    parms = list(r1 = r1, r2 = r2, alpha11 = alpha11, alpha12 = alpha12, 
                                 alpha21 = alpha21, alpha22 = alpha22)
    )
    out_dataframe <- as.data.frame(out_solu)
    out_dataframe <- out_dataframe[out_dataframe$time %in% unique(biculture_data$Time), ]
    names(out_dataframe) <- c("time", "ldensity_B4_pred", "ldensity_B5_pred")
    
    dd <- merge(ddata, out_dataframe, by.x = c("Time"), by.y = c("time"))
    Sigma <- matrix(c(sigma_b4^2, sigma_b4*sigma_b5*rho, sigma_b4*sigma_b5*rho, sigma_b5^2),
                    nrow = 2, ncol = 2, byrow = T)
     ll <- numeric(nrow(dd))
     for(i in 1:nrow(dd)) {
       
       ll[i] <- dmvnorm(x = dd[i, c("ldensity_B4", "ldensity_B5")], 
                        mean = c(dd$ldensity_B4_pred[i], dd$ldensity_B5_pred[i]),  
                        sigma = Sigma, log = T) 
     }
    
    return(-sum(ll))
  }
  
  
  parms_start <- c(r1 = 0.7, r2 = 0.4, alpha11 = 0.06, 
                   alpha12 = 0.02, alpha21 = 0.02, alpha22 = 0.06, 
                   sigma_b4 = 0.70, 
                   sigma_b5 = 0.45, 
                   rho = 0.5, 
                   N0_1 = 7, 
                   N0_2 = 6)
  lve_freq <- optim(par = parms_start, fn = SSR, 
                    method = "Nelder-Mead",
                     control = list(maxit = 100), 
                    lower = c(r1 = 0, r2 = 0, 
                              alpha11 = 0, 
                              alpha12 = 0, alpha21 = 0, 
                              alpha22 = 0, sigma_b4 = 0, 
                              sigma_b5 = 0, rho = -1, 
                              N0_1 = 0, N0_2 = 0),
                    upper= c(r1 = Inf, r2 = Inf, 
                             alpha11 = Inf, 
                             alpha12 = Inf, alpha21 = 0, 
                             alpha22 = Inf, sigma_b4 = Inf, 
                             sigma_b5 = Inf, rho = 1, 
                             N0_1 = Inf, N0_2 = Inf)
                    )
  
  lve_sum <- summary(lve_freq)
  UCI <- lve_sum$coefficients[, 1] + lve_sum$coefficients[, 2]
  LCI <- lve_sum$coefficients[, 1] - lve_sum$coefficients[, 2]
  Parameters <- c("r1", "r2", "alpha11", "alpha12", "alpha21", "alpha22")
  data.frame(Parameters = Parameters, Estimates = lve_sum$coefficients[, 1], LCI = LCI, UCI = UCI)
  
})

#### Bayesian Analysis

biculture_data_wide_nona <- biculture_data_wide[!is.na(biculture_data_wide$ldensity_B4), ]
results_backup <- vector("list", 3)
j <- 1
results_treatment <- lapply(c(18, 20, 22), function(i) {
  
  print(i)
  
  bi_datas <- biculture_data_wide_nona[biculture_data_wide_nona$Treatment == i, ]
  bi_data1 <- list(
    N_obs = nrow(bi_datas),
    T = length(ts),
    N = as.matrix(bi_datas[, c("ldensity_B4", "ldensity_B5")]),
    t0 = 0,
    ts = ts,
    time_obs = bi_datas$Time,
    nsp = 2
  )
  if( i == 18) {
    
    N01 <- c(3.6, 3.5) 
    r <- c(0.43, 0.35)
    alpha11 <- 0.07; alpha22 <- 0.07
    alpha12 <- 0.02; alpha21 <- 0.02
    
  } else if (i == 20) {
    
    N01 <- c(6.5, 5.6) 
    r <- c(0.20, 0.21)
    alpha11 <- 0.06; alpha22 <- 0.06
    alpha12 <- 0.02; alpha21 <- 0.02
    
  } else if (i == 22) {
    
    
    N01 <- c(6.1, 5.3) 
    r <- c(0.16, 0.16)
    alpha11 <- 0.06; alpha22 <- 0.06
    alpha12 <- 0.02; alpha21 <- 0.01
    
  }
  
  
  st_mod <- stan_mod_func(mod = "stan_models/generalisedLVE.stan", mod_data = bi_data1, 
                       init = list( 
                         list(r = r, alpha11 = alpha11, alpha22 = alpha22, 
                              alpha12 = alpha12, alpha21 = alpha21, N0 = N01), 
                         list(r = r*1.5, alpha11 = alpha11*1.5, alpha22 = alpha22*1.5, 
                              alpha12 = alpha12*1.5, alpha21 = alpha21*1.5, 
                              N0 = N01*1.5))
  )
  
  results_backup[[j]] <- st_mod
  j <<- j + 1
  print(par_postsums2(st_mod, c("r", "N0", "alpha11", "alpha22", "alpha12", "alpha21", "cor", "sigma")))
  return(st_mod)
  
})

results_18 <- par_postsums2(results_treatment[[1]], 
                            c("r", "N0", "alpha11", "alpha22", "alpha12", "alpha21", "cor")) %>%
  mutate(Parameters = c("r18B4", "r18B5", "N0B4", "N0B5", "alpha11", "alpha22", "alpha12", "alpha21", "rho"), 
         Type = rep(c("growth rate", "N0", "alphas", "rho"), times = c(2, 2, 4, 1) ), 
         Temperature = rep(18, times = 9))
results_20 <- par_postsums2(results_treatment[[2]], 
                            c("r", "N0", "alpha11", "alpha22", "alpha12", "alpha21", "cor")) %>%
  mutate(Parameters = c("r18B4", "r18B5", "N0B4", "N0B5", "alpha11", "alpha22", "alpha12", "alpha21", "rho"), 
         Type = rep(c("growth rate", "N0", "alphas", "rho"), times = c(2, 2, 4, 1) ),
         Temperature = rep(20, times = 9))
results_22 <- par_postsums2(results_treatment[[3]], 
                            c("r", "N0", "alpha11", "alpha22", "alpha12", "alpha21", "cor")) %>%
  mutate(Parameters = c("r18B4", "r18B5", "N0B4", "N0B5", "alpha11", "alpha22", "alpha12", "alpha21", "rho"), 
         Type = rep(c("growth rate", "N0", "alphas", "rho"), times = c(2, 2, 4, 1) ), 
         Temperature = rep(22, times = 9))
results_all <- rbind(results_18, results_20, results_22)
results_all$Temperature <- factor(results_all$Temperature)
pd <- position_dodge(1.0)
p4_list <- vector("list", 3)
j <- 1
for(i in c("growth rate", "N0", "alphas", "rho")) {
  
  if(j == 4) {
    
    p4_list[[j]] <- results_all %>% 
      filter(Type == i) %>%
      ggplot(aes(x = Parameters, y = mean, group = Temperature, colour = Temperature)) + 
      geom_errorbar(aes(ymin = LCI, ymax = UCI), colour = "black", 
                    width=.1, position = pd, size = 1.5) + 
      geom_point(position = pd, size = 9) + 
      theme_minimal() + 
      theme(axis.line = element_line(size = 0.5), 
            legend.position = "none", 
            legend.direction = "vertical") + 
      labs(y = "Posterior Mean +- Posterior SD")
    
  } else {
    
    p4_list[[j]] <- results_all %>% 
      filter(Type == i) %>%
      ggplot(aes(x = Parameters, y = mean, group = Temperature, colour = Temperature)) + 
      geom_errorbar(aes(ymin = LCI, ymax = UCI), colour = "black", 
                    width=.1, position = pd, size = 1.5) + 
      geom_point(position = pd, size = 9) + 
      theme_minimal() + 
      theme(legend.position = "none", 
            axis.line = element_line(size = 0.5)
      ) + 
      labs(y = "Posterior Mean +- Posterior SD")
  }
  
  j <<- j+1
}

p4s <- subplot(p4_list, nrows = 2)
withr::with_dir("html", htmlwidgets::saveWidget(p4s, file = "est_bayLVE.html"))



############### model comparison

stan_gompertz <- stan_mod_func(mod = "stan_models/gompertz.stan", mod_data = verhulst_data, 
                               init = list( 
                                 list(r = c(0.22, 0.23, 0.35)/2, c = c(0.72, 0.50, 1.99)/2), 
                                 list(r = c(0.22, 0.23, 0.35)*1.5, c = c(1.72, 3.50, 3.99)*1.5 )
                               ))

save.image("RImage/demonstration.RData")

