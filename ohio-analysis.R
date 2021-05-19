args = commandArgs(trailingOnly=TRUE)

main_folder <- "./ohio"
prop_file <- "fast.properties"

if (0 < length(args)) {
  prop_file <- file.path(main_folder, args[1])
}

# =============================================================================
# Libraries

library(doRNG)
library(foreach)
library(doParallel)
library(boot)
library(pracma)
library(properties)

# tidyr has to be imported before magrittr so extract() from the latter is used
library(tidyr)

library(dplyr)
library(plyr)
library(reshape2)
library(magrittr)

library(ggplot2)
theme_set(theme_bw())

library(pomp)
stopifnot(packageVersion("pomp")>="2")

# =============================================================================
# PARAMETERS

prop <- read.properties(prop_file)

output_name <- prop$output_name
output_folder <- file.path(main_folder, output_name)
if (!dir.exists(output_folder)) dir.create(output_folder)
cooking_folder <- file.path(main_folder, output_name, "bake")
if (!dir.exists(cooking_folder)) dir.create(cooking_folder)
plotting_folder <- file.path(main_folder, output_name, "plots")
if (!dir.exists(plotting_folder)) dir.create(plotting_folder)
code_folder <- file.path(main_folder, "code")
if (!dir.exists(code_folder)) dir.create(code_folder)
file_name <- "snippets"

pop_size <- as.integer(as.numeric(prop$pop_size) * as.numeric(prop$susc_frac))
exp0 <- as.integer(prop$exp0)
inf0 <- as.integer(prop$inf0)
rec0 <- as.integer(prop$rec0)

first_week <- as.integer(prop$first_week)
last_week <- as.integer(prop$last_week)
start_aware_week <- as.integer(prop$start_aware_week)

time_step <-as.numeric(prop$time_step)
num_sims <- as.integer(prop$num_sims)

free_param_names <- c("beta", "rho", "psi", "w")
free_param_box <- rbind(
  beta = c(as.numeric(prop$bounds_beta_min), as.numeric(prop$bounds_beta_max)),
  rho = c(as.numeric(prop$bounds_rho_min), as.numeric(prop$bounds_rho_max)),
  psi = c(as.numeric(prop$bounds_psi_min), as.numeric(prop$bounds_psi_max)),
  w = c(as.numeric(prop$bounds_w_min), as.numeric(prop$bounds_w_max))
)

log_trans_params <- c("beta", "psi")
logit_trans_params <- c("rho", "w")

fixed_param_names <- c("pop", "S_0", "E_0", "I_0", "R_0", "sigma", "gamma", "intervention")
fixed_param_values <- c(pop=pop_size, 
                        S_0=1-(exp0+inf0+rec0)/pop_size, E_0=exp0/pop_size, I_0=inf0/pop_size, R_0=rec0/pop_size, 
                        sigma=as.numeric(prop$sigma), gamma=as.numeric(prop$gamma), intervention=start_aware_week)

all_param_names <- c(free_param_names, fixed_param_names)

# Random seeds, keep unchanged to ensure reproducibilty of results
test_mle_seed <- as.integer(prop$test_mle_seed)
full_mle_seed <- as.integer(prop$full_mle_seed)
global_search_seed <- as.integer(prop$global_search_seed)
test_sim_seed <- as.integer(prop$test_sim_seed)
full_sim_seed <- as.integer(prop$full_sim_seed)
mcap_plik_seed <- as.integer(prop$mcap_plik_seed)

# =============================================================================
# OHIO DATASET

us_data <- read.csv(file.path(main_folder, "mumps_data_per_state.csv"))
ohio_data <- us_data[us_data$Reporting.Area == 'OHIO', ]
ohio_data[is.na(ohio_data)] <- 0

ohio_data

ggplot(data=ohio_data, aes(x=MMWR.Week, y=Mumps..Current.week, group=1)) +
  geom_line()

ohio_data <- data.frame('time' = ohio_data$MMWR.Week, 
                        'cases' = ohio_data$Mumps..Current.week)

ohio_data <- subset(ohio_data, first_week <= time & time <= last_week)

ggplot(data=ohio_data, aes(x=time, y=cases, group=1)) +
  geom_line() + ggtitle('Outbreak at Ohio State University')

ggsave(file.path(plotting_folder, "0-weekly-cases.pdf"))  

# =============================================================================
# COMPARTMENTAL MODEL

# Csnippets defining the SEIR model

rproc <- Csnippet("
  double dSE, dEI, dIR;
  double betaf;

  if (intervention <= t) {
    betaf = w * beta;
  } else {
    betaf = beta;
  }
  
  double rate[3], trans[3];

  rate[0] = betaf * I/pop;
  rate[1] = sigma;
  rate[2] = gamma;

  reulermultinom(1, S, &rate[0], dt, &trans[0]);
  reulermultinom(1, E, &rate[1], dt, &trans[1]);
  reulermultinom(1, I, &rate[2], dt, &trans[2]);

  dSE = trans[0];
  dEI = trans[1];
  dIR = trans[2];

  S += -dSE;
  E += dSE - dEI;
  I += dEI - dIR;
  R += dIR;

  C += dIR;
")

initlz <- Csnippet("
  double m = pop/(S_0 + E_0 + I_0 + R_0);

  S = nearbyint(m*S_0);
  E = nearbyint(m*E_0);
  I = nearbyint(m*I_0);
  R = nearbyint(m*R_0);

  C = 0;
")

dmeas <- Csnippet("
  double m = rho * C;
  double v = m * (1.0 - rho + psi * psi * m);
  double tol = 1.0e-18;
  if (cases > 0.0) {
    lik = pnorm(cases + 0.5, m, sqrt(v) + tol, 1, 0) - pnorm(cases - 0.5, m, sqrt(v)  + tol, 1, 0) + tol;
  } else {
    lik = pnorm(cases + 0.5, m, sqrt(v) + tol, 1, 0) + tol;
  }
  if (give_log) lik = log(lik);
")

rmeas <- Csnippet("
  double m = rho * C;
  double v = m * (1.0 - rho + psi * psi * m);
  double tol = 1.0e-18;
  cases = rnorm(m, sqrt(v) + tol);
  if (cases > 0.0) {
    cases = nearbyint(cases);
  } else {
    cases = 0.0;
  }
")

# POMP model
ohio_data %>% 
  pomp(t0=ohio_data$time[1],
       time="time",
       rprocess=euler(rproc, delta.t=time_step),
       rinit=initlz,
       rmeasure=rmeas,
       dmeasure=dmeas,
       cdir = code_folder,
       cfile = file_name,       
       accumvars=c("C"),
       statenames=c("S", "E", "I", "R", "C"),
       partrans=parameter_trans(
         log=log_trans_params,
         logit=logit_trans_params),
       paramnames=all_param_names,
       verbose = TRUE
  ) -> mdl_int

mdl_int %>% as.data.frame() %>%
  melt(id="time") %>%
  ggplot(aes(x=time, y=value))+
  geom_line()+
  facet_grid(variable~.,scales="free_y")

# =============================================================================
# GLOBAL MLE

num_test_runs <- as.integer(prop$num_test_runs)

num_guesses <- as.integer(prop$num_guesses)   
num_filter_iter <- as.integer(prop$num_filter_iter)
num_particles <- as.integer(prop$num_particles)
num_replicates <- as.integer(prop$num_replicates)

perturb_sizes <- list(beta=as.numeric(prop$perturb_size_beta),
                      rho=as.numeric(prop$perturb_size_rho),
                      psi=as.numeric(prop$perturb_size_psi),
                      w=as.numeric(prop$perturb_size_w))

cool_frac <- as.numeric(prop$cool_frac)
cool_type <- prop$cool_type

# Variables to use in the scatterplot matrix showing the result of the IF search
pair_vars <- ~loglik+beta+rho+psi+w

# Test run from single starting point in parameter space and no replicates

registerDoParallel()
#set.seed(test_mle_seed, kind="L'Ecuyer")

guess <- apply(free_param_box, 1, function(x) runif(1, x[1], x[2]))

bake(file=file.path(cooking_folder, "box_search_local.rds"), {
  foreach(i=1:num_test_runs,
          .packages='pomp', .combine=c, .options.multicore=list(set.seed=TRUE)
  ) %dopar%  
  {
    mif2(
      mdl_int,
      params=c(guess, fixed_param_values),
      Np=num_particles,
      Nmif=num_filter_iter,
      cooling.type=cool_type,
      cooling.fraction.50=cool_frac,
      rw.sd=do.call(rw.sd, perturb_sizes)
    )
  }
}) -> mifs_test

ggplot(data=melt(traces(mifs_test)),
       aes(x=iteration, y=value, group=L1, color=factor(L1))) +
  geom_line() +
  guides(color=FALSE) +
  facet_wrap(~variable, scales="free_y") +
  theme_bw()

ggsave(file.path(plotting_folder, "1-mle_local_search.pdf"))

# Full MLE with multiple starting points for the free parameters

registerDoParallel()
set.seed(full_mle_seed, kind="L'Ecuyer")

stew(file=file.path(cooking_folder, "box_search_global.rda"), {
  param_guesses <- as.data.frame(apply(free_param_box, 1, function(x) runif(num_guesses, x[1], x[2])))
  
  workers <- getDoParWorkers()
  systime <- system.time({
    res_global <- foreach(guess=iter(param_guesses,"row"), 
                         .packages='pomp', .combine=rbind, .options.multicore=list(set.seed=TRUE)
    ) %dopar% {
      mifs_local <- mif2(
        mdl_int,
        params=c(unlist(guess), fixed_param_values),
        Np=num_particles,
        Nmif=num_filter_iter,
        cooling.type=cool_type,
        cooling.fraction.50=cool_frac,
        rw.sd=do.call(rw.sd, perturb_sizes)
        )
      ll <- logmeanexp(replicate(num_replicates, logLik(pfilter(mifs_local, Np=num_particles))), se = TRUE)
      data.frame(as.list(coef(mifs_local)), loglik=ll[1], loglik.se=ll[2])
    }
  })
}, seed=global_search_seed, kind="L'Ecuyer")

res_global <- as.data.frame(res_global)

all <- ldply(list(guess=param_guesses, result=subset(res_global, loglik > max(loglik)-50)), .id="type")

pdf(file=file.path(plotting_folder, "2-mle_global_search.pdf"))
pairs(pair_vars, data=all, col=ifelse(all$type=="guess", grey(0.5), "red"), pch=16)
dev.off()

# Getting the "optimal" parameters (corresponding to highest loglik)

mle_params <- arrange(rbind(res_global), -loglik)
write.csv(mle_params, file=file.path(output_folder, "param_mle_global_search.csv"), row.names=FALSE, na="")

log_idx <- length(mle_params) - 1
mle_global <- mle_params[which.max(mle_params[,log_idx] ),] 
mle_global %>% extract(all_param_names) %>% unlist() -> theta_opt
write.csv(theta_opt, file=file.path(output_folder, "param_opt_estimates.csv"), row.names=TRUE, na="")
theta_opt

# =============================================================================
# BOOTSTRAP PARAMETER ESTIMATES AND CONFIDENCE INTERVALS

# Bootstrap parameters

bootstrap_confidence <- as.numeric(prop$bootstrap_confidence)
num_quantiles <- as.numeric(prop$num_quantiles)
num_bootstrap_replicates <- as.integer(prop$num_bootstrap_replicates)
use_bootstrap_estimates <- as.logical(prop$use_bootstrap_estimates)

# Getting the mean and CI for the parameters (from bootstrapping the top quantile of the loglik)

# Calculate quantiles
quantiles <- quantile(mle_params$loglik, prob = seq(0, 1, length = num_quantiles+1), type = 5)
top_quantile <- quantiles[num_quantiles][[1]]

# Remove outliers
mle_params_top <- mle_params[top_quantile < mle_params$loglik,]

theta_boot <- theta_opt
mle_mean_ci <- NULL
for (p in free_param_names) {
  p_boot <- boot(mle_params_top[[p]], function(x, i) mean(x[i]), R=num_bootstrap_replicates)

  # Confidence Intervals from the bootstrap:
  p_boot.ci <- boot.ci(p_boot, conf=bootstrap_confidence, type=c("norm", "basic" ,"perc", "bca"))

  p_mean <- mean(p_boot$t[,1])
  p_ci_lo <- p_boot.ci$bca[,4]
  p_ci_hi <- p_boot.ci$bca[,5]

  pdf(file=file.path(plotting_folder, paste("2-", p, "-bootstrap_values.pdf", sep="")))
  hist(p_boot$t[,1], col = "darkgray")
  dev.off()

  print(sprintf("%s %0.2f %0.2f %0.2f", p, p_mean, p_ci_lo, p_ci_hi))
  
  if (is.null(mle_mean_ci)) {
    mle_mean_ci <- data.frame("name" = c(p), "mean" = c(p_mean), "ci_lo" = c(p_ci_lo), "ci_hi" = c(p_ci_hi), stringsAsFactors = FALSE)  
  } else {
    mle_mean_ci <- rbind(mle_mean_ci, c(p, p_mean, p_ci_lo, p_ci_hi))  
  }
  
  theta_boot[p] <- p_mean
}

write.csv(theta_boot, file=file.path(output_folder, "param_bootstrap_estimates.csv"), row.names=TRUE, na="")
write.csv(mle_mean_ci, file=file.path(output_folder, "param_bootstrap_conf_intervals.csv"), row.names=FALSE, na="")

if (use_bootstrap_estimates) {
  theta <- theta_boot
} else {
  theta <- theta_opt
}

# =============================================================================
# SIMULATION W/ BEST PARAMETERS

# Some utility functions to calculate cumulative case numbers

cumulative_curve <- function(dat, t0, t1) {
  total_sum <- 0
  daily_sum <- c()
  for (i in 1:(t1-t0+1)) {
    total_sum <- total_sum + dat$cases[i]
    daily_sum <- c(daily_sum, total_sum)
  }
  return(daily_sum)
}  
  
median_simulation <- function(sdat, n) {
  all_totals <- c()
  for (i in 1:n) {
    sim <- subset(sdat, .id == i)
    tot <- sum(sim$cases)
    all_totals <- c(all_totals, tot)
  }

  # Taking the median
  n2 <- 0.5 * n
  median_idx <- order(all_totals)[n2]
  median_sim <- subset(sdat, .id == median_idx)
  
  return(median_sim)
}

calc_cumulative_counts <- function(sdat, t0, t1) {
  all_csum <- c()
  csum <- NULL
  for (t in t0:t1) {
    dt <- subset(sdat, .id != "data" & time == t)
    if (is.null(csum)) {
      csum <- dt$cases
    } else {
      csum <- csum + dt$cases
    }
    
    all_csum <- c(all_csum, csum)
  }
  return(all_csum)
}

#set.seed(test_sim_seed)

mdl_int  %>%
  simulate(params=theta, nsim=9, format="data.frame", include.data=TRUE) -> sim_data

sim_data %>%
  ggplot(aes(x=time, y=cases, group=.id, color=(.id=="data"))) +
  guides(color=FALSE) +
  geom_line() + facet_wrap(~.id, ncol=2)

ggsave(file.path(plotting_folder, "3-simulations_int.pdf"))

# =============================================================================
# PLOT OHIO OUTBREAK WITHOUT INTERVENTION

theta_noint <- theta
theta_noint["w"] <- 1
theta_noint

#set.seed(test_sim_seed)

mdl_int  %>%
  simulate(params=theta_noint, nsim=9, format="data.frame", include.data=TRUE) -> sim_data_noint

sim_data_noint %>%
  ggplot(aes(x=time, y=cases, group=.id, color=(.id=="data"))) +
  guides(color=FALSE) +
  geom_line() + facet_wrap(~.id, ncol=2)

ggsave(file.path(plotting_folder, "3-simulations_noint.pdf"))

# Running a large number of simulations

set.seed(full_sim_seed)

mdl_int %>% 
  simulate(params=theta, nsim=num_sims, format = "data.frame", include.data=TRUE) -> sim_data_int

mdl_int %>% 
  simulate(params=theta_noint, nsim=num_sims, format = "data.frame", include.data=TRUE) -> sim_data_noint

# Adding cumulative counts
cobs <- cumulative_curve(ohio_data, first_week, last_week)

all_csum <- calc_cumulative_counts(sim_data_int, first_week, last_week)
all_csum <- c(cobs, all_csum)
sim_data_int <- cbind(sim_data_int, cumulative=all_csum)

all_csum <- calc_cumulative_counts(sim_data_noint, first_week, last_week)
all_csum <- c(cobs, all_csum)
sim_data_noint <- cbind(sim_data_noint, cumulative=all_csum)

# Compare the observed data with the simulated data between the 5th and 95th percentiles
case_area_plot <- function(sdat, fname) {
  sdat %>%
    select(time, .id, cases) %>%
    mutate(data=.id=="data") %>%
    ddply(~time+data, plyr::summarize,
      p=c(0.05, 0.5, 0.95), q=quantile(cases, prob=p, names=FALSE)) %>%
    mutate(
      p=mapvalues(p, from=c(0.05, 0.5, 0.95), to=c("lo", "med", "hi")),
      data=mapvalues(data, from=c(TRUE, FALSE), to=c("data", "simulation"))
    ) %>%
    spread(p, q) %>%
    ggplot(aes(x=time, y=med, color=data, fill=data, ymin=lo, ymax=hi)) +
           geom_ribbon(alpha=0.2)
  ggsave(file.path(plotting_folder, fname))  
}

case_area_plot(sim_data_int, "4-simulated_percentiles-with_int.pdf")
case_area_plot(sim_data_noint, "4-simulated_percentiles-wout_int.pdf")

# Getting the median cumulative curve for simulated data with intervention
num_sims_int <- length(unique(sim_data_int$.id)) - 1
sim_int <- median_simulation(sim_data_int, num_sims_int)
csim_int <- cumulative_curve(sim_int, 1, last_week - first_week + 1)

# Getting the median cumulative curve for the simulations w/out intervention
num_sims_noint <- length(unique(sim_data_noint$.id)) - 1
sim_noint <- median_simulation(sim_data_noint, num_sims_noint)
csim_noint <- cumulative_curve(sim_noint, 1, last_week - first_week + 1)

df <- data.frame('time' = seq(first_week, last_week),
                 'obs_data' = cobs,
                 'sim_data_int' = csim_int,
                 'sim_data_noint' = csim_noint)

ggplot(df, aes(time)) + 
  geom_line(aes(y = obs_data, colour = "Real Data")) + 
  geom_line(aes(y = sim_data_int, colour = "Simulation with Intervention")) +
  geom_line(aes(y = sim_data_noint, colour = "Simulation without Intervention")) +
  ylab('Cumulative Number of Cases')

ggsave(file.path(plotting_folder, "4-cumulative_cases_median-comparison.pdf"))

# Compare the observed cumulative counts with the simulated counts between the 5th and 95th percentiles 
cumulative_area_plot <- function(sdat, ymax, fname) {
  sdat %>%
    select(time, .id, cumulative) %>%
    mutate(data=.id=="data") %>%
    ddply(~time+data, plyr::summarize,
      p=c(0.05, 0.5, 0.95), q=quantile(cumulative, prob=p, names=FALSE)) %>%
    mutate(
      p=mapvalues(p, from=c(0.05, 0.5, 0.95), to=c("lo", "med", "hi")),
      data=mapvalues(data, from=c(TRUE, FALSE), to=c("data", "simulation"))
    ) %>%
    spread(p, q) %>%
    ggplot(aes(x=time, y=med, color=data, fill=data, ymin=lo, ymax=hi)) + ylim(0, ymax) +
           geom_ribbon(alpha=0.2)
  ggsave(file.path(plotting_folder, fname))
}

cmax_int <- max(sim_data_int$cumulative)
cmax_noint <- max(sim_data_noint$cumulative)
cumulative_area_plot(sim_data_int, cmax_int, "4-cumulative_cases_percentiles-comparison-with_int.pdf")
cumulative_area_plot(sim_data_noint, cmax_noint, "4-cumulative_cases_percentiles-comparison-wout_int.pdf")

# =============================================================================
# SHOW HOW OUTBREAK SIZE CHANGES (comparison to the simulation size AND actual size)

outbreak_size <- function(intervention_week) {
  theta_global_int_week <- theta
  theta_global_int_week["intervention"] <- intervention_week
  
  mdl_int %>% 
    simulate(params=c(theta_global_int_week), 
             nsim=num_sims, format = "data.frame", include.data=TRUE) -> simulation_data

  sizes <- c()
  sizes_aweek <- c()
  for (i in 1:num_sims) {
    temp <- subset(simulation_data, .id == i)
    temp_aweek <- subset(data_aweek, .id == i)
    size <- sum(temp$cases)
    sizes <- c(sizes, size)
    sizes_aweek <- c(sizes_aweek, sum(temp_aweek$cases))
  }

  actual_size <- sum(ohio_data$cases)

  final_size <- mean(sizes)
  final_size_aweek <- mean(sizes_aweek)

  percentage <- (final_size_aweek - final_size) / final_size_aweek
  percentage_2 <- (actual_size - final_size) / actual_size

  return(c(final_size, percentage * 100, percentage_2 * 100))
}

outbreak_size <- function(intervention_week) {
  theta_global_int_week <- theta
  theta_global_int_week["intervention"] <- intervention_week
  
  mdl_int %>% 
    simulate(params=c(theta_global_int_week), 
             nsim=num_sims, format = "data.frame", include.data=TRUE) -> simulation_data

  sizes <- c()
  sizes_aweek <- c()
  for (i in 1:num_sims) {
    temp <- subset(simulation_data, .id == i)
    temp_aweek <- subset(data_aweek, .id == i)
    size <- sum(temp$cases)
    sizes <- c(sizes, size)
    sizes_aweek <- c(sizes_aweek, sum(temp_aweek$cases))
  }

  actual_size <- sum(ohio_data$cases)

  final_size <- mean(sizes)
  final_size_aweek <- mean(sizes_aweek)

  percentage <- (final_size_aweek - final_size) / final_size_aweek
  percentage_2 <- (actual_size - final_size) / actual_size

  return(c(final_size, percentage * 100, percentage_2 * 100))
}

red_size <- reduction_size
#red_size <- reduction_size_actual

ggplot(data.frame('time' = seq(1, start_aware_week - 1, length = start_aware_week - 1), 'reduction' = red_size),
       aes(x = time, y = red_size)) + geom_point() + 
       stat_smooth(aes(x = time, y = red_size), method = "lm", formula = y ~ sigmoid(x, a = -1, b = 12), se = TRUE) + 
       ylab('Reduction (%)') + xlab('Week of Intervention') + 
       ggtitle('Reduction in Outbreak Size')

reduction_line = lm(red_size ~ time, 
                    data = data.frame('time' = seq(1, start_aware_week - 1, length=start_aware_week - 1),
                                      'reduction' = red_size))

sink(file.path(output_folder, "reduction_fit_sigmoid.txt"))
summary(reduction_line)
sink()

ggsave(file.path(plotting_folder, "5-outbreak_reduction.pdf"))

# =============================================================================
# CONFIDENCE INTERVALS FOR GLOBAL PARAMS VIA PROFILE LIKELIHOOD

# MCAP settings

calculate_mcap_cis <- as.logical(prop$calculate_mcap_cis)

mcap_confidence <- as.numeric(prop$mcap_confidence)

mcap_profdes_len <- as.integer(prop$mcap_profdes_len)
mcap_nprof <- as.integer(prop$mcap_nprof)
mcap_replicates <- as.integer(prop$mcap_replicates)
mcap_num_particles_pfilt <- as.integer(prop$mcap_num_particles_pfilt)

mcap_num_particles <- as.integer(prop$mcap_num_particles)
mcap_num_filter_iter <- as.integer(prop$mcap_num_filter_iter)

mcap_cool_type <- prop$mcap_cool_type
mcap_cool_frac <- as.numeric(prop$mcap_cool_frac)
mcap_cool_frac_lastif <- as.numeric(prop$mcap_cool_frac_lastif)

mcap_lambda <- list(beta=as.numeric(prop$mcap_lambda_beta),
                    rho=as.numeric(prop$mcap_lambda_rho),
                    psi=as.numeric(prop$mcap_lambda_psi),
                    w=as.numeric(prop$mcap_lambda_w))

mcap_ngrid <- list(beta=as.numeric(prop$mcap_ngrid_beta),
                   rho=as.numeric(prop$mcap_ngrid_rho),
                   psi=as.numeric(prop$mcap_ngrid_psi),
                   w=as.numeric(prop$mcap_ngrid_w))

# Function to calculate the MCAP CIs 

mcap <- function(lp, parameter, confidence, lambda, Ngrid) {
  smooth_fit <- loess(lp ~ parameter, span=lambda)
  parameter_grid <- seq(min(parameter), max(parameter), length.out = Ngrid)
  smoothed_loglik <- predict(smooth_fit, newdata=parameter_grid)
  smooth_arg_max <- parameter_grid[which.max(smoothed_loglik)]
  dist <- abs(parameter - smooth_arg_max)
  included <- dist < sort(dist)[trunc(lambda*length(dist))]
  maxdist <- max(dist[included])
  weight <- rep(0, length(parameter))
  weight[included] <- (1 - (dist[included]/maxdist)^3)^3
  quadratic_fit <- lm(lp ~ a + b, weight=weight,
                      data = data.frame(lp=lp, b=parameter, a=-parameter^2))
  b <- unname(coef(quadratic_fit)["b"] )
  a <- unname(coef(quadratic_fit)["a"] )
  m <- vcov(quadratic_fit)
  
  var_b <- m["b", "b"]
  var_a <- m["a", "a"]
  cov_ab <- m["a", "b"]
  
  se_mc_squared <- (1 / (4 * a^2)) * (var_b - (2 * b/a) * cov_ab + (b^2 / a^2) * var_a)
  se_stat_squared <- 1/(2*a)
  se_total_squared <- se_mc_squared + se_stat_squared
  
  delta <- qchisq(confidence, df=1) * ( a * se_mc_squared + 0.5)
  loglik_diff <- max(smoothed_loglik) - smoothed_loglik
  ci <- range(parameter_grid[loglik_diff < delta])
  list(lp=lp, parameter=parameter, confidence=confidence,
       quadratic_fit=quadratic_fit, quadratic_max=b/(2*a),
       smooth_fit=smooth_fit,
       fit=data.frame(parameter=parameter_grid,
                      smoothed=smoothed_loglik,
                      quadratic=predict(quadratic_fit, list(b = parameter_grid, 
                                                            a = -parameter_grid^2))),
       mle=smooth_arg_max, ci=ci, delta=delta,
       se_stat=sqrt(se_stat_squared), se_mc=sqrt(se_mc_squared), 
       se=sqrt(se_total_squared)
  )
}

# Generates the model for the profile likelihood calculation

plik <- function(pname, pval0, pval1) {
  registerDoParallel()
  
  bake(file=file.path(cooking_folder, paste(pname, "_mcap.rds", sep="")), {  
    
    desel <- which(names(mle_params) %in% c("loglik", "loglik.se", pname))
    mle_params %>% 
      subset(
        loglik > max(loglik,na.rm=TRUE) - 20,
        select=-desel
      ) %>% 
      melt(id=NULL) %>% 
      daply(~variable, function(x)range(x$value)) -> box
    
    starts <- profile_design(pname=seq(pval0, pval1, length=mcap_profdes_len),
                             lower=box[,1], upper=box[,2],
                             nprof=mcap_nprof)
    names(starts) <- sub("pname", pname,names(starts))
    
    psizes <- perturb_sizes[names(perturb_sizes) != pname]
    
    foreach(params=iter(starts, "row"),
            .combine=rbind,
            .packages="pomp",
            .options.multicore=list(set.seed=TRUE),
            .options.mpi=list(seed=mcap_plik_seed, chunkSize=1)
    ) %dopar% {
      mf <- mif2(mdl_int,
                 params=unlist(params),
                 Np=mcap_num_particles,
                 Nmif=mcap_num_filter_iter,
                 cooling.type=mcap_cool_type,
                 cooling.fraction.50=mcap_cool_frac,
                 rw.sd=do.call(rw.sd, psizes)
      )
      mf <- mif2(mf, 
                 Np=mcap_num_particles,
                 Nmif=mcap_num_filter_iter,
                 cooling.fraction.50=mcap_cool_frac_lastif)
      ll <- logmeanexp(replicate(mcap_replicates, 
                                 logLik(pfilter(mf, Np=mcap_num_particles_pfilt))), se=TRUE)
      data.frame(as.list(coef(mf)), loglik=ll[1], loglik.se=ll[2])
    }
  }) -> pmodel
  
  return(pmodel)
}

# Iterate over all the free parameters in the model to calculate their CIs... may take a while!

if (calculate_mcap_cis) {
  par_names <- row.names(free_param_box)
  for (i in 1:nrow(free_param_box)) {
    name <- par_names[i]  
    print(sprintf("Calculating CI for %s...", name))
    
    row <- free_param_box[i,]
    mdl <- plik(pname=name, pval0=row[1], pval1=row[2])
    
    par_range <- seq(row[1], row[2], length=mcap_profdes_len)
    log_likelihoods <- c()
    for (val in par_range) {
      likelihoods <- subset(mdl, abs(mdl[[name]]-val)<1)$loglik
      if (length(likelihoods) == 0) next
      log_likelihoods <- c(log_likelihoods, max(likelihoods))
    }
    
    x <- mcap(log_likelihoods, par_range, mcap_confidence, mcap_lambda[[i]], mcap_ngrid[[i]])
    if (i == 1) {
      cis <- data.frame("name" = c(name), "x0" = c(x$ci[1]), "x1" = c(x$ci[2]), stringsAsFactors = FALSE)  
    } else {
      cis <- rbind(cis, c(name, x$ci[1], x$ci[2]))  
    }
    print(sprintf("%s %0.2f %0.2f", name, x$ci[1], x$ci[2]))    
    
    ggplot(x$fit, aes(parameter, quadratic)) + geom_line() + 
      geom_vline(xintercept=c(x$ci[1], x$ci[2]), linetype=4, colour='red') +
      geom_point(data = data.frame('parameters'=par_range, 'loglik'=log_likelihoods), 
                 aes(parameters, log_likelihoods)) 
    ggsave(file.path(plotting_folder, paste("6-", name, "_ci.pdf", sep="")))    
  }
  write.csv(cis, file=file.path(output_folder, "param_confidence_intervals.csv"), row.names=FALSE, na="")  
}
