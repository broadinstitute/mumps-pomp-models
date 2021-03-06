args = commandArgs(trailingOnly=TRUE)

main_folder <- "./harvard"
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
# DATA

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

harv_data <- read.csv(file.path(main_folder, "cases.csv"))
harv_data <- data.frame(time = harv_data$Time, cases = harv_data$Case.Number)
head(harv_data)

# =============================================================================
# PARAMETERS

pop_size <- as.integer(as.numeric(prop$pop_size) * as.numeric(prop$susc_frac))
exp0 <- as.integer(prop$exp0)
inf0 <- as.integer(prop$inf0)
rec0 <- as.integer(prop$rec0)

time_step <-as.numeric(prop$time_step)
num_sims <- as.integer(prop$num_sims)

total_days <- as.integer(prop$total_days)
total_weeks <- as.integer(prop$total_weeks)
new_diag_day <- as.integer(prop$new_diag_day)
spring_day0 <- as.integer(prop$spring_day0)
spring_day1 <- as.integer(prop$spring_day1)
summer_day0 <- as.integer(prop$summer_day0)

free_param_names <- c("beta", "gamma", "rho", "psi", "p", "q")
free_param_box <- rbind(
  beta = c(as.numeric(prop$bounds_beta_min), as.numeric(prop$bounds_beta_max)),
  gamma = c(as.numeric(prop$bounds_gamma_min), as.numeric(prop$bounds_gamma_max)),
  rho = c(as.numeric(prop$bounds_rho_min), as.numeric(prop$bounds_rho_max)),
  psi = c(as.numeric(prop$bounds_psi_min), as.numeric(prop$bounds_psi_max)),
  p = c(as.numeric(prop$bounds_p_min), as.numeric(prop$bounds_p_max)),
  q = c(as.numeric(prop$bounds_q_min), as.numeric(prop$bounds_q_max))
)

log_trans_params <- c("beta", "psi", "q")
logit_trans_params <- c("gamma", "rho", "p")

fixed_param_names <- c("pop", "S_0", "E_0", "I_0", "R_0", "sigma", "intervention", "spring0", "spring1", "summer0")
fixed_param_values <- c(pop=pop_size, 
                        S_0=1-(exp0+inf0+rec0)/pop_size, E_0=exp0/pop_size, I_0=inf0/pop_size, R_0=rec0/pop_size, 
                        sigma=as.numeric(prop$sigma), intervention=new_diag_day, 
                        spring0=spring_day0, spring1=spring_day1, summer0=summer_day0)

all_param_names <- c(free_param_names, fixed_param_names)

# Random seeds, keep unchanged to ensure reproducibilty of results
test_mle_seed <- as.integer(prop$test_mle_seed)
full_mle_seed <- as.integer(prop$full_mle_seed)
global_search_seed <- as.integer(prop$global_search_seed)
test_sim_seed <- as.integer(prop$test_sim_seed)
full_sim_seed <- as.integer(prop$full_sim_seed)
mcap_plik_seed <- as.integer(prop$mcap_plik_seed)

# =============================================================================
# CLEAN DATA

# convert daily data to weekly data 
harv_data$week = floor(harv_data$time/7)

weekly_data <- data.frame(week = harv_data$week, cases = harv_data$cases)
weekly_data <- group_by(weekly_data, week)

sum_cases = c()
for (i in unique(weekly_data$week)){
  temp = subset(weekly_data, week == i)
  sum_cases = c(sum_cases, sum(temp$cases))
}

weekly_data <- data.frame(week = unique(weekly_data$week), cases = sum_cases)

other_data = data.frame('time' = seq(1, 600), 'cases' = rep(0, 600))

# keep the relevant columns 
harv_data = harv_data[,1:2]

#FINAL data
temp = merge(x = other_data, y = harv_data, by = "time", all.x = TRUE)
temp[is.na(temp)] <- 0
data = data.frame('time' = temp$time, 'cases' = temp$cases.x + temp$cases.y)
data

# plot both daily and weekly curves 
ggplot(data=subset(data, time <= total_days), aes(x=time, y=cases, group=1)) + geom_line()
ggsave(file.path(plotting_folder, "0-daily-cases.pdf"))

ggplot(data=subset(weekly_data, week <= total_weeks), aes(x=week, y=cases, group=1)) + geom_line()
ggsave(file.path(plotting_folder, "0-weekly-cases.pdf"))

data=subset(data,time <= total_days)

# =============================================================================
# COMPARTMENTAL MODEL

# C snippets

rproc <- Csnippet("
  double betaf, gammaf;
  double rate[3], trans[3];

  if (intervention <= t) {
      gammaf = q * gamma;
  } else {
      gammaf = gamma;
  }

  if ((spring0 <= t && t <= spring1) || summer0 <= t) {
    betaf = p * beta;
  } else {
    betaf = beta;
  }

  rate[0] = betaf * I/pop;
  rate[1] = sigma;
  rate[2] = gammaf;

  // transitions between classes
  reulermultinom(1, S, &rate[0], dt, &trans[0]);
  reulermultinom(1, E, &rate[1], dt, &trans[1]);
  reulermultinom(1, I, &rate[2], dt, &trans[2]);

  S += -trans[0];
  E += trans[0] - trans[1];
  I += trans[1] - trans[2];
  R += trans[2];

  // Assigning the right number to the accumulation variable that's used
  // in the observation model is absolutely critical!!!!
  C += trans[2]; // We are observing the number of infectious cases that get quarantined when identified
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
data %>% 
  pomp(t0=data$time[1],
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
                      gamma=as.numeric(prop$perturb_size_gamma),
                      rho=as.numeric(prop$perturb_size_rho),
                      psi=as.numeric(prop$perturb_size_psi),
                      p=as.numeric(prop$perturb_size_p),
                      q=as.numeric(prop$perturb_size_q))

cool_frac <- as.numeric(prop$cool_frac)
cool_type <- prop$cool_type

# Variables to use in the scatterplot matrix showing the result of the IF search
pair_vars <- ~loglik+beta+gamma+rho+psi+p+q

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
# SIMULATION W/ BEST PARAMETERS

# Some utility functions to calculate cumulative case numbers

rem_low_count_simulations <- function(sdat, n) {
  wlim <- 2
  wday <- 7
  
  all_totals <- c()
  for (i in 1:n) {
    sim <- subset(sdat, .id == i)
    tot <- sum(sim$cases)
    sim$week <- floor(sim$time/wday)
    wdat <- data.frame(week = sim$week, cases = sim$cases)
    wdat <- group_by(wdat, week)

    wcases <- c()
    for (w in unique(wdat$week)){
      temp <- subset(wdat, week == w)
      wcases <- c(wcases, sum(temp$cases))
    }
    wdat <- data.frame(week = unique(wdat$week), cases = wcases)
    if (any(wcases[1:(length(wcases)-wlim)] == 0)) {
      sdat <- subset(sdat, .id != i)
    } else{
      all_totals <-  c(all_totals, tot)
    }
  }
  #print(mean(all_totals))
  
  # Make ids consecutive
  uniq <- unique(sdat$.id)
  uid <- 1
  for (u in uniq) {
    if (u == 'data') next
    sdat$.id[sdat$.id == u] <- uid
    uid <- uid + 1  
  }
  
  return(sdat)
}

cumulative_curve <- function(dat, len) {
  total_sum <- 0
  daily_sum <- c()
  for (i in 1:len) {
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

calc_cumulative_counts <- function(sdat, len) {
  all_csum <- c()
  csum <- NULL
  for (t in 1:len) {
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
  simulate(params=theta_opt, nsim=9, format="data.frame", include.data=TRUE) -> sim_data

sim_data %>%
  ggplot(aes(x=time, y=cases, group=.id, color=(.id=="data"))) +
  guides(color=FALSE) +
  geom_line() + facet_wrap(~.id, ncol=2)

ggsave(file.path(plotting_folder, "3-simulations_int.pdf"))

# =============================================================================
# PLOT HARVARD OUTBREAK WITHOUT INTERVENTION

theta_noq <- theta_opt
theta_noq["q"] <- 1

#set.seed(test_sim_seed)

mdl_int  %>%
  simulate(params=theta_noq, nsim=9, format="data.frame", include.data=TRUE) -> sim_data_noq

sim_data_noq %>%
  ggplot(aes(x=time, y=cases, group=.id, color=(.id=="data"))) +
  guides(color=FALSE) +
  geom_line() + facet_wrap(~.id, ncol=2)

ggsave(file.path(plotting_folder, "3-simulations_noint.pdf"))

# =============================================================================
# COMPARE CUMULATIVE NUMBERS OF OUTBREAK W/ AND W/O INTERVENTION

# Running a large number of simulations with and without intervention

set.seed(full_sim_seed)

mdl_int %>% 
  simulate(params=theta_opt, nsim=num_sims, format = "data.frame", include.data=TRUE) %>% 
  rem_low_count_simulations(num_sims) -> sim_data_int

mdl_int %>% 
  simulate(params=theta_noq, nsim=num_sims, format = "data.frame", include.data=TRUE) %>%  
  rem_low_count_simulations(num_sims) -> sim_data_noint

# Adding cumulative counts
cobs <- cumulative_curve(data, total_days)

all_csum <- calc_cumulative_counts(sim_data_int, total_days)
all_csum <- c(cobs, all_csum)
sim_data_int <- cbind(sim_data_int, cumulative=all_csum)

all_csum <- calc_cumulative_counts(sim_data_noint, total_days)
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
csim_int <- cumulative_curve(sim_int, total_days)

# Getting the median cumulative curve for the simulations w/out interventionsim_noint <- median_simulation(sim_data_noint, num_sims)
num_sims_noint <- length(unique(sim_data_noint$.id)) - 1
sim_noint <- median_simulation(sim_data_noint, num_sims_noint)
csim_noint <- cumulative_curve(sim_noint, total_days)

df <- data.frame('time' = seq(1, total_days), 
                 'obs_data' = cobs,
                 'sim_data_int' = csim_int,
                 'sim_data_noint' = csim_noint)

ggplot(df, aes(time)) + 
  geom_line(aes(y = obs_data, colour = "Real Data")) + 
  geom_line(aes(y = sim_data_int, colour = "Simulation with Intervention")) +
  geom_line(aes(y = sim_data_noint, colour = "Simulation without Intervention")) +
  geom_vline(xintercept = new_diag_day, colour = 'black', linetype = 3)  + 
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

cmax <- 0.8 * max(sim_data_int$cumulative, sim_data_noint$cumulative)
cumulative_area_plot(sim_data_int, cmax, "4-cumulative_cases_percentiles-comparison-with_int.pdf")
cumulative_area_plot(sim_data_noint, cmax, "4-cumulative_cases_percentiles-comparison-wout_int.pdf")

# =============================================================================
# SHOW HOW OUTBREAK SIZE CHANGES (comparison to the simulation size AND actual size)

outbreak_size <- function(intervention_day) {
  theta_global_int_day <- theta_opt
  theta_global_int_day["intervention"] <- intervention_day
  
  mdl_int %>% 
    simulate(params=c(theta_global_int_day), 
             nsim=num_sims, format = "data.frame", include.data=TRUE) -> simulation_data

  sizes <- c()
  sizes_dday <- c()
  for (i in 1:num_sims) {
    temp <- subset(simulation_data, .id == i)
    temp_dday <- subset(data_dday, .id == i)
    size <- sum(temp$cases)
    sizes <- c(sizes, size)
    sizes_dday <- c(sizes_dday, sum(temp_dday$cases))
  }
  
  actual_size <- sum(data$cases)

  final_size <- mean(sizes)
  final_size_dday <- mean(sizes_dday)

  percentage <- (final_size_dday - final_size) / final_size_dday
  percentage_2 <- (actual_size - final_size) / actual_size

  return(c(final_size, percentage * 100, percentage_2 * 100))
}

theta_dday <- theta_opt
theta_dday["intervention"] <- new_diag_day

mdl_int %>% 
   simulate(params=theta_dday, nsim=num_sims, format = "data.frame", include.data=TRUE) -> data_dday

reduction_size = c()
reduction_size_actual = c()
for (i in seq(1, new_diag_day - 1, length=new_diag_day - 1)) {
    size = outbreak_size(floor(i))
    reduction_size = c(reduction_size, size[2])
    reduction_size_actual = c(reduction_size_actual, size[3])
}

red_size <- reduction_size
#red_size <- reduction_size_actual

ggplot(data.frame('time' = seq(1, new_diag_day - 1, length = new_diag_day - 1), 'reduction' = red_size),
       aes(x = time, y = red_size)) + geom_point() + 
       stat_smooth(aes(x = time, y = red_size), method = "lm", formula = y ~ x, se = TRUE) + 
       ylab('Reduction (%)') + xlab('Day of Intervention') + 
       ggtitle('Reduction in Outbreak Size')

reduction_line = lm(red_size ~ time, 
                    data = data.frame('time' = seq(1, new_diag_day - 1, length=new_diag_day - 1),
                                      'reduction' = red_size))

sink(file.path(output_folder, "reduction_fit_line.txt"))
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
                    gamma=as.numeric(prop$mcap_lambda_gamma),
                    rho=as.numeric(prop$mcap_lambda_rho),
                    psi=as.numeric(prop$mcap_lambda_psi),
                    p=as.numeric(prop$mcap_lambda_p),
                    q=as.numeric(prop$mcap_lambda_q))


mcap_ngrid <- list(beta=as.numeric(prop$mcap_ngrid_beta),
                   gamma=as.numeric(prop$mcap_ngrid_gamma),
                   rho=as.numeric(prop$mcap_ngrid_rho),
                   psi=as.numeric(prop$mcap_ngrid_psi),
                   p=as.numeric(prop$mcap_ngrid_p),
                   q=as.numeric(prop$mcap_ngrid_q))

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

logitTransform <- function(p) {
  p <- min(0.99999, max(0.00001, p))
  return(log(p/(1-p)))
}
sigmoidTransform <- function(p) { 1/(1+exp(-p)) }

logTransform <- function(p) { 
  p <- max(0.00001, p)
  return(log(p))
}
expTransform <- function(p) { exp(p) }

# Iterate over all the free parameters in the model to calculate their CIs... may take a while!

if (calculate_mcap_cis) {
  par_names <- row.names(free_param_box)
  for (i in 1:nrow(free_param_box)) {
    name <- par_names[i]  
    print(sprintf("Calculating CI for %s...", name))
    
    row <- free_param_box[i,]
    mdl <- plik(pname=name, pval0=row[1], pval1=row[2])
    
    par_range <- seq(row[1], row[2], length=mcap_profdes_len)
    par_range_delta <- (row[2] - row[1]) / (2 * (mcap_profdes_len-1))
    log_likelihoods <- c()
    par_values <- c()
    for (val in par_range) {
      likelihoods <- subset(mdl, abs(mdl[[name]]-val)<par_range_delta)$loglik
      if (length(likelihoods) == 0) next
      log_likelihoods <- c(log_likelihoods, max(likelihoods))
      
      if (is.element(name, logit_trans_params)) {
        logitVal <- logitTransform(val)
        par_values <- c(par_values, logitVal)
      } else if (is.element(name, log_trans_params)) {
        logVal = logTransform(val)
        par_values <- c(par_values, logVal)
      } else {
        par_values <- c(par_values, val)
      }
    }
    
    # Discarding outlier values that are too low
    par_prof <- data.frame(par_values, log_likelihoods)
    Q <- quantile(par_prof$log_likelihoods, probs=c(.25, .75), na.rm = FALSE)
    iqr <- IQR(par_prof$log_likelihoods)
    par_prof <- subset(par_prof, par_prof$log_likelihoods > (Q[1] - 1.5*iqr))
    
    # Getting the CI from the MCAP function and transforming back to the original scale if needed
    x <- mcap(par_prof$log_likelihoods, par_prof$par_values, mcap_confidence, mcap_lambda[[i]], mcap_ngrid[[i]])
    tci <- x$ci
    if (is.element(name, logit_trans_params)) {
      tci[1] <- sigmoidTransform(x$ci[1])
      tci[2] <- sigmoidTransform(x$ci[2])
      var_label <- paste("logit(", name, ")")
    } else if (is.element(name, log_trans_params)) {
      tci[1] <- expTransform(x$ci[1])
      tci[2] <- expTransform(x$ci[2])
      var_label <- paste("log(", name, ")")
    } else {
      var_label <- name
    }
    print(sprintf("%s %0.2f %0.2f", name, tci[1], tci[2]))
    
    ggplot(x$fit, aes(parameter, quadratic)) + geom_line() + 
      geom_vline(xintercept=c(x$ci[1], x$ci[2]), linetype=4, colour='red') +
      geom_point(data = data.frame('parameters'=par_prof$par_values, 'loglik'=par_prof$log_likelihoods), 
                 aes(parameters, par_prof$log_likelihoods)) + xlab(var_label) + ylab("-loglik")
    ggsave(file.path(plotting_folder, paste("6-", name, "_ci.pdf", sep="")))
    
    if (i == 1) {
      mcap_cis <- data.frame("name" = c(name), "x0" = c(tci[1]), "x1" = c(tci[2]), stringsAsFactors = FALSE)  
    } else {
      mcap_cis <- rbind(mcap_cis, c(name, tci[1], tci[2]))
    }
  }
  write.csv(mcap_cis, file=file.path(output_folder, "param_mcap_confidence_intervals.csv"), row.names=FALSE, na="")  
}