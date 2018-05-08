#' Simulate population according to linearised SIS model
#' 
#' @param lambda The force of infection
#' @param gamma The rate of recovery (inverse of infectious period)
#' @param max_time Maximum length of time (days) for which to simulate individuals for
#' @param num_individuals Number of individuals to seperately simulate
#' 
#' @return List where each element contains the time series information of an individual.
#' 
#' @examples 
#' simulations <- simulateLinearisedSIS(lambda=1/60, gamma=1/20, num_individuals=400)
simulateLinearedSIS <- function(lambda, gamma, max_time=365, num_individuals=1)
{
  simulations <- lapply(1:num_individuals, function(i) { 
    l <- runLinearisedSISModel(lambda, gamma, max_time)
    l$id <- i
    return (l)})
  
  return (simulations)
}

#' Convert model simulations to panel data
#' 
#' @param simulations A list where each element contains the time series information for a single individual, and also the ID of that individual.
#' @param sampling_times A list where each element contain the times to sample that specific individual.
#' 
#' @return A data frame containing the panel data for the entire population.
#' 
#' @examples
#' simulations <- simulateLinearedSIS(lambda=1/60, gamma=1/20, num_individuals=100)
#' sampling_times <- lapply(1:num_individuals, function(i) { rep(x=1, times=365)})
#' panel_data <- convertSimulationsToPanelData(simulations, sampling_times)
convertSimulationsToPanelData <- function(simulations, sampling_times)
{
  num_individuals <- length(simulations)
  max_time <- max(sapply(1:num_individuals, function(i) { max(simulations[[i]]$t)}))
  panel_data <- lapply(1:num_individuals, function(i) {
    t <- 0
    counter <- 1
    simulation <- simulations[[i]]
    sampling_time <- sampling_times[[i]]
    observation_times <- c()
    state <- c()
    while (t < max_time & counter < length(sampling_time))
    {
      observation_times <- c(observation_times, t)
      idx <- tail(which(simulation$t <= t), n=1)
      if (simulation$S[idx] == 1)
        state <- c(state, 0)
      else
        state <- c(state, 1)
      
      t <- t + sampling_time[counter]
      counter <- counter+1
    }
    
    return (data.frame("state"=state, "t"=observation_times, "id"=i))
  })
  
  panel_data <- do.call(rbind, panel_data)
  return (panel_data)
}

#' Plot marginal posteriors from fitted data.
#' 
#' Plots the marginal posteriors for both lambda and gamma, as well as 95% CIs (solid lines) and 99% CIs (dashed lines)
#' 
#' @param stan_fit The fitted result from stan
#' @param lambda True value of lambda (-1 if true value is unknown)
#' @param gamma True value of gamma (-1 if true value is unknown)
plot_fits <- function(stan_fit, lambda_ = -1, gamma_ = -1)
{
  df <- as.data.frame(stan_fit) %>% select(lambda, gamma) 
  df$R0 <- (df$lambda+df$gamma)/df$gamma
  df <- df %>% melt(id.vars=c())
  cis <- group_by(df, variable) %>% summarise(lower95 = quantile(value, 0.025), upper95 = quantile(value, 0.975),
                                              lower99 = quantile(value, 0.005), upper99 = quantile(value, 0.995))
  
  int_frame <- data.frame("variable"=c("lambda","gamma","R0"), "value"=c(lambda_, gamma_, (lambda_+gamma_)/(gamma_)))
  
  p <- ggplot() + geom_density(aes(x=value), data=df) + facet_wrap(~variable, scales = "free") + 
    labs(x="Parameter Value", y="Density") +
    geom_vline(aes(xintercept=lower95), data=cis, colour="red") + geom_vline(aes(xintercept=upper95), data=cis, colour="red") +
    geom_vline(aes(xintercept=lower99), data=cis, colour="red", linetype="dashed") + geom_vline(aes(xintercept=upper99), data=cis, colour="red", linetype="dashed")
  if (lambda_ == -1 | gamma_ == -1)
  {
    return (p)
    
  }
  
  return (p + geom_vline(aes(xintercept=value), data = int_frame, colour="blue")) 
}

extract_means <- function(stan_fit)
{
  df <- as.data.frame(stan_fit) %>% select(lambda,gamma)
  
  lambda_mean <- mean(df$lambda)
  gamma_mean <- mean(df$gamma)
  r0_mean <- mean((df$lambda+df$gamma)/df$gamma)
  lambda_cis <- quantile(df$lambda, c(0.025, 0.975))
  gamma_cis <- quantile(df$gamma, c(0.025, 0.975))
  r0_cis <- quantile((df$lambda+df$gamma)/df$gamma, c(0.025, 0.975))
  
  lret <- list("lambda"=c("mean"=lambda_mean, "lower95"=as.numeric(lambda_cis[1]), "upper95"=as.numeric(lambda_cis[2])),
               "gamma"=c("mean"=gamma_mean, "lower95"=as.numeric(gamma_cis[1]), "upper95"=as.numeric(gamma_cis[2])),
               "R0"=c("mean"=r0_mean, "lower95"=as.numeric(r0_cis[1]), "upper95"=as.numeric(r0_cis[2])))
  
  class(lret) <- "tmiestimates"
  return (lret)
}

print.tmiestimates <- function(data)
{
  cat("Parameter Estimates from TMI Fit:\n")
  cat("Lambda:\n")
  cat("  Mean: ",data$lambda["mean"], " (1/", round(1/data$lambda["mean"], digits = 2), ")\n", sep="")
  cat("  95% CI: [", data$lambda["lower95"], ", ", data$lambda["upper95"], "] ([1/", round(1/data$lambda["lower95"], digits = 2), ", 1/",
      round(1/data$lambda["upper95"], digits = 2), "])\n", sep="")
  
  cat("Gamma:\n")
  cat("  Mean: ",data$gamma["mean"], " (1/ ", round(1/data$gamma["mean"], digits = 2), ")\n", sep="")
  cat("  95% CI: [", data$gamma["lower95"], ", ", data$gamma["upper95"], "] ([1/", round(1/data$gamma["lower95"], digits = 2), ", 1/",
      round(1/data$gamma["upper95"], digits = 2), "])\n", sep="")
  
  cat("R0:\n")
  cat("  Mean: ",data$R0["mean"], "\n", sep="")
  cat("  95% CI: [", data$R0["lower95"], ", ", data$R0["upper95"], "]", sep="")
  
}

#' Set up panel data for fitting
#' 
#' @param panel_data The input panel data
#' @param formula A formula giving the state/time dependence. The left hand side should be the state variable, the right hand side the time. The state variable should
#' have the format where 0 indicates not infected, and 1 indicates infected.
#' @param id The column of the panel_data that refers to the (unique) individual ID. Should be provided unquoted.
#' 
#' @return A list containing the data for use for fitting, of class "tmidata".
#' 
#' @examples 
#' data <- setupData(panel_data = convertSimulationsToPanelData(simulations), formula = state~t, id = id)
setupData <- function(panel_data, formula, id=id)
{
  quoid <- enquo(id)
  unique_ids <- (select(panel_data, !!quoid) %>% unique())[,1]
  num_individuals <- length(unique_ids)

  num_obs <- sapply(unique_ids, function(i) { nrow(filter(panel_data, (!!quoid) == i))})
  
  while (sum(num_obs == 1) > 0) 
  {
    #Silently drop those with 1 observation.
    keep <- which(num_obs > 1)
    unique_ids <- unique_ids[keep]
    num_individuals <- length(unique_ids)
    num_obs <- sapply(unique_ids, function(i) { nrow(filter(panel_data, (!!quoid) == i))})
  }
  
  times <- lapply(unique_ids, function(i) {
    t <- filter(panel_data, (!!quoid) == i)[[all.vars(formula)[2]]]
    padding <- max(num_obs) - length(t)
    return (c(t, rep(x=0, times = padding)))
  })
  times <- do.call(rbind, times)
  S <- lapply(unique_ids, function(i) { 
    f <- as.numeric(filter(panel_data, (!!quoid) == i)[[all.vars(formula)[1]]])
    padding <- max(num_obs) - length(f)
    return (c(f, rep(x=0, times=padding)))
    
  })
  S <- do.call(rbind, S)
  
  stan_inputs = list("S"=S, "times"=times, "num_individuals"=num_individuals, "num_obs"=num_obs)
  
  class(stan_inputs) <- "tmidata"
  
  return (stan_inputs)
}

#' Fit model using stan
#' 
#' @param data Panel data to fit for, of class "tmidata" (such as that generated by the setupData function)
#' @param iterations Number of iterations to run each chain for
#' @param model Model code for stan. The default model is that for the two-state Markov chain.
#' @param ... Other parameters which are passed to `rstan::sampling`.
#' 
#' @return The stan fit.
fit_stan <- function(data, iterations = 10000, model = tmi_model, ...)
{
  if (class(data) != "tmidata")
    stop("Argument data must have class tmidata")
  
  fit <- sampling(model, data = data, iter=iterations, pars=c("lambda","gamma"), ...)
  
  return (fit)
}

#' Get the per-individual sore durations
#' 
#' @param panel_data Panel data to fit for
#' @param formula A formula of the form state~time, where state is the state variable (1=infectious) and time is the observation time
#' @param id The ID variable which identifies each individual (provide unquoted)
#' 
#' @return A vector of the empricial durations of each skin sore infection, assuming continuous observation.
getSoreDurations <- function(panel_data, formula, id=ID)
{
  idquo <- enquo(id)
  unique_ids <- (panel_data %>% select(!!idquo) %>% unique())[,1]
  l <- lapply(unique_ids, function(i) {
    df <- filter(panel_data, (!!idquo) == i)
    #Get infectious presentations:
    state <- df[,all.vars(formula)[1]]
    times <- df[,all.vars(formula)[2]]
    currently_infectious <- state[1]
    if (currently_infectious == 1) {
      infection_time <- times[1]
    }
    row <- 1
    sore_durations <- c()
    while (row < length(state))
    {
      if (currently_infectious == 1)
      {
        #Find next 0:
        next_susc <- which(state == 0)
        next_susc <- next_susc[which(next_susc>row)][1]
        if (is.na(next_susc))
          break
        sore_durations <- c(sore_durations, times[next_susc]-infection_time)
        currently_infectious <- 0
        row <- next_susc
      }
      
      else
      {
        next_inf <- which(state == 1)
        next_inf <- next_inf[which(next_inf>row)][1]
        if (is.na(next_inf))
          break
        infection_time <- times[next_inf]
        row <- next_inf
        currently_infectious <- 1
      }
    }
    
    return (sore_durations)
  })
  
  ret <- do.call(c,l)
  class(ret) <- "tmiempirical"
  return (ret)
}

print.tmiempirical <- function(means)
{
  cat("Total number of skin sore episodes:",length(means),"\n")
  cat("Empirical mean duration:",mean(means, na.rm=T),"\n")
  cat("95% CIs: [",quantile(means, probs=0.025, na.rm=T),", ",quantile(means,probs=0.975, na.rm=T),"]\n", sep="")
}

#' Generates sampling times to be similar to that from already collected panel data
#' 
#' @param panel_data Panel data from which to draw the observation distribution from
#' @param time_variable The variable corresponding to time (unquoted)
#' @param id The ID variable for each individual
#' @param num_individuals Number of individuals to generate sampling times for
#' @param max_time Maximum time to observe each individual
#' 
#' @return Returns a list of sampling times for each individual.
generateSamplingTimesFromPanelData <- function(panel_data, time_variable, id=id, num_individuals, max_time=365)
{
  quoid <- enquo(id)
  quotime <- enquo(time_variable)
  unique_ids <- (select(panel_data, !!quoid) %>% unique())[,1]
  ages_seen <- sapply(unique_ids, function(i) {
    d <- filter(panel_data, (!!quoid) == i)
    times <- select(d, !!quotime)[,1]
    return (diff(c(0, times)))
  })
  ages_seen <- do.call(c, ages_seen)
  time_dist <- density(ages_seen)
  bw <- time_dist$bw
  
  sampling_times <- lapply(1:num_individuals, function(i) {
    t <- 0
    times <- c()
    while (t < max_time) {
      s <- sample(ages_seen, size=1)
      next_time <-  s + rnorm(n=1, mean=0, sd=sqrt(bw))
      if (next_time < 1) 
        next_time <- 1
      t <- t + next_time
      times <- c(times, t)
    }
    while (length(times) <= 1)
    {
      t <- 0
      times <- c()
      while (t < max_time) {
        s <- sample(ages_seen, size=1)
        next_time <-  s + rnorm(n=1, mean=0, sd=sqrt(bw))
        if (next_time < 1) 
          next_time <- 1
        t <- t + next_time
        times <- c(times, t)
      }
    }
    return (diff(c(0, times)))
  })
  
  return (sampling_times)
}

