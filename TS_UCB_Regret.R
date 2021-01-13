#### This code compares the regret over two algortihms: UCB and Thompson Sampling ####

library(purrr)
library(ggplot2)
library(dplyr)

#######################################################################################
#### Thompson Sampling with no previous knwoledge ####
#######################################################################################

# sets the seed
set.seed(12345)

# let us have three arms, drawn from a bernoulli distribution
arm1_initial <- NULL
arm2_initial <- NULL
arm3_initial <- NULL
arms <- list(arm1_initial, arm2_initial, arm3_initial)
arms_initial <- list(arm1_initial, arm2_initial, arm3_initial)

# true unknown means for the three arms
mean1 <- 0.50
mean2 <- 0.55
mean3 <- 0.45
true_means <- c(mean1, mean2, mean3)
max_true_mean <- max(true_means)

# set of rounds
T2_TS <- 10000

# number of simulations per time step
N_TS <- 150

## Thompson Sampling ##

# prior parameters no knowledge
a1_prior <- 2
b1_prior <- 2
a2_prior <- 2
b2_prior <- 2
a3_prior <- 2
b3_prior <- 2
prior_matrix <- matrix(c(a1_prior, b1_prior, a2_prior, 
                         b2_prior, a3_prior, b3_prior), ncol = 3, byrow = F)


# initial posterior parameters
posterior_matrix <- prior_matrix

# regret matrix: colums represent different simulations, number of colums = number of simulations performed
regretTS_T2 <- matrix(rep(NA, T2_TS * N_TS), ncol = N_TS, byrow = F)

# first for loop iterates over the simulations
for (k in 1:N_TS) {
  # resets the posterior parameters for the next simulation
  posterior_matrix <- prior_matrix
  # resets the arms
  arms <- arms_initial
  # this for loop  creates a simulation
  for (i in  1:T2_TS) {
    samplingTS_T2 <- rep(NA, length(true_means))
    # samples from the posteriors
    for (j in 1:length(true_means)){
      samplingTS_T2[j] <- rbeta(1, posterior_matrix[1, j], posterior_matrix[2, j])
    }
    # identifies which arm was chosen
    index <- which.max(samplingTS_T2)
    
    #computes the regret
    if (i > 1) {
      regretTS_T2[i, k] <- regretTS_T2[i - 1, k] + max_true_mean - true_means[index]
    } else {
      regretTS_T2[i, k] <- max_true_mean - true_means[index]
    }
    
    # updates the arms
    arms[[index]] <- c(arms[[index]], as.integer(rbernoulli(1, true_means[index])))
    
    # updates the posterior parameters
    posterior_matrix[1, index] <- posterior_matrix[1, index] + last(arms[[index]])
    posterior_matrix[2, index] <- posterior_matrix[2, index] + 1 - last(arms[[index]])
  }
}

# calculates the average regrets over the simulations
average_regret_no_knowledgeTS_T2 <- rowMeans(regretTS_T2)

#######################################################################################
#### UCB ####
#######################################################################################
# let us have three arms, drawn from a bernoulli distribution
arms <- arms_initial

# set of rounds
T2_UCB <- 10000

# number of simulations
N_UCB <- 150

# creates an empty vector to store accumulated regret values
regretUCB_T2 <- matrix(rep(NA, T2_UCB * N_UCB), ncol = N_UCB, byrow = F)

## UCB algorithm ##

# makes multiple simulations
for (k in 1:N_UCB) {
  # resets the arms
  arms <- arms_initial
  # creates one simulation
  for (i in 1:T2_UCB) {
    # creates an empty vector to store the sampled values and resets it every iteration
    samplingUCB_T2 <- rep(NA, length(true_means))
    # executes the first step of the UCB algorithm, i.e. chooses each arm once
    if (i <= length(true_means)) {
      arms[[i]] <- c(arms[[i]], 1)
      if (i > 1) {
        regretUCB_T2[i, k] <- regretUCB_T2[i - 1, k] + max_true_mean - true_means[i]
      } else {
        regretUCB_T2[i, k] <- max_true_mean - true_means[i]
      }
    } else {
      # sampling
      for (j in 1:length(true_means)) {
        samplingUCB_T2[j] <- (sum(as.integer(rbernoulli(length(arms[[j]]), true_means[j])) * arms[[j]]) / sum(arms[[j]])) + sqrt(2 * log(j) / sum(arms[[j]]))
      }
      # finds which arm was chosen
      indexUCB_T2 <- which.max(samplingUCB_T2)
      
      # updates the arms
      arms[[indexUCB_T2]] <- c(arms[[indexUCB_T2]], 1)
      regretUCB_T2[i, k] <- regretUCB_T2[i - 1, k] + max_true_mean - true_means[indexUCB_T2]
    }
  }
}

# calculates the average regrets over the simulations
average_regretUCB_T2 <- rowMeans(regretUCB_T2)

## ploting the cumulative regret for the two algorithms
plot_regret_T2 <- tbl_df(data.frame(1:T2_TS, average_regret_no_knowledgeTS_T2, average_regretUCB_T2))

plot_regret_T2 %>%
  ggplot(aes(x = 1:T2_TS, y = average_regret_no_knowledgeTS_T2)) +
  geom_line(color = "red") +
  geom_line(y = average_regretUCB_T2, color = "blue") +
  coord_cartesian(ylim = c(0, 70)) +
  ggtitle("Thompson Sampling (equal priors) vs  UCB") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "Time", y = "Cumulative regret") +
  geom_text(x = 3000, y = 50, label = "UCB", col = "blue") +
  geom_text(x = 7500, y = 42, label = "Thompson Sampling", col = "red")

# plotting the priors for Thompson sampling
# sets the number of data points that will be used
N_priors <- 1000
# creates a vector of values for the x-axis
x_priors <- seq(from = 0, to = 1, length.out = N_priors)
y_prior1 <- dbeta(x_priors, shape1 = 2, shape2 = 2)
y_prior2 <- dbeta(x_priors, shape1 = 2, shape2 = 2)
y_prior3 <- dbeta(x_priors, shape1 = 2, shape2 = 2)

priors_plot <- tbl_df(data.frame(x_priors, y_prior1, y_prior2, y_prior3))

priors_plot %>% 
  ggplot(aes(x = x_priors, y = y_prior1)) +
  geom_line(color = "black") +
  geom_line(y = y_prior2, color = "green4") +
  geom_line(y = y_prior3, color = "darkorange3") +
  ggtitle("Prior distributions when no previous knowledge is available") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "", y = "Beta priors") +
  geom_vline(xintercept = 0.5, linetype="dotted", color = "black", size=1.5) +
  geom_vline(xintercept = 0.55, linetype="dotted", color = "green4", size=1.5) +
  geom_vline(xintercept = 0.45, linetype="dotted", color = "darkorange3", size=1.5)

#######################################################################################
#### Thompson Sampling Bad Priors ####
#######################################################################################
arms <- arms_initial

# prior parameters no knowledge
a1_prior <- 10
b1_prior <- 2
a2_prior <- 2
b2_prior <- 10
a3_prior <- 11
b3_prior <- 2
prior_matrix <- matrix(c(a1_prior, b1_prior, a2_prior, 
                         b2_prior, a3_prior, b3_prior), ncol = 3, byrow = F)


# regret matrix: colums represent different simulations, number of colums = number of simulations performed
regretTS_T2 <- matrix(rep(NA, T2_TS * N_TS), ncol = N_TS, byrow = F)

# first for loop iterates over the simulations
for (k in 1:N_TS) {
  # resets the posterior parameters for the next simulation
  posterior_matrix <- prior_matrix
  # resets the arms
  arms <- arms_initial
  # this for loop  creates a simulation
  for (i in  1:T2_TS) {
    samplingTS_T2 <- rep(NA, length(true_means))
    # samples from the posteriors
    for (j in 1:length(true_means)){
      samplingTS_T2[j] <- rbeta(1, posterior_matrix[1, j], posterior_matrix[2, j])
    }
    # identifies which arm was chosen
    index <- which.max(samplingTS_T2)
    
    #computes the regret
    if (i > 1) {
      regretTS_T2[i, k] <- regretTS_T2[i - 1, k] + max_true_mean - true_means[index]
    } else {
      regretTS_T2[i, k] <- max_true_mean - true_means[index]
    }
    
    # updates the arms
    arms[[index]] <- c(arms[[index]], as.integer(rbernoulli(1, true_means[index])))
    
    # updates the posterior parameters
    posterior_matrix[1, index] <- posterior_matrix[1, index] + last(arms[[index]])
    posterior_matrix[2, index] <- posterior_matrix[2, index] + 1 - last(arms[[index]])
  }
}

# calculates the average regrets over the simulations
average_regret_bad_knowledgeTS_T2 <- rowMeans(regretTS_T2)

# plotting the priors for Thompson sampling
# sets the number of data points that will be used
N_priors <- 1000
# creates a vector of values for the x-axis
x_priors <- seq(from = 0, to = 1, length.out = N_priors)
y_prior1 <- dbeta(x_priors, shape1 = 10, shape2 = 2)
y_prior2 <- dbeta(x_priors, shape1 = 2, shape2 = 10)
y_prior3 <- dbeta(x_priors, shape1 = 11, shape2 = 2)

priors_plot <- tbl_df(data.frame(x_priors, y_prior1, y_prior2, y_prior3))

priors_plot %>% 
  ggplot(aes(x = x_priors, y = y_prior1)) +
  geom_line(color = "black") +
  geom_line(y = y_prior2, color = "green4") +
  geom_line(y = y_prior3, color = "darkorange3") +
  coord_cartesian(ylim = c(0, 4.5)) +
  ggtitle("Prior distributions when bad knowledge is available") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "", y = "Beta priors") +
  geom_vline(xintercept = 0.5, linetype="dotted", color = "black", size=1.5) +
  geom_vline(xintercept = 0.55, linetype="dotted", color = "green4", size=1.5) +
  geom_vline(xintercept = 0.45, linetype="dotted", color = "darkorange3", size=1.5)


#######################################################################################
#### Thompson Sampling with good priors ####
#######################################################################################
arms <- arms_initial

# prior parameters no knowledge
a1_prior <- 2
b1_prior <- 10
a2_prior <- 10
b2_prior <- 2
a3_prior <- 2
b3_prior <- 11
prior_matrix <- matrix(c(a1_prior, b1_prior, a2_prior, 
                         b2_prior, a3_prior, b3_prior), ncol = 3, byrow = F)


# regret matrix: colums represent different simulations, number of colums = number of simulations performed
regretTS_T2 <- matrix(rep(NA, T2_TS * N_TS), ncol = N_TS, byrow = F)

# first for loop iterates over the simulations
for (k in 1:N_TS) {
  # resets the posterior parameters for the next simulation
  posterior_matrix <- prior_matrix
  # resets the arms
  arms <- arms_initial
  # this for loop  creates a simulation
  for (i in  1:T2_TS) {
    samplingTS_T2 <- rep(NA, length(true_means))
    # samples from the posteriors
    for (j in 1:length(true_means)){
      samplingTS_T2[j] <- rbeta(1, posterior_matrix[1, j], posterior_matrix[2, j])
    }
    # identifies which arm was chosen
    index <- which.max(samplingTS_T2)
    
    #computes the regret
    if (i > 1) {
      regretTS_T2[i, k] <- regretTS_T2[i - 1, k] + max_true_mean - true_means[index]
    } else {
      regretTS_T2[i, k] <- max_true_mean - true_means[index]
    }
    
    # updates the arms
    arms[[index]] <- c(arms[[index]], as.integer(rbernoulli(1, true_means[index])))
    
    # updates the posterior parameters
    posterior_matrix[1, index] <- posterior_matrix[1, index] + last(arms[[index]])
    posterior_matrix[2, index] <- posterior_matrix[2, index] + 1 - last(arms[[index]])
  }
}

# calculates the average regrets over the simulations
average_regret_good_knowledgeTS_T2 <- rowMeans(regretTS_T2)

# plotting the priors for Thompson sampling
# sets the number of data points that will be used
N_priors <- 1000
# creates a vector of values for the x-axis
x_priors <- seq(from = 0, to = 1, length.out = N_priors)
y_prior1 <- dbeta(x_priors, shape1 = 2, shape2 = 10)
y_prior2 <- dbeta(x_priors, shape1 = 10, shape2 = 2)
y_prior3 <- dbeta(x_priors, shape1 = 2, shape2 = 11)

priors_plot <- tbl_df(data.frame(x_priors, y_prior1, y_prior2, y_prior3))

priors_plot %>% 
  ggplot(aes(x = x_priors, y = y_prior1)) +
  geom_line(color = "black") +
  geom_line(y = y_prior2, color = "green4") +
  geom_line(y = y_prior3, color = "darkorange3") +
  coord_cartesian(ylim = c(0, 4.5)) +
  ggtitle("Prior distributions when good knowledge is available") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "", y = "Beta priors") +
  geom_vline(xintercept = 0.5, linetype="dotted", color = "black", size=1.5) +
  geom_vline(xintercept = 0.55, linetype="dotted", color = "green4", size=1.5) +
  geom_vline(xintercept = 0.45, linetype="dotted", color = "darkorange3", size=1.5)

#######################################################################################
#### Comparing all scenarios of Thompson Sampling with UCB ####
#######################################################################################
## ploting the cumulative regret for the two algorithms
plot_all_regret_T2 <- tbl_df(data.frame(1:T2_TS, average_regret_no_knowledgeTS_T2, average_regretUCB_T2, average_regret_bad_knowledgeTS_T2,
                                        average_regret_good_knowledgeTS_T2))

plot_regret_T2 %>%
  ggplot(aes(x = 1:T2_TS, y = average_regret_no_knowledgeTS_T2)) +
  geom_line(color = "red") +
  geom_line(y = average_regretUCB_T2, color = "blue") +
  geom_line(y = average_regret_bad_knowledgeTS_T2, color = "black") +
  geom_line(y = average_regret_good_knowledgeTS_T2, color = "green4") +
  coord_cartesian(ylim = c(0, 165)) +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Thompson Sampling with different priors vs UCB") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "Time", y = "Cumulative regret") +
  geom_text(x = 3000, y = 140, label = "Bad Priors") +
  geom_text(x = 5000, y = 20, label = "Good Priors", col = "green4") +
  geom_text(x = 5000, y = 65, label = "UCB", col = "blue") +
  geom_text(x = 8000, y = 40, label = "Equal Priors", col = "red")

