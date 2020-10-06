############################################# Problem 1 ########################################

#1c)

com_cold_sim <- function(lambda, mu) {
  tTimes <- c(0)
  infected = FALSE #The person starts off being healthy
  while (tail(tTimes, 1) < 1000*365){
    #Add a number of days to the vector being either sick or healthy. Days are added increasingly
    if (infected) {
      tTimes = c(tTimes, tail(tTimes,1) + rexp(1, mu))
    }
    else {
      tTimes = c(tTimes, tail(tTimes,1) + rexp(1, lambda))
    }
    infected = !infected # The state of the person switches
  }
  return(tTimes)
}

tTimes <- com_cold_sim(1/100, 1/7) 

xVals = rep(c(1, 2), times = length((tTimes))/2) #A list that alternate between 1 and 2, only to make plots

plot(NULL, NULL, xlim = c(0, 5*365), ylim = c(0.9, 2.1), xlab = "Days", ylab = "State", cex.lab = 1.5, cex.axis = 1.5)
for(i in 1:(length(xVals)-1)){
  lines(tTimes[i:(i+1)], rep(xVals[i], 2), lwd = 2)
}

sim <- function(lambda, mu) {
  days <- 0
  days_infected <- 0
  infected = FALSE #The person starts off being healthy
  while (days + days_infected < 1000*365){
    if (infected) {
      days_infected = days_infected + rexp(1, mu)
    }
    else {
      days = days + rexp(1, lambda)
    }
    infected = !infected
  }
  return(c(days + days_infected, days_infected))
}

result <- sim(1/100, 1/7)
cat("Long run ratio of days infected: ", result[2]/result[1])

###########################################Problem 2#################################################

#2a)

#We first made the correlation function definded in the exercise:
corr <- function(t1, t2){
  (1 + 15*abs(t1-t2))*exp(-15*abs(t1 - t2))
}

#If we multiply this by the standard deviations of each of Y(theta), 
#we arrive at the covariance function:
cov <- function(t1, t2){
  0.25 * corr(t1, t2)
}
#We now make the function that takes in the evaluation points 
#and returns the conditional distribution of the unevaluated 
#thetas. 
#We will use this function both in a) and in c). 

predict <- function(t, Yt){
  n <- length(t)
  
  t_vec <- rep(0, 51)
  for (i in 1:51){
    t_vec[i] <- round(0.25 + 0.005*(i-1), 3)#all points/thetas on the grid 
  }
  tA_vec <- t_vec[!t_vec %in% t] #thetas except the ones we want to condition on
  
  #We now make the different covariance matrices used in 
  #one of the theorems from the lectures:
  cov_BB <- matrix(0, n, n)
  for (i in 1:n){
    for (j in 1:n){
      cov_BB[i, j] <- cov(t[i], t[j])
    }
  }
  
  cov_AB = matrix(0, 51 - n, n)
  
  for (i in 1:(51 - n)){
    for (j in 1:n){
      cov_AB[i, j] <- cov(tA_vec[i], t[j])
    }
  }
  
  cov_BA <- t(cov_AB)
  
  cov_AA <- matrix(0, 51-n, 51-n)
  for (i in 1:(51-n)){
    for (j in 1:(51-n)){
      cov_AA[i, j] <- cov(tA_vec[i], tA_vec[j])
    }
  }
   
  #Then we predict the 51 - n unknown values and find their covariance matrix
  
  mu_A <- rep(0.5, 51-n)
  mu_B <- rep(0.5, n)
  
  xA_vec <- mu_A + cov_AB %*% solve(cov_BB) %*% as.matrix(Yt - mu_B)
  cov <- cov_AA - cov_AB %*% solve(cov_BB) %*% cov_BA
  
  #We return the nonevaluated thetas, the mean vector and the covariance 
  #matrix for the conditional distribution.
  list(as.vector(tA_vec), as.vector(xA_vec), cov)
}
#We then make a function that uses the previous function to make 90% 
#confidence intervals for each nonevaluated theta. It then plots the 
#predictions (in blue) along with the upper adn lower limits of 
#the confidence intervals for each of the thetas. 

plot_results <- function(t, Yt){
  n <- length(t)
  results <- predict(t, Yt)
  thetas <- results[[1]]
  pred_thetas <- results[[2]]
  cov_mat <- results[[3]]
  
  #Plot of predictions
  plot(thetas, pred_thetas, type = "l", col = "blue", ylim = c(0.2, 0.8))
  
  #Making the confidence intervals
  pred_low <- rep(0, 51-n)
  for (i in 1:(51-n)){
    pred_low[i] <- pred_thetas[i] - 1.645 * sqrt(cov_mat[i, i])
  }
  
  pred_high <- rep(0, 51-n)
  for (i in 1:(51-n)){
    pred_high[i] <- pred_thetas[i] + 1.645 * sqrt(cov_mat[i, i])
  }
  #There is a 90% chance that the value lies between the two red lines. 
  lines(thetas, pred_low, col = "red")
  lines(thetas, pred_high, col = "red")
}
#Let us now test the functions for the evaluation points in a)

known_thetas_1 <- round(c(0.300, 0.350, 0.390, 0.410, 0.450), 3)
known_values_1 <- c(0.5, 0.32, 0.40, 0.35, 0.60)

plot_results(known_thetas_1, known_values_1)

#2b)

#The following function uses predict to calculate the probabilities for 
#each of the normally distributed Y(theta)s to be less than a given value val.
less_than <- function(t, Yt, val){
  n <- length(t)
  results <- predict(t, Yt)
  thetas <- results[[1]]
  pred_thetas <- results[[2]]
  cov_mat <- results[[3]]
  
  prob <- rep(0, 51-n) 
  for (i in 1:(51-n)){
    prob[i] <- pnorm(val, pred_thetas, sqrt(cov_mat[i, i]))
  }
  plot(thetas, prob)
  prob
}
#In this exercise, val = 0.3 and we use the same evaluation points as in a):
less_than(known_thetas_1, known_values_1, 0.3)
#From this plot, we see that the chance of Y(theta) is far greater for the 
#first and last thetas than the ones in the midle. 

#2c)
#We now have an extra evaluation point, so we use this new set of points as 
#input in predict.First we visualize the predictions and confidence intervals as in a), 
#Then we plot the probabilities of the values being less than 0.30 as in b). 

known_thetas_2 <- round(c(0.300, 0.330, 0.350, 0.390, 0.410, 0.450), 3)
known_values_2 <- c(0.5, 0.40, 0.32, 0.40, 0.35, 0.60)

plot_results(known_thetas_2, known_values_2)

less_than(known_thetas_2, known_values_2, 0.3)
#The highest probability is stil close to the edges, that is, for the biggest and the 
#smallest thetas. Since the highest probability of being less than 0.3 is at 0.25, 
#we should evaluate there. 

