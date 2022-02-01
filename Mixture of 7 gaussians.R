#' Expectation Step of the EM Algorithm
#'
#' Calculate the posterior probabilities (soft labels) that each component
#' has to each data point.
#'
#' @param sd.vector Vector containing the standard deviations of each component
#' @param mu.vector Vector containing the mean of each component
#' @param alpha.vector Vector containing the mixing weights  of each component
#' @return Named list containing the loglik and posterior.df
e_step <- function(x, mu.vector, sd.vector, alpha.vector) {
  comp1.prod <- dnorm(x, mu.vector[1], sd.vector[1]) * alpha.vector[1]
  comp2.prod <- dnorm(x, mu.vector[2], sd.vector[2]) * alpha.vector[2]
  comp3.prod <- dnorm(x, mu.vector[3], sd.vector[3]) * alpha.vector[3]
  comp4.prod <- dnorm(x, mu.vector[4], sd.vector[4]) * alpha.vector[4]
  comp5.prod <- dnorm(x, mu.vector[5], sd.vector[5]) * alpha.vector[5]
  comp6.prod <- dnorm(x, mu.vector[6], sd.vector[6]) * alpha.vector[6]
  comp7.prod <- dnorm(x, mu.vector[7], sd.vector[7]) * alpha.vector[7]
  sum.of.comps <- comp1.prod + comp2.prod +comp3.prod +comp4.prod +
                  comp5.prod +comp6.prod + comp7.prod
  comp1.post <- comp1.prod / sum.of.comps
  comp2.post <- comp2.prod / sum.of.comps
  comp3.post <- comp3.prod / sum.of.comps
  comp4.post <- comp4.prod / sum.of.comps
  comp5.post <- comp5.prod / sum.of.comps
  comp6.post <- comp6.prod / sum.of.comps
  comp7.post <- comp7.prod / sum.of.comps
  
  sum.of.comps.ln <- log(sum.of.comps, base = exp(1))
  sum.of.comps.ln.sum <- sum(sum.of.comps.ln)
  
  list("loglik" = sum.of.comps.ln.sum,
       "posterior.df" = cbind(comp1.post, comp2.post,
                              comp3.post, comp4.post,
                              comp5.post, comp6.post,
                              comp7.post))
}

#' Maximization Step of the EM Algorithm
#'
#' Update the Component Parameters
#'
#' @param x Input data.
#' @param posterior.df Posterior probability data.frame.
#' @return Named list containing the mean (mu), variance (var), and mixing
#'   weights (alpha) for each component.
m_step <- function(x, posterior.df) {
  comp1.n <- sum(posterior.df[, 1])
  comp2.n <- sum(posterior.df[, 2])
  comp3.n <- sum(posterior.df[, 3])
  comp4.n <- sum(posterior.df[, 4])
  comp5.n <- sum(posterior.df[, 5])
  comp6.n <- sum(posterior.df[, 6])
  comp7.n <- sum(posterior.df[, 7])
  
  comp1.mu <- 1/comp1.n * sum(posterior.df[, 1] * x)
  comp2.mu <- 1/comp2.n * sum(posterior.df[, 2] * x)
  comp3.mu <- 1/comp3.n * sum(posterior.df[, 3] * x)
  comp4.mu <- 1/comp4.n * sum(posterior.df[, 4] * x)
  comp5.mu <- 1/comp5.n * sum(posterior.df[, 5] * x)
  comp6.mu <- 1/comp6.n * sum(posterior.df[, 6] * x)
  comp7.mu <- 1/comp7.n * sum(posterior.df[, 7] * x)
  
  comp1.var <- sum(posterior.df[, 1] * (x - comp1.mu)^2) * 1/comp1.n
  comp2.var <- sum(posterior.df[, 2] * (x - comp2.mu)^2) * 1/comp2.n
  comp3.var <- sum(posterior.df[, 3] * (x - comp3.mu)^2) * 1/comp3.n
  comp4.var <- sum(posterior.df[, 4] * (x - comp4.mu)^2) * 1/comp4.n
  comp5.var <- sum(posterior.df[, 5] * (x - comp5.mu)^2) * 1/comp5.n
  comp6.var <- sum(posterior.df[, 6] * (x - comp6.mu)^2) * 1/comp6.n
  comp7.var <- sum(posterior.df[, 7] * (x - comp7.mu)^2) * 1/comp7.n
  
  comp1.alpha <- comp1.n / length(x)
  comp2.alpha <- comp2.n / length(x)
  comp3.alpha <- comp3.n / length(x)
  comp4.alpha <- comp4.n / length(x)
  comp5.alpha <- comp5.n / length(x)
  comp6.alpha <- comp6.n / length(x)
  comp7.alpha <- comp7.n / length(x)
  
  list("mu" = c(comp1.mu, comp2.mu,
                comp3.mu, comp4.mu,
                comp5.mu,comp6.mu,comp7.mu),
       
       "var" = c(comp1.var, comp2.var,
                comp3.var, comp4.var,
                comp5.var, comp6.var, comp7.var),
       
       "alpha" = c(comp1.alpha, comp2.alpha,
                   comp3.alpha, comp4.alpha,
                   comp5.alpha, comp6.alpha, comp7.alpha))
}


library("dplyr")



sim = log(rnorm(10000)^2) 
# or sim = log(rchisq(100000,1))

sim.kmeans <- kmeans(sim, 7)


sim.kmeans["size"]

sim.kmeans.cluster <- sim.kmeans$cluster

sim.df <- data.frame(x = sim, cluster = sim.kmeans.cluster)

sim.summary.df <- sim.df %>%
  group_by(cluster) %>%
  summarize(mu = mean(x), variance = var(x), std = sd(x), size = n())

sim.summary.df %>%
  select(cluster, mu, variance, std)

sim.summary.df <- sim.summary.df %>%
        mutate(alpha = size / sum(size))


for (i in 1:1000) {
  if (i == 1) {
    # Initialization
    e.step <- e_step(sim, sim.summary.df[["mu"]], sim.summary.df[["std"]],
                     sim.summary.df[["alpha"]])
    m.step <- m_step(sim, e.step[["posterior.df"]])
    cur.loglik <- e.step[["loglik"]]
    loglik.vector <- e.step[["loglik"]]
  } else {
    # Repeat E and M steps till convergence
    e.step <- e_step(sim, m.step[["mu"]], sqrt(m.step[["var"]]), 
                     m.step[["alpha"]])
    m.step <- m_step(sim, e.step[["posterior.df"]])
    loglik.vector <- c(loglik.vector, e.step[["loglik"]])
    
    loglik.diff <- abs((cur.loglik - e.step[["loglik"]]))
    if(loglik.diff < 1e-6) {
      break
    } else {
      cur.loglik <- e.step[["loglik"]]
    }
  }
}


tail(loglik.vector)
m.step

#Kim et.al 1998 mixture
#mixture
q = c(	0.04395,	0.24566,	0.34001,	0.25750,	0.10556,	0.00002,	0.00730  )
b = c(  2.77786,  1.79518,  0.61942, -1.08819, -3.97281, -8.56686, -10.12999 )-1.2704
w = c(	0.16735,	0.34023,	0.64009,	1.26261,	2.61369,	5.17950,	5.79596  )

sort(q)
sort(c(m.step$alpha))

#or via LaplacesDemon Package
library(LaplacesDemon)
n = 10000
plot(density(log(rnorm(n,0,1)^2)),col = 'red', lwd=2.5, lty=2)
lines(density(rnormm(n, q , b , sqrt(w))) , col ="green")
lines(density(rnormm(n, m.step$alpha , m.step$mu , sqrt(m.step$var))), col = "blue")
lines(density(rnorm(1000000,-1.2704 ,sqrt(4.93) )), col = 'black', lwd=3, lty=2)

legend("topleft",                                     
       legend = c("chisq(1) distr",
                  "mixture of 7 normals variables used by Kim et.al 1998",
                  "mixture of 7 normals variables estimated by this code",
                  " normal with -1.2704 and 4.93 ") ,
       col = c("red", 'green','blue','black'),
       lwd=3, ncol = 1,  cex = 0.72)

mean(log(rnorm(n,0,1)^2))
var(log(rnorm(n,0,1)^2))
mean(rnormm(n, q , b , sqrt(w)))
var(rnormm(n, q , b , sqrt(w)))
mean(rnormm(n, m.step$alpha , m.step$mu , sqrt(m.step$var)))
var(rnormm(n, m.step$alpha , m.step$mu , sqrt(m.step$var)))

