# install.packages("pracma")
library(pracma)

## define kappa in Theorem 1 of the paper 
kappa1 <- function(theta1,mu,s) {- mu + theta1/2 - sqrt( (mu - theta1/2)^2 +2*s) };  
kappa2 <- function(theta1,mu,s) {- mu + theta1/2 + sqrt( (mu - theta1/2)^2 +2*s) };
Fcomp <-function(theta,theta1, mu,s)
  { kappa2(theta1,mu, s) * exp(kappa1(theta1,mu,s)*theta) - kappa1(theta1,mu,s) * exp(kappa2(theta1,mu,s)*theta)  }


## True parameters in the model
theta0 = 0.1;
theta1 = 1;
theta2 = 5;

## (1) Under the null of no stimulus, i.e., mu=0
## (1a) False alarm rate or CDF 
mu = 0;
F0 <- function(s) { Fcomp(theta0, theta1,mu,s) / (s* Fcomp(theta2,theta1,mu,s)) }

f0_vals <-  invlap(F0, 10,2000,1000)
plot(f0_vals$x, f0_vals$y,
     type = "l",
     col = "black",
     lwd = 3,
     main = expression(paste("Error Probability under ", mu == 0)),
     xlab = "t",
     ylab = expression('P'(T <= t)))
grid()  

# (1) PDF under the null of no stimulus
mu = 0;
Fs0 <- function(s) { Fcomp(theta0, theta1,mu,s) / Fcomp(theta2,theta1,mu,s) }

# Inverse Laplace numerically at some t values
f0_vals <-  invlap(Fs0, 10,2000,1000)
plot(f0_vals$x, f0_vals$y,
     type = "l",
     col = "black",
     lwd = 3,
     main = expression(paste("PDF under ", mu == 0)),
     xlab = "t",
     ylab = expression('f'(t)))
grid()  

# (2) PDF under the critical stimulus of mu = 0.5
mu = 0.5;
Fs2 <- function(s) { Fcomp(theta0, theta1,mu,s) / Fcomp(theta2,theta1,mu,s) }
f2_vals <-   invlap(Fs2,  1, 100, 1000)

plot(f2_vals$x, f2_vals$y,
     type = "l",
     col = "blue",
     lwd = 3,
     main = expression(paste("PDF under ", mu == 0.5)),
     xlab = "t",
     ylab = expression(f(t)))
grid()


# (3) PDF under the stimulus of mu=1
mu = 1;
Fs1 <- function(s) { Fcomp(theta0, theta1,mu,s) / Fcomp(theta2,theta1,mu,s) }
f1_vals <-   invlap(Fs1,  1, 30, 1000)

plot(f1_vals$x, f1_vals$y,
     type = "l",
     col = "red",
     lwd = 3,
     main = expression(paste("PDF under ", mu == 1)),
     xlab = "t",
     ylab = expression('f'(t)))
grid()  




