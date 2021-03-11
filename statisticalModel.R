# Load Packages and Data
# Required packages
library(tidyverse)
library(qqplotr)
library(GGally)
library(coda)
library(kableExtra)

# Import data
X<-read.csv("D:/statistics/X.csv",header=FALSE)
Y<-read.csv("D:/statistics/Y.csv",header=FALSE)
T<-read.csv("D:/statistics/Time.csv",header=FALSE)
eeg_data<-bind_cols(T,Y,X)
names(eeg_data)<-c("t","y","x1","x2","x3","x4")
eeg_data %>% head (n=10) %>%
kbl() %>%
  kable_paper() %>%
  scroll_box(width = "900px", height = "250px")

#----------------- Task 1: Preliminary data analysis  -----------------------
# Time series plots (input and output EEG signals)

out_plot <- ggplot(eeg_data, aes(x=t, y=y)) +
            geom_line(color = "darkred") + 
            xlab("Time")+ ylab("y")+ ggtitle( "Output EEG signals")
out_plot

inpt_df <- eeg_data %>%
           select(-y) %>%
           gather(key = "variable", value = "value", -t)

inpt_plot <-ggplot(inpt_df, aes(x = t, y = value)) + 
            geom_line(aes(color = variable)) + 
            xlab("Time")+ ggtitle( "Input EEG signals")
inpt_plot


# Distribution for each EEG signal

sgnl_data<-eeg_data %>%
           select(-t)
sgnl_hist <- ggplot(gather(sgnl_data), aes(value)) 
sgnl_hist + geom_histogram(bins = 10) + facet_wrap(~key, scales = "free_x")


# Correlation and scatter plots (between different input EEG signals and the output EEG)

ggpairs(sgnl_data, title="Correlation and scatterplots among EEG signals") 


#----- Task 2: Regression- modelling the relationship between EEG signals -----

# predict data using the given models
sgnl_data %>%
mutate(xx1=x4,xx2=x1^2,xx3=x1^3,xx4=x3^4,xx0=1) %>%
select(xx1,xx2,xx3,xx4,xx0)-> M1
M1 %>% head (n=10) %>%
kbl() %>%
  kable_paper() %>%
  scroll_box(width = "900px", height = "250px")

sgnl_data %>%
mutate(xx1=x3^3,xx2=x3^4,xx0=1) %>%
select(xx1,xx2,xx0)-> M2
M2 %>% head (n=10) %>%
kbl() %>%
  kable_paper() %>%
  scroll_box(width = "900px", height = "250px")

sgnl_data %>%
mutate(xx1=x2,xx2=x1^3,xx3=x3^4,xx0=1) %>%
select(xx1,xx2,xx3,xx0)-> M3
M3 %>% head (n=10) %>%
kbl() %>%
  kable_paper() %>%
  scroll_box(width = "900px", height = "250px")

sgnl_data %>%
mutate(xx1=x4,xx2=x1^3,xx3=x3^4,xx0=1) %>%
select(xx1,xx2,xx3,xx0)-> M4
M4 %>% head (n=10) %>%
kbl() %>%
  kable_paper() %>%
  scroll_box(width = "900px", height = "250px")

sgnl_data %>%
mutate(xx1=x4,xx2=x1^2,xx3=x1^3,xx4=x3^4,xx5=x1^4,xx0=1) %>%
select(xx1,xx2,xx3,xx4,xx5,xx0)-> M5
M5 %>% head (n=10) %>%
kbl() %>%
  kable_paper() %>%
  scroll_box(width = "900px", height = "250px")

# Task 2.1

ols<-function(X,y)
     {
       X<-as.matrix(X)
       Y<-as.matrix(y)
       X.X<-solve(t(X) %*% X)       
       b<-b<-X.X %*% (t(X) %*% Y)
       colnames(b)<- "OLS_parm" 
       return(b)}

theta_M1<-ols(M1,Y)
theta_M2<-ols(M2,Y)
theta_M3<-ols(M3,Y)
theta_M4<-ols(M4,Y)
theta_M5<-ols(M5,Y)
# Theta for Model-1
theta_M1
# Theta for Model-2
theta_M2
# Theta for Model-3
theta_M3
# Theta for Model-4
theta_M4
# Theta for Model-5
theta_M5

# Task 2.2

RSS<- function(X,y,theta)
      { 
       X<-as.matrix(X)
       Y<-as.matrix(y)
       B<-as.matrix(theta)
       Yhat<- X %*% B
       E<-Y-Yhat
       rss<-as.numeric(sum(E^2))
       return(rss)}

rss1<-RSS(M1,Y,theta_M1)
rss2<-RSS(M2,Y,theta_M2)
rss3<-RSS(M3,Y,theta_M3)
rss4<-RSS(M4,Y,theta_M4)
rss5<-RSS(M5,Y,theta_M5)
rss_all<-data.frame(Model=paste("Model",1:5, sep="_"),RSS=c(rss1,rss2,rss3,rss4,rss5))
rss_all %>% 
kbl() %>%
  kable_paper() %>%
  scroll_box(width = "900px", height = "200px")

# Task 2.3

log_L<-function(rss,n)
       {
        sigma2<-rss/(n-1)
        L<- -(n/2)*log(2*pi) - (n/2)*log(sigma2) -(1/2*sigma2)*rss
        return(L)}

L1<-log_L(rss1,nrow(sgnl_data))
L2<-log_L(rss2,nrow(sgnl_data))
L3<-log_L(rss3,nrow(sgnl_data))
L4<-log_L(rss4,nrow(sgnl_data))
L5<-log_L(rss5,nrow(sgnl_data))
L_all<-data.frame(Model=paste("Model",1:5, sep="_"),Likelihood=c(L1,L2,L3,L4,L5))
L_all %>% 
kbl() %>%
  kable_paper() %>%
  scroll_box(width = "900px", height = "200px")

# Task 2.4 

AIC_BIC<-function(logL,theta,n)
         {
          k<-length(as.vector(theta))
          aic<- 2*k - 2*logL
          bic<- k*log(n) - 2*logL
          return(c(aic,bic))}

ab1<-AIC_BIC(L1,theta_M1,nrow(sgnl_data))
ab2<-AIC_BIC(L2,theta_M2,nrow(sgnl_data))
ab3<-AIC_BIC(L3,theta_M3,nrow(sgnl_data))
ab4<-AIC_BIC(L4,theta_M4,nrow(sgnl_data))
ab5<-AIC_BIC(L5,theta_M5,nrow(sgnl_data))
aic_bic_all<-data.frame(Model=paste("Model",1:5, sep="_"),AIC=c(ab1[1],ab2[1],ab3[1],ab4[1],ab5[1]),BIC=c(ab1[2],ab2[2],ab3[2],ab4[2],ab5[2]))
aic_bic_all %>% 
kbl() %>%
  kable_paper() %>%
  scroll_box(width = "900px", height = "200px")

# Task 2.5

errors<- function(X,y,theta)
      { 
       X<-as.matrix(X)
       Y<-as.matrix(y)
       B<-as.matrix(theta)
       Yhat<- X %*% B
       E<-Y-Yhat
       return(as.numeric(E))}

mod_res<-bind_cols(errors(M1,Y,theta_M1),
                        errors(M2,Y,theta_M2),
                        errors(M3,Y,theta_M3),
                        errors(M4,Y,theta_M4),
                        errors(M5,Y,theta_M5))
names(mod_res)<-paste("Error_M",1:ncol(mod_res),sep="")
mod_res %>% head (n=10) %>%
kbl() %>%
  kable_paper() %>%
  scroll_box(width = "900px", height = "250px")

# Error distribution
error_hist <- ggplot(gather(mod_res), aes(value)) 
error_hist + geom_histogram(bins = 10) + facet_wrap(~key, scales = "free_x")

# qq-plot for all model residual
qq_plot <- ggplot(data = gather(mod_res), mapping = aes(sample = value, color = key, fill = key)) +
           stat_qq_band(alpha=0.5) +
           stat_qq_line() +
           stat_qq_point() +
           facet_wrap(~ key) +
           labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
qq_plot

# Task 2.6

# Model 3 is better model because its have smaller AIC and BIC values and the residual of this model are normal. 

# Task 2.7

           eeg_data %>%
           mutate(xx1=x2,xx2=x1^3,xx3=x3^4) %>%
           select(t,y,xx1,xx2,xx3)-> bst_data

dt <-sort(sample(nrow(bst_data), nrow(bst_data)*.7))
train<-bst_data[dt,]
test<-bst_data[-dt,]

bst_mod<- lm(y ~ xx1+xx2+xx3 , data=train)
fit_CI <- predict(bst_mod, newdata=test,interval="confidence",level = 0.95)
se <- predict(bst_mod, newdata=test, se.fit=TRUE)$se.fit

fit_df <- cbind(test, fit_CI,se)
fit_df %>% head (n=10) %>%
kbl() %>%
  kable_paper() %>%
  scroll_box(width = "900px", height = "250px")

# Plot for prediction
ggplot(fit_df, aes(t, y))+
    geom_point() +
    geom_line(aes(y=fit), color = "blue", linetype = "solid")+
    geom_line(aes(y=lwr), color = "red", linetype = "dashed")+
    geom_line(aes(y=upr), color = "red", linetype = "dashed")+
    geom_errorbar(aes(ymax = fit+se, ymin = fit-se), width = 0.009)+
    xlab("Time")

# Task 3

# (1)

# Parameters with largest absolute values in your least squares estimation
larg_par<-order(abs(theta_M3),decreasing =TRUE)[1:2]
two_par<-theta_M3[larg_par,] # theta_bias & theta_1

# (2)

# Uniform prior
theta_bs<- function () {
  return (runif(1, min=-0.5, max=0.5))
}
theta_1 <- function () {
  return (runif(1, min=-0.5, max=0.5))
}

# The function to simulate a data set 
sim_data <- function (X,theta_bs, theta_1,B) { 
         larg_par<-order(abs(B),decreasing =TRUE)[1:2]
         B[larg_par]<-c(theta_bs, theta_1)
         Yhat<- as.matrix(X) %*% as.matrix(B)
         return(as.vector(Yhat))}

# Function to compute the quantiles, We choose to use 3 quantiles.
  comp_qntls <- function(data) {
  return (quantile(data, probs=c(0.1, 0.5, 0.9)))
}

# Distance to compare a simulated sample to the observed data
com_data <- function (true, simulated) {
  distance <- sqrt(sum(mapply(function(x,y) (x-y)^2, true, simulated)))
  return(distance)
}

# Accept or reject based on the threshold values
decision <- function (true, simulated, tslod) {
  distance<- com_data(comp_qntls(true), comp_qntls(simulated))
  if((distance < tslod) ) {
    return(TRUE) 
   }else return(FALSE)
}
# Function for rejection ABC rule
rejection <- function (obs_data,X,B, n_iter, tslod, dec_function) {
  obs_data<- as.vector(unlist(obs_data))
  n<- length(obs_data)
  acp_rej <- vector(length = n_iter)
  draw_theta0 <- vector(length = n_iter, mode = "numeric")
  draw_theta1 <- vector (length = n_iter, mode = "numeric")
  for (i in 1:n_iter){
    theta0 <- theta_bs()
    theta1 <- theta_1()
    parameters <- list("theta0"=theta0, "theta1"=theta1 )
    syn_data <- sim_data(X,theta0, theta1,B)
    acp_rej[i] <- dec_function(obs_data, syn_data, tslod)
    draw_theta0[i] <- theta0
    draw_theta1[i] <- theta1
  }
  return(data.frame(cbind("acp_rej" = acp_rej, "draw_theta0" = draw_theta0, "draw_theta1" = draw_theta1)))
}


# (3)

res<-rejection(Y,M3,theta_M3, 10000, 0.5, decision )
res %>% head (n=10) %>%
kbl() %>%
  kable_paper() %>%
  scroll_box(width = "900px", height = "250px")

# How many samples have been accepted out of 10000?
sum(res$acp_rej)

# (4) 

post_mcmc <- mcmc(res[which(res$acp_rej==1),c(2,3)])
summary(post_mcmc)
# Plotting the posterior distributions
plot(post_mcmc)