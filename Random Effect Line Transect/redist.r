# 21 Nov 2012; J. Laake
# Script to simulate data from random scale 
# half-normal detection function for line transect data and
# fit random and fixed scale detection functions to the data.
# simdata creates simulated data 
# fitdata fits the model to the data and returns model;
# if beta_eps set to NULL the fixed scale model is fitted.
# if wrong=TRUE it shows results with incorrect likelihood
# plotfit uses histline from mrds package to show fit.
#
library(mrds)
# Half-normal detection function
gx=function(x,sigma)
{
  exp(-.5*(x/sigma)^2)
}
# Half-normal * density(eps)
gx.eps=function(x,eps,par)
{
    sigma=exp(par[1]+eps*exp(par[2]))
    gx(x,sigma)*dnorm(eps,0,1)
}
# Average half-normal integrated over distribution of eps
avg_gx=function(x,par,weps=5)
{
    gx=vector("numeric",length=length(x))
    geps.x=function(eps,x,par)gx.eps(x,eps,par)
    for(i in 1:length(x))
     gx[i]=integrate(geps.x,-weps,weps,x=x[i],par=par)$value
    gx
}
# Average detection probability * density eps
mu=function(eps,w,par)
{
  mu=vector("numeric",length=length(eps))
  for(i in 1:length(eps))
  {
    sigma=exp(par[1]+eps[i]*exp(par[2]))
    mu[i]=sigma*sqrt(2*pi)*(pnorm(w,0,sigma)-.5)*dnorm(eps[i],0,1)
  }
  mu
}
# Average detection probability - integrated over eps
avg_mu=function(w,par,weps=5)
	integrate(mu,-weps,weps,w=w,par=par)$value
# Negative log-likelihood function 
flnl=function(par,x,w,weps=5,wrong=FALSE)
{
  lnl=0
  for(i in 1:length(x))
  {
     if(length(par)>1)
		 if(wrong)
			 lnl=lnl-log(integrate(fx.eps,-weps,weps,x=x[i],w=w,par=par)$value)
		 else
             lnl=lnl-log(avg_gx(x[i],par)/avg_mu(w,par,weps))
	 else
	 {
		 sigma=exp(par[1])
		 mu=sigma*sqrt(2*pi)*(pnorm(w,0,sigma)-.5)
		 lnl=lnl-log(gx(x[i],sigma=sigma))+log(mu)
	 }
  }
  lnl
}
# pdf - for incorrect likelihood 
fx=function(x,sigma,w)
{
	mu=sigma*sqrt(2*pi)*(pnorm(w,0,sigma)-.5)
	exp(-.5*(x/sigma)^2)/mu
}
# pdf * density(eps) - for incorrect likelihood 
fx.eps=function(eps,x,par,w)
{
	sigma=exp(par[1]+eps*exp(par[2]))
	fx(x,sigma,w)*dnorm(eps,0,1)
}

# n        - sample size
# beta     - scale intercept
# beta_eps - beta for random effect sigma
# w        - width of transect
# reject   - if TRUE uses brute force rejection sampling to 
#            generate distances;if w=Inf, it uses a large w 
#            to avoid any appreciable truncation
# b        - number of deviates generated in a batch for 
#             rejection sampling
#         
# If reject==FALSE, it generates distances directly with 
# rnorm but if this approach is used, the estimated beta with the
# correct likelihood will not match the generating beta but 
# this is not bias in the usual sense. It will match the beta
# from the incorrect likelihood.
simdata=function(n=500,beta=2,beta_eps=-1,w=Inf,reject=TRUE,b=10000)
{
	if(reject)
	{
		if(w==Inf)
			w=3*exp(beta+3*exp(beta_eps))
		x=NULL
        while(length(x)<n)
		{
			u=runif(b,0,w)
			sigma=exp(beta+rnorm(b,0,1)*exp(beta_eps))
			seen=gx(u,sigma)>runif(b,0,1)
			x=c(x,u[seen])
		}
		return(x[1:n])
	}else
	{
		sigma=exp(beta+rnorm(n,0,1)*exp(beta_eps))
		x=abs(rnorm(n,0,sigma))
		return(x)
	}
}

plotfit=function(result,nclass=NULL,weps=5,main=NULL)
{
# Create plot of fit; uses histline from mrds
par=result$model$par
if(length(par)<2)
{
	# compute fixed mu
	avg_mu_est=exp(par)*sqrt(2*pi)*(pnorm(result$model$w,0,exp(par))-.5)
}else
{
	# Compute average_mu
	avg_mu_est=avg_mu(w=result$model$w,par=par,weps=weps)
}
# Compute Nhat if W<Inf
if(result$model$w!=Inf)
{
	Nhat=length(x)/(avg_mu_est/result$model$w)
}else
	Nhat=NULL
result$Nhat=Nhat
x=result$model$x
max_x=max(x)
if(is.null(nclass))
	nints=ceiling(sqrt(length(x)))
else
	nints=nclass
int_width=max_x/nints
breaks=int_width*(0:nints)
ints=(0:100)*max_x/100
hh=hist(x,plot=F,breaks=breaks)
bars=(hh$counts/(sum(hh$counts)*int_width))*avg_mu_est
if(length(par)<2)
{
	gx=gx(ints,sigma=exp(par))
}else
	gx=avg_gx(ints,par=par)
mrds:::histline(bars,breaks,ylim=c(0,max(c(1,bars))),
		xlab="Distance",ylab="Detection probability",main=main)
lines(ints,gx)
return(result)
}
# examples
set.seed(123)
x=simdata(n=500,w=25,beta=1,beta_eps=-.5)
par(mfrow=c(1,2))
results_random=fitdata(x,w=Inf,beta_eps=-.5)
results_random=plotfit(results_random,nclass=30)
results_random_wrong=fitdata(x,w=Inf,beta_eps=-.5,wrong=TRUE)
results_random_wrong=plotfit(results_random_wrong,nclass=30)

x=simdata(n=500,w=Inf,beta_eps=-.5,reject=FALSE)
par(mfrow=c(1,2)) 
results_random=fitdata(x,w=Inf,beta_eps=-.5)
results_random=plotfit(results_random,nclass=30)
results_random_wrong=fitdata(x,w=Inf,beta_eps=-.5,wrong=TRUE)
results_random_wrong=plotfit(results_random_wrong,nclass=30)


results_fixed=fitdata(x,w=Inf,beta=1,beta_eps=NULL)
win.graph()
results_fixed=plotfit(results_fixed,main="Fixed scale R fit")
#
# ADMB section --doesn't work at present
#
prepare_admb()
compile_admb("hnre",re=TRUE)
# write out data file for untruncated analysis with ADMB
con=file("hnre.dat",open="wt")
write(length(x),con,append=FALSE)
write(-1,con,append=TRUE)
write(x,con,ncol=1,append=TRUE)
close(con)
# run model and get results
run_admb("hnre")
admb_model_untrunc=read_admb("hnre")
admb_model_untrunc$model$par=coef(admb_model_untrunc)
admb_model_untrunc$model$par[2]=log(admb_model_untrunc$model$par[2])
admb_model_untrunc$model$w=Inf
admb_model_untrunc$model$x=x
win.graph()
admb_model_untrunc=plotfit(admb_model_untrunc,main="Random scale ADMB fit")

# write out data file for untruncated analysis with ADMB
truncx=x[x<=25]
con=file("hnre.dat",open="wt")
write(length(truncx),con,append=FALSE)
write(25,con,append=TRUE)
write(truncx,con,ncol=1,append=TRUE)
close(con)
run_admb("hnre")
admb_model_trunc=read_admb("hnre")
admb_model_trunc=read_admb("hnre")
admb_model_trunc$model$par=coef(admb_model_trunc)
admb_model_trunc$model$par[2]=log(admb_model_trunc$model$par[2])
admb_model_trunc$model$w=25
admb_model_trunc$model$x=truncx
win.graph()
admb_model_trunc=plotfit(admb_model_trunc,main="Random scale ADMB w=25 fit")



