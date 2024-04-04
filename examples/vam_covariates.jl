using VirtualAgeModels
using DataFrames
using Distributions
using RCall
using Optim

R"""
	require(VAM)

    dataDF<-data.frame(System=c(rep(1,5),rep(2,4),rep(3,1),rep(4,5)),Time=c(0.800,2.646,3.190,3.916,4.109,0.910,1.127,1.245, 1.349, 0.541,1.397,1.463,2.406,2.506,3.159),Type=rep(-1,15))
    dataCov<-data.frame(cov1=c(4.336,5.615,4.770,4.655),cov2=c(0,0,0,1))
    mle<-mle.vam(System&Time&Type~(ARAInf(0.4) | Weibull(0.001,1.5| 1*cov1 + -2*cov2)),data=dataDF,data.covariates = dataCov)
    theta<-c(0.005,2.0,0.5,1.0,-2.0)
	res <- list()
	res$lnL <- logLik(mle,theta,TRUE,FALSE,FALSE)
	res$dlnL<- logLik(mle,theta,FALSE,TRUE,FALSE)
	res$d2lnL <- logLik(mle,theta,FALSE,FALSE,TRUE)
	res$C <- contrast(mle,theta,TRUE,FALSE,FALSE)
	res$dC <- contrast(mle,theta,FALSE,TRUE,FALSE)
	res$d2C <- contrast(mle,theta,FALSE,FALSE,TRUE)
	"""
res = @rget res
dataCov = @rget dataCov
dataDF= @rget dataDF
θ =  @rget theta

m = @vam( System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,1.5| 1*cov1 + -2*cov2)), data=dataDF, datacov=dataCov)

lnL = res[:lnL]
dlnL = res[:dlnL]
d2lnL = res[:d2lnL]
contrast(m, θ, profile=false)-lnL
gradient(m, θ, profile=false)-dlnL
hessian(m, θ, profile=false)-d2lnL
C = res[:C]
dC = res[:dC]
d2C = res[:d2C]
contrast(m,θ)-C
gradient(m, θ)[2:end]-dC
hessian(m, θ)[2:end,2:end]-d2C


contrast(m, θ, profile=false)
gradient(m, θ, profile=false)





m=5
dfcov=DataFrame(cov1=rand(Uniform(),m), cov2=rand(Uniform(),m), cov3=rand(Uniform(),m))

m = @vam( Time & Type ~ (ARAInf(0.4) | Weibull(0.001,1.5| 1*cov1 + -2*cov2 + 3*cov3)), datacov=dfcov)
df = rand(m, 29)

ml = mle(m,df,dfcov)
contrast(ml)
data(ml)
VirtualAgeModels.covariates(ml.model)
contrast(ml,[0.001,0.5,0.4,1.0,-2.0,3.0])
gradient(ml)
params(ml)
m.params_cov
m.vars_cov
m.datacov

R"""
require(VAM)
dataDF<-data.frame(System=c(rep(1,5),rep(2,4),rep(3,1),rep(4,5)),Time=c(0.800,2.646,3.190,3.916,4.109,0.910,1.127,1.245, 1.349, 0.541,1.397,1.463,2.406,2.506,3.159),Type=rep(-1,15))
dataCov<-data.frame(cov1=c(4.336,5.615,4.770,4.655),cov2=c(0,0,0,1))
mle<-mle.vam(System&Time&Type~(ARAInf(0.8)|Weibull(0.15,2.3|0.6*cov1-0.9*cov2)),data=dataDF,data.covariates = dataCov)
theta<-c(0.15,2.3,0.8,0.6,-0.9)
res = VAM:::contrast.mle.vam(mle,par0=theta)
cres <- coef(mle)
"""
@rget dataCov
@rget dataDF
@rget res
@rget cres

mbis = @vam(System & Time & Type ~ (ARAInf(0.8)|Weibull(0.15,2.3|0.6cov1 + -0.9cov2)))
ml = mle(mbis,dataDF,dataCov;method=NelderMead())
mlbis = mle(mbis,dataDF,dataCov)

params(ml)
sum(params(ml) .- cres)
contrast(ml)

m2 = @vam( time & type ~ (ARAInf(0.4) | Weibull(0.001,1.5)) )
df2 = rand(m2, 29)
contrast(m2, df2)
gradient(m2, df2)