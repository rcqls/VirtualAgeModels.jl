using VirtualAgeModels
using DataFrames
#using RData
using RCall

R"""
require(VAM)
data("AMC_Amb")
AMC_mle<-mle.vam(Time & Type ~ (ARAInf(0.6)|Weibull(1,3)),data=AMC_Amb)
res <- coef(AMC_mle)
"""

@rget AMC_Amb
@rget res

m = @vam(Time & Type ~ (ARA∞(0.6) | Weibull(1.0,3.0)))
ml = mle(m, AMC_Amb)
ml.model
data(m)
data(ml)
params(m)
params(ml)
sum(abs.(res .- params(ml)))

m2 = @vam Time & Type ~ (ARA∞(0.6) | Weibull(1.0,3.0)) data=AMC_Amb
data(m2)
ml2 = mle(m2)
data(m2)
params(m2)
params(ml2)
sum(abs.(res .- params(ml2)))