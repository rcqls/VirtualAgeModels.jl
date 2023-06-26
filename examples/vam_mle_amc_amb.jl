using VAM
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
mle(m, AMC_Amb)
VAM.params(m)
sum(abs.(res .- params(m)))