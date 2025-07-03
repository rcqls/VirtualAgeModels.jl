using VirtualAgeModels
using DataFrames
#using RData
using RCall
#using Plots

R"""
require(VAM)
data("AMC_Amb")
AMC_mle<-mle.vam(Time & Type ~ (ARAInf(0.6)|Weibull(1,3)),data=AMC_Amb)
res <- coef(AMC_mle)
"""

@rget AMC_Amb
@rget res

@macroexpand  @vam(Time & Type ~ (ARA∞(0.6) | Weibull(1.0,3.0)))
m = @vam Time & Type ~ (ARA∞(0.6) | Weibull(1.0,3.0)) # data=AMC_Amb
m = @vam Time & Type ~ (ARA∞(0.6) | Weibull(1.0,3.0)) data=AMC_Amb
data!(m, AMC_Amb)
data(m)
ml = mle(m, AMC_Amb)
ml.optim
params(ml)

m2 = @vam Time & Type ~ (ARA∞(0.6) | Weibull(1.0,3.0)) data=AMC_Amb
ml2 = mle(m2)
params(ml2)
# plot(m, VirtualAgeModels.VirtualAgePlot)
# plot(m, VirtualAgeModels.IntensityPlot)
# plot(m, VirtualAgeModels.CummulativeIntensityPlot)
# plot(m, VirtualAgeModels.ConditionalCummulativeDensityFunctionPlot)
# plot(m, VirtualAgeModels.ConditionalSurvivalPlot)
# plot(m, VirtualAgeModels.ConditionalProbabilityDensityFunctionPlot)

# AMC_Amb
# ex = @macroexpand @vam(Time & Type ~ (AGAN() | Weibull(1.0,3.0)), data = AMC_Amb)
# eval(ex)
# m2 = @vam(Time & Type ~ (AGAN() | Weibull(1.0,3.0)), data = AMC_Amb)
# plot(m2)

# VirtualAgeModels.symbol2typeplot
plot(m, :v)
plot!(m, :i)
plot(m, :i)
plot(m, :I)
plot(m, :F)
plot(m, :f)
plot(m, :S)