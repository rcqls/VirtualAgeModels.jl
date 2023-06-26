
## RCpp: sim.vam(  ~ (ARA1(.9) | Weibull(.001,2.5)) & (ARAInf(.4) | AtIntensity(0.2)))
m = @vam(Temps & Type ~ (ARA1(.9) | Weibull(0.001,2.5)) & (ARAInf(.7) | AtIntensity(0.2)))
p = params(m)

# Random.seed!(3)
# df = simulate(m, @stop(s < 10))
# freqtable(df.type)

# df
#@run mle(m, params(m), df; method=Newton())
# if false
#     res = mle(m, params(m), df; fixed=[1, 3])
#     res = mle(m, params(m), df)
#     res.params
#     # res = mle(m, params(m), df, method=LBFGS())
#     contrast(m, params(m), df)
#     gradient(m, params(m), df)
#     hessian(m, params(m), df)
# end

R"""
require(VAM)
#set.seed(3)
simCMPM<-sim.vam( time & type ~ (ARA1(.9) | Weibull(.001,2.5)) & (ARAInf(.7) | AtIntensity(0.2)))
df<-simulate(simCMPM, Size>=100)
"""

@rget df

R"""
mleCMPM <- mle.vam( time & type ~ (ARA1(.9) | Weibull(.001,2.5)) & (ARAInf(.7) | AtIntensity(0.2)), data=df)
theta<-c(.001,2.5,.9,.7)
list(logLik(mleCMPM,theta,TRUE,FALSE,FALSE),
logLik(mleCMPM,theta,FALSE,TRUE,FALSE),
contrast(mleCMPM,theta,TRUE,FALSE,FALSE),
contrast(mleCMPM,theta,FALSE,TRUE,FALSE))
"""

res = mle(m, p, df)
res.params
contrast(m, p, df; profile=false)
gradient(m, p, df; profile=false)
hessian(m, p, df; profile=false)
contrast(m, p, df)
gradient(m, p, df)
hessian(m, p, df)

