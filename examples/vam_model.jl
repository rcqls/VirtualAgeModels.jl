using VAM
using Distributions

vam = @vam(time & type ~ (ARAInf(0.4) | Weibull(0.001,2.5)))
VAM.init!(vam)

ex_f = :(time & type ~ (ARAInf(0.4) | Weibull(0.001,2.5)) )
ex_f1 =:(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)) & (AGAN()))
ex_f2 = :(Temps & Type ~ (ARA1(.5) | Weibull(0.01,2.5)) & (ARAInf(.7)+ARAInf(.3)+ ABAO()|Periodic(12,[0.6,0.4]) * AtIntensity(1.2)))
ex_f.head
ex_f.args
ex_f.args[3].head
ex_f.args[3].args
ex_f2.args[3].args
VAM.parse_model(ex_f)
ex_f1.head
ex_f1.args
ex_f1.args[2].head
Meta.isexpr(ex_f.args[2].args[2], :call)
ex_f.args[2]
Meta.isexpr(ex_f1.args[2].args[2], :call)
ex_f1.args[2].args[2].args

vars = if Meta.isexpr(ex_f.args[2].args[2], :call)
    vcat(ex_f.args[2].args[2].args[2:3], ex_f.args[2].args[3])
else
    ex_f.args[2].args[2:3]
end

vars1 = if Meta.isexpr(ex_f1.args[2].args[2], :call)
    vcat(ex_f1.args[2].args[2].args[2:3], ex_f1.args[2].args[3])
else
    ex_f1.args[2].args[2:3]
end
VAM.parse_model(ex_f1)
VAM.isbayesian(VAM.parse_model(ex_f2))

ex_f_b = :(Time & Type ~ (ARAInf(~Uniform()) | Weibull(~Uniform(),~Uniform(2,4))))
m = VAM.parse_model(ex_f_b)
m.models[1].priors[1]
VAM.isbayesian(m)


res = VAM.parse_covariates(:(Weibull(0.001,2.5)))
length(res)
res[1] 

ex_f_cov = :(time & type ~ (ARAInf(0.4) | Weibull(0.001,2.5, .3| 1*cov1)))
#VAM.parse_model(ex_f_cov)
ex_fm = ex_f_cov.args[3].args[3]
res = VAM.parse_covariates(ex_fm)
length(res)
res[1]
res[2]

ex_fm.args
index = findall(e -> Meta.isexpr(e,:call) && e.args[1] in [:|,:+], ex_fm.args)[1]
vcat(ex_fm.args[1:index-1],ex_fm.args[index].args[2],Expr(:call,:+,ex_fm.args[index].args[3]))
 
ex_f_cov = :(time & type ~ (ARAInf(0.4) | Weibull(0.001,2.5| 1*cov1 + -2cov2 + 3cov3)) )
ex_fm = ex_f_cov.args[3].args[3]
res = VAM.parse_covariates(ex_fm)
length(res)
res[1]
res[2]

VAM.parse_covariates(:(Weibull(0.15,2.3|0.6*cov1-0.9*cov2)))