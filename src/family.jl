abstract type FamilyModel end

## TODO: FamilyCompute NOT USED from now, used in the non yet implemented  

mutable struct FamilyCompute
    dHR::Vector{Float64}
    dHL::Vector{Float64}
    dhR::Vector{Float64}
    dhL::Vector{Float64}
    d2HR::Vector{Float64}
    d2HL::Vector{Float64}
    d2h::Vector{Float64}
    dhd::Vector{Float64}
    FamilyCompute() = new()
end

function init!(fc::FamilyCompute, fm::FamilyModel)
    nbd = nbparams(fm) -1
    fc.dHR = zeros(nbd)
    fc.dHL = zeros(nbd)
    fc.dhR = zeros(nbd)
    fc.dhL = zeros(nbd)
    fc.dhd = zeros(nbd)
    nbd2 = (nbd + 1) * nbd ÷ 2
    fc.d2HR = zeros(nbd2)
    fc.d2HL = zeros(nbd2)
    fc.d2h = zeros(nbd2)
end
mutable struct WeibullFamilyModel <: FamilyModel
    α::Parameter #Float64
    β::Parameter #Float64
    comp::FamilyCompute
    priors::Priors
    WeibullFamilyModel() = new()
end
function WeibullFamilyModel(α::Parameter, β::Parameter)
    fm = WeibullFamilyModel()
    fm.α, fm.β =  α, β
    fm.priors = [nothing, nothing]
    fm.comp = FamilyCompute()
    init!(fm.comp, fm)
    return fm
end
params(fm::WeibullFamilyModel)::Parameters = [fm.α,fm.β]
params!(fm::WeibullFamilyModel, p::Parameters) = begin; fm.α,fm.β = p; nothing; end
nbparams(fm::WeibullFamilyModel) = 2
function WeibullFamilyModel(α::Prior, β::Prior)
    fm = WeibullFamilyModel(0.0, 0.0)
    fm.priors = [α, β]
    return fm
end

hazard_rate(f::WeibullFamilyModel, x::Float64)::Float64 = (x<=0 ? 0 : f.α * f.β * x^(f.β - 1) )
inverse_hazard_rate(f::WeibullFamilyModel, x::Float64)::Float64 = (x<=0 ? 0 : (x/f.α/f.β)^(1/(f.β-1)))
cumulative_hazard_rate(f::WeibullFamilyModel, x::Float64)::Float64 = (x<=0 ? 0 : f.α * x^f.β)
inverse_cumulative_hazard_rate(f::WeibullFamilyModel, x::Float64)::Float64 = (x<=0 ? 0 : (x/f.α)^(1/f.β))
hazard_rate_derivative(f::WeibullFamilyModel, x::Float64)::Float64 = (x<=0 ? 0 : f.α * f.β * (f.β - 1) * x^(f.β - 2) )

function hazard_rate_param_derivative(f::WeibullFamilyModel, x::Float64,right::Bool)::Vector{Float64}
    dh = right ? f.comp.dhR : f.comp.dhL
    dh[1] = x<=0 ? 0 : f.α * (1 + f.β * log(x)) * x^(f.β-1)
    return dh
end
function cumulative_hazard_rate_param_derivative(f::WeibullFamilyModel, x::Float64,right::Bool)::Vector{Float64}
    dH = right ? f.comp.dHR : f.comp.dHL
    dH[1] = x<=0 ? 0 : f.α * log(x) * x^f.β
    return dH
end
function hazard_rate_derivative_param_derivative(f::WeibullFamilyModel, x::Float64)::Vector{Float64}
    f.comp.dhd[1] = x<=0 ? 0 : f.α * (2 * f.β - 1 + f.β * (f.β - 1) * log(x)) * x^(f.β - 2)
    return f.comp.dhd
end
function hazard_rate_2derivative(f::WeibullFamilyModel, x::Float64)
    return x<=0 ? 0 : f.α * f.β * (f.β - 1) * (f.β - 2) * x^(f.β - 3)
end
function hazard_rate_param_2derivative(f::WeibullFamilyModel, x::Float64)::Vector{Float64}
    f.comp.d2h[1] = x<=0 ? 0 : f.α * (2 + f.β * log(x)) * log(x) * x^(f.β - 1)
    return f.comp.d2h
end
function cumulative_hazard_rate_param_2derivative(f::WeibullFamilyModel, x::Float64, right::Bool)::Vector{Float64}
    d2H = right ? f.comp.d2HR : f.comp.d2HL
    d2H[1]= x<=0 ? 0 : f.α * log(x)^2 * x^f.β
    return d2H
end

const LDorder = 5
const bxLim = 0.000001
mutable struct LogLinearFamilyModel <:  FamilyModel
    α::Parameter #Float64
    β::Parameter #Float64
    comp::FamilyCompute
    priors::Priors
    LogLinearFamilyModel() = new()
end
function  LogLinearFamilyModel(α::Parameter, β::Parameter)
    fm =  LogLinearFamilyModel()
    fm.α, fm.β =  α, β
    fm.priors = [nothing, nothing] 
    fm.comp = FamilyCompute()
    init!(fm.comp, fm)
    return fm
end
params(fm::LogLinearFamilyModel)::Parameters = [fm.α,fm.β]
params!(fm::LogLinearFamilyModel, p::Parameters) = begin;fm.α,fm.β = p; nothing; end
nbparams(fm::LogLinearFamilyModel) = 2
function LogLinearFamilyModel(α::Prior, β::Prior)
    fm = LogLinearFamilyModel(0.0, 0.0)
    fm.priors = [α, β]
    return fm
end

hazard_rate(f::LogLinearFamilyModel, x::Float64)::Float64 = (x<=0 ? 0 : f.α * exp(f.β * x) )
inverse_hazard_rate(f::LogLinearFamilyModel, x::Float64)::Float64 = log(x/f.α)/f.β
function cumulative_hazard_rate(f::LogLinearFamilyModel, x::Float64)::Float64
    res = 0.0
    if abs(f.β*x) < bxLim 
      prec = f.β * x/2
      res = 1 + prec
      for i in 1:(LDorder - 1)
        prec = prec * f.β*x / (i+1)
        res += prec
      end
      res = f.α * res * x
    else
      res=f.α * (exp(f.β * x) - 1)/f.β
    end
    return res
end
inverse_cumulative_hazard_rate(f::LogLinearFamilyModel, x::Float64)::Float64 = log(1 + x * f.β / f.α) / f.β
hazard_rate_derivative(f::LogLinearFamilyModel, x::Float64)::Float64 = (x<0 ? 0 : f.α * f.β * exp(f.β * x))

function hazard_rate_param_derivative(f::LogLinearFamilyModel, x::Float64,right::Bool)::Vector{Float64}
    dh = right ? f.comp.dhR : f.comp.dhL
    dh[1]= f.α * x * exp(f.β * x)
    return dh
end
function cumulative_hazard_rate_param_derivative(f::LogLinearFamilyModel, x::Float64,right::Bool)::Vector{Float64}
    prec = 0.0
    dH = right ? f.comp.dHR : f.comp.dHL
    if abs(f.β * x) < bxLim
        prec = f.β * x / 6
        res = 0.5 + 2 * prec
        for i in 1:(LDorder - 1)
          prec *=  f.β * x / (i+3)
          res += (i + 2) * prec
        end
        res *= f.α * x^2
    else 
        res = f.α * (x * exp(x * f.β) / f.β - (exp(f.β * x) - 1) / f.β^2)
    end
  
    dH[1] = res
    return dH
end
function hazard_rate_derivative_param_derivative(f::LogLinearFamilyModel, x::Float64)::Vector{Float64}
    f.comp.dhd[1]= f.α * exp(f.β * x) * (1 + f.β*x)
    return f.comp.dhd
end
function hazard_rate_2derivative(f::LogLinearFamilyModel, x::Float64)::Float64
    return x<=0 ? 0 : f.α * f.β^2 * exp(f.β*x)
end
function hazard_rate_param_2derivative(f::LogLinearFamilyModel, x::Float64)::Vector{Float64}
    f.comp.d2h[1] = f.α * x^2 * exp(f.β*x)
    return f.comp.d2h
end
function cumulative_hazard_rate_param_2derivative(f::LogLinearFamilyModel, x::Float64, right::Bool)::Vector{Float64}
    prec = 0.0
    d2H = right ? f.comp.d2HR : f.comp.d2HL
    if abs(f.β * x) < bxLim
        prec = f.β * x /24
        res = (2.0/3) + 6 * prec
        for i=1:(LDorder - 1)
          prec *= f.β * x / (i + 4)
          res += (i + 2) * (i + 3) * prec
        end
        res *= f.α * x^3
    else 
        res = f.α * (x^2 * exp(x * f.β) / f.β - 2 * x * exp(x * f.β) / f.β^2 + 2 * (exp(f.β*x)-1)/f.β^3)
    end
  
    d2H[1]= res
    return d2H
end

mutable struct Weibull3FamilyModel <: FamilyModel
    α::Parameter #Float64
    β::Parameter #Float64
    δ::Parameter #Float64
    comp::FamilyCompute
    priors::Priors
    Weibull3FamilyModel() = new()
end
function Weibull3FamilyModel(α::Parameter, β::Parameter, δ::Parameter)
    fm = Weibull3FamilyModel()
    fm.α, fm.β, fm.δ =  α, β, δ
    fm.priors = [nothing, nothing, nothing]
    fm.comp = FamilyCompute()
    init!(fm.comp, fm)
    return fm
end
params(fm::Weibull3FamilyModel)::Parameters = [fm.α,fm.β,fm.δ]
params!(fm::Weibull3FamilyModel, p::Parameters) = begin; fm.α,fm.β,fm.δ = p; nothing; end
nbparams(fm::Weibull3FamilyModel) = 3
function Weibull3FamilyModel(α::Prior, β::Prior, δ::Prior)
    fm = Weibull3FamilyModel(0.0, 0.0, 0.0)
    fm.priors = [α, β, δ]
    return fm
end

hazard_rate(f::Weibull3FamilyModel, x::Float64)::Float64 = x<0 ? 0 : f.α * f.β * (x + f.δ)^(f.β - 1)
inverse_hazard_rate(f::Weibull3FamilyModel, x::Float64)::Float64 = x<=0 ? 0 : (x/f.α/f.β)^(1/(f.β-1)) - f.δ
cumulative_hazard_rate(f::Weibull3FamilyModel, x::Float64)::Float64 = x<=0 ? 0 : f.α * ((x + f.δ)^f.β - f.δ^f.β)
inverse_cumulative_hazard_rate(f::Weibull3FamilyModel, x::Float64)::Float64 = x<=0 ? 0 : (f.δ^f.β + x/f.α)^(1/f.β) - f.δ
hazard_rate_derivative(f::Weibull3FamilyModel, x::Float64)::Float64 = x<=0 ? 0 : f.α * f.β * (f.β - 1) * (x + f.δ)^(f.β - 2) 

function hazard_rate_param_derivative(f::Weibull3FamilyModel, x::Float64,right::Bool)::Vector{Float64}
    dh = right ? f.comp.dhR : f.comp.dhL
    dh[1] = x==0 ? 0 : f.α * (1 + f.β * log(x + f.δ)) * (x + f.δ)^(f.β - 1)
    dh[2] = x<=0 ? 0 : f.α * f.β * (f.β - 1) * (x+f.δ)^(f.β - 2)
    return dh
end
function cumulative_hazard_rate_param_derivative(f::Weibull3FamilyModel, x::Float64,right::Bool)::Vector{Float64}
    dH = right ? f.comp.dHR : f.comp.dHL
    dH[1] = x==0 ? 0 : f.α * (log(x + f.δ) * (x+f.δ)^f.β - log(f.δ) * f.δ^f.β)
    dH[2] = x<=0 ? 0 : f.α * f.β * ((x + f.δ)^(f.β - 1) - f.δ^(f.β - 1))
    return dH
end
function hazard_rate_derivative_param_derivative(f::Weibull3FamilyModel, x::Float64)::Vector{Float64}
    f.comp.dhd[1] = x==0 ? 0 : f.α * (2 * f.β - 1 + f.β * (f.β - 1) * log(x + f.δ)) * (x+f.δ)^(f.β-2)
    f.comp.dhd[2] = x<=0 ? 0 : f.α * f.β * (f.β - 1) * (f.β - 2) * (x + f.δ)^(f.β - 3)
    return f.comp.dhd
end
function hazard_rate_2derivative(f::Weibull3FamilyModel, x::Float64)::Float64
    return x <= 0 ? 0 : f.α * f.β * (f.β - 1) * (f.β - 2) * (x + f.δ)^(f.β - 3)
end
function hazard_rate_param_2derivative(f::Weibull3FamilyModel, x::Float64)::Vector{Float64}
    f.comp.d2h[1] = x==0 ? 0 : f.α * (2 + f.β * log(x + f.δ)) * log(x + f.δ) * (x + f.δ)^(f.β - 1)
    f.comp.d2h[2] = x==0 ? 0 : f.α * (x + f.δ)^(f.β - 2) * (2 * f.β - 1 + f.β * (f.β - 1) * log(x + f.δ))
    f.comp.d2h[3] = x<=0 ? 0 : f.α * f.β * (f.β - 1) * (f.β - 2) * (x + f.δ)^(f.β - 3)
    return f.comp.d2h;
end
function cumulative_hazard_rate_param_2derivative(f::Weibull3FamilyModel, x::Float64, right::Bool)::Vector{Float64}
    d2H = right ? f.comp.d2HR : f.comp.d2HL
    d2H[1] = x==0 ? 0 : f.α * (log(x + f.δ)^2 * (x + f.δ)^f.β - log(f.δ)^2 * (f.δ^f.β))
    d2H[2] = x==0 ? 0 : f.α * ((x + f.δ)^(f.β - 1) * (f.β * log(x + f.δ) + 1) - f.δ^(f.β-1) * (f.β * log(f.δ) + 1))
    d2H[3] = x<=0 ? 0 : f.α * f.β * (f.β-1) * ((x + f.δ)^(f.β - 2) - f.δ^(f.β-2))
    return d2H
end

