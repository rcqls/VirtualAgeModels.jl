@enum ProposalMode normal=1 lognormal=2

mutable struct Bayesian
    mle::MLE
    model::Model
    estimates::DataFrame
    proposal::Function # see proposal!, rand and pdf methods below
    σ::Float64
    mode::ProposalMode
    Bayesian() = new()
end

function bayesian(model::Model, data::DataFrame; mode::ProposalMode = normal, method = Newton())
    if isbayesian(model)
        b = Bayesian()
        b.model = model
        b.mode = mode
        b.σ = 0.1
        b.mle = MLE(model, data)
        return b
    else
        return nothing
    end
end

function proposal!(b::Bayesian, mode::ProposalMode, σ::Float64)
    b.mode = mode
    b.σ = σ
    if b.mode == normal
        b.proposal = μ::Parameter -> Normal(μ,b.σ)
    elseif b.mode == lognormal
        tmp = 1+b.σ^2/μ^2
        b.proposal = μ::Parameter -> LogNormal(log(μ/sqrt(tmp)),sqrt(log(tmp)))
    end 
end

import Base.rand
rand(b::Bayesian, μ::Parameter) = rand(b.proposal(μ))
import Distributions.pdf
pdf(b::Bayesian, μ::Parameter, x::Float64) = pdf(b.proposal(μ),x)

function mcmc(b::Bayesian, θ::Parameters; nb::Int=10000, burn::Int=1000, profile::Bool=true, mode = ProposalMode = normal, σ=0.1)
    nbparams = length(θ)
    proposal!(b, mode, σ)
    priors_ = priors(b.model)
    curθ, oldθ =copy(θ), copy(θ)
    oldL = contrast(b.mle,oldθ,profile=profile)
    curL = oldL
    ind, θhat, αhat = Int[], Parameter[], Parameter[]
    for i in 1:nb
        for j in 1+profile:nbparams
            curθ[j]=rand(b,oldθ[j])
            curL=contrast(b.mle, curθ, profile=profile)
            r = exp(curL-oldL) * pdf(priors_[j], curθ[j])/ pdf(priors_[j], oldθ[j])
            if b.mode == lognormal
                r *= pdf(b, curθ[j], oldθ[j]) / pdf(b, oldθ[j], curθ[j])
            end
            if r > rand()
                if i >= burn
                    push!(ind, j)
                    push!(θhat,curθ[j])
                    if profile
                        push!(αhat, αEst(b.mle,curθ))
                    end
                end
                oldθ[j] = curθ[j]
                oldL = curL
            else
                curθ[j] = oldθ[j]
            end
        end
    end
    res = DataFrame(ind=ind, θ=θhat)
    if profile
        res[!,:α] = αhat
    end
    return res
end

