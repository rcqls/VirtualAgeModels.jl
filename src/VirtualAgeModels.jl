module VirtualAgeModels

using Random, DataFrames, Optim, Distributions, Plots
import Base.first
import Plots.plot


export Simulator
export parse_model, @vam, @stop, params, params!, nbparams, data, data!, mle, bayesian
export contrast, gradient, hessian
export plot


abstract type AbstractModel end

const Prior = Union{Nothing,Distribution}
const Priors = Vector{Prior}
const Parameter = Float64
const Parameters = Vector{Float64}
# const debugprint = println
const debugprint(x) = ()

include("tool.jl")
include("compute.jl")
include("formula.jl")
include("family.jl")
include("maintenance_model.jl")
include("maintenance_policy.jl")
include("model.jl")
include("plots.jl")
include("simulate.jl")
include("mle.jl")
include("bayesian.jl")

end
