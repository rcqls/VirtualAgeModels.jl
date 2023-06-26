mutable struct Compute
    S1::Float64
    S2::Float64
    S0::Float64
    S3::Float64
    S4::Float64

    dS1::Vector{Float64}
    dS2::Vector{Float64}
    dS3::Vector{Float64}
    dS4::Vector{Float64} 

    d2S1::Vector{Float64}
    d2S2::Vector{Float64}
    d2S3::Vector{Float64}
    # dims
    nbm::Int #maintenance
    nb2m::Int #maintenance deriv order 2
    nbc::Int #covariate
    nbd::Int #params for deriv
    nb2d::Int #deriv order 2
    Compute() = new()
end

function Compute(m::AbstractModel)
    comp = Compute()
    init_dims!(comp,m)
    init!(comp,deriv=true)
    return comp
end

function init_dims!(c::Compute, m::AbstractModel)
    c.nbm = m.nb_params_maintenance
    c.nbd = m.nb_params_maintenance + m.nb_params_family - 1 + m.nb_params_cov
    c.nb2m = m.nb_params_maintenance * (m.nb_params_maintenance + 1) ÷ 2
    c.nb2d = c.nbd * (c.nbd + 1) ÷ 2
    c.nbc = m.nb_params_cov
end

function init!(c::Compute; deriv::Bool = false)
    c.S0, c.S1, c.S2, c.S3, c.S4 = zeros(5)
    if deriv
        c.dS1, c.dS2, c.dS3=zeros(c.nbd + c.nbc), zeros(c.nbd), zeros(c.nbm)
        if c.nbc > 0
            c.dS4 = zeros(c.nbc)
        end
        c.d2S1, c.d2S2 = zeros(c.nb2d), zeros(c.nb2d)  #inferior part of the hessian matrice by lines
        c.d2S3 = zeros(c.nb2m)
    end
end