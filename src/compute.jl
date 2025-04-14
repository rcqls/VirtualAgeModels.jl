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
    #d2S4::Vector{Float64}#LD: second order derivatives of S4 are null !
    # dims
    nbm::Int #number of maintenance params
    nb2m::Int #maintenance params for deriv order 2
    nbc::Int #number of covariate params
    nb2c::Int#covariate params deriv for deriv order 2 #LD
    nbd::Int #total number of params
    nb2d::Int #deriv order 2
    nbd_ssm::Int #total number of params except maintenance ones params#LD
    nb2d_ssm::Int #deriv order 2
    Compute() = new()
end

function Compute(m::AbstractVirtualAgeModel)
    comp = Compute()
    init_dims!(comp,m)
    init!(comp,deriv=true)
    return comp
end

function init_dims!(c::Compute, m::AbstractVirtualAgeModel)
    c.nbm = m.nb_params_maintenance
    c.nbd = m.nb_params_maintenance + m.nb_params_family - 1 + m.nb_params_cov
    c.nbc = m.nb_params_cov
    c.nbd_ssm = m.nb_params_maintenance + m.nb_params_family - 1
    c.nb2m = ind_nb(c.nbm)
    c.nb2d = ind_nb(c.nbd)
    c.nb2c = ind_nb(c.nbc)
    c.nb2d_ssm = ind_nb(c.nbd_ssm)
end

function init!(c::Compute; deriv::Bool = false)
    c.S0, c.S1, c.S2, c.S3, c.S4 = zeros(5)
    if deriv
        c.dS1, c.dS2, c.dS3=zeros(c.nbd), zeros(c.nbd_ssm), zeros(c.nbm)
        c.d2S1, c.d2S2, c.d2S3=zeros(c.nb2d), zeros(c.nb2d_ssm), zeros(c.nb2m)
        if c.nbc > 0
            c.dS4 = zeros(c.nbc)
            #c.d2S4 = zeros(c.nb2c)
        end
    end
end