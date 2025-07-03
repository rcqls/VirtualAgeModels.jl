function mle(model::VirtualAgeModel, θ::Vector{Float64},  data::DataFrame=DataFrame(), datacov::DataFrame=DataFrame(); profile::Bool = true, fixed::Union{Vector{Int},Vector{Bool}} = Bool[], method = Newton())
    mle = MLE(model, data, datacov)
    # TODO: check boundary for fixed
    if fixed isa Vector{Bool}
        res = Int[]
        for (i,v) in enumerate(fixed)
            if v
                push!(res, i)
            end
        end
        fixed = res
    end
    # Apply profile likelihood only when α (at index 1) is not fixed
    profile = !(1 in fixed)
    unfixed = setdiff(1:length(θ),fixed)
    p = θ[unfixed]
    function f(θ′)
        θ[unfixed] = θ′
        -contrast(mle, θ, profile = profile)
    end
    function g!(storage, θ′)
        θ[unfixed] = θ′
        dlnL = gradient(mle, θ, profile=profile)
        storage .= -dlnL[unfixed]
    end
    res = nothing
    if method isa Optim.ZerothOrderOptimizer
        res = optimize(f, p, method=method)
    elseif method isa Optim.FirstOrderOptimizer
        res = optimize(f, g!, p, method=method)
    elseif method isa Optim.SecondOrderOptimizer
        function h!(storage, θ′)
            θ[unfixed] = θ′
            storage .= -hessian(mle, θ, profile=profile)[unfixed, unfixed]
        end
        res = optimize(f, g!, h!, p, method=method)
    end
    p = θ
    p[unfixed] = Optim.minimizer(res)
    mle.optim.result = res
    mle.optim.fixed = fixed
    mle.optim.profile = profile
    params!(mle, p)
    return mle
end

function mle(model::VirtualAgeModel, data::DataFrame=DataFrame(), datacov::DataFrame=DataFrame(); fixed::Union{Vector{Int},Vector{Bool}} = Bool[], method = Newton())
    θ = params(model)
    mle(model, θ, data, datacov; fixed=fixed, method=method)
end

function contrast(model::VirtualAgeModel, θ::Vector{Float64}, data::DataFrame=DataFrame(), datacov::DataFrame=DataFrame(); profile::Bool=true)::Float64
    m = MLE(model, data, datacov)
    return contrast(m, θ, profile = profile)
end

contrast(model::VirtualAgeModel, data::DataFrame=DataFrame(), datacov::DataFrame=DataFrame(); profile::Bool=true) = contrast(model, params(model), data; profile = profile, datacov = datacov)

function gradient(model::VirtualAgeModel, θ::Vector{Float64}, data::DataFrame=DataFrame(), datacov::DataFrame=DataFrame(); profile::Bool=true)::Vector{Float64}
    m = MLE(model, data, datacov)
    return gradient(m, θ, profile = profile)
end

gradient(model::VirtualAgeModel, data::DataFrame=DataFrame(), datacov::DataFrame=DataFrame(); profile::Bool=true) = gradient(model, params(model), data; profile = profile, datacov = datacov)

function hessian(model::VirtualAgeModel, θ::Vector{Float64}, data::DataFrame=DataFrame(), datacov::DataFrame=DataFrame(); profile::Bool=true)::Matrix{Float64}
    m = MLE(model, data, datacov)
    return hessian(m, θ, profile = profile)
end

hessian(model::VirtualAgeModel, data::DataFrame=DataFrame(), datacov::DataFrame=DataFrame(); profile::Bool=true) = hessian(model, params(model), data; profile = profile, datacov = datacov)

mutable struct MLEOptim
    profile::Bool # !(1 in fixed)
    fixed::Vector{Int}
    result
    MLEOptim() = new()
end

mutable struct MLE
    model::VirtualAgeModel
    left_censors::Vector{Int} #CAREFUL: this is a vector of indices!
    left_censor::Int #left_censor for current system
    optim::MLEOptim
    comp::Compute
    MLE() = new()
end

function MLE(model::VirtualAgeModel)::MLE
    mle = MLE()
    mle.model = model
    init!(mle.model)
    mle.comp = Compute(mle.model)
    mle.optim = MLEOptim()
    left_censors!(mle, Int[])
    return mle
end

function MLE(model::VirtualAgeModel, data::DataFrame, datacov::DataFrame)::MLE
    mle = MLE()
    mle.model = model
    init!(mle.model)

    if !isempty(data)
        data!(mle.model, data, datacov)
    end
    mle.comp = Compute(mle.model)
    mle.optim = MLEOptim()
    left_censors!(mle, Int[])
    return mle
end

function show(io::IO, mle::MLE)
    println(io, mle.model)
    print(io, mle.optim)
end

function params(mle::MLE)::Vector{Float64}
    if mle.optim.profile
        [αEst(mle, params(mle.model)); params(mle.model)[2:end]]
    else
        params(mle.model)
    end
end
params!(m::MLE, θ::Vector{Float64}) = params!(m.model, θ)

## TODO: deal with left_censors
function left_censors!(m::MLE, left_censors::Vector{Int})
    m.left_censors = left_censors
    m.left_censor = 0
end

function select_left_censor(mle::MLE, i::Int)
    if length(mle.left_censors) >= i 
        mle.left_censor=left_censors[i]
    end
end

# Rcpp -> init_mle_vam_for_current_system
# deriv means gradient and hessian
function init_mle(mle::MLE; deriv::Bool=false)
    for mm in mle.model.models
        init!(mm)
    end

    mle.model.Vright = 0
    mle.model.Vright = 0
    mle.model.A = 1
    mle.model.k = 1
    mle.model.id_mod = 0 #id of current model
    mle.model.id_params = 1
    init!(mle.model.comp, deriv=deriv)

    for type in mle.model.type
        if type < 0 
            mle.model.comp.S0 += 1
        end
    end
    if deriv
        # mle.model.dVright = zeros(mle.model.comp.nbd)
        # mle.model.dA = zeros(mle.model.comp.nbd)
        # nb2m = mle.model.nb_params_maintenance * (mle.model.nb_params_maintenance + 1) ÷ 2
        # mle.model.d2Vright = zeros(nb2m)
        # mle.model.d2A = zeros(nb2m)
        fill!(mle.model.dVright, 0.0)
        fill!(mle.model.dA, 0.0)
        fill!(mle.model.d2Vright, 0.0)
        fill!(mle.model.d2A, 0.0)
    end
end

function contrast(mle::MLE, θ::Vector{Float64}; profile::Bool=true)::Float64
    res = 0
    α = θ[1] #;//save current value of alpha

    θ[1] = 1 #//Rmk: alpha replaces θ[0] => a bit weird!

    init!(mle.comp)

    params!(mle.model,θ);
    # //printf("System %d\n",1);
    data!(mle.model, 1)
    mle.model.nb_params_cov > 0 && select_current_system(mle.model, 1, true)
    select_left_censor(mle, 1)
    debugprint("Contrast System 1")
    contrast_current(mle)
    # //only if multi-system
    if mle.model.nb_system > 1
        for i in 2:mle.model.nb_system 
            #
            debugprint("Contrast System $i");
            data!(mle.model, i)
            mle.model.nb_params_cov > 0 && select_current_system(mle.model, i, true)
            select_left_censor(mle, i)
            contrast_current(mle)
        end
    end

    # //DEBUG: println("alpha=$alpha,S0=$(mle.comp.S0),S1=$(mle.comp.S1),S2=$(mle.comp.S2),S3=$(mle.comp.S3),S4=$(mle.comp.S4)")
    # //printf("params=(%lf,%lf)\n",model.params_cov[0],model.params_cov[1]);
    # // log-likelihood (with constant +S0*(log(S0)-1))
    if profile
        θ[1] = mle.comp.S0 / mle.comp.S1
        res = -log(mle.comp.S1) * mle.comp.S0 + mle.comp.S2 +  mle.comp.S0 * (log( mle.comp.S0) - 1) +  mle.comp.S3;
        params!(mle.model, θ) #//also memorize the current value for alpha which is not 1 in fact
    else
        θ[1] = α
        res = log(α) * mle.comp.S0 + mle.comp.S2 - α * mle.comp.S1 + mle.comp.S3
        params!(mle.model, θ) #//also memorize the current value for alpha which is not 1 in fact
    end
    mle.model.nb_params_cov > 0 && (res += mle.comp.S4)

    θ[1] = α #//LD:changed for bayesian
    return res
end

contrast(mle::MLE; profile::Bool=true) = contrast(mle, params(mle); profile=profile)

function contrast_current(mle::MLE)
    init_mle(mle)
    n = length(mle.model.time)
    while mle.model.k < n
        # //printf("  Time=%f, Type=%d\n",model.time[model.k+1],model.type[model.k+1]);
        contrast_update_current(mle)
        #// previous model for the next step
        type = mle.model.type[mle.model.k + 1]
        if type < 0 
            type = 0
        end
        #//model.indMode = (type < 0 ? 0 : type);
        update_maintenance!(mle.model, type)
    end
    contrast_update_S(mle)
end

function contrast_update_current(mle::MLE; gradient::Bool=false, hessian::Bool=false)
    update_Vleft!(mle.model, gradient=gradient, hessian=hessian)
    mle.model.hVleft = hazard_rate(mle.model.family, mle.model.Vleft)
    debugprint("jl ($gradient, $hessian): hVleft=$(mle.model.hVleft) $(params(mle))")
    mle.model.indType = (mle.model.type[mle.model.k + 1] < 0 ? 1.0 : 0.0)
    if mle.model.k >= mle.left_censor 
        mle.model.comp.S1 += cumulative_hazard_rate(mle.model.family, mle.model.Vleft) - cumulative_hazard_rate(mle.model.family, mle.model.Vright)
    end
    mle.model.comp.S2 += log(mle.model.hVleft) * mle.model.indType
    mle.model.comp.S3 += log(mle.model.A) * mle.model.indType
end

function contrast_update_S(mle::MLE)
    #//model updated for current system: mle.comp.S1,S2,S0,dS1,dS2
    tmp = mle.model.comp.S1
    mle.comp.S2 += mle.model.comp.S2
    mle.comp.S0 += mle.model.comp.S0
    mle.comp.S3 += mle.model.comp.S3
    if mle.model.nb_params_cov > 0
        compute_covariates(mle.model) #//initialize model.sum_cov
        tmp *= exp(mle.model.sum_cov)
        mle.comp.S4 += mle.model.comp.S0 * mle.model.sum_cov
        #//printf("(S0=%lf) * (sum_cov=%lf) = (S4 =%lf)\n",model.comp.S0, mle.model.sum_cov,model.comp.S0 * mle.model.sum_cov);
    end
    mle.comp.S1 += tmp
    #//printf("Conclusion : mle.comp.S1=%f, S2=%f, S0=%f, S4=%f\n",model.comp.S1,model.comp.S2,model.comp.S0,model.comp.S4);
end

function gradient(mle::MLE, θ::Vector{Float64}; profile::Bool=true)::Vector{Float64}
    res = zeros(mle.model.nb_params_family + mle.model.nb_params_maintenance + mle.model.nb_params_cov)
    α=θ[1] #save current value of alpha

    θ[1]=1
    init!(mle.comp, deriv = true)
    params!(mle.model, θ)
    data!(mle.model, 1)
    debugprint("gradient System 1")
    mle.model.nb_params_cov > 0 && select_current_system(mle.model, 1, true)
    select_left_censor(mle, 1)
    gradient_current(mle)
    # //only if multi-system
    if mle.model.nb_system > 1
        for i in 2:mle.model.nb_system 
            data!(mle.model, i)
            debugprint("gradient System $i")
            mle.model.nb_params_cov > 0 && select_current_system(mle.model, i, true)
            select_left_censor(mle, i)
            gradient_current(mle)
        end
    end

        #verif dSommes
        if false
            print("jl: dS1:")
            for i in 1:(mle.model.nb_params_family + mle.model.nb_params_maintenance + mle.model.nb_params_cov -1)
                print("[$i]=$(mle.comp.dS1[i]), ")
            end
            print("\n")
            print("jl: dS2:")
            for i in 1:(mle.model.nb_params_family + mle.model.nb_params_maintenance -1)
                print("[$i]=$(mle.comp.dS2[i]), ")
            end
            print("\n")
            print("jl: d23:")
            for i in 1:(mle.model.nb_params_maintenance)
                print("[$i]=$(mle.comp.dS3[i]), ")
            end
            print("\n")
            print("nb param: family=$(mle.model.nb_params_family), maintenance=$(mle.model.nb_params_maintenance), cov=$(mle.model.nb_params_cov)\n")
            end


    # compute gradient
    θ[1] = profile ? mle.comp.S0 / mle.comp.S1 : α

    params!(mle.model, θ) # also memorize the current value for alpha which is not 1 in fact

    res[1] = profile ? 0.0 : mle.comp.S0/α - mle.comp.S1
    
    np = 1
    for i in 1:(mle.model.nb_params_family - 1)
        res[i + np] = -mle.comp.dS1[i] * θ[1] + mle.comp.dS2[i]
    end
    np += mle.model.nb_params_family - 1
    if mle.model.nb_params_maintenance > 0
        for i in 1:mle.model.nb_params_maintenance
            res[i + np] = -mle.comp.dS1[i + np - 1] * θ[1] + mle.comp.dS2[i + np - 1] + mle.comp.dS3[i]
        end
    end
    np += mle.model.nb_params_maintenance
    if mle.model.nb_params_cov > 0
        for i in 1:mle.model.nb_params_cov
            #LD
            res[i + np] = -mle.comp.dS1[i + np - 1] * θ[1] + mle.comp.dS4[i]
            #res[i + np] = -mle.comp.dS1[i + np] * θ[1] + mle.comp.dS4[i]
            # println("jl: res[$(i + np)] = $(-mle.comp.dS1[i + np]) * $(θ[1]) + $(mle.comp.dS4[i])")
        end
    end
    θ[1] = α ## BIZARRE!
    return res
end

gradient(mle::MLE; profile::Bool=true) = gradient(mle, params(mle); profile=profile)

function gradient_current(mle::MLE)
    init_mle(mle, deriv = true)
    n = length(mle.model.time)
    while mle.model.k < n
        gradient_update_current(mle)
        type = mle.model.type[mle.model.k + 1]
        if type < 0 
            type = 0
        end
        # //model.indMode = (type < 0 ? 0 : type)
        update_maintenance!(mle.model, type, gradient = true)
    end
    contrast_update_S(mle)
    #//precomputation of covariate term to multiply (in fact just exp)
    for i = 1:(mle.model.nb_params_family - 1)
        gradient_update_dS_family(mle, i)
    end
    np = mle.model.nb_params_family - 1
    for i in 1:mle.model.nb_params_maintenance
        gradient_update_dS_maintenance(mle, i + np,i)
    end
    if mle.model.nb_params_cov >= 1 #LD: > 1
        np += mle.model.nb_params_maintenance #LD: mle.model.nb_params_cov
        for i in 1:mle.model.nb_params_cov
            gradient_update_dS_covariate(mle, i + np, i)
        end
    end
end

function gradient_update_current(mle::MLE)
    contrast_update_current(mle, gradient =true)

    cumhVright_param_derivative = cumulative_hazard_rate_param_derivative(mle.model.family, mle.model.Vright, true)
    cumhVleft_param_derivative=cumulative_hazard_rate_param_derivative(mle.model.family, mle.model.Vleft, false)
    hVleft_param_derivative=hazard_rate_param_derivative(mle.model.family, mle.model.Vleft, false)
    for i in 1:(mle.model.nb_params_family - 1)
        if mle.model.k >= mle.left_censor 
            mle.model.comp.dS1[i] +=  cumhVleft_param_derivative[i]-cumhVright_param_derivative[i]
        end
        mle.model.comp.dS2[i] += hVleft_param_derivative[i] / mle.model.hVleft * mle.model.indType
    end
    hVright=hazard_rate(mle.model.family, mle.model.Vright)
    dhVleft=hazard_rate_derivative(mle.model.family, mle.model.Vleft)
    # printf("k:%d,hVright:%lf,dhVleft:%lf,indType:%lf\n",model.k,hVright,dhVleft,model.indType);
    np = mle.model.nb_params_family - 1
    for i in 1:mle.model.nb_params_maintenance
        if mle.model.k >= mle.left_censor 
            mle.model.comp.dS1[i + np] += mle.model.hVleft * mle.model.dVleft[i] - hVright * mle.model.dVright[i]
        end
        # printf("dS1[%d]=(%lf,%lf,%lf),%lf,",ii+1,model.hVleft,model.dVleft[ii],model.dVright[ii],model.dS1[ii+1]);
        mle.model.comp.dS2[i + np] +=  dhVleft * mle.model.dVleft[i] / mle.model.hVleft * mle.model.indType
        #//printf("dS2[%d]=%lf,",ii+1,model.dS2[ii+1]);
        mle.model.comp.dS3[i] +=  mle.model.dA[i] / mle.model.A * mle.model.indType
    end
    #//printf("\n");
end

function gradient_update_dS_maintenance(mle::MLE,  i::Int, ii::Int) 
    mle.comp.dS1[i] += mle.model.comp.dS1[i] * (mle.model.nb_params_cov > 0 ? exp(mle.model.sum_cov) : 1.0)
    mle.comp.dS2[i] += mle.model.comp.dS2[i]
    mle.comp.dS3[ii] += mle.model.comp.dS3[ii]
end

function gradient_update_dS_family(mle::MLE, i::Int)
    mle.comp.dS1[i] += mle.model.comp.dS1[i] * (mle.model.nb_params_cov > 0 ? exp(mle.model.sum_cov) : 1.0)
    mle.comp.dS2[i] += mle.model.comp.dS2[i]
end

function gradient_update_dS_covariate(mle::MLE, i::Int, ii::Int)
    #//nb_params_cov > 0 necessarily
    cov=covariate(mle.model, ii)
    mle.comp.dS1[i] += mle.model.comp.S1 * cov * exp(mle.model.sum_cov)
    #//dS2[i]=0
    mle.comp.dS4[ii] += mle.model.comp.S0 * cov
    # println("jl: dS4[$ii]=$(mle.comp.dS4[ii])")
end

function hessian(mle::MLE, θ::Vector{Float64}; profile::Bool=true)::Matrix{Float64}
    res = zeros(mle.model.nb_params_family + mle.model.nb_params_maintenance + mle.model.nb_params_cov, mle.model.nb_params_family + mle.model.nb_params_maintenance + mle.model.nb_params_cov)
    α=θ[1] #save current value of alpha

    θ[1]=1
    init!(mle.comp, deriv = true)
    params!(mle.model, θ)
    data!(mle.model, 1)
    debugprint("hessian System 1")
    mle.model.nb_params_cov > 0 && select_current_system(mle.model, 1, true)
    select_left_censor(mle, 1)
    hessian_current(mle)

    # //only if multi-system
    if mle.model.nb_system > 1
        for i in 2:mle.model.nb_system 
            data!(mle.model, i)
            debugprint("hessian System $i")
            mle.model.nb_params_cov > 0 && select_current_system(mle.model, i, true)
            select_left_censor(mle, i)
            hessian_current(mle)
        end
    end

    #verif d2Sommes
    if false
    print("jl: d2S1:")
    for i in 1:(mle.model.nb_params_family + mle.model.nb_params_maintenance + mle.model.nb_params_cov -1)
        for j in 1:i
            ij = ind_ij(i, j)
            print("[$ij]=$(mle.comp.d2S1[ij]), ")
        end
    end
    print("\n")
    print("jl: d2S2:")
    for i in 1:(mle.model.nb_params_family + mle.model.nb_params_maintenance -1)
        for j in 1:i
            ij = ind_ij(i, j)
            print("[$ij]=$(mle.comp.d2S2[ij]), ")
        end
    end
    print("\n")
    print("jl: d2S3:")
    for i in 1:(mle.model.nb_params_maintenance)
        for j in 1:i
            ij = ind_ij(i, j)
            print("[$ij]=$(mle.comp.d2S3[ij]), ")
        end
    end
    print("\n")
    print("nb param: family=$(mle.model.nb_params_family), maintenance=$(mle.model.nb_params_maintenance), cov=$(mle.model.nb_params_cov)\n")
    end
    
    # LD: compute hessian
    θ[1] = profile ? mle.comp.S0 / mle.comp.S1 : α

    params!(mle.model, θ) # also memorize the current value for alpha which is not 1 in fact

    if !profile
        res[1, 1] = - mle.comp.S0 / α^2
        for i in 1:(mle.model.nb_params_family + mle.model.nb_params_maintenance + mle.model.nb_params_cov -1)
            res[1,i + 1] = - mle.comp.dS1[i]
            res[i + 1,1] = res[1,i + 1]
        end
    else
        for i in 1:(mle.model.nb_params_family + mle.model.nb_params_maintenance + mle.model.nb_params_cov - 1)
            for j in 1:i
                res[i + 1,j + 1] = mle.comp.S0 / mle.comp.S1^2 * mle.comp.dS1[i] * mle.comp.dS1[j]
                res[j + 1,i + 1] = res[i + 1,j + 1]
            end
        end
    end

    for i in 1:(mle.model.nb_params_family - 1)
        for j in 1:i
            ij = ind_ij(i, j)
            res[i + 1,j + 1] += mle.comp.d2S2[ij] - θ[1] * mle.comp.d2S1[ij]
            res[j+ 1,i + 1] = res[i + 1,j + 1]
        end
    end
    np = mle.model.nb_params_family
    if mle.model.nb_params_maintenance > 0
        for i in 1:mle.model.nb_params_maintenance
            for j in 1:(mle.model.nb_params_family - 1)
                ij = ind_ij(i + np - 1, j)
                res[i + np,j + 1] += mle.comp.d2S2[ij] - θ[1] * mle.comp.d2S1[ij]
                res[j + 1,i + np] = res[i + np,j + 1]
            end
            for j in 1:i #Warning: j<=i
                ij = ind_ij(i + np - 1,j + np - 1)
                ij_2 = ind_ij(i ,j)
                res[i + np,j + np] += mle.comp.d2S2[ij] - θ[1] * mle.comp.d2S1[ij] + mle.comp.d2S3[ij_2]
                res[j + np,i + np] = res[i + np,j + np]
            end
        end
    end
    np += mle.model.nb_params_maintenance
    if mle.model.nb_params_cov > 0
        #LD: second order derivatives for S4 are null !
        #consequently only the ssecond order of S1 contributes to the second order derivatives of the likelihood
        for i in (np + 1):(np + mle.model.nb_params_cov)
            for j in 1:(i - 1) #Warning: j<=i and S1 does not depend on alpha
                ij = ind_ij(i - 1, j)
                res[i ,j + 1] += - θ[1] * mle.comp.d2S1[ij]
                res[j+1 ,i] = res[i ,j + 1]
            end
        end
        #for i in 1:mle.model.nb_params_cov
        #    for j in 1:(mle.model.nb_params_family + mle.model.nb_params_maintenance - 1)
        #        ij = ind_ij(i + np - 1, j)
        #        res[i + np,j + 1] += - θ[1] * mle.comp.d2S1[ij]
        #        res[j + 1,i + np] = res[i + np,j + 1]
        #    end
        #    for j in 1:i #Warning: j<=i
        #        ij = ind_ij(i + np - 1,j + np - 1)
        #        ij_2 = ind_ij(i ,j)
        #        res[i + np,j + np] += - θ[1] * mle.comp.d2S1[ij] + mle.comp.d2S4[ij_2]
        #        res[j + np,i + np] = res[i + np,j + np]
        #    end
        #end
    end
# End LD

#    # println("dS1=$(mle.comp.dS1) dS2=$(mle.comp.dS2) ")
#    # println("d2S1=$(mle.comp.d2S1) d2S2=$(mle.comp.d2S3) d2S2=$(mle.comp.d2S3) ")
#    # //compute hessian
#    if profile
#        res[1,1] = 0
#        for i in 1:(mle.model.nb_params_family - 1)
#            ii = ind_ij(i, i)
#            res[1, i + 1] = 0
#            res[i + 1, 1] = 0
#            res[i + 1, i + 1] = mle.comp.dS1[i]^2 / mle.comp.S1^2 * mle.comp.S0 - mle.comp.d2S1[ii]/mle.comp.S1 * mle.comp.S0 + mle.comp.d2S2[ii]
#            for j in 1:i # ?? or for j in 0:(i - 1)
#                ij = ind_ij(i, j)
#                #//i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
#                res[i + 1, j + 1] = mle.comp.dS1[i] * mle.comp.dS1[j] / mle.comp.S1^2 * mle.comp.S0 - mle.comp.d2S1[ij] / mle.comp.S1 * mle.comp.S0 + mle.comp.d2S2[ij]
#                res[j + 1, i + 1] = res[i + 1,j + 1]
#            end
#        end
#        for i in mle.model.nb_params_family:(mle.model.nb_params_maintenance + mle.model.nb_params_family - 1)#avt LD : -2)
#            ii = ind_ij(i, i)
#            res[1, i + 1] = 0
#            res[i + 1, 1] = 0
#            #res[i + 1, i + 1] = mle.comp.dS1[i]^2 / mle.comp.S1^2 * mle.comp.S0 - mle.comp.d2S1[ii] / mle.comp.S1 * mle.comp.S0 + mle.comp.d2S2[ii] + mle.comp.d2S3[(i - 1 -(mle.model.nb_params_family - 1)) * (i -(mle.model.nb_params_family - 1)) ÷ 2 + i - (mle.model.nb_params_family-1)]
#            res[i + 1, i + 1] = mle.comp.dS1[i]^2 / mle.comp.S1^2 * mle.comp.S0 - mle.comp.d2S1[ii] / mle.comp.S1 * mle.comp.S0 + mle.comp.d2S2[ii] + mle.comp.d2S3[ind_ij(i-mle.model.nb_params_maintenance+1, i-mle.model.nb_params_maintenance+1)]#LD
#            for j in 1:(mle.model.nb_params_family-1)
#                ij = ind_ij(i, j)
#                #//i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
#                res[i + 1, j + 1] = mle.comp.dS1[i] * mle.comp.dS1[j] / mle.comp.S1^2 * mle.comp.S0 - mle.comp.d2S1[ij] / mle.comp.S1 * mle.comp.S0 + mle.comp.d2S2[ij]
#                res[j + 1, i + 1] = res[i + 1, j + 1]
#            end
#            for j in mle.model.nb_params_family:i
#                ij = ind_ij(i, j)
#                #//i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
#                res[i + 1, j + 1] = mle.comp.dS1[i] * mle.comp.dS1[j] / mle.comp.S1^2 * mle.comp.S0 - mle.comp.d2S1[ij] / mle.comp.S1 * mle.comp.S0 + mle.comp.d2S2[ij] + mle.comp.d2S3[(i - (mle.model.nb_params_family-1)) * (i - (mle.model.nb_params_family - 1) + 1) ÷ 2 + j - (mle.model.nb_params_family - 1)]
#                res[j + 1, i + 1] = res[i + 1, j + 1]
#            end
#        end
#        for i in (mle.model.nb_params_maintenance + mle.model.nb_params_family):(mle.model.nb_params_maintenance + mle.model.nb_params_family + mle.model.nb_params_cov - 1)
#            ii = ind_ij(i, i)
#            res[1, i + 1] = 0
#            res[i + 1, 1] = 0
#            res[i + 1, i + 1] = mle.comp.dS1[i] ^2 / mle.comp.S1^2 * mle.comp.S0 - mle.comp.d2S1[ii] / mle.comp.S1 * mle.comp.S0
#            for j in 1:i
#                ij = ind_ij(i, j)
#                res[i + 1,j + 1] = mle.comp.dS1[i] * mle.comp.dS1[j] / mle.comp.S1^2 * mle.comp.S0 - mle.comp.d2S1[ij] / mle.comp.S1 * mle.comp.S0
#                res[j + 1,i + 1] = res[i + 1,j + 1]
#            end
#        end
#        θ[1] = mle.comp.S0 / mle.comp.S1
#        params!(mle.model, θ) #;//also memorize the current value for alpha which is not 1 in fact
#    else
#
#        res[1, 1] = -mle.comp.S0 / α^2
#        for i in 1:(mle.model.nb_params_family-1)
#            ii = ind_ij(i, i)
#            res[1,i + 1] = -mle.comp.dS1[i]
#            res[i + 1,1] = -mle.comp.dS1[i]
#            res[i + 1,i + 1] = mle.comp.d2S2[ii] - α * mle.comp.d2S1[ii]
#            for j in 1:i
#                ij = ind_ij(i, j)
#                #//i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
#                res[i + 1,j + 1] = mle.comp.d2S2[ij] - α * mle.comp.d2S1[ij]
#                res[j + 1,i + 1] = res[i + 1,j + 1]
#            end
#        end
#        for i in mle.model.nb_params_family:(mle.model.nb_params_maintenance + mle.model.nb_params_family - 1)
#            ii = ind_ij(i, i)
#            res[1, i + 1] = -mle.comp.dS1[i]
#            res[i + 1, 1] = -mle.comp.dS1[i]
#            res[i + 1, i + 1] = mle.comp.d2S2[ii] - α * mle.comp.d2S1[ii] + mle.comp.d2S3[ind_ij(i-(mle.model.nb_params_family-1), i - (mle.model.nb_params_family - 1))]
#            for j in 1:(mle.model.nb_params_family - 1)
#                ij = ind_ij(i, j)
#                #//i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
#                res[i + 1, j + 1] = mle.comp.d2S2[ij] - α * mle.comp.d2S1[ij]
#                res[j + 1, i + 1] = res[i + 1, j + 1]
#            end
#            for j in mle.model.nb_params_family:(i - 1)
#                ij = ind_ij(i, j)
#                #//i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
#                res[i + 1, j + 1] = mle.comp.d2S2[ij] - α * mle.comp.d2S1[ij] + mle.comp.d2S3[ind_ij(i - (mle.model.nb_params_family - 1), j - (mle.model.nb_params_family - 1))]
#                res[j + 1, i + 1] = res[i + 1, j + 1]
#            end
#        end
#        for i in (mle.model.nb_params_maintenance+mle.model.nb_params_family):(mle.model.nb_params_maintenance+mle.model.nb_params_family + mle.model.nb_params_cov - 1)
#            ii = ind_ij(i, i)
#            res[1, i + 1] = -mle.comp.dS1[i]
#            res[i + 1, 1] = -mle.comp.dS1[i]
#            res[i + 1, i + 1] = -α * mle.comp.d2S1[ii]
#            for j in 1:i
#                ij = ind_ij(i, j)
#                res[i + 1, j + 1] = -α * mle.comp.d2S1[ij]
#                res[j + 1, i + 1] = res[i + 1 , j + 1]
#            end
#        end
#        θ[1] = α
#        params!(mle.model, θ) #also memorize the current value for alpha which is not 1 in fact
#
#    end
    θ[1] = α # LD:changed for bayesian
    return res
end

hessian(mle::MLE; profile::Bool=true) = hessian(mle, params(mle); profile=profile)

function hessian_current(mle::MLE)
    npf, npm, npc = mle.model.nb_params_family - 1, mle.model.nb_params_maintenance, mle.model.nb_params_cov
    init_mle(mle, deriv = true)
    n = length(mle.model.time)
    while mle.model.k < n
        hessian_update_current(mle)
        type = mle.model.type[mle.model.k + 1]
        if type < 0 
            type = 0
        end
        # //model.indMode = (type < 0 ? 0 : type)
        update_maintenance!(mle.model, type, hessian = true)
    end
    contrast_update_S(mle)
    for i = 1:npf
        gradient_update_dS_family(mle, i)
        for j in 1:i
            #//i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
            ij = ind_ij(i, j)
           mle.comp.d2S1[ij] += mle.model.comp.d2S1[ij] * (mle.model.nb_params_cov > 0 ? exp(mle.model.sum_cov) : 1.0)
           mle.comp.d2S2[ij] += mle.model.comp.d2S2[ij]
        end
    end
    np = npf
    for i in 1:npm
        gradient_update_dS_maintenance(mle, i + np,i)
        for j in 1:i
            #//i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
            ij = ind_ij(i + np, j)
           mle.comp.d2S1[ij] += mle.model.comp.d2S1[ij] * (npc > 0 ? exp(mle.model.sum_cov) : 1.0)
           mle.comp.d2S2[ij] += mle.model.comp.d2S2[ij]
            #//ii and j(<=ii) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
            ij = ind_ij(i, j)
           mle.comp.d2S3[ij] += mle.model.comp.d2S3[ij]
        end
        for j=(i + 1):(i + np)
            # //i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
            ij = ind_ij(i + np, j)
           mle.comp.d2S1[ij] += mle.model.comp.d2S1[ij] * (npc > 0 ? exp(mle.model.sum_cov) : 1.0)
           mle.comp.d2S2[ij] += mle.model.comp.d2S2[ij]
        end
    end
    np += npm #LD: npc
    if npc > 0
        for i in 1:npc
            gradient_update_dS_covariate(mle, i + np, i)
            for j in 1:np
                ij = ind_ij(i + np, j)
                mle.comp.d2S1[ij] += covariate(mle.model, i) * exp(mle.model.sum_cov) * mle.model.comp.dS1[j]
            end
            for j in (np + 1):(np+i)#LD:i
                ij = ind_ij(i + np, j)
                mle.comp.d2S1[ij] += covariate(mle.model, i) * covariate(mle.model, j - np) * exp(mle.model.sum_cov) * mle.model.comp.S1
            end
        end
    end
end

function hessian_update_current(mle::MLE)
    npf = mle.model.nb_params_family - 1
    npm = mle.model.nb_params_maintenance
    contrast_update_current(mle, gradient =true, hessian = true)

    cumhVright_param_derivative = cumulative_hazard_rate_param_derivative(mle.model.family, mle.model.Vright, true)
    cumhVleft_param_derivative = cumulative_hazard_rate_param_derivative(mle.model.family, mle.model.Vleft, false)
    hVleft_param_derivative = hazard_rate_param_derivative(mle.model.family, mle.model.Vleft, false)
    cumhVright_param_2derivative = cumulative_hazard_rate_param_2derivative(mle.model.family, mle.model.Vright,true)
    cumhVleft_param_2derivative = cumulative_hazard_rate_param_2derivative(mle.model.family, mle.model.Vleft,false)
    hVleft_param_2derivative = hazard_rate_param_2derivative(mle.model.family, mle.model.Vleft)
    for i in 1:npf
        if mle.model.k >= mle.left_censor 
            mle.model.comp.dS1[i] +=  cumhVleft_param_derivative[i] - cumhVright_param_derivative[i]
            #print("j0: $(mle.model.comp.dS1[i])\n")##LDprint
        end
        mle.model.comp.dS2[i] += hVleft_param_derivative[i] / mle.model.hVleft * mle.model.indType ;
        #print("j1: $(mle.model.comp.dS2[i])\n")##LDprint
        for j in 1:i
            ij = ind_ij(i, j)
            if mle.model.k >= mle.left_censor 
                mle.model.comp.d2S1[ij] += cumhVleft_param_2derivative[ij] - cumhVright_param_2derivative[ij]
                #print("j2: $(mle.model.comp.d2S1[ij])\n")##LDprint
            end
            mle.model.comp.d2S2[ij] += (hVleft_param_2derivative[ij] / mle.model.hVleft - hVleft_param_derivative[i] * hVleft_param_derivative[j] / mle.model.hVleft^2) * mle.model.indType
            #print("j3: $(mle.model.comp.d2S2[ij])\n")##LDprint
        end
    end
    hVright = hazard_rate(mle.model.family, mle.model.Vright)
    dhVleft = hazard_rate_derivative(mle.model.family, mle.model.Vleft)
    dhVright = hazard_rate_derivative(mle.model.family, mle.model.Vright)
    hVright_param_derivative = hazard_rate_param_derivative(mle.model.family, mle.model.Vright,true)
    dhVleft_param_derivative = hazard_rate_derivative_param_derivative(mle.model.family, mle.model.Vleft)
    d2hVleft = hazard_rate_2derivative(mle.model.family, mle.model.Vleft)
    # //printf("k:%d,hVright:%lf,dhVleft:%lf,indType:%lf\n",model.k,hVright,dhVleft,model.indType);
    for i in 1:npm
        if mle.model.k >= mle.left_censor
            mle.model.comp.dS1[i + npf] += mle.model.hVleft * mle.model.dVleft[i] - hVright * mle.model.dVright[i];
            #print("j4: $(mle.model.comp.dS1[i + npf])\n")##LDprint
        end
        # //printf("dS1[%d]=(%lf,%lf,%lf),%lf,",i+1,model.hVleft,model.dVleft[i],model.dVright[i],model.dS1[i+1]);
        mle.model.comp.dS2[i + npf] +=  dhVleft * mle.model.dVleft[i] / mle.model.hVleft * mle.model.indType
        #print("j5: $(mle.model.comp.dS2[i + npf])\n")##LDprint
        #//printf("dS2[%d]=%lf,",i+1,model.dS2[i+1]);
        #//column 0 and i+1 corresponds to the line indice of (inferior diagonal part of) the hessian matrice
        mle.model.comp.dS3[i] +=  mle.model.dA[i] / mle.model.A * mle.model.indType;
        for j in 1:npf
            ij = ind_ij(i + npf, j)
            if mle.model.k >= mle.left_censor
                mle.model.comp.d2S1[ij] += hVleft_param_derivative[j] * mle.model.dVleft[i] - hVright_param_derivative[j] * mle.model.dVright[i]
                #print("j6: $(mle.model.comp.d2S1[ij])\n")##LDprint
            end
            mle.model.comp.d2S2[ij] +=  dhVleft_param_derivative[j] * mle.model.dVleft[i] / mle.model.hVleft * mle.model.indType - hVleft_param_derivative[j] * dhVleft * mle.model.dVleft[i] / mle.model.hVleft^2 * mle.model.indType
            #print("j7: $(mle.model.comp.d2S2[ij])\n")##LDprint
        end
        for j=1:i
            #//i+1 and j+1(<=i+1) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
            ijp = ind_ij(i + npf, j + npf)
            ij = ind_ij(i, j)
            if mle.model.k >= mle.left_censor
                #print("j-8: $(mle.model.comp.d2S1[ijp])\n")##LDprint
                #print("jj8: 0:$ijp,$ij, 1:$(dhVleft),  2:$(mle.model.dVleft[i]), 3:$(mle.model.dVleft[j]), 4:$(mle.model.hVleft), 5:$(mle.model.d2Vleft[ij]), 6:$(dhVright), 7:$(mle.model.dVright[i]), 8:$(mle.model.dVright[j]), 9:$(hVright), 10:$(mle.model.d2Vright[ij])\n")
                mle.model.comp.d2S1[ijp] += dhVleft * mle.model.dVleft[i] * mle.model.dVleft[j] + mle.model.hVleft * mle.model.d2Vleft[ij] - dhVright * mle.model.dVright[i] * mle.model.dVright[j] - hVright * mle.model.d2Vright[ij]
                #print("j8:  $(mle.model.comp.d2S1[ijp])\n")##LDprint
            end
            #print("j-9:  $(mle.model.comp.d2S2[ijp])\n")##LDprint
            mle.model.comp.d2S2[ijp] += (mle.model.dVleft[i] * mle.model.dVleft[j] * (d2hVleft / mle.model.hVleft - (dhVleft / mle.model.hVleft)^2) + dhVleft * mle.model.d2Vleft[ij] / mle.model.hVleft) * mle.model.indType
            #print("j9:  $(mle.model.comp.d2S2[ijp])\n")##LDprint
            mle.model.comp.d2S3[ij] += (mle.model.d2A[ij] / mle.model.A - mle.model.dA[i] * mle.model.dA[j] / mle.model.A^2) * mle.model.indType
        end
    end
end

function αEst(mle::MLE, param::Vector{Float64})
    #contrast(mle, param) #//To compute mle.comp.S1 and S0
    return mle.comp.S0 / mle.comp.S1
end

data!(mle::MLE, i::Int) = data!(mle.model, i)
data(mle::MLE) = data(mle.model)

#     //delegate from model cache!
#     List get_virtual_age_infos(double by,double from, double to) {
#         return mle.model.get_virtual_age_infos(by,from,to);
#     }

#     DataFrame get_selected_data(int i) {
#         return mle.model.get_selected_data(i);
#     }

#     void set_covariates(List covariates_) {
#         mle.model.set_covariates(covariates_);
#     }

