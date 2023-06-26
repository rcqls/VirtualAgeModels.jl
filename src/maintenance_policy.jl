abstract type AbstractMaintenancePolicy end
function first(mp::AbstractMaintenancePolicy); end

abstract type MaintenancePolicyWithExternalModel <: AbstractMaintenancePolicy end
function update_external_model(mp::MaintenancePolicyWithExternalModel, model::AbstractModel)
    mod = mp.model # external model attached to MP

    if isnothing(mod)
        return model
    else

        
        # //update everything needed to compute the update step (see end of simulation step)
        mod.k = model.k
        mod.time = model.time  
        mod.type = model.type
        if mod.k > 0 # Not at the init step!
            ## VERY VERY IMPORTANT: Since we want to update at previous step
            mod.k -= 1

            # And then as what it is applied at the end of the simulation step but for mod instead of model!
            update_Vleft(mod) #, false,false) #; //mod->idMod is not yet updated!
            mod.idMod = model.idMod #//mod->idMod is then updated for the next task!
            update!(mod.models[mod.idMod], mod) #,false,false)
        end
        if model.nb_paramsCov > 0
            mod.datacov=model.datacov
            mod.params_cov=model.params_cov
            mod.nb_paramsCov=model.nb_params_cov
        end
    end
    
end
mutable struct PeriodicMaintenancePolicy <: AbstractMaintenancePolicy
    from::Float64
    by::Float64
    prob::Vector{Float64}
end

PeriodicMaintenancePolicy(by::Real, prob::Vector{Float64}) = PeriodicMaintenancePolicy(0.0, by, prob)
PeriodicMaintenancePolicy(from::Real, by::Real) = PeriodicMaintenancePolicy(from, by, [1])
PeriodicMaintenancePolicy(by::Real) = PeriodicMaintenancePolicy(0, by)

type_size(mp::PeriodicMaintenancePolicy)::Int = length(mp.prob)

function first(mp::PeriodicMaintenancePolicy); end

function update(mp::PeriodicMaintenancePolicy, model::AbstractModel)::NamedTuple{(:time, :type), Tuple{Float64, Int64}}
    current=model.time[model.k]
	time = mp.from + (floor((current - mp.from)/mp.by) + 1) * mp.by
 
	r=rand(1)[1]
    t = 1
	for p in mp.prob[1:end-1]
        if r < p
            break
        else 
            t += 1
            r -= p
        end
    end
	# //printf("from=%d,t=%d, %lf\n",get_from_type(),t,prob[0]);
	return (time=time, type=t)
end


mutable struct AtTimesMaintenancePolicy <:  AbstractMaintenancePolicy
    times::Vector{Float64}
    i::Int
    k::Int
    cycle::Bool
    differentTypeIfCM::Bool
end

type_size(mp::AtTimesMaintenancePolicy)::Int = mp.differentTypeIfCM ? 2 : 1

mutable struct AtIntensityMaintenancePolicy <: MaintenancePolicyWithExternalModel
    level::Float64
    model::Union{Nothing, AbstractModel}
end
AtIntensityMaintenancePolicy(level::Float64) = AtIntensityMaintenancePolicy(level, nothing)

type_size(mp::AtIntensityMaintenancePolicy)::Int = 1
function update(mp::AtIntensityMaintenancePolicy, model::AbstractModel)::NamedTuple{(:time, :type), Tuple{Float64, Int64}}
    
    mod=update_external_model(mp, model)
	u = mod.A
	if mod.nb_params_cov > 0 
        # u *= exp(compute_covariates(mod))
    end
	# 	//ToRemove: double next_time2=model->virtual_age_inverse(model->family->inverse_hazardRate(level[0]));
	time=virtual_age_inverse(mod,inverse_hazard_rate(mod.family, mp.level/u))
    return (time=time, type=1)
end

mutable struct AtVirtualAgeMaintenancePolicy <:  MaintenancePolicyWithExternalModel
    level::Float64
    external_model::AbstractModel
end

type_size(mp::AtVirtualAgeMaintenancePolicy)::Int = 1


mutable struct AtFailureProbabilityMaintenancePolicy <: AbstractMaintenancePolicy
    level::Vector{Float64}
    external_model::AbstractModel
end

type_size(mp::AtFailureProbabilityMaintenancePolicy)::Int = 1
 

struct MaintenancePolicyList <: AbstractMaintenancePolicy
    policies::Vector{AbstractMaintenancePolicy}
    from_type::Vector{Int}
end

function MaintenancePolicyList(policies::Vector{AbstractMaintenancePolicy})
    from_type = [0]
    from = 0
    for policy in policies[1:end - 1]
        from += type_size(policy)
        push!(from_type, from)
    end
    MaintenancePolicyList(policies, from_type)
end

function type_size(mp::MaintenancePolicyList)::Int
    s = 0
    for policy in mp.policies
        s += type_size(policy)
    end
    return s;
end

function first(mp::MaintenancePolicyList)
	for policy in mp.policies
        first(policy)
    end
end


function update(mp::MaintenancePolicyList,model::AbstractModel)
    time, type = update(mp.policies[1], model)
    ts = type_size(mp.policies[1])
    for policy in mp.policies[2:end] 
    	time2, type2 = update(policy, model)
    	##println("$time2 < $time ? ")
        if time2 < time 
            time, type = time2, type2
            type += ts
        end
        ts += type_size(policy)
    end
    return (time=time, type=type)
end