import Base.rand

function rand(model::Model, stop::Union{Int, Vector{Any}}; system::Int=1, datacov::DataFrame=DataFrame())::DataFrame
    stop_policy_ = stop_policy(stop)
    return rand(model, stop_policy_, system = system, datacov = datacov)
end

rand(model::Model; system::Int = 1, datacov::DataFrame=DataFrame()) = rand(model, 100, system = system, datacov = datacov)

function rand(model::Model, stop_policy::Expr; system::Int=1, datacov::DataFrame=DataFrame())::DataFrame
    if has_maintenance_policy(model)
        first(model.maintenance_policy)
    end
    if !isempty(datacov)
        covariates!(model, datacov)
    end
    if model.nb_params_cov > 0
        system = size(model.datacov)[1]
    end

    data = DataFrame()
    for syst in 1:system
        init_sim!(model)
        run = true
        while run
            u = log(model.rand())::Float64
            if model.nb_params_cov > 0
            #   u *= compute_covariates(sim) #;//set_current_system launched in R for simulation
            end
            timeCM = virtual_age_inverse(model, inverse_cumulative_hazard_rate(model.family, cumulative_hazard_rate(model.family, virtual_age(model,model.time[model.k]))-u))
            #   TODO: submodels
            id_mod = 0
            if has_maintenance_policy(model)
                timePM, typePM = update(model.maintenance_policy, model) # //# Peut-être ajout Vright comme argument de update
                if timePM < timeCM && timePM < model.time[model.k]
                    #print("Warning: PM ignored since next_time(=%lf)<current_time(=%lf) at rank %d.\n",timePM,model->time[model->k],model->k);
                    print("warning")
                end
            end
            if !has_maintenance_policy(model) || timeCM < timePM || timePM < model.time[model.k]
                push!(model.time,timeCM)
                push!(model.type, -1)
                id_mod=0
            else
                push!(model.time, timePM)
                #//DEBUG[distrib type1]: typeCptAP++;if(typePM==1) type1CptAP++;printf("typePM=%d\n",typePM);
                push!(model.type, typePM)
                id_mod=typePM
            end
            update_Vleft!(model) #, false,false)
            update_maintenance!(model, id_mod) #false,false)
            run = ok(model, stop_policy)
            ## TODO work on stop later
        end
        data = vcat(data,DataFrame(system=syst, time=model.time, type=model.type)[2:length(model.time),:]) #LD: vcat(data,DataFrame(system=syst, time=sim.model.time, type=sim.model.type))
    end
    ## println(data)
    if system == 1
        data = data[:,[:time, :type]]
    end
    df = data #LD: df = data[2:size(data,1),:]
    if (system > 1) && (length(model.varnames)==2) #LD: (system > 1)
        rename!(df, vcat(["system"], model.varnames) )
    else
        rename!(df, model.varnames)
    end
    data!(model, df)
    df
end

function init_sim!(model::Model)
    #// Almost everything in the 5 following lines are defined in model->init_computation_values() (but this last one initializes more than this 5 lines)
    model.Vright=0
    model.A=1
    model.k=1
    for mm in model.models
        init!(mm)
    end
    model.id_mod=0 #// Since no maintenance is possible!
    model.time = [0.0]
    model.type = [-1]
end

function ok(model::Model, stop_policy::Expr)::Bool
    s = length(model.time) - 1 # 1st is 0 time to be removed when returned
    t = model.time[model.k]
    eval(:(s=$s))
    eval(:(t=$t))
    eval(stop_policy)
end

function stop_policy(stop::Union{Nothing, Int, Vector{Any}})::Expr
    return if stop isa Int
        Expr(:call, :<, :s,  stop)
    elseif isnothing(stop)
        Expr(:call, :<, :s,  100)
    else  
        formula_translate(Expr(:call, stop...))
    end
end

### Simulator is no more than model with stop_policy embedded together
mutable struct Simulator
    model::Model
    stop_policy::Expr
end

function Simulator(model::Model, stop::Union{Nothing, Int, Vector{Any}})::Simulator
    sim = Simulator(model, stop_policy(stop))
    init_sim!(model)
    return sim
end

rand(sim::Simulator; system::Int=1, datacov::DataFrame=DataFrame())::DataFrame = rand(sim.model, sim.stop_policy; system = system, datacov = datacov)

stop_policy!(sim::Simulator, stop::Union{Nothing, Int,Vector{Any}}) = sim.stop_policy = stop_policy(stop)