# Time & Type ~ (ABAO()|Weibull(1,3)
# Time & Type ~ (ARA1(~Beta(1.08,0.108)) | Weibull(~NonInform(),~Gamma(32,0.097)))

# ~ (ARA1(.5) | Weibull(0.01,2.5)) & (ARAInf(.7)+ARAInf(.3)|Periodic(12,prob=c(0.6,0.4))))

# Systeme & Temps & Type ~ (ARA1(.5) | Weibull(0.01,2.5)) & (ARAInf(.7)+ARAInf(.3)|Periodic(12,prob=c(0.6,0.4)))

function parse_model(ex_f::Expr)
    m = Model()
    m.formula = ex_f
    if Meta.isexpr(ex_f, :call)
        ex_m = ex_f
        varnames = ["Time", "Type"] # default names
        ## model detection
        if ex_f.args[1] == :&
            ## No names given on the left side of :~ => def_names used
            ex_m.args[2] = ex_m.args[2].args[2] # remove unused the tilde
        elseif ex_f.args[1] == :~
            ## Left part is names and right part is the model
            vars = if Meta.isexpr(ex_f.args[2].args[2], :call) # System as first variable
                vcat(ex_f.args[2].args[2].args[2:3], ex_f.args[2].args[3])
            else
                ex_f.args[2].args[2:3]
            end
            varnames = string.(vars)
            ex_m = ex_f.args[3]
        end
        m.varnames = varnames 
        ## parsing model (ex_m)
        m.models = AbstractMaintenanceModel[]
        #print(ex_m.args)
        if ex_m.args[2].args[1] == :|
            ex_cm = ex_m.args[2]
            parse_cm!(m, ex_cm)
            ## PMs (Preventive Maintenances) and MPs (Maintenance Policies)
            ex_pm = ex_m.args[3]
            if ex_pm.args[1] == :|
                # PMs
                ex_pms = ex_pm.args[2]
                if ex_pms.args[1] == :+
                    # several PMs
                    for pm in ex_pms.args[2:end]    
                        #push!(m.models,eval(pm))
                        add_maintenance_model!(m, pm)
                    end
                else
                    # only 1 PM
                    #push!(m.models,eval(ex_pms))
                    add_maintenance_model!(m, ex_pms)
                end
                # Maintenance policies
                ex_mps = ex_pm.args[3]
                if ex_mps.args[1] == :*
                    # several MPs
                    maintenance_policies = AbstractMaintenancePolicy[]
                    for mp in ex_mps.args[2:end]
                        push!(maintenance_policies,eval(complete_name!(mp,1,"MaintenancePolicy")))
                    end
                    m.maintenance_policy = MaintenancePolicyList(maintenance_policies)
                else
                    # only 1 MP
                    m.maintenance_policy = eval(complete_name!(ex_mps,1,"MaintenancePolicy"))
                end
            else
                # No MP (maintenance policy) only PMs

                if ex_pm.args[1] == :+
                    # several PMs
                    for pm in ex_pm.args[2:end]
                        ## push!(m.models,eval(pm))
                        add_maintenance_model!(m, pm)
                    end
                else
                    # only 1 PM
                    #push!(m.models,eval(ex_pm))
                    add_maintenance_model!(m, ex_pm)
                end 

            end
        else
            ##No PM
            parse_cm!(m, ex_m)
        end
        return m
    end
end

function parse_cm!(m::AbstractModel,ex_cm::Expr)
    if ex_cm.args[1] == :|
        #push!(m.models,eval(ex_cm.args[2])) # CM (Corrective Maintenance)
        add_maintenance_model!(m,ex_cm.args[2])
        fm = ex_cm.args[3]
        #fm.args[1] = complete_name(fm.args[1], "FamilyModel")
        if isa(fm.args[2], Expr) && fm.args[2].head == :parameters
            # covariates
            fm2=Expr(:call, fm.args[1], fm.args[3:end]..., fm.args[2].args[1].args)
            #m.family = eval(complete_name!(fm2, 1, "FamilyModel"))
            add_family_model!(m,fm2)
        else
            #m.family = eval(complete_name!(fm, 1, "FamilyModel"))
            add_family_model!(m,fm)
        end
    end
end

function parse_bayesian_parameters!(ex::Expr)
    ex.args = map(e -> Meta.isexpr(e, :call) && e.args[1] == :~  ? e.args[2] : e , ex.args)
end

function parse_covariates(ex_fm::Expr)
    ex = copy(ex_fm)
    tmp = findall(e -> Meta.isexpr(e,:call) && e.args[1] in [:|,:+, :-], ex_fm.args)
    if !isempty(tmp)
        index = tmp[1]
        op = ex.args[index].args[1]
        if op == :|
            # Weibull(0.001,2.5| 1*cov1)
            ex.args = vcat(ex_fm.args[1:index-1],ex_fm.args[index].args[2])
            return (ex,Expr(:call,:+,ex_fm.args[index].args[3]))
        elseif op in  [:+, :-]
            # Weibull(0.001,2.5| 1*cov1 + -2cov2 + 3cov3)
            ex.args = vcat(ex_fm.args[1:index-1],ex_fm.args[index].args[2].args[2])
            return (ex, Expr(:call, op, ex_fm.args[index].args[2].args[3],ex_fm.args[index].args[3:end]...))
        end
    end
    return (ex,)
end

function add_family_model!(m::AbstractModel,ex_fm::Expr)
    res = parse_covariates(ex_fm)
    if length(res) == 1
        parse_bayesian_parameters!(ex_fm)
        m.family = eval(complete_name!(ex_fm, 1, "FamilyModel"))
    elseif length(res) == 2
        parse_bayesian_parameters!(res[1])
        m.family = eval(complete_name!(res[1], 1, "FamilyModel"))
        covariates!(m, res[2])
    end
end

function add_maintenance_model!(m::AbstractModel,ex_mm::Expr)
    pipe_index = findall(e -> Meta.isexpr(e,:call) && e.args[1] == :|, ex_mm.args)
    if !isempty(pipe_index) && length(pipe_index) == 1
    # OLD: if Meta.isexpr(ex_mm.args[end],:call) && ex_mm.args[end].args[1] == :|
        # OLD: ex: GQR(...|f) -> GQR(...,f) but not GQR_ARAm(...|3,f) -> GQR_ARAm(...,3,f)
        ex = copy(ex_mm) # Do not change original ex_mm
        # OLD: ex.args = vcat(ex.args[1:end-1],ex.args[end].args[2:end])
        i = pipe_index[1]
        ex.args = vcat(ex.args[1:i-1],ex.args[i].args[2:end],ex.args[i+1:end])
        parse_bayesian_parameters!(ex)
        push!(m.models,eval(ex))
    else
        parse_bayesian_parameters!(ex_mm)
        push!(m.models,eval(ex_mm))
    end
end

macro vam(ex_f)
    return parse_model(ex_f)
end

macro stop(ex_s)
    return formula_translate(ex_s).args
end

function formula_translate(ex_f::Expr)
    Meta.parse(replace(string(ex_f), "size" => "s", "time" => "t"))
end

function complete_name!(ex::Expr, i::Int, append::String)::Expr
    s = string(ex.args[i])
    if (length(s) <= length(append)) || !(s[end-length(append)+1:end] == append )
        ex.args[i] = Symbol(string(ex.args[i]) * append)
    end
    return ex
end

