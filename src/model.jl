mutable struct Model <: AbstractModel
	# Main info
    k::Int # current index
	Δt::Float64
	time::Vector{Float64}
	type::Vector{Int}
	current_system::Int #current system

	# variable names
	varnames::Vector{String}

	# number
	nb_system::Int #number of system
	nbPM::Int
	nb_params_maintenance::Int
    nb_params_family::Int
    nb_params_cov::Int

	# index
	id_mod::Int # current maintenance model
	id_params::Int # current index in maintenance model parameters
	id_params_list::Vector{Int}
	mu::Int # max memory of the model

	# Model specification
	family::FamilyModel
	models::Vector{AbstractMaintenanceModel}
	maintenance_policy::AbstractMaintenancePolicy
	data::Vector{DataFrame}

	# For computation
	indType::Float64
    Vleft::Float64
    Vright::Float64
    hVleft::Float64
	dVleft::Vector{Float64}
    dVright::Vector{Float64}
	d2Vleft::Vector{Float64}
    d2Vright::Vector{Float64}
	A::Float64
	dA::Vector{Float64}
	d2A::Vector{Float64}
	VR_prec::Vector{Float64}
	dVR_prec::Vector{Float64}
    d2VR_prec::Vector{Float64}
	comp::Compute

	# Additional Covariates stuff
	datacov::DataFrame
	vars_cov::Vector{Symbol}
	params_cov::Vector{Float64}
	sum_cov::Float64 #to save the computation
	expr_cov::Union{Nothing, Expr}

	# Formula
	formula::Expr

	Model() = begin
		m = new()
		# init!(m)
		init_covariates!(m)
		return m
	end
end

function init!(m::Model)
		m.k = 1  # current position in data
		m.Δt = 0
		m.current_system = 1 # current system
		m.nb_system = 1 #number of system
		if isdefined(m, :models) && length(m.models) > 0
			m.nbPM = length(m.models) - 1
		else
			m.nbPM = 0
		end
		m.id_mod = 0
		m.id_params = 1
		m.nb_params_maintenance=0
		cur_id_params = 1
		m.id_params_list = Int[]
		for (id, mm) in enumerate(m.models)
			push!(m.id_params_list, cur_id_params)
			m.nb_params_maintenance += nbparams(mm)
			cur_id_params += nbparams(mm)
		end
		m.nb_params_family = nbparams(m.family)
		m.mu = max_memory(m)

		m.data=DataFrame[]

		## internal
		m.time = Float64[]
		m.type = Int[]

		m.indType = 0

		m.Vleft = 0
		m.Vright = 0
		m.hVleft = 0


		m.dVleft = zeros(m.nb_params_maintenance)
		m.dVright = zeros(m.nb_params_maintenance)
		
		nb2d = m.nb_params_maintenance * (m.nb_params_maintenance + 1) ÷ 2

		m.d2Vleft = zeros(nb2d)
		m.d2Vright = zeros(nb2d)

		m.A = 0
		m.dA = zeros(m.nb_params_maintenance)
		m.d2A = zeros(nb2d)
		if m.mu > 0
			m.VR_prec = zeros(m.mu)
			m.dVR_prec = zeros(m.mu * m.nb_params_maintenance)
			m.d2VR_prec = zeros(m.mu * nb2d)
		end
		m.comp = Compute(m)
		return nothing
end

function init_covariates!(m::Model)
	# Additional Covariates stuff
	m.nb_params_cov = 0
	m.datacov = DataFrame()
	m.params_cov = Parameter[]
	m.sum_cov = 0 #to save the computation
	m.expr_cov = nothing
end

function inc!(m::Model)
	m.k += 1
	m.Δt = m.time[m.k] - m.time[m.k - 1]
end

nbparams(m::Model)::Int = m.nb_params_family + m.nb_params_maintenance + m.nb_params_cov

params(m::Model)::Parameters = cat(params(m.family),(map(m.models) do mm;params(mm); end)...,m.params_cov,dims=1)

function params!(m::Model, θ::Vector{Float64})
	from, to = 1, nbparams(m.family)
	params!(m.family,θ[from:to])
	for mm in m.models
		if nbparams(mm) > 0
			from = to + 1
			to = from + nbparams(mm) - 1
			params!(mm, θ[from:to])
		end
	end
	if m.nb_params_cov > 0
		for i in 1:m.nb_params_cov
			m.params_cov[i] = θ[to + i]
		end
	end
end

priors(m::Model)::Priors = cat(m.family.priors,(map(m.models) do mm;mm.priors; end)...,dims=1)

function init_compute!(m::Model)
	init!(m.comp)
	for mm in m.models
		init!(mm)
	end
end

virtual_age(m::Model, x::Float64)::Float64 = m.Vright + (x  - m.time[m.k]) * m.A
virtual_age_inverse(m::Model, x::Float64) = (x - m.Vright) / m.A + m.time[m.k]

function update_Vleft!(m::Model; gradient::Bool=false, hessian::Bool=false)
	# /*if(model->k < 10) printf("Vleft:%lf\n", model->Vleft);*/
	m.Vleft = virtual_age(m, m.time[m.k + 1])
	# //printf("Vleft:%lf\n", model->Vleft);
	if hessian
		if m.nb_params_maintenance > 0
			for i in 1:m.nb_params_maintenance
				m.dVleft[i] = m.dVright[i] + (m.time[m.k+1] - m.time[m.k]) * m.dA[i]
				for j in 1:i
					ij = ind_ij(i, j)
					m.d2Vleft[ij]= m.d2Vright[ij] + (m.time[m.k+1]  - m.time[m.k]) * m.d2A[ij]
				end
			end
		end
	elseif gradient
		if m.nb_params_maintenance > 0
			for i in 1:m.nb_params_maintenance
				m.dVleft[i]= m.dVright[i] + (m.time[m.k+1]  - m.time[m.k]) * m.dA[i]
			end
		end
	end
end

function update_maintenance!(m::Model, id_mod::Int; gradient::Bool=false, hessian::Bool=false)
	# id_params is used inside update! for gradient and hessian
	m.id_params = m.id_params_list[1 + id_mod]
	# println("up maint $id_mod $(m.id_params)")
	update!(m.models[1 + id_mod], m; gradient=gradient, hessian=hessian)
	m.id_mod = id_mod
end


function data!(m::Model,data::Union{DataFrame,Vector{DataFrame}})
	if data isa Vector{DataFrame}
		# TODO: considering maybe sub-dataframe form m.varnames
		m.data = data
		m.nb_system = length(data)
	else
		m.data=DataFrame[]
		if m.varnames ∩ names(data) == m.varnames
			data2 = data[:, m.varnames]
			if size(data2, 2) == 2
				prepend!(data2[!,1], 0)
				prepend!(data2[!,2], 1)
				push!(m.data, data2)
				m.nb_system = 1
			else # multi-system with size(data2, 2) ==3
				systs = sort(unique(data2[!,1]))
				for i=systs
					data3 = data2[data2[!,1] .== i, 2:3]
					prepend!(data3[!,1], 0)
					prepend!(data3[!,2], 1)
					push!(m.data, data3)
				end
				m.nb_system = length(m.data)
			end
		elseif size(data,2) == 2
			df = vcat(DataFrame(time=0, type=1),data)
			# the 2 column are renamed 
			rename!(df,m.varnames)
			push!(m.data, df)
			
			m.nb_system = 1
		elseif size(data,2) == 3
			m.nb_system=maximum(data[!,1])
			for syst in sort(unique(data[!,1]))
				data3 = data[data[!,1].==syst,2:3]
				prepend!(data3[!,1], 0)
				prepend!(data3[!,2], 1)
				rename!(data3,length(m.varnames) == 3 ? m.varnames[2:end] : m.varnames)
				push!(m.data, data3)
			end
		end
	end
	data!(m,1) #;//default when only one system no need to
end

function data!(m::Model, i::Int)
	if length(m.data) >= i
		data=m.data[i]
		m.time = data[!,1]
		m.type = data[!,2]
		m.current_system = i
	end
end

function data(m::Model, i::Int)::DataFrame
	data!(m, i) #;//Skipped if data is unset (see above)
	return DataFrame(time=m.time,type=m.type)
end

function init_virtual_age_infos(m::Model)
		# int i;
    	# k=0;
    	# idMod=0; //id of current model
    	# S1 = 0;
    	# Vright=0;
    	# A=1;
    	# for(i=0;i<nbPM + 1;i++) models->at(i)->init();
end

function get_virtual_age_info(m::Model, from::Float64, to::Float64, by::Float64, expCov::Float64)
# 	double s=ceil((to-from)/by);
# 	int n=static_cast<int>(s);
# //printf("ici=%d,%lf (%lf,%lf,%lf)\n",n,s,to,from,by);
# 	std::vector<double> t(n+1);
# 	std::vector<double> v(n+1);
# 	std::vector<double> h(n+1); //i as intensity
# 	std::vector<double> H(n+1); //I for cumulative intensity
# 	std::vector<double> F(n+1); //F for conditional cumulative distribution function
# 	std::vector<double> S(n+1); //S for conditional survival function
# 	std::vector<double> f(n+1); //S for conditional survival function

# 	t[0]=from;t[n]=to;
# 	v[0]=virtual_age(from);v[n]=virtual_age(to);
# 	h[0]=expCov*A*family->hazardRate(v[0]);h[n]=expCov*A*family->hazardRate(v[n]);
# 	H[0]=S1;H[n]=S1+expCov*(family->cumulative_hazardRate(v[n])-family->cumulative_hazardRate(v[0]));
# 	F[0]=0;F[n]=1-exp(-expCov*(family->cumulative_hazardRate(v[n])-family->cumulative_hazardRate(v[0])));
# 	S[0]=1;S[n]=exp(-expCov*(family->cumulative_hazardRate(v[n])-family->cumulative_hazardRate(v[0])));
# 	f[0]=expCov*A*family->hazardRate(v[0]);f[n]=expCov*A*family->hazardRate(v[n])*exp(-expCov*(family->cumulative_hazardRate(v[n])-family->cumulative_hazardRate(v[0])));
# 	double by_t=(t[n]-t[0])/s;
# 	double by_v=(v[n]-v[0])/s;

# 	for(int i=1;i<n;i++) {
# 		t[i]=t[i-1]+by_t;//printf("t[%d]=%lf\n",i,t[i]);
# 		v[i]=v[i-1]+by_v;
# 		h[i]=expCov*A*family->hazardRate(v[i]);
# 		H[i]=S1+expCov*(family->cumulative_hazardRate(v[i])-family->cumulative_hazardRate(v[0]));
# 		F[i]=1-exp(-expCov*(family->cumulative_hazardRate(v[i])-family->cumulative_hazardRate(v[0])));
# 		S[i]=exp(-expCov*(family->cumulative_hazardRate(v[i])-family->cumulative_hazardRate(v[0])));
# 		f[i]=expCov*A*family->hazardRate(v[i])*exp(-expCov*(family->cumulative_hazardRate(v[i])-family->cumulative_hazardRate(v[0])));
# 	}

# 	return DataFrame::create(
# 		_["t"]=NumericVector(t.begin(),t.end()),
# 		_["v"]=NumericVector(v.begin(),v.end()),
# 		_["i"]=NumericVector(h.begin(),h.end()),
# 		_["I"]=NumericVector(H.begin(),H.end()),
# 		_["F"]=NumericVector(F.begin(),F.end()),
# 		_["S"]=NumericVector(S.begin(),S.end()),
# 		_["f"]=NumericVector(f.begin(),f.end())
# 	);
end

function get_virtual_age_infos(m::Model, from::Float64, to::Float64, by::Float64)

	# // Only one system first!
	# init_virtual_age_infos();
	# double expCov=1;
	# if (nb_paramsCov>0) expCov=exp(compute_covariates());
	# int n=time.size() - 1;
	# List res(n);
	# while(k < n) {
	# 	//printf("k=%d/n=%d,(%lf,%lf)\n",k,n,time[k],time[k+1]);
	# 	update_Vleft(false,false);
	# 	if(from > time[k] || time[k+1] > to ) res[k] = R_NilValue;
	# 	else res[k]=get_virtual_age_info(time[k],time[k+1],by,expCov);
	# 	S1 += expCov*(family->cumulative_hazardRate(Vleft) - family->cumulative_hazardRate(Vright));
	# 	//gradient_update_for_current_system();
	# 	int type2=type[k + 1];
	# 	if(type2 < 0) type2=0;
	# 	models->at(type2)->update(false,false);
	# }
	# return res;
end

function select_current_system(m::Model, i::Int, compute::Bool)
	#//Covariates related
	m.current_system = i
	#//simulation: compute=false since only computation in c++ and set_current_system in R
	#//mle: compute=true since both computation and select_current_system in c++
	if compute
		compute_covariates(m)
	end
end

# //Covariates related
function covariates!(m::Model,formula::Expr)
	if !isnothing(formula)
		m.expr_cov = formula
		m.params_cov = Parameter[]
		m.vars_cov = Symbol[]
		for (p, v) in map(x -> x.args[2:end], formula.args[2:end])
			push!(m.params_cov, p)
			push!(m.vars_cov, v)
		end
		m.nb_params_cov = length(m.params_cov)
	end
end

covariates!(m::Model) = covariates!(m, m.expr_cov)

function covariates!(m::Model, data::DataFrame) 
	m.datacov = data[!,m.vars_cov]
	m.nb_params_cov = size(m.datacov)[2]
end

function compute_covariates(m::Model)
	m.sum_cov = 0.0
	# println(m.params_cov)
	# println(m.datacov)
	# println( (m.current_system, m.nb_params_cov, m.params_cov) )
	for j in 1:m.nb_params_cov
		m.sum_cov += m.params_cov[j] * m.datacov[m.current_system, m.vars_cov[j]]
		# //printf("syst=%d,j=%d,th=%lf,params_cov=%lf\n",current_system,j,params_cov[j],var[current_system]);
	end
	return m.sum_cov
end

covariate(m::Model, j::Int) = m.datacov[m.current_system, j]

has_maintenance_policy(m::Model)::Bool = isdefined(m,:maintenance_policy) #|| !isnothing(m.maintenance_policy)

function max_memory(m::Model)::Int
	maxmem = 1
	for mm in m.models
		if isdefined(mm,:m)
			if mm.m > maxmem
				maxmem = mm.m
			end
		end
	end
	return maxmem
end

function isbayesian(m::Model)::Bool
	return all(map(isbayesian, m.models))
end

# Used inside ModelTest do guess the r formula and RData
function rterms(m::Model, data::DataFrame, datacov::DataFrame=DataFrame())
	f = m.formula
	df_rexpr = string("data.frame(", 
		join(map(names(data)) do var
			string(var,"=c(",join(data[!,var],","),")")
		end,","),
	")")
	vars_str = replace(string(f.args[2]),"(" => "", ")" => "")
	model_str = ""
	if f.args[3].args[1] == :&
		# CM + PM
		cm_str = string("(",f.args[3].args[2],")")
		pm_str = string("(",f.args[3].args[3],")")
		model_str = cm_str * " & " * pm_str
	elseif f.args[3].args[1] == :|
		# CM only
		model_str = string("(",f.args[3],")")
	end
	vam_repxr = replace(string(vars_str, " ~ ", replace( model_str  ,"FamilyModel" => "")),"∞" => "Inf")
	dfcov_rexpr = if isempty(datacov)
		"NULL"
	else
		string("data.frame(", 
			join(map(names(datacov)) do var
				string(var,"=c(",join(datacov[!,var],","),")")
			end,","),
		")")
	end
	[df_rexpr, vam_repxr, dfcov_rexpr]
end