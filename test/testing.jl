using OrderedCollections
import Base.empty!

const ModelDict = Dict{Symbol, Any}
function empty!(md::ModelDict)
	for (k, _) in md
		delete!(md, k)
	end
end

mutable struct ModelTest
	models::ModelDict
	results::ModelDict
	ModelTest() = new(ModelDict(),ModelDict())
end

import Base.insert!
function insert!(modtest::ModelTest, models::Vararg{Pair{Symbol, ModelDict}})
	for model in models
		m = model.second
		if !(:r in  keys(m))
			m[:r] = VAM.rterms(m[:vam], m[:data], (:datacov in keys(m)) ? m[:datacov] : DataFrame())
			if :rform in keys(m)
				m[:r][2] = m[:rform]
			end

		end
		modtest.models[model.first] = m
	end
end

function update!(modtest::ModelTest, key::Symbol)
	model = modtest.models[key]
	# println("update")
	# println(model[:vam])
	θ=model[:θ]
	result = Dict()
	m = if :datacov in keys(model)
		VAM.MLE(model[:vam], model[:data], model[:datacov])
	else
		VAM.MLE(model[:vam], model[:data])
	end
	# println("update2")
	# println(m.model)
	R"""
	require(VAM)
	dataDF <- eval(parse(text=$(model[:r][1])))
	form <-  eval(parse(text=$(model[:r][2])))
	covDF <- eval(parse(text=$(model[:r][3])))
	mle <- mle.vam(form,data=dataDF, data.covariates=covDF)
	theta <- $(model[:θ])
	res <- list()
	res$lnL <- logLik(mle,theta,TRUE,FALSE,FALSE)
	res$dlnL<- logLik(mle,theta,FALSE,TRUE,FALSE)
	res$d2lnL <- logLik(mle,theta,FALSE,FALSE,TRUE)
	res$C <- contrast(mle,theta,TRUE,FALSE,FALSE)
	res$dC <- contrast(mle,theta,FALSE,TRUE,FALSE)
	res$d2C <- contrast(mle,theta,FALSE,FALSE,TRUE)
	"""
	result[:r] = @rget res
	result[:jl] = OrderedDict(
		:lnL => contrast(m, θ, profile=false),
		:dlnL => gradient(m, θ, profile=false),
		:d2lnL => hessian(m, θ, profile=false),
		:C => contrast(m, θ),
		:dC => gradient(m, θ)[2:end],
		:d2C => hessian(m, θ)[2:end,2:end]
	)
	modtest.results[key]=result
end

function update!(modtest::ModelTest)
	for (key, _) in modtest.models
		update!(modtest, key)
	end
end

function test(modtest::ModelTest;atol=0.0000000001)
	for (key, result) in modtest.results
		@testset verbose = true "model $key" begin
			@testset "result $k" for k in [:lnL, :dlnL, :C, :dC]
				@test result[:jl][k] ≈ result[:r][k] atol=atol
			end
		end
	end
end

function test(modtest::ModelTest, key::Symbol;atol=0.0000000001)
	result = modtest.results[key]
	@testset verbose = true "model $key" begin
		@testset "result $k" for k in [:lnL, :dlnL, :C, :dC]
			@test result[:jl][k] ≈ result[:r][k] atol=atol
		end
	end
end

function empty!(modtest::ModelTest; mode=:all)
	if mode ∈ [:all, :models]
		empty!(modtest.models)
	end
	if mode ∈ [:all, :results]	
		empty!(modtest.results)
	end
end