using VirtualAgeModels
using DataFrames
using RCall

##### Comparison between VAM and VirtualAgeModels.jl
randr() = convert(Float64,R"runif(1)")

randr()

function simulcompare()
	## 1) VirtualAgeModels.jl simulation with randr (using runif(1) from R thanks to RCall) 

	R"set.seed(12)"
	m = @vam(Time & Type ~ (ARAInf(.4) | Weibull(.001,2.5)))
	m.rand = randr
	dfJl = rand(m, 30)
	## 2) VAM simulation
	R"""
	require(VAM)
	set.seed(12)
	dfR = simulate(sim.vam( ~ (ARAInf(.4) | Weibull(.001,2.5))),30)
	"""
	@rget dfR
	## 3) comparison

	dfR == dfJl
end

simulcompare()

function simulcompare_multi()
	## 1) VirtualAgeModels.jl simulation with randr (using runif(1) from R thanks to RCall) 

	R"set.seed(12)"
	m = @vam(System & Time & Type ~ (ARAInf(.4) | Weibull(.001,2.5)))
	m.rand = randr
	dfJl = rand(m, 10; system = 5)
	## 2) VAM simulation
	R"""
	require(VAM)
	set.seed(12)
	dfR = simulate(sim.vam( ~ (ARAInf(.4) | Weibull(.001,2.5))), 10, nb.syst=5)
	"""
	@rget dfR
	## 3) comparison

	dfR == dfJl
end

simulcompare_multi()

### Rest of examples
m = @vam(Time & Type ~ (ARAInf(.4) | Weibull(.001,2.5)))
rand(m)

@stop (size < 30)

s = Simulator(
	@vam(Temps & Type ~ (ARAInf(.4) | Weibull(.001,2.5))),
	@stop (size < 30)
)

df = rand(s)

s = Simulator(
	@vam(Temps & Type ~ (ARA1(.5) | Weibull(0.01,2.5)) & (ARAInf(.7)+ARAInf(.3)+ AGAN()|Periodic(12,[0.6,0.4]) * AtIntensity(0.2))),
  	@stop (size < 30) & (time < 30)
)

s = Simulator(
	#@vam(Temps & Type ~ (ARA1(.5) | Weibull(0.01,2.5)) & (ARAInf(.7)+ARAInf(.3)+ ABAO()|Periodic(12,[0.6,0.4]) * AtIntensity(0.2))),
	@vam(Temps & Type ~ (ABAO() | Weibull(.001,2.7)) & (ARAInf(.6)+ARAInf(-.2)+AGAN() | Periodic(12,[0.6,0.4]) * AtIntensity(0.6))),
	@stop (size < 30) & (time < 30)
)



df = rand(s)

m = @vam(Temps & Type ~ (ARA1(.5) | Weibull(0.01,2.5)) & (ARAInf(.7)+ARAInf(.3)+ ABAO()|Periodic(12,[0.6,0.4]) * AtIntensity(.2)))
df = rand(m, @stop(size < 30))
df = rand(m,30)

params(s.model)

df2 = rand(s, system=30)

df = rand(
    @vam(Temps & Type ~ (ARA1(.5) | Weibull(0.01,2.5)) & (ARAInf(.7)+ARAInf(.3)+ ABAO()|Periodic(12,[0.6,0.4]) * AtIntensity(1.2))),
  	@stop (size < 30) & (time < 30)
)

vam = @vam(Temps & Type ~ (ARA1(.5) | Weibull(0.01,2.5)) & (ARAInf(.7)+ARAInf(.3)+ ABAO()|Periodic(12,[0.6,0.4]) * AtIntensity(1.2)))
rand(vam, @stop (size < 30) & (time < 30))
rand(vam, @stop(size < 30), system=30)
