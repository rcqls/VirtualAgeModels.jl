using VAM
using DataFrames

s = simulator(
	@vam(Temps & Type ~ (ARA1(.5) | Weibull(0.01,2.5)) & (ARAInf(.7)+ARAInf(.3)+ ABAO()|Periodic(12,[0.6,0.4]) * AtIntensity(0.2))),
  	@stop (size < 30) & (time < 30)
)

df = simulate(s)

m = @vam(Temps & Type ~ (ARA1(.5) | Weibull(0.01,2.5)) & (ARAInf(.7)+ARAInf(.3)+ ABAO()|Periodic(12,[0.6,0.4]) * AtIntensity(.2)))
df = simulate(m, @stop(size < 30))


params(s.model)

df2 = simulate(s,system=30)

df = simulate(
    @vam(Temps & Type ~ (ARA1(.5) | Weibull(0.01,2.5)) & (ARAInf(.7)+ARAInf(.3)+ ABAO()|Periodic(12,[0.6,0.4]) * AtIntensity(1.2))),
  	@stop (size < 30) & (time < 30)
)

vam = @vam(Temps & Type ~ (ARA1(.5) | Weibull(0.01,2.5)) & (ARAInf(.7)+ARAInf(.3)+ ABAO()|Periodic(12,[0.6,0.4]) * AtIntensity(1.2)))
simulate(vam, @stop (size < 30) & (time < 30))
simulate(vam, @stop(size < 30), system=30)
