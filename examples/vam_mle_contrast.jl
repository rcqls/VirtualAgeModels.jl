using VirtualAgeModels
using DataFrames
using RCall

function testLD()

	(((contrast(m, θ, profile=false)-lnL).^2).^0.5,
	maximum(((gradient(m, θ, profile=false)-dlnL).^2).^0.5),
	maximum(((hessian(m, θ, profile=false)-d2lnL).^2).^0.5),
	maximum(((contrast(m, θ)-C).^2).^0.5),
	maximum(((gradient(m, θ)[2:end]-dC).^2).^0.5),
	maximum(((hessian(m, θ)[2:end,2:end]-d2C).^2).^0.5))
end


df = DataFrame(time=[3.36],type=[-1])
m = @vam time & type ~ (ARAInf(0.4) | Weibull(0.001,2.5)) data=df
θ = [0.3,1.8,0.6]
lnL = -2.30449245951301
dlnL = [-5.52619699756427,-1.45367181592636,0]
d2lnL = [
	-11.1111111111111 -10.7372278181901 0;
	-10.7372278181901 -4.21250787723964 0;
	0 0 0
]
contrast(m, θ, profile=false)
gradient(m, θ, profile=false)
hessian(m, θ, profile=false)
C = -1.62415430907299
dC = [0.555555555555556,0]
d2C = [-0.308641975308642 0; 0 0]
contrast(m,θ)
gradient(m, θ)
hessian(m, θ)
# @run contrast(m,θ)
# @run gradient(m, θ)
# @run hessian(m, θ)
maximum(testLD())

df = DataFrame(time=[3.36],type=[-1])
m = @vam Time & Type ~ (ARAInf(0.4) | LogLinear(0.001,2.5)) data=df
θ = [0.3,0.8,0.6]
lnL = -3.65431355894635
dlnL = [-13.7944691820681, -8.74189899224908, 0]
d2lnL = [-11.1111111111111 -40.3396633074969 0;
-40.3396633074969 -31.9886643027399  0;
0 0 0]
contrast(m, θ, profile=false)
gradient(m, θ, profile=false)
hessian(m, θ, profile=false)
C = -1.15270302129339
dC = [1.00478465516967,0]
d2C = reshape([-0.678445876029819,0,0,0], 2, 2)
contrast(m,θ)
gradient(m, θ)
hessian(m, θ)
maximum(testLD())

df = DataFrame(time=[3.36],type=[-1])
m =  @vam Time & Type ~ (ARAInf(0.4) | Weibull3(0.001,2.5,5.0)) data=df
θ = [0.3,1.8,4,0.6]
lnL =  -6.28349650594271
dlnL = [-20.8805277090384,-14.1662361334174,-0.920550636729322,0]
d2lnL = hcat([-11.1111111111111,-55.726172072379,-3.43082096301078,0],[-55.726172072379,-36.7534932861936,-3.48854152770058,0],[-3.43082096301078,-3.48854152770058,0.0228198059801746,0],[0,0,0,0])
contrast(m, θ, profile=false)
gradient(m, θ, profile=false)
hessian(m, θ, profile=false)
C = -2.00229062844351
dC = [0.250199288093325,-0.0329926542472398,0]
d2C = hcat([-0.0292037862530474,-0.0369910685383343,0],[-0.0369910685383343,0.01048162449677,0],[0,0,0])
contrast(m,θ)
gradient(m, θ)
hessian(m, θ)
maximum(testLD())



df = DataFrame(time=[3.36,4.04,4.97,5.16], type=[-1,-1,-1,-1])
m = @vam Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)) data = df
θ = [0.3,0.8,0.6]
lnL = -7.37830963135462
dlnL = [9.33348796771948,5.77076155284033,1.19923836457015]
d2lnL = hcat([-44.4444444444444,-5.26230480903023,-0.723044448779292],[-5.26230480903023,-7.7471885008781,-6.30726171420585],[-0.723044448779292,-6.30726171420585,0.435684125802398])
contrast(m, θ, profile=false)
gradient(m, θ, profile=false)
hessian(m, θ, profile=false)
C = -5.36231016699152
dC = [2.08694474533596,0.693079297460398]
d2C = hcat([-4.31732300351524,-3.55104259186756],[-3.55104259186756,-0.351517905180765])
contrast(m,θ)
gradient(m, θ)
hessian(m, θ)
testLD()
maximum(testLD())

df = DataFrame(time=[3.36,4.04,4.97,5.16], type=[-1,-1,-1,-1])
m = @vam Time & Type ~ (ARA1(0.4) | Weibull(0.001,2.5)) data=df 
m = @vam Time & Type ~ (ARAm(0.4,1) | Weibull(0.001,2.5)) data=df 
θ = [0.3,0.8,0.6]
lnL = -7.60531410020218
dlnL = [9.43046924122556,6.95299388770745,0.950641033424141]
d2lnL = hcat([-44.4444444444444,-5.58984454112956,-0.51705781427391],[-5.58984454112956,-8.18607969766148,-5.04449078731948],[-0.51705781427391,-5.04449078731948,1.70955451742894])
contrast(m, θ, profile=false)
gradient(m, θ, profile=false)
hessian(m, θ, profile=false)
C = -5.5202288755881
dC = [2.90098061507367,0.575831839583211]
d2C = hcat([-4.65895384634265,-3.11529381117442],[-3.11529381117442,0.845549840945982])
contrast(m,θ)
gradient(m, θ)
hessian(m, θ)
testLD()
maximum(testLD())

df = DataFrame(time=[3.36,4.04,4.97,5.16], type=[-1,-1,-1,-1])
m = @vam Time & Type ~ (AGAN() | Weibull(0.001,2.5)) data=df 
θ = [0.3,0.8]
lnL = -6.90098251895707
dlnL = [8.75359407434948,3.3717767832448]
d2lnL = hcat([-44.4444444444444,-2.40399936753948],[-2.40399936753948,-7.6652713149061])
contrast(m, θ, profile=false)
gradient(m, θ, profile=false)
hessian(m, θ, profile=false)
C = -5.25256034408912
dC = [1.99329429144225]
d2C = hcat(-9.26821752088532)
contrast(m,θ)
gradient(m, θ)
hessian(m, θ)
maximum(testLD())


df = DataFrame(time=[3.36,4.04,4.97,5.16], type=[-1,-1,-1,-1])
m = @vam Time & Type ~ (ABAO() | LogLinear(0.001,2.5)) data=df
θ = [0.3,0.8]
lnL = -13.6870254780068
dlnL = [-62.9837808690102,-73.9249749593491]
d2lnL = hcat([-44.4444444444444,-304.849916531164],[-304.849916531164,-390.943849373404])
contrast(m, θ, profile=false)
gradient(m, θ, profile=false)
hessian(m, θ, profile=false)
C = -1.77041141396439
dC = [1.55193690275558]
d2C = hcat(-4.47702275378921)
contrast(m,θ)
gradient(m, θ)
hessian(m, θ)
maximum(testLD())

df = DataFrame(
	time=[18.09,52.07,95.71,145.75,198.7,220.9],
	type=[-1,-1,-1,-1,-1,-1]
)
m = @vam Time & Type ~ (ARAm(0.5,2) | Weibull(0.001,2.5)) data=df
θ = [0.03,2.4,0.7]
lnL = -2600.46156675627
dlnL = [-87032.7099834185,-10989.6661601188,6745.47891744838]
d2lnL = [
	-6666.66666666667 -367174.735520771 225153.425274492;
	-367174.735520771 -46307.9946188132 32956.7724964216;
	225153.425274492 32956.7724964216 -16880.4479240297
]
contrast(m, θ, profile=false)
gradient(m, θ, profile=false)
hessian(m, θ, profile=false)
C = -25.948370978711
dC = [0.32106229767103,6.36255820206048]
d2C = [
	-0.909264674311954 3.87406673594902;
	3.87406673594902 -0.394244003389716
]
contrast(m,θ)
#@run 
gradient(m, θ)
hessian(m, θ)
testLD()
maximum(testLD())

df = DataFrame(
	system=vcat([repeat([x],i) for (x, i) in [(1,4), (2,4), (3,4), (4,4), (5,7)]]...),
	time=[3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98],
	type=[1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0]
)
m = @vam System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)) & (AGAN()) data=df
m.nb_system
θ = [0.3,1.8,0.3]
lnL = -20.373593780469
dlnL = [-25.6662791768442,-8.48131134572237,5.63552763127959]
d2lnL = [
	-100 -52.7375595999022 27.2348357573035;
	-52.7375595999022 -19.8545917834264 17.1554852531126;
	27.2348357573035 17.1554852531126 -5.23495386861318
	]
contrast(m, θ, profile=false)
gradient(m, θ, profile=false)
hessian(m, θ, profile=false)
C = -18.237304657921
dC = [-1.18653460936235,1.86834447361728]
d2C = [
	-3.90302067583932 3.61294880573467;
	3.61294880573467 -1.23797415150246
	]
contrast(m,θ)
gradient(m, θ)
hessian(m, θ)
testLD()
maximum(testLD())

df =  DataFrame(
	System=vcat(repeat(1:4,inner=4),repeat([5],7)),
	Time=[3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98],
	Type=[1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0]
)
m = @vam System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)) & (ARAInf(0.5)) data=df
m.nb_system = 5
θ = [0.3,1.8,0.3,0.8]
lnL = -20.7592907965932
dlnL = [-31.3904767850125,-9.15468712792062,6.37367814316912,3.41532811689743]
d2lnL = [
	-100 -60.5303216502559 30.8177749831935 29.6544563881339;
	-60.5303216502559 -20.7750026445059 19.227056494657 8.24534883396676;
	30.8177749831935 19.227056494657 -6.01682216912879 -4.05851834664032;
	29.6544563881339 8.24534883396676 -4.05851834664032 -11.086407735654
		]
contrast(m, θ, profile=false)
gradient(m, θ, profile=false)
hessian(m, θ, profile=false)
C = -17.7866638232513
dC = [0.130510297108522,1.64630982399109,-1.13359008787765]
d2C = [
	-2.82300923948606 3.10568724954822 -3.76042913146761;
	3.10568724954822 -1.35443592315717 1.21888246564082;
	-3.76042913146761 1.21888246564082 -6.87069738841726
	]
contrast(m,θ)
gradient(m, θ)
hessian(m, θ)
testLD()
maximum(testLD())


df =  DataFrame(
	Time=[18.09,52.07,95.71,145.75],
	Type=[-1,-1,-1,-1]
)
m = @vam Time & Type ~ (GQR_ARA1(0.7,-1.2|log) | Weibull(0.001,2.5)) data=df
θ =  [0.03,2.4,0.7,-1.2]
lnL = -5227.28210215578
dlnL = [-174523.718461965,-27932.4679482722,-15588.9239900853,2925.83535839578]
d2lnL = [
	-4444.44444444444 -931711.328605459 -519891.725970003 97574.7515774708;
	-931711.328605459 -149124.678752068 -89338.3200091043 17794.0729848851;
	-519891.725970003 -89338.3200091043 -29458.0754936798 8538.61396441627;
	97574.7515774708 17794.0729848851 8538.61396441627 -474.24229613857
		]
contrast(m, θ, profile=false)
gradient(m, θ, profile=false)
hessian(m, θ, profile=false)
C = -20.2814577676388
dC = [-2.46616536000733,-4.07878379460467,0.827470378552931]
d2C = [
	-0.707449731980772 -2.33853540176398 0.658809838099249;
	-2.33853540176398 3.22719216928062 0.213056341848732;
	0.658809838099249 0.213056341848732 0.406344931337982
	]
contrast(m,θ)
gradient(m, θ)
hessian(m, θ)
testLD()
maximum(testLD())

df =  DataFrame(
	Time=[18.09,52.07,95.71,145.75],
	Type=[-1,-1,-1,-1]
)
m = @vam Time & Type ~ (GQR(0.7|sqrt) | Weibull(0.001,2.5)) data=df
θ =  [0.03,2.4,0.7]
lnL = -244.558483645129
dlnL = [-8207.95806347584,-788.351238674561,-1050.82462756074]
d2lnL = [
	-4444.44444444444 -26754.9654997246 -35501.3463704942;
	-26754.9654997246 -2579.29345017032 -3912.38114820971;
	-35501.3463704942 -3912.38114820971 -3897.23540002875
		]
contrast(m, θ, profile=false)
gradient(m, θ, profile=false)
hessian(m, θ, profile=false)
C = -14.8642260155259
dC = [1.46759523701393,-2.80862496401469]
d2C = [
	-0.759600585249061 -2.10353871070821;
	-2.10353871070821 -9.8224362844605
	]
contrast(m,θ)
gradient(m, θ)
hessian(m, θ)
testLD()
maximum(testLD())

df =  DataFrame(
	Time=[18.09,52.07,95.71,145.75,198.7,220.9,230],
	Type=[-1,-1,-1,-1,-1,-1,0]
)
m = @vam Time & Type ~ (GQR_ARAm(1.2,0.5|3) | Weibull(0.001,2.5)) data=df
θ =  [0.03,2.4,1.3,0.7]
lnL = -19524.0996600096
dlnL = [-651440.845101926,-98931.83785181,-130451.573242292,44948.9374958108]
d2lnL = [
	-6666.66666666667 -3298691.66087083 -4349260.39370564 1498505.11901772;
	-3298691.66087083 -501269.384415686 -717595.854835935 254816.142137408;
	-4349260.39370564 -717595.854835935 -837390.260758198 292604.491724178;
	1498505.11901772 254816.142137408 292604.491724178 -115028.435473886
		]
contrast(m, θ, profile=false)
gradient(m, θ, profile=false)
hessian(m, θ, profile=false)
C = -29.4078957809114
dC = [-1.46081479729706,-13.8073591192265,7.58145006849591]
d2C = [
	-1.13865627103432 -7.02814050027214 3.92396572758071;
	-7.02814050027214 -8.65037145925458 2.91612651943195;
	3.92396572758071 2.91612651943195 -2.45231214454718
	]
contrast(m,θ)
gradient(m, θ)
hessian(m, θ)
testLD()
maximum(testLD())

df =  DataFrame(
	System=vcat(repeat(1:2,inner=4),repeat([3],14)),
	Time=[3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78],
	Type=[1,1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,1,1,-1,1,-1,1,-1,0]
)
m = @vam System & Time & Type ~ (ARAInf(0.3) | Weibull(0.001,2.5)) & (ARAm(0.6|3)) data=df
m2 = @vam System & Time & Type ~ (ARAm(0.3|22) | Weibull(0.001,2.5)) & (ARAm(0.6|3)) data=df
θ =  [0.03,2.4,0.3,0.7]
contrast(m, θ, profile=false)-contrast(m2, θ, profile=false)
gradient(m, θ, profile=false)-gradient(m2, θ, profile=false)
hessian(m, θ, profile=false)-hessian(m2, θ, profile=false)
contrast(m,θ)-contrast(m2,θ)
gradient(m, θ)-gradient(m2, θ)
hessian(m, θ)-hessian(m2, θ)


df =  DataFrame(
	System=vcat(repeat(1:2,inner=4),repeat([3],14)),
	Time=[3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78],
	Type=[1,1,1,1, -1,1,1,0, 1,-1,1,-1,-1,1,1,1,-1,1,-1,1,-1,0]
)
m = @vam System & Time & Type ~ (GQR(1.) | Weibull(0.001,2.5)) & (GQR_ARAm(0.9,0.6|3,log)) data=df
θ =  [0.03,2.4,1.1,0.9,0.7]
R"""
	require(VAM)
	dataDF <- data.frame(System=c(rep(1,4),rep(2,4),rep(3,14)),Time=c(3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78),Type=c(1,1,1,1, -1,1,1,0, 1,-1,1,-1,-1,1,1,1,-1,1,-1,1,-1,0),row.names=1:22)
    theta<-c(0.03,2.4,1.1,0.9,0.7)
	mle <- mle.vam(System & Time & Type ~ (GQR(1.) | Weibull(0.001,2.5)) & (GQR_ARAm(0.9,0.6|3,log)),data=dataDF)
	res <- list()
	res$lnL <- logLik(mle,theta,TRUE,FALSE,FALSE)
	res$dlnL<- logLik(mle,theta,FALSE,TRUE,FALSE)
	res$d2lnL <- logLik(mle,theta,FALSE,FALSE,TRUE)
	res$C <- contrast(mle,theta,TRUE,FALSE,FALSE)
	res$dC <- contrast(mle,theta,FALSE,TRUE,FALSE)
	res$d2C <- contrast(mle,theta,FALSE,FALSE,TRUE)
	"""
	res = @rget res
lnL = res[:lnL]
dlnL = res[:dlnL]
d2lnL = res[:d2lnL]
contrast(m, θ, profile=false)-lnL
gradient(m, θ, profile=false)-dlnL
hessian(m, θ, profile=false)-d2lnL
C = res[:C]
dC = res[:dC]
d2C = res[:d2C]
contrast(m,θ)-C
gradient(m, θ)[2:end]-dC
hessian(m, θ)[2:end,2:end]-d2C
testLD()
maximum(testLD())


df = DataFrame(System=[1,1,1,1,2,2,2,3],Time=[3.36,4.04,4.97,5.16,2.34,3.46,5.02,4],Type=[-1,-1,-1,0,-1,-1,-1,0])
m = @vam System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)) data = df
θ = [0.3,0.8,0.6]
R"""
	require(VAM)
	dataDF <- data.frame(System=c(1,1,1,1,2,2,2,3),Time=c(3.36,4.04,4.97,5.16,2.34,3.46,5.02,4),Type=c(-1,-1,-1,0,-1,-1,-1,0))
    theta<-c(0.3,0.8,0.6)
	mle <- mle.vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)),data=dataDF)
	res <- list()
	res$lnL <- logLik(mle,theta,TRUE,FALSE,FALSE)
	res$dlnL<- logLik(mle,theta,FALSE,TRUE,FALSE)
	res$d2lnL <- logLik(mle,theta,FALSE,FALSE,TRUE)
	res$C <- contrast(mle,theta,TRUE,FALSE,FALSE)
	res$dC <- contrast(mle,theta,FALSE,TRUE,FALSE)
	res$d2C <- contrast(mle,theta,FALSE,FALSE,TRUE)
	"""
	res = @rget res
lnL = res[:lnL]
dlnL = res[:dlnL]
d2lnL = res[:d2lnL]
(contrast(m, θ, profile=false)-lnL)/lnL
maximum(abs.((gradient(m, θ, profile=false)-dlnL)/dlnL))<0.00001
(hessian(m, θ, profile=false)-d2lnL)/d2lnL
C = res[:C]
dC = res[:dC]
d2C = res[:d2C]
(contrast(m,θ)-C)/C
(gradient(m, θ)[2:end]-dC)/dC
(hessian(m, θ)[2:end,2:end]-d2C)/d2C


