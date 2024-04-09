using Test
using VirtualAgeModels
using DataFrames
using RCall

include("testing.jl")

modtest = ModelTest()

insert!(modtest, 

    :W => Dict(
        :θ =>  [0.3,1.8,0.6],
        :data => DataFrame(Temps=[3.361],Type=[-1]),
        :vam => @vam(Temps & Type ~ (ARA∞(0.4) | Weibull(0.001,2.5)))
    ),
    :LL => Dict(
        :θ =>  [0.3,1.8,0.6],
        :data => DataFrame(Temps=[3.361],Type=[-1]),
        :vam => @vam(Temps & Type ~ (ARA∞(0.4) | LogLinear(0.001,2.5)))
    ),

    :W_AGAN => Dict(
		:θ =>  [0.3,0.8],
		:data => DataFrame(Temps=[3.36,4.04,4.97,5.16],Type=[-1,-1,-1, -1]),
		:vam => @vam(Temps & Type ~ (AGAN() | Weibull(0.001,2.5))) #,
	),
    :LL_ABAO => Dict(
		:θ =>  [0.3,0.8],
		:data => DataFrame(Temps=[3.36,4.04,4.97,5.16],Type=[-1,-1,-1, -1]),
		:vam => @vam(Temps & Type ~ (ABAO() | LogLinear(0.001,2.5))) #,
	),


	:W_ARA∞ => Dict(
		:θ =>  [0.3,1.8,0.6],
		:data => DataFrame(Temps=[3.36, 4.1],Type=[-1, 0]),
		:vam => @vam(Temps & Type ~ (ARA∞(0.4) | Weibull(0.001,2.5))) #,
		# :r => [
		# 	"data.frame(Time=c(3.36, 4.1),Type=c(-1, 0))",
		# 	"Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5))"
		# 	]
	),
	:W_ARA∞_ARA∞ => Dict(
		:θ =>  [0.3,1.8,0.3,0.8],
		:data => DataFrame(Temps=[3.36, 4.1, 6.1, 7.2],Type=[-1, 1, 1, -1]),
		:vam => @vam(Temps & Type ~ (ARA∞(0.4) | Weibull(0.001,2.5)) & (ARA∞(0.5))) #,
		# :r => [
		# 	"data.frame(Time=c(3.36, 4.1),Type=c(-1, 0))",
		# 	"Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5))"
		# 	]
	), 
	:W_ARA1 => Dict(
		:θ =>  [0.3,1.8,0.6],
		:data => DataFrame(Temps=[3.36, 4.1],Type=[-1, 0]),
		:vam => @vam(Temps & Type ~ (ARA1(0.4) | Weibull(0.001,2.5))) #,
		# :r => [
		# 	"data.frame(Time=c(3.36, 4.1),Type=c(-1, 0),row.names=1:2)",
		# 	"Time & Type ~ (ARA1(0.4) | Weibull(0.001,2.5))"
		# ] 
	),
    :W_ARA1_Cens => Dict(
		:θ =>  [0.3,1.8,0.6],
		:data => DataFrame(Temps=[3.36,4.04,4.97,5.16],Type=[-1,-1,-1, 0]),
		:vam => @vam(Temps & Type ~ (ARA1(0.4) | Weibull(0.001,2.5))) #,
		# :r => [
		# 	"data.frame(Time=c(3.36, 4.1),Type=c(-1, 0),row.names=1:2)",
		# 	"Time & Type ~ (ARA1(0.4) | Weibull(0.001,2.5))"
		# ] 
	),
    #TODO: Left censored
    #:W_ABAO_LeftCens => Dict(
	#	:θ =>  [0.3,2.5],
	#	:data => DataFrame(Temps=[2,3.36,4.04,4.97,5.16],Type=[0,-1,-1,-1, -1]),
	#	:vam => @vam(Temps & Type ~ (ABAO() | Weibull(0.001,2.5))) 
	#),
	:W_ARA∞bis => Dict(
		:θ => [0.3,0.8,0.6],
		:data => DataFrame(Time=[3.36,4.04,4.97,5.16],Type=[-1,-1,-1,-1]),
		:vam => @vam(Time & Type ~ (ARA∞(0.4) | Weibull(0.001,2.5)))
	),
	:W_ARA1bis => Dict(
		:θ => [0.3,0.8,0.6],
		:data => DataFrame(Time=[3.36,4.04,4.97,5.16],Type=[-1,-1,-1,-1]),
		:vam => @vam(Time & Type ~ (ARA1(0.4) | Weibull(0.001,2.5)))
	),
    :W_ARAm2 => Dict(
        :θ => [0.03,2.4,0.7],
        :data => DataFrame(Time=[18.09,52.07,95.71,145.75,198.7,220.9],Type=[-1,-1,-1,-1,-1,-1]),
        :vam => @vam(Time & Type ~ (ARAm(0.5|2) | Weibull(0.001,2.5)))
    ),
    :W_ARAm3 => Dict(
        :θ => [0.03,2.4,0.7],
        :data => DataFrame(Time=[18.09,52.07,95.71,145.75,198.7,220.9],Type=[-1,-1,-1,-1,-1,-1]),
        :vam => @vam(Time & Type ~ (ARAm(0.5|3) | Weibull(0.001,2.5)))
    ),
    :W_ARAm4 => Dict(
        :θ => [0.03,2.4,0.7],
        :data => DataFrame(Time=[18.09,52.07,95.71,145.75,198.7,220.9],Type=[-1,-1,-1,-1,-1,-1]),
        :vam => @vam(Time & Type ~ (ARAm(0.5|4) | Weibull(0.001,2.5)))
    ),
    :W_ARAm1 => Dict(
        :θ => [0.03,2.4,0.7],
        :data => DataFrame(Time=[18.09,52.07,95.71,145.75,198.7,220.9],Type=[-1,-1,-1,-1,-1,-1]),
        :vam => @vam(Time & Type ~ (ARAm(0.5|1) | Weibull(0.001,2.5)))
    ),
    :W_ARAm1_vs_ARA1 => Dict(
        :θ => [0.03,2.4,0.7],
        :data => DataFrame(Time=[18.09,52.07,95.71,145.75,198.7,220.9],Type=[-1,-1,-1,-1,-1,-1]),
        :vam => @vam(Time & Type ~ (ARA1(0.5) | Weibull(0.001,2.5)))
    ),
    :W_ARAm∞ => Dict(
        :θ => [0.03,2.4,0.7],
        :data => DataFrame(Time=[18.09,52.07,95.71,145.75,198.7,220.9],Type=[-1,-1,-1,-1,-1,-1]),
        :vam => @vam(Time & Type ~ (ARAm(0.5|8) | Weibull(0.001,2.5)))
    ),
    :W_ARAm∞_vs_ARA∞ => Dict(
        :θ => [0.03,2.4,0.7],
        :data => DataFrame(Time=[18.09,52.07,95.71,145.75,198.7,220.9],Type=[-1,-1,-1,-1,-1,-1]),
        :vam => @vam(Time & Type ~ (ARAInf(0.5) | Weibull(0.001,2.5)))
    ),
    :W_ARAm2_ARAm4_MS => Dict(
        :θ  => [0.3,1.8,0.3,0.7],
        :data => DataFrame(System=vcat(repeat(1:2,inner=4),repeat([3],14)),Time=[3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78],Type=[1,1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,1,1,-1,1,-1,1,-1,0]),
        :vam => @vam(System & Time & Type ~ (ARAm(0.3|2) | Weibull(0.001,2.5)) & (ARAm(0.6|4)))
    ),
    :W_ARAm4_ARAm2_MS => Dict(
        :θ  => [0.3,1.8,0.3,0.7],
        :data => DataFrame(System=vcat(repeat(1:2,inner=4),repeat([3],14)),Time=[3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78],Type=[1,1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,1,1,-1,1,-1,1,-1,0]),
        :vam => @vam(System & Time & Type ~ (ARAm(0.3|4) | Weibull(0.001,2.5)) & (ARAm(0.6|2)))
    ),
    :W_QR_ARAm2_MS => Dict(
        :θ  => [0.3,1.8,1.3,0.7],
        :data => DataFrame(System=vcat(repeat(1:2,inner=4),repeat([3],14)),Time=[3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78],Type=[1,1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,1,1,-1,1,-1,1,-1,0]),
        :vam => @vam(System & Time & Type ~ (QR(0.3) | Weibull(0.001,2.5)) & (ARAm(0.6|2)))
    ),
    :W_ARAInf_ARAm2_MS => Dict(
        :θ  => [0.3,1.8,0.3,0.7],
        :data => DataFrame(System=vcat(repeat(1:2,inner=4),repeat([3],14)),Time=[3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78],Type=[1,1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,1,1,-1,1,-1,1,-1,0]),
        :vam => @vam(System & Time & Type ~ (ARAInf(0.3) | Weibull(0.001,2.5)) & (ARAm(0.6|2)))
    ),
    :W_ARAInf_ARAm3_MS => Dict(
        :θ  => [0.3,1.8,0.3,0.7],
        :data => DataFrame(System=vcat(repeat(1:2,inner=4),repeat([3],14)),Time=[3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78],Type=[1,1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,1,1,-1,1,-1,1,-1,0]),
        :vam => @vam(System & Time & Type ~ (ARAInf(0.3) | Weibull(0.001,2.5)) & (ARAm(0.6|3)))
    ),
    :W_ARA1_ARAm2_MS => Dict(
        :θ  => [0.3,1.8,0.3,0.7],
        :data => DataFrame(System=vcat(repeat(1:2,inner=4),repeat([3],14)),Time=[3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78],Type=[1,1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,1,1,-1,1,-1,1,-1,0]),
        :vam => @vam(System & Time & Type ~ (ARA1(0.3) | Weibull(0.001,2.5)) & (ARAm(0.6|2)))
    ),
    :W_ARA1_ARAm3_MS => Dict(
        :θ  => [0.3,1.8,0.3,0.7],
        :data => DataFrame(System=vcat(repeat(1:2,inner=4),repeat([3],14)),Time=[3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78],Type=[1,1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,1,1,-1,1,-1,1,-1,0]),
        :vam => @vam(System & Time & Type ~ (ARA1(0.3) | Weibull(0.001,2.5)) & (ARAm(0.6|3)))
    ),
    :W_ARAm3_QR_MS => Dict(
        :θ  => [0.3,1.8,0.3,0.7],
        :data => DataFrame(System=vcat(repeat(1:2,inner=4),repeat([3],14)),Time=[3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78],Type=[1,1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,1,1,-1,1,-1,1,-1,0]),
        :vam => @vam(System & Time & Type ~ (ARAm(0.3|3) | Weibull(0.001,2.5)) & (QR(0.6)))
    ),
	:LL_ARA∞ => Dict(
		:θ =>  [0.3,1.8,0.6],
		:data => DataFrame(Temps=[3.36, 4.1],Type=[-1, 0]),
		:vam => @vam(Temps & Type ~ (ARA∞(0.4) | LogLinear(0.001,2.5))) #,
		# :r => [
		# 	"data.frame(Time=c(3.36, 4.1),Type=c(-1, 0))",
		# 	"Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5))"
		# 	]
	), 
	:LL_ARA1 => Dict(
		:θ =>  [0.3,1.8,0.6],
		:data => DataFrame(Temps=[3.36, 4.1],Type=[-1, 0]),
		:vam => @vam(Temps & Type ~ (ARA1(0.4) | LogLinear(0.001,2.5))) #,
		# :r => [
		# 	"data.frame(Time=c(3.36, 4.1),Type=c(-1, 0),row.names=1:2)",
		# 	"Time & Type ~ (ARA1(0.4) | Weibull(0.001,2.5))"
		# ] 
	),
	:LL_ARA∞bis => Dict(
		:θ => [0.3,0.8,0.6],
		:data => DataFrame(Time=[3.36,4.04,4.97,5.16],Type=[-1,-1,-1,-1]),
		:vam => @vam(Time & Type ~ (ARA∞(0.4) | LogLinear(0.001,2.5)))
	),
	:LL_ARA1bis => Dict(
		:θ => [0.3,0.8,0.6],
		:data => DataFrame(Time=[3.36,4.04,4.97,5.16],Type=[-1,-1,-1,-1]),
		:vam => @vam(Time & Type ~ (ARA1(0.4) | LogLinear(0.001,2.5)))
	),
	:W_ARA∞_ARA∞_ARA1 => Dict(
		:θ => [0.3,1.8,0.3,0.4,0.7],
		:data => DataFrame(Time=[3.36,4.04,4.97,5.16,7.16],Type=[1,2,2,1,-1]),
		:vam => @vam(Time & Type ~ (ARA∞(0.4) | Weibull(0.001,2.5)) & (ARA∞(0.5)+ARA1(0.7)))
	),
	:W_ARA∞_AGAN_ARA1 => Dict(
		:θ => [0.3,1.8,0.3,0.7],
		:data => DataFrame(Time=[3.36,4.04,4.97,5.16,7.16],Type=[1,2,2,1,-1]),
		:vam => @vam(Time & Type ~ (ARA∞(0.4) | Weibull(0.001,2.5)) & (AGAN()+ARA1(0.7)))
	),
    :W_MS => Dict(
        :θ => [0.3,0.8,0.6],
        :data => DataFrame(System=[1,2],Time=[3.36,2.34],Type=[-1,-1]),
        :vam => @vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5))) 
    ),
    :W_ARA∞_MS => Dict(
        :θ => [0.3,0.8,0.6],
        :data => DataFrame(System=[1,1,1,1,2,2,2,3],Time=[3.36,4.04,4.97,5.16,2.34,3.46,5.02,4],Type=[-1,-1,-1,0,-1,-1,-1,0]),
        :vam => @vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5))) 
    ),
    :W_ARA∞_ARA∞_MS => Dict(
        :θ => [0.3,1.8,0.3,0.8],
        :data => DataFrame(System=vcat(repeat(1:4,inner=4),repeat([5],7)),Time=[3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98],Type=[1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0]),
        :vam => @vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)) & (ARAInf(0.5))) 
    ),
    :W_ARAInf_AGAN_MS => Dict(
        :θ => [0.3,1.8,0.3],
        :data => DataFrame(System=vcat(repeat(1:4,inner=4),repeat([5],7)),Time=[3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98],Type=[1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0]),
        :vam => @vam(System & Time & Type ~ (ARA1(0.4) | Weibull(0.001,2.5)) & (AGAN())) 
    ),
    :W_ARAInf_ABAO_MS => Dict(
        :θ => [0.3,1.8,0.3],
        :data => DataFrame(System=vcat(repeat(1:4,inner=4),repeat([5],7)),Time=[3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98],Type=[1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0]),
        :vam => @vam(System & Time & Type ~ (ARA1(0.4) | Weibull(0.001,2.5)) & (ABAO())) 
    ),
    :W_AGAN_ARAInf_MS => Dict(
        :θ => [0.3,1.8,0.3],
        :data => DataFrame(System=vcat(repeat(1:4,inner=4),repeat([5],7)),Time=[3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98],Type=[1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0]),
        :vam => @vam(System & Time & Type ~ (AGAN() | Weibull(0.001,2.5)) & (ARAInf(0.5))) 
    ),
    :W_ABAO_ARAInf_MS => Dict(
        :θ => [0.3,1.8,0.3],
        :data => DataFrame(System=vcat(repeat(1:4,inner=4),repeat([5],7)),Time=[3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98],Type=[1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0]),
        :vam => @vam(System & Time & Type ~ (ABAO() | Weibull(0.001,2.5)) & (ARAInf(0.5))) 
    ),
    :W_ARA1_ARA1_MS => Dict(
        :θ => [0.3,1.8,0.3,0.8],
        :data => DataFrame(System=vcat(repeat(1:4,inner=4),repeat([5],7)),Time=[3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98],Type=[1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0]),
        :vam => @vam(System & Time & Type ~ (ARA1(0.4) | Weibull(0.001,2.5)) & (ARA1(0.5))) 
    ),
    :W_ARA1_ARA∞_MS => Dict(
        :θ => [0.3,1.8,0.3,0.8],
        :data => DataFrame(System=vcat(repeat(1:4,inner=4),repeat([5],7)),Time=[3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98],Type=[1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0]),
        :vam => @vam(System & Time & Type ~ (ARA1(0.4) | Weibull(0.001,2.5)) & (ARAInf(0.5))) 
    ),
    :W_ARA∞_ARA1_MS => Dict(
        :θ => [0.3,1.8,0.3,0.8],
        :data => DataFrame(System=vcat(repeat(1:4,inner=4),repeat([5],7)),Time=[3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98],Type=[1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0]),
        :vam => @vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)) & (ARA1(0.5))) 
    ),
    :W_ARA∞_QR_MS => Dict(
        :θ => [0.3,2.2,0.3,0.4],
        :data => DataFrame(System=vcat(repeat(1:4,inner=4),repeat([5],7)),Time=[3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98],Type=[1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0]),
        :vam => @vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)) & (QR(0.7))) 
    ),
    :W_ARA∞_GQR_MS => Dict(
        :θ => [0.3,2.2,0.3,0.4],
        :data => DataFrame(System=vcat(repeat(1:4,inner=4),repeat([5],7)),Time=[3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98],Type=[1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0]),
        :vam => @vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)) & (GQR(0.7))) 
    ),
    :W_QR_ARA∞_MS => Dict(
        :θ => [0.3,2.2,0.7,0.4],
        :data => DataFrame(System=vcat(repeat(1:4,inner=4),repeat([5],7)),Time=[3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98],Type=[1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0]),
        :vam => @vam(System & Time & Type ~ (QR(0.2) | Weibull(0.001,2.5)) & (ARAInf(0.7))) 
    ),
    :W_GQR_ARA∞_MS => Dict(
        :θ => [0.3,2.2,0.7,0.4],
        :data => DataFrame(System=vcat(repeat(1:4,inner=4),repeat([5],7)),Time=[3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98],Type=[1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0]),
        :vam => @vam(System & Time & Type ~ (GQR(0.2) | Weibull(0.001,2.5)) & (ARAInf(0.7))) 
    ),
    :W_ARA1_QR_MS => Dict(
        :θ => [0.3,2.2,0.3,0.4],
        :data => DataFrame(System=vcat(repeat(1:4,inner=4),repeat([5],7)),Time=[3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98],Type=[1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0]),
        :vam => @vam(System & Time & Type ~ (ARA1(0.4) | Weibull(0.001,2.5)) & (QR(0.7))) 
    ),
    :W_ARA1_GQR_MS => Dict(
        :θ => [0.3,2.2,0.3,0.4],
        :data => DataFrame(System=vcat(repeat(1:4,inner=4),repeat([5],7)),Time=[3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98],Type=[1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0]),
        :vam => @vam(System & Time & Type ~ (ARA1(0.4) | Weibull(0.001,2.5)) & (GQR(0.7))) 
    ),
    :W_QR_ARA1_MS => Dict(
        :θ => [0.3,2.2,0.7,0.4],
        :data => DataFrame(System=vcat(repeat(1:4,inner=4),repeat([5],7)),Time=[3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98],Type=[1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0]),
        :vam => @vam(System & Time & Type ~ (QR(0.2) | Weibull(0.001,2.5)) & (ARA1(0.7))) 
    ),
    :W_GQR_ARA1_MS => Dict(
        :θ => [0.3,2.2,0.7,0.4],
        :data => DataFrame(System=vcat(repeat(1:4,inner=4),repeat([5],7)),Time=[3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98],Type=[1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0]),
        :vam => @vam(System & Time & Type ~ (GQR(0.2) | Weibull(0.001,2.5)) & (ARA1(0.7))) 
    ),
    :W_ARA1_GQRARA1log_MS => Dict(
        :θ => [0.3,2.2,0.7,0.4,-0.9],
        :data => DataFrame(System=vcat(repeat(1:4,inner=4),repeat([5],7)),Time=[3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98],Type=[1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0]),
        :vam => @vam(System & Time & Type ~ (ARA1(0.2) | Weibull(0.001,2.5)) & (GQR_ARA1(0.7,-1.3|log))) 
    ),
    :W_AGAN_GQRARA1log_MS => Dict(
        :θ => [0.3,2.2,0.4,-0.9],
        :data => DataFrame(System=vcat(repeat(1:4,inner=4),repeat([5],7)),Time=[3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98],Type=[1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0]),
        :vam => @vam(System & Time & Type ~ (AGAN() | Weibull(0.001,2.5)) & (GQR_ARA1(0.7,-1.3|log))) 
    ),
    :W_ARA1_GQRARA∞log_MS => Dict(
        :θ => [0.3,2.2,0.7,0.4,-0.9],
        :data => DataFrame(System=vcat(repeat(1:4,inner=4),repeat([5],7)),Time=[3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98],Type=[1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0]),
        :vam => @vam(System & Time & Type ~ (ARA1(0.2) | Weibull(0.001,2.5)) & (GQR_ARAInf(0.7,-1.3|log))) 
    ),
    :W_ARA∞_ARA∞_ARA1_MS => Dict(
        :θ => [0.3,1.8,0.3,0.8,-1],
        :data => DataFrame(System=vcat(repeat(1:4,inner=4),repeat([5],10)),Time=[3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43],Type=[1,2,2,1, -1,-1,-1,0, 1,-1,-1,2, -1,2,2,0, 1,-1,1,-1,-1,2,1,2,-1,0]),
        :vam => @vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)) & (ARAInf(0.5)+ARA1(-1.))) 
    ),
    :W_ARA∞_AGAN_ARA1_MS => Dict(
        :θ => [0.3,1.8,0.3,-1],
        :data => DataFrame(System=vcat(repeat(1:4,inner=4),repeat([5],10)),Time=[3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43],Type=[1,2,2,1, -1,-1,-1,0, 1,-1,-1,2, -1,2,2,0, 1,-1,1,-1,-1,2,1,2,-1,0]),
        :vam => @vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)) & (AGAN()+ARA1(-1.))) 
    ),
    :W_ARA∞_AGAN_QR_MS => Dict(
        :θ => [0.3,1.8,0.3,0.7],
        :data => DataFrame(System=vcat(repeat(1:2,inner=4),repeat([3],14)),Time=[3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78],Type=[1,2,2,1, -1,2,2,0, 1,-1,1,-1,-1,2,1,2,-1,1,2,2,1,0]),
        :vam => @vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)) & (AGAN()+QR(0.7))) 
    ),
    :W_ARA∞_ABAO_QR_MS => Dict(
        :θ => [0.3,1.8,0.3,0.7],
        :data => DataFrame(System=vcat(repeat(1:2,inner=4),repeat([3],14)),Time=[3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78],Type=[1,2,2,1, -1,2,2,0, 1,-1,1,-1,-1,2,1,2,-1,1,2,2,1,0]),
        :vam => @vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)) & (ABAO()+QR(0.7))) 
    ),
    :W_ARA∞_AGAP_QR_MS => Dict(
        :θ => [0.3,1.8,0.3,0.7],
        :data => DataFrame(System=vcat(repeat(1:2,inner=4),repeat([3],14)),Time=[3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78],Type=[1,2,2,1, -1,2,2,0, 1,-1,1,-1,-1,2,1,2,-1,1,2,2,1,0]),
        :vam => @vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)) & (AGAP()+QR(0.7))) 
    ),
    :W_ARA∞_QAGAN_QR_MS => Dict(
        :θ => [0.3,1.8,0.3,0.7],
        :data => DataFrame(System=vcat(repeat(1:2,inner=4),repeat([3],14)),Time=[3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78],Type=[1,2,2,1, -1,2,2,0, 1,-1,1,-1,-1,2,1,2,-1,1,2,2,1,0]),
        :vam => @vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)) & (QAGAN()+QR(0.7))) 
    ),
    :W_ARA∞_QR_MS => Dict(
        :θ => [0.3,2.2,0.3,0.4],
        :data => DataFrame(System=vcat(repeat(1:4,inner=4),repeat([5],7)),Time=[3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98],Type=[1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0]),
        :vam => @vam(System & Time & Type ~ (ARAInf(0.2) | Weibull(0.001,2.5)) & (QR(0.7))) 
    ),
    :W_ARA1_GQR_log_MS => Dict(
        :θ => [0.3,2.2,0.7,0.4],
        :data => DataFrame(System=vcat(repeat(1:4,inner=4),repeat([5],7)),Time=[3.36,4.04,4.97,5.16, 2.34,3.46,5.02,5.45, 1.18,2.22,3.14,4.83, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98],Type=[1,1,1,1, -1,-1,-1,0, 1,-1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,0]),
        :vam => @vam(System & Time & Type ~ (ARA1(0.2) | Weibull(0.001,2.5)) & (GQR(0.7|log))) 
    ),
    :W_QR_MS => Dict(
        :θ => [0.03,2.4,0.7],
        :data => DataFrame(System=[1,1,1,1,2,2,2,3],Time=[18.09,52.07,95.71,145.75,15.02,45.1,82,20.1],Type=[-1,-1,-1,0,-1,-1,-1,-1]),
        :vam => @vam(System & Time & Type ~ (QR(0.7) | Weibull(0.001,2.5))) 
    ),
    :W_GQR_MS => Dict(
        :θ => [0.03,2.4,0.7],
        :data => DataFrame(System=[1,1,1,1,2,2,2,3],Time=[18.09,52.07,95.71,145.75,15.02,45.1,82,20.1],Type=[-1,-1,-1,0,-1,-1,-1,-1]),
        :vam => @vam(System & Time & Type ~ (GQR(0.7) | Weibull(0.001,2.5))) 
    ),
    :W_ARA∞_AGAP_QR_MS => Dict(
        :θ => [0.3,1.8,0.3,0.7],
        :data => DataFrame(System=vcat(repeat(1:2,inner=4),repeat([3],14)),Time=[3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78],Type=[1,2,2,1, -1,2,2,0, 1,-1,1,-1,-1,2,1,2,-1,1,2,2,1,0]),
        :vam => @vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)) & (AGAP()+QR(0.7)))
    ),
    :W_ARA∞_QAGAN_QR_MS => Dict(
        :θ => [0.3,1.8,0.3,0.7],
        :data => DataFrame(System=vcat(repeat(1:2,inner=4),repeat([3],14)),Time=[3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78],Type=[1,2,2,1, -1,2,2,0, 1,-1,1,-1,-1,2,1,2,-1,1,2,2,1,0]),
        :vam => @vam(System & Time & Type ~ (ARAInf(0.4) | Weibull(0.001,2.5)) & (QAGAN()+QR(0.7)))
    ),
    :W_QR => Dict(
        :θ => [0.03,2.4,0.7],
        :data => DataFrame(Time=[18.09,52.07,95.71,145.75],Type=[-1,-1,-1,-1]),
        :vam => @vam(Time & Type ~ (QR(0.7) | Weibull(0.001,2.5))) 
    ),
    :W_GQR => Dict(
        :θ => [0.03,2.4,0.7],
        :data => DataFrame(Time=[18.09,52.07,95.71,145.75],Type=[-1,-1,-1,-1]),
        :vam => @vam(Time & Type ~ (GQR(0.7) | Weibull(0.001,2.5))) 
    ),
    :W_GQRLog => Dict(
        :θ => [0.03,2.4,0.7],
        :data => DataFrame(Time=[18.09,52.07,95.71,145.75],Type=[-1,-1,-1,-1]),
        :vam => @vam(Time & Type ~ (GQR(0.7|log) | Weibull(0.001,2.5))) 
    ),
    :W_GQRSqrt => Dict(
        :θ => [0.03,2.4,0.7],
        :data => DataFrame(Time=[18.09,52.07,95.71,145.75],Type=[-1,-1,-1,-1]),
        :vam => @vam(Time & Type ~ (GQR(0.7|sqrt) | Weibull(0.001,2.5))) 
    ),
    :W_GQRARA1Log => Dict(
        :θ => [0.03,2.4,0.7,-1.2],
        :data => DataFrame(Time=[18.09,52.07,95.71,145.75],Type=[-1,-1,-1,-1]),
        :vam => @vam(Time & Type ~ (GQR_ARA1(0.7,-1.2|log) | Weibull(0.001,2.5))) 
    ),
    :W_GQRARA1 => Dict(
        :θ => [0.03,2.4,0.7,-1.2],
        :data => DataFrame(Time=[18.09,52.07,95.71,145.75],Type=[-1,-1,-1,-1]),
        :vam => @vam(Time & Type ~ (GQR_ARA1(0.7,-1.2) | Weibull(0.001,2.5))) 
    ),
    :W_GQRARA∞Log => Dict(
        :θ => [0.03,2.4,0.7,-1.2],
        :data => DataFrame(Time=[18.09,52.07,95.71,145.75],Type=[-1,-1,-1,-1]),
        :vam => @vam(Time & Type ~ (GQR_ARA∞(0.1,0.5|log) | Weibull(0.001,2.5))) 
    ),
    :W_GQRARA∞ => Dict(
        :θ => [0.03,2.4,0.7,-1.2],
        :data => DataFrame(Time=[18.09,52.07,95.71,145.75],Type=[-1,-1,-1,-1]),
        :vam => @vam(Time & Type ~ (GQR_ARA∞(0.1,0.5) | Weibull(0.001,2.5))) 
    ),
    :W_GQRARAm3Log => Dict(
         :θ => [0.03,2.4,1.3,0.7],
         :data => DataFrame(Time=[18.09,52.07,95.71,145.75,198.7,220.9,230],Type=[-1,-1,-1,-1,-1,-1,0]),
         :vam => @vam(Time & Type ~ (GQR_ARAm(1.2,0.5|3,log) | Weibull(0.001,2.5))) 
     ),
    :W_GQRARAm3 => Dict(
          :θ => [0.03,2.4,1.3,0.7],
          :data => DataFrame(Time=[18.09,52.07,95.71,145.75,198.7,220.9,230],Type=[-1,-1,-1,-1,-1,-1,0]),
          :vam => @vam(Time & Type ~ (GQR_ARAm(1.2,0.5|3) | Weibull(0.001,2.5))) 
      ),
    :W_GQRARAm3_ARA1_MS => Dict(
        :θ => [0.03,2.4,1.1,0.3,0.7],
        :data => DataFrame(System=vcat(repeat(1:2,inner=4),repeat([3],14)),Time=[3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78],Type=[1,1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,1,1,-1,1,-1,1,-1,0]),
        :vam => @vam(System & Time & Type ~ (GQR_ARAm(1.,0.3|3) | Weibull(0.001,2.5)) & (ARA1(0.6))) 
    ),
    :W_GQRARAm3log_ARA1_MS => Dict(
        :θ => [0.03,2.4,1.1,0.3,0.7],
        :data => DataFrame(System=vcat(repeat(1:2,inner=4),repeat([3],14)),Time=[3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78],Type=[1,1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,1,1,-1,1,-1,1,-1,0]),
        :vam => @vam(System & Time & Type ~ (GQR_ARAm(1.,0.3|3,log) | Weibull(0.001,2.5)) & (ARA1(0.6))) 
    ),
    :W_GQRARAm3_ARAInf_MS => Dict(
        :θ => [0.03,2.4,1.1,0.3,0.7],
        :data => DataFrame(System=vcat(repeat(1:2,inner=4),repeat([3],14)),Time=[3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78],Type=[1,1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,1,1,-1,1,-1,1,-1,0]),
        :vam => @vam(System & Time & Type ~ (GQR_ARAm(1.,0.3|3) | Weibull(0.001,2.5)) & (ARAInf(0.6))) 
    ),
    :W_GQRARAm3log_ARAInf_MS => Dict(
        :θ => [0.03,2.4,1.1,0.3,0.7],
        :data => DataFrame(System=vcat(repeat(1:2,inner=4),repeat([3],14)),Time=[3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78],Type=[1,1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,1,1,-1,1,-1,1,-1,0]),
        :vam => @vam(System & Time & Type ~ (GQR_ARAm(1.,0.3|3,log) | Weibull(0.001,2.5)) & (ARAInf(0.6))) 
    ),
    :W_GQRARAm1_ARAInf_MS => Dict(
        :θ => [0.03,2.4,1.1,0.3,0.7],
        :data => DataFrame(System=vcat(repeat(1:2,inner=4),repeat([3],14)),Time=[3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78],Type=[1,1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,1,1,-1,1,-1,1,-1,0]),
        :vam => @vam(System & Time & Type ~ (GQR_ARAm(1.,0.3|1) | Weibull(0.001,2.5)) & (ARAInf(0.6))) 
    ),
    :W_GQRARAm1log_ARAInf_MS => Dict(
        :θ => [0.03,2.4,1.1,0.3,0.7],
        :data => DataFrame(System=vcat(repeat(1:2,inner=4),repeat([3],14)),Time=[3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78],Type=[1,1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,1,1,-1,1,-1,1,-1,0]),
        :vam => @vam(System & Time & Type ~ (GQR_ARAm(1.,0.3|1,log) | Weibull(0.001,2.5)) & (ARAInf(0.6))) 
    ),
    :W_GQR_GQRARAm3_MS => Dict(
        :θ => [0.03,2.4,1.1,0.9,0.7],
        :data => DataFrame(System=vcat(repeat(1:2,inner=4),repeat([3],14)),Time=[3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78],Type=[1,1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,1,1,-1,1,-1,1,-1,0]),
        :vam => @vam(System & Time & Type ~ (GQR(1.) | Weibull(0.001,2.5)) & (GQR_ARAm(0.9,0.6|3))) 
    ),
    :W_GQR_GQRARAm3log_MS => Dict(
        :θ => [0.03,2.4,1.1,0.9,0.7],
        :data => DataFrame(System=vcat(repeat(1:2,inner=4),repeat([3],14)),Time=[3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78],Type=[1,1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,1,1,-1,1,-1,1,-1,0]),
        :vam => @vam(System & Time & Type ~ (GQR(1.) | Weibull(0.001,2.5)) & (GQR_ARAm(0.9,0.6|3,log))) 
    ),
    :W_GQRlog_GQRARAm3sqrt_MS => Dict(
        :θ => [0.03,2.4,1.1,0.9,0.7],
        :data => DataFrame(System=vcat(repeat(1:2,inner=4),repeat([3],14)),Time=[3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78],Type=[1,1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,1,1,-1,1,-1,1,-1,0]),
        :vam => @vam(System & Time & Type ~ (GQR(1.0|log) | Weibull(0.001,2.5)) & (GQR_ARAm(0.9,0.6|3,sqrt))) 
    ),
    :W_QR_GQRARAm3log_MS => Dict(
        :θ => [0.03,2.4,1.1,0.9,0.7],
        :data => DataFrame(System=vcat(repeat(1:2,inner=4),repeat([3],14)),Time=[3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78],Type=[1,1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,1,1,-1,1,-1,1,-1,0]),
        :vam => @vam(System & Time & Type ~ (QR(1.) | Weibull(0.001,2.5)) & (GQR_ARAm(0.9,0.6|3,log))) 
    ),
    :W_GQQARAm4_GQRARAm3_MS => Dict(
        :θ => [0.03,2.4,1.1,0.8,0.9,0.7],
        :data => DataFrame(System=vcat(repeat(1:2,inner=4),repeat([3],14)),Time=[3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78],Type=[1,1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,1,1,-1,1,-1,1,-1,0]),
        :vam => @vam(System & Time & Type ~ (GQR_ARAm(1.,0.7|4) | Weibull(0.001,2.5)) & (GQR_ARAm(0.9,0.6|3))) 
    ),
    :W_GQRARAm4sqrt_GQRARAm3log_MS => Dict(
        :θ => [0.03,2.4,1.1,0.8,0.9,0.7],
        :data => DataFrame(System=vcat(repeat(1:2,inner=4),repeat([3],14)),Time=[3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78],Type=[1,1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,1,1,-1,1,-1,1,-1,0]),
        :vam => @vam(System & Time & Type ~ (GQR_ARAm(1.,0.7|4,sqrt) | Weibull(0.001,2.5)) & (GQR_ARAm(0.9,0.6|3,log))) 
    ),
    :W_GQRARAm4sqrt_ARAm3_MS => Dict(
        :θ => [0.03,2.4,1.1,0.8,0.9],
        :data => DataFrame(System=vcat(repeat(1:2,inner=4),repeat([3],14)),Time=[3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78],Type=[1,1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,1,1,-1,1,-1,1,-1,0]),
        :vam => @vam(System & Time & Type ~ (GQR_ARAm(1.,0.7|4,sqrt) | Weibull(0.001,2.5)) & (ARAm(0.9|1))) 
    ),
    :W_GQRsqrt_GQRARA4log_AGAN => Dict(
        :θ => [0.3,2.3,0.9,1.3,0.7],
        :data => DataFrame(System=vcat(repeat(1:2,inner=4),repeat([3],15)),Time=[3.36,4.04,4.97,5.16, 0.78,2.36,4.05,4.97, 2.45,2.78,3.56,4.23,5.32,6.43,6.98,7.51,8.02,9.43,10.2,11.5,12,13.78,15.2],Type=[1,1,-1,1, -1,1,1,0, 1,-1,1,-1,-1,1,1,1,-1,2,1,-1,1,-1,0]),
        :vam => @vam(System & Time & Type ~ (GQR(0.8|sqrt) | Weibull(0.001,2.5)) & (GQR_ARAm(1.2,0.9|2,log) + AGAN())) 
    ),
    :W_ARA∞_2cov => Dict(
        :θ => [0.15,2.3,0.8,0.6,-0.9],
        :data => DataFrame(System=vcat(repeat([1],5),repeat([2],4),repeat([3],1),repeat([4],5)),Time=[0.800,2.646,3.190,3.916,4.109,0.910,1.127,1.245, 1.349, 0.541,1.397,1.463,2.406,2.506,3.159],Type=repeat([-1],15)),
        :datacov => DataFrame(cov1=[4.336,5.615,4.770,4.655],cov2=[0,0,0,1]),
        :vam => @vam(System & Time & Type ~ (ARAInf(0.8) | Weibull(0.15,2.3 | 0.6 * cov1 + -0.9 * cov2))),
        :rform => "System & Time & Type ~ (ARAInf(0.8) | Weibull(0.15,2.3 | 0.6 * cov1 - 0.9 * cov2))"
    ),
    :W_ARA1_1cov => Dict(
        :θ => [0.15,2.3,0.8,0.6],
        :data => DataFrame(System=repeat([1],5),Time=[0.800,2.646,3.190,3.916,4.109],Type=repeat([-1],5)),
        :datacov => DataFrame(cov1=[4.336]),
        :vam => @vam(System & Time & Type ~ (ARA1(0.8) | Weibull(0.15,2.3 | 0.6 * cov1))),
        :rform => "System & Time & Type ~ (ARA1(0.8) | Weibull(0.15,2.3 | 0.6 * cov1))"
    ),
    :L_GQRARAm3_log_1cov => Dict(
        :θ => [10,0.8,0.9,0.6,0.05],
        :data => DataFrame(System=vcat(repeat([1],5),repeat([2],4),repeat([3],1),repeat([4],5)),Time=[0.800,2.646,3.190,3.916,4.109,0.910,1.127,1.245, 1.349, 0.541,1.397,1.463,2.406,2.506,3.159],Type=repeat([-1],15)),
        :datacov => DataFrame(cov1=[4.336,5.615,4.770,4.655]),
        :vam => @vam(System&Time&Type~(GQR_ARAm(0.9,0.6|3,log)|LogLinear(0.01,0.8|0.6*cov1))),
        :rform => "System&Time&Type~(GQR_ARAm(0.9,0.6|3,log)|LogLinear(0.01,0.8|0.6*cov1))"
    ),
    :W_GQRARAm3_AGAN_2cov => Dict(
        :θ => [0.1,2.5,0.9,0.3,0.6,-0.9],
        :data => DataFrame(System=vcat(repeat([1],5),repeat([2],4),repeat([3],1),repeat([4],5)),Time=[0.800,2.646,3.190,3.916,4.109,0.910,1.127,1.245, 1.349, 0.541,1.397,1.463,2.406,2.506,3.159],Type=[1,-1,1,1,0, -1,-1,1,-1, 0, 1,-1,-1,1,1]),
        :datacov => DataFrame(cov1=[4.336,5.615,4.770,4.655],cov2=[1,0,0,1]),
        :vam => @vam(System&Time&Type~(GQR_ARAm(0.9,0.6|3,log)|Weibull(0.1,2.1|0.6*cov1+0.7*cov2))&(AGAN())),
        :rform => "System&Time&Type~(GQR_ARAm(0.9,0.6|3,log)|Weibull(0.1,2.1|0.6*cov1+0.7*cov2))&(AGAN())"
    ),
    :W_ABAO_AGAN_2cov => Dict(
        :θ => [0.1,2.5,0.6,-0.9],
        :data => DataFrame(System=vcat(repeat([1],5),repeat([2],4),repeat([3],1),repeat([4],5)),Time=[0.800,2.646,3.190,3.916,4.109,0.910,1.127,1.245, 1.349, 0.541,1.397,1.463,2.406,2.506,3.159],Type=[1,-1,1,1,0, -1,-1,1,-1, 0, 1,-1,-1,1,1]),
        :datacov => DataFrame(cov1=[4.336,5.615,4.770,4.655],cov2=[1,0,0,1]),
        :vam => @vam(System&Time&Type~(ABAO()|Weibull(0.1,2.1|0.6*cov1+0.7*cov2))&(AGAN())),
        :rform => "System&Time&Type~(ABAO()|Weibull(0.1,2.1|0.6*cov1+0.7*cov2))&(AGAN())"
    ),
    :W_GQRARAm_ARA1_GQRARAm_2cov => Dict(
        :θ => [0.11,2.5,1.15,0.62,0.92,1.25,0.75,0.65,-0.56],
        :data => DataFrame(System=vcat(repeat([1],10),repeat([2],8),repeat([3],9),repeat([4],10)),Time=[0.79,1.583,1.761,2,2.652,2.705,3.445,4,4.57,4.901,0.069,0.477,0.682,1.004,1.327,1.704,1.83,1.848,0.735,1.004,1.408,1.874,2,2.466,2.514,2.919,3.376,1.02,1.401,2,2.216,2.752,3.426,3.592,4,4.35,4.82],Type=[-1,-1,-1,2,1,-1,1,2,-1,-1,-1,-1,-1,-1,1,-1,-1,-1,1,-1,-1,1,2,-1,-1,1,1,-1,-1,2,-1,1,1,-1,2,-1,1]),
        :datacov => DataFrame(cov1=[4.336,5.615,4.770,4.655],cov2=[1,0,0,1]),
        :vam => @vam(System&Time&Type~(GQR_ARAm(1.1,0.6|3,sqrt)|Weibull(0.1,2.1|0.6*cov1+(-0.7)*cov2))&(ARA1(0.9)+GQR_ARAm(1.3,0.8|2,log))),
        :rform => "System&Time&Type~(GQR_ARAm(1.1,0.6|3,sqrt)|Weibull(0.1,2.1|0.6*cov1+(-0.7)*cov2))&(ARA1(0.9)+GQR_ARAm(1.3,0.8|2,log))"
    ),
)

# insert!(modtest,
#     :W_ARA∞_2cov => Dict(
#         :θ => [0.15,2.3,0.8,0.6,-0.9],
#         :data => DataFrame(System=vcat(repeat([1],5),repeat([2],4),3,repeat([4],5)),Time=[0.800,2.646,3.190,3.916,4.109,0.910,1.127,1.245, 1.349, 0.541,1.397,1.463,2.406,2.506,3.159],Type=repeat([-1],15)),
#         :datacov => DataFrame(cov1=[4.336,5.615,4.770,4.655],cov2=[0.0,0,0,1]),
#         :vam => @vam(System&Time&Type~(ARAInf(0.8)|Weibull(0.15,2.3|0.6*cov1 + -0.9*cov2))),
#         :rform => "System & Time & Type ~ (ARAInf(0.8)|Weibull(0.15,2.3|0.6*cov1 - 0.9*cov2))"
#     )
# )





# key = :W_ARA∞_QAGAN_QR_MS
#key = :W_GQRLog
#key = :W_GQRSqrt
# key = :W_ARAm2
# key = :W_ARAm3 
# key = :W_ARAm4
# key = :W_ARAm2_ARAm4_MS
#key = :W_GQRARA1Log
# key = :W_GQRARAm3Log


# key = :W_ARA∞_2cov
# modtest.models[key][:r]
# update!(modtest, key)
# test(modtest, key)
# modtest.results[key][:r]
# modtest.results[key][:jl]

#update!(modtest, :W_ARA∞_MS)
#test(modtest, :W_ARA∞_MS)
#update!(modtest, :W_ARAm1)
#test(modtest, :W_ARAm1)
#update!(modtest, :W_ARAm1_vs_ARA1)
#test(modtest, :W_ARAm1_vs_ARA1)
#modtest.results[:W_ARAm1][:r]
#modtest.results[:W_ARAm1][:jl]
#modtest.results[:W_ARAm1_vs_ARA1][:r]
#modtest.results[:W_ARAm1_vs_ARA1][:jl]

#update!(modtest, :W_ARAm∞)
#test(modtest, :W_ARAm∞)
#update!(modtest, :W_ARAm∞_vs_ARA∞)
#test(modtest, :W_ARAm∞_vs_ARA∞)
#modtest.results[:W_ARAm∞][:r]
#modtest.results[:W_ARAm∞][:jl]
#modtest.results[:W_ARAm∞_vs_ARA∞][:r]
#modtest.results[:W_ARAm∞_vs_ARA∞][:jl]

#empty!(modtest)

update!(modtest)
test(modtest)

VirtualAgeModels.isbayesian(modtest.models[:W_ARAm∞][:vam].models[1])