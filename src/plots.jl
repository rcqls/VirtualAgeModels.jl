abstract type VirtualAgePlot end
abstract type IntensityPlot end
abstract type CummulativeIntensityPlot end
abstract type ConditionalCummulativeDensityFunctionPlot end
abstract type ConditionalSurvivalPlot end
abstract type ConditionalProbabilityDensityFunctionPlot end

const symbol2typeplot = Dict(
    :v => VirtualAgePlot, :virtual_age => VirtualAgePlot,
    :i => IntensityPlot, :intensity => IntensityPlot,
    :I => CummulativeIntensityPlot, :cummulative_intensity => CummulativeIntensityPlot,
    :F => ConditionalCummulativeDensityFunctionPlot, :conditional_cummulative_density_function => ConditionalCummulativeDensityFunctionPlot,
    :S => ConditionalSurvivalPlot, :conditional_survival => ConditionalSurvivalPlot,
    :f => ConditionalProbabilityDensityFunctionPlot, :conditional_probability_density_function => ConditionalProbabilityDensityFunctionPlot
)

import Plots.plot
import Plots.plot!

function plot(m::Model, type::Symbol=:v)
    plot(m, symbol2typeplot[type])
end

function plot!(m::Model, type::Symbol)
    plot!(m, symbol2typeplot[type])
end

function plot(m::Model, ::Type{VirtualAgePlot})
    infos = virtual_age_infos(m, m.time[1], m.time[end])
    # println(infos)
    plot(infos.x, infos.y, legend=nothing, color=:blue)
end

function plot!(m::Model, ::Type{VirtualAgePlot})
    infos = virtual_age_infos(m, m.time[1], m.time[end])
    # println(infos)
    plot!(infos.x, infos.y, legend=nothing, color=:blue)
end

function plot(m::Model, ::Type{IntensityPlot})
    infos = virtual_age_infos(m, m.time[1], m.time[end], type=:i)
    # println(infos)
    plot(infos.x, infos.y, legend=nothing, color=:blue)
end

function plot!(m::Model, ::Type{IntensityPlot})
    infos = virtual_age_infos(m, m.time[1], m.time[end], type=:i)
    # println(infos)
    plot!(infos.x, infos.y, legend=nothing, color=:blue)
end

function plot(m::Model, ::Type{CummulativeIntensityPlot})
    infos = virtual_age_infos(m, m.time[1], m.time[end], type=:I)
    # println(infos)
    plot(infos.x, infos.y, legend=nothing, color=:blue)
end

function plot!(m::Model, ::Type{CummulativeIntensityPlot})
    infos = virtual_age_infos(m, m.time[1], m.time[end], type=:I)
    # println(infos)
    plot!(infos.x, infos.y, legend=nothing, color=:blue)
end

function plot(m::Model, ::Type{ConditionalCummulativeDensityFunctionPlot})
    infos = virtual_age_infos(m, m.time[1], m.time[end], type=:F)
    # println(infos)
    plot(infos.x, infos.y, legend=nothing, color=:blue)
end

function plot!(m::Model, ::Type{ConditionalCummulativeDensityFunctionPlot})
    infos = virtual_age_infos(m, m.time[1], m.time[end], type=:F)
    # println(infos)
    plot!(infos.x, infos.y, legend=nothing, color=:blue)
end

function plot(m::Model, ::Type{ConditionalSurvivalPlot})
    infos = virtual_age_infos(m, m.time[1], m.time[end], type=:S)
    # println(infos)
    plot(infos.x, infos.y, legend=nothing, color=:blue)
end

function plot!(m::Model, ::Type{ConditionalSurvivalPlot})
    infos = virtual_age_infos(m, m.time[1], m.time[end], type=:S)
    # println(infos)
    plot!(infos.x, infos.y, legend=nothing, color=:blue)
end

function plot(m::Model, ::Type{ConditionalProbabilityDensityFunctionPlot})
    infos = virtual_age_infos(m, m.time[1], m.time[end], type=:f)
    # println(infos)
    plot(infos.x, infos.y, legend=nothing, color=:blue)
end

function plot!(m::Model, ::Type{ConditionalProbabilityDensityFunctionPlot})
    infos = virtual_age_infos(m, m.time[1], m.time[end], type=:f)
    # println(infos)
    plot!(infos.x, infos.y, legend=nothing, color=:blue)
end