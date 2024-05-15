using VirtualAgeModels
using DataFrames
using RCall

m = @vam(System & Time & Type ~ (ARA∞(0.6) | Weibull(1.0,3.0)))

R"""
require(VAM)
data("OREngines_CM")
"""

@rget OREngines_CM


## data! with DataFrame
data!(m, OREngines_CM)

# m.data[1]
data(m)

m.data

data2 = copy(m.data)


## data! with DataFrame[]

empty(m.data)

data2 = convert.(DataFrame,data2)

data!(m, data2)

m.data
