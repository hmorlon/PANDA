print("loading needed packages")
using RCall
using StatsBase
using Distributions
using Random
using BenchmarkTools
using JLD2
println(" done")

print("loading Tree")
include("Tree.jl")
println(" done")

print("loading sim_ClaDS")
include("sim_ClaDS.jl")
println(" done")

print("loading gibbs_functions")
include("gibbs_functions.jl")
println(" done")

print("loading ClaDS2_enhance")
include("ClaDS2_enhance.jl")
println(" done")

print("loading mcmc_ClaDS2")
include("mcmc_ClaDS2.jl")
println(" done")

print("loading LTT")
include("LTT.jl")
println(" done")

print("loading update_rates")
include("update_rates.jl")
println(" done")

print("loading update_edges")
include("update_edges.jl")
println(" done")

print("loading ClaDS2_enhance")
include("ClaDS2_enhance_ETR.jl")
println(" done")

print("loading prune")
include("prune.jl")
println(" done")
