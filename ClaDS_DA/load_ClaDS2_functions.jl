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
include("update_edges_quad.jl")
println(" done")

print("loading ClaDS2_enhance")
include("ClaDS2_enhance_ETR.jl")
println(" done")

print("loading prune")
include("prune.jl")
println(" done")

print("loading ETR2 functions")
include("ETR2.jl")
include("get_parent_rates_ETR2.jl")
include("change_edges_ETR2.jl")
include("gibbs_ETR2.jl")
include("update_edges_ETR2.jl")
include("mcmc_ETR2.jl")
include("extract_from_ETR2.jl")
include("LTT_ETR2.jl")
include("gelman.jl")
println(" done")

print("loading mean_rates")
include("rates_through_time.jl")
println(" done")

print("loading mcmc_ClaDS2")
include("mcmc_ClaDS2.jl")
println(" done")

print("loading LTT")
include("LTT.jl")
println(" done")

include("load_ETR2.jl")
