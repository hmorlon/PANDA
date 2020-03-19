include("load_ClaDS2_functions.jl")

R"library(RPANDA)"
R"data(Cetacea)"
@rget Cetacea
tree = ape2Tree(Cetacea)

plot_ClaDS(tree)

sampling_fraction = 87/89

burn = 0.25
thin = 10

Random.seed!(813)
n_iter = 50

sampler = run_ClaDS2(tree,          # the phylogeny, see above for how to load it
    n_iter,                         # the number of iterations between 2 gelman statistics computations
    print_state = 25,               # when to print the state of the chain, put to 0 not to print it
    thin = thin,                    # the thinning parameter
    burn = burn,                    # the burning parameter (for computation of the gelman)
    plot_chain = true,              # are the chains plotted during the run ? (every n_iter iterations)
    plot_tree= 50,                  # when to plot the augmented tree, put to 0 not to plot it
    f=sampling_fraction,            # the sampling fraction
    end_it = Inf, ltt_steps = 5)


#=
    The function also saves the number of alive species in the augmented tree at given time steps.
    The number of time steps is given by ltt_steps

    if youdon't want to plot anything during the run of the function,
     plot_tree should be set to 0 and plot_chain to false (the default values)
=#


result = (tree, sampler, sampling_fraction)
@save "Example_Ceatacea.jld2" result

#plot the result in Julia
@load "Example_Ceatacea.jld2" result
sampler = result[2]
tree = result[1]
MAPS = sampler[2]
npar = tree.n_nodes + 3


R"par(mfrow = c(1,2))"
# plot the tree with painted branches :
plot_ClaDS(tree, map(exp, MAPS[5:npar]))

# plot the number of lineages through time for the augmented tree :
plot_LTT_chain(sampler, tree, npar, n_ltt = 100)
#= the black line is the number of lineages in the reconstructed tree (the LTT),
The thin blue lines are the result for on given iteration in the MCMC, with the 95%
interval displayed with the darker blue lines, the dotted blue lines summarize the result in the chain,
the blue one being the Maximum A Posteriory at each time point and the green one the mean
To increase the number of time points, change the ltt_steps parameter in run_ClaDS2
=#

#save the results as a Rdata file
@load "Example_Ceatacea.jld2" result
sampler = result[2]
tree = result[1]
sample_fraction = result[3]
sampler_to_Rdata(tree, sampler, "ClaDS2_Cetacea.Rdata" ; sample_fraction = sample_fraction)

# you can then access and plot the results in R as shown in "plot_result.R"
