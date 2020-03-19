include("load_ClaDS2_functions.jl")

R"library(RPANDA)"

# load the tree in R
R"data(Cetacea)"
@rget Cetacea
tree = ape2Tree(Cetacea)

#= more generally, the code would be

reval("""
    Rtree = read.tree("phylo.tre")               # here the code you would use to load the tree in R
""")
@rget Rtree                                      # import the tree in Julia
tree tree = ape2Tree(Rtree)                      # and make it a Julia object

=#


plot_ClaDS(tree)

Random.seed!(813)
n_iter = 400

sample_fraction = 87/89
sampler = run_ClaDS2(tree, n_iter,plot_tree = Int64(n_iter/4),
    print_state = Int64(n_iter/4), plot_chain = true,
    f=87/89, ltt_steps = 50)
#=
The function also saves the number of alive species in the augmented tree at given time steps.
The number of time steps is given by ltt_steps

plot_tree should be set to 0 and plot_chain to false (the default values) so that nothing is plotted
during the run of the function
=#


result = (tree, sampler, sample_fraction)
@save "Example_Ceatacea.jld2" result


#plot the result in Julia
@load "Example_Ceatacea.jld2" result
sampler = result[2]
tree = result[1]
MAPS = sampler[2]
npar = tree.n_nodes + 3

R"par(mfrow = c(1,2))"
plot_ClaDS(tree, map(exp, MAPS[5:npar]))
plot_LTT_chain(sampler, tree, npar, n_ltt = 100)


#save the results as a Rdata file
@load "Example_Ceatacea.jld2" result
sampler = result[2]
tree = result[1]
sample_fraction = result[3]
sampler_to_Rdata(tree, sampler, "ClaDS2_Cetacea.Rdata" ; sample_fraction = sample_fraction)
