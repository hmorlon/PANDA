function branch_time_id(tree::Tree, times::Array{Float64,1})
    nb = node_depths_base(tree)
    ne = node_depths(tree)

    rd= maximum(ne)
    ids = Array{Array{Int64,1},1}(undef,0)

    for t in times
        if t <= rd
            idt = Array{Int64,1}(undef,0)

            for i in 1:length(nb)
                if nb[i] <= t < ne[i]
                    push!(idt,i)
                end
            end
            push!(ids,idt)
        end
    end

    return ids
end

function branch_time_id(tree::Tree, edge_trees::Array{EdgeTreeRates2,1}, times::Array{Float64,1})
    nb = node_depths_base(tree, edge_trees)

    ne = node_depths(tree,edge_trees)

    rd= maximum(ne)
    ids = Array{Array{Int64,1},1}(undef,0)

    for t in times
        if t <= rd
            idt = Array{Int64,1}(undef,0)

            for i in 1:length(nb)
                if nb[i] <= t < ne[i]
                    push!(idt,i+1)
                end
            end
            push!(ids,idt)
        end
    end

    return ids
end


function time_rates(ids::Array{Array{Int64,1},1}, row::Array{Float64,1})
    mean_rates = Array{Float64,1}(undef,0)

    for id in ids
        mr = 0.
        n = 0.
        for i in id
            n += 1.
            mr += row[i]
        end
        push!(mean_rates, mr/n)
    end
    return mean_rates
end

function time_rates(tree::Tree, row::Array{Float64,1}, times::Array{Float64,1})
    ids = branch_time_id(tree, times)

    time_rates(ids, row)
end

function time_rates(tree::Tree, times::Array{Float64,1})
    ids = branch_time_id(tree, times)
    row = extract_rates(tree)
    time_rates(ids, row)
end

function time_rates(tree::Tree, edge_trees::Array{EdgeTreeRates2,1}, times::Array{Float64,1})
    ids = branch_time_id(tree, edge_trees, times)
    row = extract_rates(tree, edge_trees)
    time_rates(ids, row)
end

function time_rates(tree::Tree, row::Array{Float64,1}, times::Array{Float64,1})
    ids = branch_time_id(tree, times)
    time_rates(ids, row)
end

function plot_time_rates(tree::Tree, times::Array{Float64,1})
    tr = time_rates(tree, times)
    t = times[1:length(tr)]

    @rput t
    @rput tr

    reval("""
        library(scales)
        plot(100000, axes = F, xlim = range(t), ylim = log(range(tr)), xlab = "time", ylab = "mean_rate")
        lines(t, log(tr), col = "black", lwd = 6, lty = 1)
    """)
end

function plot_time_rates(tree::Tree, row::Array{Float64,1}, times::Array{Float64,1})
    tr = time_rates(tree, row, times)
    t = times[1:length(tr)]

    @rput t
    @rput tr

    reval("""
        library(scales)
        plot(100000, axes = F, xlim = range(t), ylim = log(range(tr)), xlab = "time", ylab = "mean_rate")
        lines(t, log(tr), col = "black", lwd = 6, lty = 1)
        axis(1)
        axis(2)
    """)
end

function plot_time_rates(tree::Tree, row::Array{Float64,1}, times::Array{Float64,1}, miny, maxy)
    tr = time_rates(tree, row, times)
    t = times[1:length(tr)]

    @rput t
    @rput tr
    @rput miny
    @rput maxy

    reval("""
        library(scales)
        plot(100000, axes = F, xlim = range(t), ylim = c(miny,maxy), xlab = "time", ylab = "mean_rate")
        lines(t, log(tr), col = "black", lwd = 6, lty = 1)
        axis(1)
        axis(2)
    """)
end

function plot_time_rates(tree::Tree, sim_rates::Array{Float64,1}, row::Array{Float64,1}, times::Array{Float64,1})
    tr = time_rates(tree, row, times)
    sim_tr = time_rates(tree, sim_rates, times)
    t = times[1:length(tr)]

    @rput t
    @rput tr
    @rput sim_tr

    reval("""
        library(scales)
        plot(100000, axes = F, xlim = range(t), ylim = log(range(c(sim_tr, tr))), xlab = "time", ylab = "mean_rate")
        lines(t, log(sim_tr), col = "orange1", lwd = 6, lty = 1)
        lines(t, log(sim_tr), col = "orange3", lwd = 4, lty = 2)
        lines(t, log(tr), col = "cadetblue1", lwd = 5, lty = 1)
        lines(t, log(tr), col = "cadetblue4", lwd = 3, lty = 6)
        axis(1)
        axis(2)
    """)
end


function plot_time_rates_mcmc(tree::Tree, sampler ; nplot = 50, miny = -1, maxy = 1, alpha_col = 0.05, burn = 0.)
    times = sampler[1][6]
    N = length(sampler[1][8][1])
    tr = time_rates(tree,  map(exp, sampler[2][4:(tree.n_nodes+3)]), times)
    t = times[1:length(tr)]
    plot_id = Array{Int64,1}(floor.(range(2+ floor(burn * N), length=nplot, stop=N)))
    mean_mr = zeros(length(tr))
    a = Array{Float64,1}(undef,0)
    map_mr = [deepcopy(a) for i in 1:length(tr)]
    @rput t
    @rput tr
    @rput miny
    @rput maxy
    @rput alpha_col

    reval("""
        library(scales)
        plot(100000, axes = F, xlim = range(t), ylim = log(range(tr))+c(miny,maxy), xlab = "time", ylab = "mean_rate")
    """)
    n = 0

    colors = ["deepskyblue2" ; "orange2" ; "orangered3"]
    for i in 1:3
        color = colors[i]
        @rput color
        for k in plot_id
            y = sampler[1][8][i][k]
            for j in 1:length(y)
                push!(map_mr[j],log(y[j]))
            end
            @rput y
            reval("""lines(t, log(y), col = alpha(color, alpha = alpha_col), lwd = 2)""")
            mean_mr .+= y
            n += 1
        end
    end

    mean_mr ./= n
    @rput mean_mr
    @rput map_mr
    reval("""
        maps = sapply(map_mr, function(x){D=density(x); return(D[[1]][which.max(D[[2]])])})
    """)
    reval("""
        lines(t, log(tr), col = "cadetblue1", lwd = 5, lty = 1)
        lines(t, log(tr), col = "cadetblue4", lwd = 3, lty = 6)
        lines(t, maps, col = "darkseagreen1", lwd = 5, lty = 1)
        lines(t, maps, col = "darkseagreen4", lwd = 3, lty = 6)
        axis(1)
        axis(2)
    """)

end

function plot_time_rates_mcmc(tree::Tree, mrChains, tr, sim_mr, times ;
    nplot = 50, miny = -1, maxy = 1, alpha_col = 0.05, burn = 0., extract_at = Array{Int64,1}(undef,0))
    N = length(mrChains[1])
    t = times[1:length(tr)]
    plot_id = Array{Int64,1}(floor.(range(2 + floor(burn * N), length=nplot, stop=N)))
    mean_mr = zeros(length(tr))
    a = Array{Float64,1}(undef,0)
    map_mr = [deepcopy(a) for i in 1:length(tr)]
    @rput t
    @rput tr
    @rput miny
    @rput maxy
    @rput alpha_col
    @rput sim_mr

    reval("""
        library(scales)
        plot(100000, axes = F, xlim = range(t), ylim = log(range(c(sim_mr,tr)))+c(miny,maxy), xlab = "time", ylab = "mean_rate")
    """)
    n = 0

    #colors = ["deepskyblue2" ; "orange2" ; "orangered3"]
    colors = ["deepskyblue2" ; "deepskyblue2" ; "deepskyblue2"]

    for i in 1:3
        color = colors[i]
        @rput color
        for k in plot_id
            y = mrChains[i][k]
            for j in 1:length(y)
                push!(map_mr[j],log(y[j]))
            end
            @rput y
            reval("""lines(t, log(y), col = alpha(color, alpha = alpha_col), lwd = 2)""")
            mean_mr .+= y
            n += 1
        end
    end

    mean_mr ./= n
    @rput mean_mr
    @rput map_mr
    reval("""
        maps = sapply(map_mr, function(x){D=density(x); return(D[[1]][which.max(D[[2]])])})
        quant = sapply(map_mr, function(x){quantile(x, probs = c(0.05,0.95))})
    """)
    reval("""
        lines(t, quant[1,], col = "deepskyblue3", lwd = 2)
        lines(t, quant[2,], col = "deepskyblue3", lwd = 2)
        lines(t, log(tr), col = "cadetblue1", lwd = 4, lty = 1)
        lines(t, log(tr), col = "cadetblue4", lwd = 2, lty = 6)
        lines(t, log(sim_mr), col = "orange1", lwd = 5, lty = 1)
        lines(t, log(sim_mr), col = "orange3", lwd = 3, lty = 6)
        lines(t, maps, col = "darkseagreen1", lwd = 4, lty = 1)
        lines(t, maps, col = "darkseagreen4", lwd = 2, lty = 6)
        axis(1)
        axis(2)
    """)

    @rget maps
    return maps[extract_at]

end


function plot_time_rates_mcmc_lpm(tree::Tree, mrChains, times, time, lambda, tr ;
    nplot = 50, miny = -1, maxy = 3, alpha_col = 0.05, burn = 0., extract_at = Array{Int64,1}(undef,0))

    rd = times[length(mrChains[1][1])]

    t_sim = deepcopy(time .* (-1.) .+ rd)
    ind = [1]
    for i in 2:length(t_sim)
        if t_sim[i]>0
            push!(ind,i)
        end
    end

    t_sim = t_sim[ind]
    push!(t_sim, 0.)

    l_sim = deepcopy(lambda)
    l_sim = l_sim[ind]
    push!(l_sim, l_sim[end])

    N = length(mrChains[1])
    t = times[1:length(mrChains[1][1])]
    plot_id = Array{Int64,1}(floor.(range(2 + floor(burn * N), length=nplot, stop=N)))
    mean_mr = zeros(length(t))
    a = Array{Float64,1}(undef,0)
    map_mr = [deepcopy(a) for i in 1:length(t)]
    @rput t
    @rput miny
    @rput maxy
    @rput alpha_col
    @rput time
    @rput lambda
    @rput t_sim
    @rput l_sim
    @rput tr

    reval("""
        library(scales)
        plot(100000, axes = F, xlim = range(t), ylim = log(range(c(lambda, tr)))+c(miny,maxy), xlab = "time", ylab = "mean_rate")
    """)
    n = 0

    #colors = ["deepskyblue2" ; "orange2" ; "orangered3"]
    colors = ["deepskyblue2" ; "deepskyblue2" ; "deepskyblue2"]

    for i in 1:3
        color = colors[i]
        @rput color
        for k in plot_id
            y = mrChains[i][k]
            for j in 1:length(y)
                push!(map_mr[j],log(y[j]))
            end
            @rput y
            reval("""lines(t, log(y), col = alpha(color, alpha = alpha_col), lwd = 2)""")
            mean_mr .+= y
            n += 1
        end
    end

    mean_mr ./= n
    @rput mean_mr
    @rput map_mr
    reval("""
        maps = sapply(map_mr, function(x){D=density(x); return(D[[1]][which.max(D[[2]])])})
        quant = sapply(map_mr, function(x){quantile(x, probs = c(0.05,0.95))})
    """)
    reval("""
        lines(t, quant[1,], col = "deepskyblue3", lwd = 2)
        lines(t, quant[2,], col = "deepskyblue3", lwd = 2)
        lines(t, log(tr), col = "cadetblue1", lwd = 4, lty = 1)
        lines(t, log(tr), col = "cadetblue4", lwd = 2, lty = 6)
        lines(t_sim, log(l_sim), col = "orange1", lwd = 5, lty = 1, type = "s")
        lines(t_sim, log(l_sim), col = "orange3", lwd = 3, lty = 6, type = "s")
        lines(t, maps, col = "darkseagreen1", lwd = 4, lty = 1)
        lines(t, maps, col = "darkseagreen4", lwd = 2, lty = 6)
        axis(1)
        axis(2)
    """)

    @rget maps
    return maps[extract_at]

end


function plot_time_rates_mcmc(tree::Tree, mrChains::Array{Any,1};
    nplot = 50, miny = -1, maxy = 1, alpha_col = 0.05, burn = 0., extract_at = Array{Int64,1}(undef,0))
    N = length(mrChains[1])
    nt = length(mrChains[1][1])
    new_tree = Tree(tree.offsprings, 0., tree.attributes)
    root_depth = maximum(node_depths(new_tree)) * 1.01
    t = [0:nt...] * root_depth/nt
    t = t[1:(end-1)]

    plot_id = Array{Int64,1}(floor.(range(2 + floor(burn * N), length=nplot, stop=N)))
    mean_mr = zeros(length(t))
    a = Array{Float64,1}(undef,0)
    map_mr = [deepcopy(a) for i in 1:length(t)]
    @rput t
    @rput miny
    @rput maxy
    @rput alpha_col

    reval("""
        library(scales)
        plot(100000, axes = F, xlim = range(t), ylim = c(miny,maxy), xlab = "time", ylab = "mean_rate")
    """)
    n = 0

    #colors = ["deepskyblue2" ; "orange2" ; "orangered3"]
    colors = ["deepskyblue2" ; "deepskyblue2" ; "deepskyblue2"]

    for i in 1:3
        color = colors[i]
        @rput color
        for k in plot_id
            y = mrChains[i][k]
            for j in 1:length(y)
                push!(map_mr[j],log(y[j]))
            end
            @rput y
            println(length(y))
            reval("""lines(t, log(y), col = alpha(color, alpha = alpha_col), lwd = 2)""")
            mean_mr .+= y
            n += 1
        end
    end

    mean_mr ./= n
    @rput mean_mr
    @rput map_mr
    reval("""
        maps = sapply(map_mr, function(x){D=density(x); return(D[[1]][which.max(D[[2]])])})
        quant = sapply(map_mr, function(x){quantile(x, probs = c(0.05,0.95))})
    """)
    reval("""
        lines(t, quant[1,], col = "deepskyblue3", lwd = 2)
        lines(t, quant[2,], col = "deepskyblue3", lwd = 2)
        lines(t, maps, col = "darkseagreen1", lwd = 4, lty = 1)
        lines(t, maps, col = "darkseagreen4", lwd = 2, lty = 6)
        axis(1)
        axis(2)
    """)

    @rget maps
    return maps[extract_at]

end
