function node_depths_base(tree::Tree)
    function aux(subtree, node_times, current_depth)
        if subtree.n_nodes < 2
            pushfirst!(node_times, current_depth)
        else
            aux(subtree.offsprings[2],node_times, current_depth + subtree.branch_length)
            aux(subtree.offsprings[1],node_times, current_depth + subtree.branch_length)
            pushfirst!(node_times,current_depth)
        end
    end

    node_times = Array{Float64,1}(undef,0)
    aux(tree, node_times, 0.)
    #pushfirst!(node_times,0.)
    return node_times
end

function node_depths_base(tree::Tree, edge_trees::Union{Array{EdgeTree,1},Array{EdgeTreeRates,1}})
    function aux(subtree, node_times, current_depth)
        if subtree.n_nodes < 2
            pushfirst!(node_times, current_depth)
        else
            aux(subtree.offsprings[2],node_times, current_depth + subtree.branch_length)
            aux(subtree.offsprings[1],node_times, current_depth + subtree.branch_length)
            pushfirst!(node_times,current_depth)
        end
    end

    node_times = node_depths_base(tree)

    node_times_full = Array{Float64,1}(undef,0)
    for i in length(edge_trees):-1:1
        aux(edge_trees[i].tree, node_times_full, node_times[i])
    end
    #pushfirst!(node_times,0.)
    return node_times_full
end

function node_depths(tree::Tree)
    function aux(subtree, node_times, current_depth)
        if subtree.n_nodes < 2
            pushfirst!(node_times, current_depth + subtree.branch_length)
        else
            aux(subtree.offsprings[2],node_times, current_depth + subtree.branch_length)
            aux(subtree.offsprings[1],node_times, current_depth + subtree.branch_length)
            pushfirst!(node_times,current_depth+ subtree.branch_length)
        end
    end

    node_times = Array{Float64,1}(undef,0)
    aux(tree, node_times, 0.)
    #pushfirst!(node_times,0.)
    return node_times
end

function node_depths(tree::Tree, edge_trees::Union{Array{EdgeTree,1},Array{EdgeTreeRates,1}})
    function aux(subtree, node_times, current_depth)
        if subtree.n_nodes < 2
            pushfirst!(node_times, current_depth + subtree.branch_length)
        else
            aux(subtree.offsprings[2],node_times, current_depth + subtree.branch_length)
            aux(subtree.offsprings[1],node_times, current_depth + subtree.branch_length)
            pushfirst!(node_times,current_depth+ subtree.branch_length)
        end
    end

    node_times = node_depths_base(tree)

    node_times_full = Array{Float64,1}(undef,0)
    for i in length(edge_trees):-1:1
        aux(edge_trees[i].tree, node_times_full, node_times[i])
    end
    #pushfirst!(node_times,0.)
    return node_times_full
end

function LTT(tree::Tree)
    function aux(subtree, node_times, events, current_depth)
        if subtree.n_nodes < 2
            pushfirst!(node_times, current_depth + subtree.branch_length)
            if subtree.extant
                pushfirst!(events, 0)
            else
                pushfirst!(events, -1)
            end
        else
            aux(subtree.offsprings[2],node_times,events, current_depth + subtree.branch_length)
            aux(subtree.offsprings[1],node_times,events, current_depth + subtree.branch_length)
            pushfirst!(node_times,current_depth+ subtree.branch_length)
            pushfirst!(events,1)
        end
    end

    node_times = Array{Float64,1}(undef,0)
    events = Array{Int64,1}(undef,0)
    aux(tree, node_times, events, 0.)
    indices = sortperm(node_times)
    events = events[indices]
    node_times = node_times[indices]
    if length(events) > 1
        while events[end-1] == 0
            pop!(events)
            pop!(node_times)
        end
    end
    events[1] += 1
    events = cumsum(events, dims = 1)


    return node_times, events
end

function LTT(tree::Tree, times::Array{Float64,1})
    n = length(times)

    function aux(subtree::Tree, from_root::Float64, id_time::Int64, ltt::Array{Int64,1})
        if id_time <= n
            if (subtree.branch_length + from_root) <= times[id_time]
                if subtree.n_nodes < 2
                    if !subtree.extant
                        ltt[id_time] += -1
                    end
                else
                    ltt[id_time] += 1
                    aux(subtree.offsprings[2],from_root+ subtree.branch_length, id_time, ltt)
                    aux(subtree.offsprings[1],from_root+ subtree.branch_length, id_time, ltt)
                end
            else
                aux(subtree,from_root, id_time+1, ltt)
            end
        end

    end

    ltt = fill(0, length(times))
    ltt[1] = 1
    aux(tree, 0., 1, ltt)
    ltt = cumsum(ltt, dims = 1)


    return times, ltt

end

function LTT(tree::Tree, edge_trees::Union{Array{EdgeTree,1},Array{EdgeTreeRates,1}}, times::Array{Float64,1}, live_nd::Array{Float64,1})
    n = length(times)

    function aux(subtree::Tree, from_root::Float64, id_time::Int64, ltt::Array{Int64,1})
        if id_time <= n
            if (subtree.branch_length + from_root) <= times[id_time]
                if subtree.n_nodes < 2
                    if !subtree.extant
                        ltt[id_time] += -1
                    end
                else
                    ltt[id_time] += 1
                    aux(subtree.offsprings[2],from_root+ subtree.branch_length, id_time, ltt)
                    aux(subtree.offsprings[1],from_root+ subtree.branch_length, id_time, ltt)
                end
            else
                aux(subtree,from_root, id_time+1, ltt)
            end
        end

    end

    ltt = fill(0, length(times))
    ltt[1] = 1
    aux(tree, 0., 1, ltt)
    #println("1 $ltt ")
    for i in 1:length(edge_trees)
        if edge_trees[i].tree.n_nodes > 1
            #print("$i $(live_nd[i+1] + maximum(node_depths(edge_trees[i][1]))) ; ")
            aux(edge_trees[i].tree, live_nd[i+1], 1, ltt)
        end
    end
    #println("2 $ltt ")
    ltt = cumsum(ltt, dims = 1)
    #println("3 $ltt ")

    return times, ltt
end

function root_tree(complete_tree::Tree)

    function aux(subtree::Tree)
        if subtree.n_nodes < 2
            return Tree(subtree.offsprings, 0., subtree.attributes)
        else
            if (subtree.offsprings[1].extant & subtree.offsprings[2].extant)
                return Tree(subtree.offsprings, 0., subtree.attributes)
            elseif subtree.offsprings[1].extant
                aux(subtree.offsprings[1])
            else
                aux(subtree.offsprings[2])
            end
        end
    end

    aux(complete_tree)
end

function plot_LTT_chain(mcmc_sampler, complete_tree::Tree, extant_tree::Tree, y_max::Int64 ; n_ltt = 100, alpha_col = 0.05, burn = 0.25, thin = 1)
    chains = mcmc_sampler[1][1]
    npar = mcmc_sampler[1][4][1].n_nodes + 3 + (mcmc_sampler[1][4][1].n_nodes+1)/2
    ltt_times = mcmc_sampler[1][6]
    maps = mcmc_sampler[2]

    pruned_tree = root_tree(complete_tree)
    sim_ltt = LTT(pruned_tree, ltt_times)[2]
    new_tree = Tree(extant_tree.offsprings, 0., extant_tree.attributes)
    live_ltt = LTT(new_tree, ltt_times)[2]

    @rput chains
    @rput maps
    @rput npar
    @rput ltt_times
    @rput n_ltt
    @rput y_max
    @rput sim_ltt
    @rput alpha_col
    @rput burn
    @rput thin
    @rput live_ltt

    reval("""
        require(coda)
        n_row = length(chains[[1]][[1]])
        ini = floor(burn * n_row)
        if(ini == 0) ini = 1
        it = seq(ini, n_row, thin)
        npar2 = length(chains[[1]])
        map_chains = mcmc.list(lapply(1:3, function(i){
            mcmc(sapply(1:npar2, function(j){
                if(j<4 || (j>npar )){#}& j<(npar +3))) {
                    chains[[i]][[j]][it]
                }else{
                    log(chains[[i]][[j]][it])
                }
        }))}))
        unlist_chains = (sapply(1:npar2, function(k){
            c(map_chains[[1]][,k], map_chains[[2]][,k], map_chains[[3]][,k])
            }))

        #maps_log=sapply(1:npar2, function(i){D=density(unlist_chains[,i]);
        #    return(D[[1]][which.max(D[[2]])])})

        means=sapply(1:npar2, function(i){if(i>= npar +3){log(mean(unlist_chains[,i]))}else{mean(unlist_chains[,i])}})
    """)

    reval("""
        library(scales)
        plot(100000, axes = F, xlim = range(ltt_times), ylim = log(c(2,y_max)), xlab = "time", ylab = "n_lineage")
        y_lab = c(2,5,10,50,100,200,500,1000,2000,5000,10000,20000,50000)
        y_lab = y_lab[y_lab<=y_max]
        lines(ltt_times, log(live_ltt), col = "black", lwd = 6, lty = 1)

        id_plot = unique(floor(seq(ini, n_row, length.out = n_ltt)))#unique(floor(seq(2,length(chains[[1]][[1]]), length.out = n_ltt)))
        id_ltt = (npar+3):length(chains[[1]])
        Ys = c()
        for (i in 1:3){
            for (j in id_plot){
                y = sapply(id_ltt, function(k){log(chains[[i]][[k]][j])})
                Ys = rbind(Ys,y)
                lines(ltt_times, y, col = alpha("deepskyblue2", alpha = alpha_col), lwd = 2)
                }
            }
        quant = sapply(1:ncol(Ys), function(i){quantile(Ys[,i], probs = c(0.05,0.95))})

        lines(ltt_times, quant[1,], col = "deepskyblue3", lwd = 2)
        lines(ltt_times, quant[2,], col = "deepskyblue3", lwd = 2)
        lines(ltt_times, log(sim_ltt), col = "orange1", lwd = 6, lty = 1)
        lines(ltt_times, log(sim_ltt), col = "orange3", lwd = 4, lty = 2)
        lines(ltt_times, log(maps[id_ltt]), col = "cadetblue1", lwd = 5, lty = 1)
        lines(ltt_times, log(maps[id_ltt]), col = "cadetblue4", lwd = 3, lty = 6)

        lines(ltt_times, means[id_ltt], col = "darkseagreen1", lwd = 5, lty = 1)
        lines(ltt_times, means[id_ltt], col = "darkseagreen4", lwd = 3, lty = 6)

        axis(1)
        axis(2, at = log(y_lab), lab = y_lab)

        #print(log(maps[id_ltt]) - maps_log[id_ltt])
    """)
end

function plot_LTT_chain_extinct(mcmc_sampler, complete_tree::Tree, extant_tree::Tree, y_max::Int64 ; n_ltt = 100, alpha_col = 0.05, burn = 0.25, thin = 1)
    chains = mcmc_sampler[1][1]
    npar = mcmc_sampler[1][4][1].n_nodes + 3 + (mcmc_sampler[1][4][1].n_nodes+1)/2
    ltt_times = mcmc_sampler[1][6]
    maps = mcmc_sampler[2]

    pruned_tree = root_tree(complete_tree)
    sim_ltt = LTT(pruned_tree, ltt_times)[2]
    new_tree = Tree(extant_tree.offsprings, 0., extant_tree.attributes)
    live_ltt = LTT(new_tree, ltt_times)[2]

    @rput chains
    @rput maps
    @rput npar
    @rput ltt_times
    @rput n_ltt
    @rput y_max
    @rput sim_ltt
    @rput alpha_col
    @rput live_ltt
    @rput burn

    reval("""
        require(coda)
        n_row = length(chains[[1]][[1]])
        ini = floor(burn * n_row)
        if(ini == 0) ini = 1
        it = seq(ini, n_row, thin)
        npar2 = length(chains[[1]])
        map_chains = mcmc.list(lapply(1:3, function(i){
            mcmc(sapply(1:npar2, function(j){
                if(j<4 || (j>npar )){#}& j<(npar +3))) {
                    chains[[i]][[j]][it]
                }else{
                    log(chains[[i]][[j]][it])
                }
        }))}))
        unlist_chains = (sapply(1:npar2, function(k){
            c(map_chains[[1]][,k], map_chains[[2]][,k], map_chains[[3]][,k])
            }))

        #maps_log=sapply(1:npar2, function(i){D=density(unlist_chains[,i]);
        #    return(D[[1]][which.max(D[[2]])])})

        means=sapply(1:npar2, function(i){if(i>= npar){(mean(unlist_chains[,i]))}else{mean(unlist_chains[,i])}})
    """)

    reval("""
        library(scales)
        plot(100000, axes = F, xlim = range(ltt_times), ylim = range(c(-1,y_max, 1.5*(sim_ltt - live_ltt))), xlab = "time", ylab = "n_lineage")
        y_lab = c(2,5,10,50,100,200,500,1000,2000,5000,10000,20000,50000)
        y_lab = y_lab[y_lab<=y_max]
        id_plot = unique(floor(seq(ini, n_row, length.out = n_ltt)))#unique(floor(seq(2,length(chains[[1]][[1]]), length.out = n_ltt)))
        id_ltt = (npar+3):length(chains[[1]])
        Ys = c()
        for (i in 1:3){
            for (j in id_plot){
                y = sapply(id_ltt, function(k){(chains[[i]][[k]][j])- live_ltt[k -npar-2]})
                Ys = rbind(Ys,y)
                lines(ltt_times, y, col = alpha("deepskyblue2", alpha = alpha_col), lwd = 2)
                }
            }
        quant = sapply(1:ncol(Ys), function(i){quantile(Ys[,i], probs = c(0.05,0.95))})
        n_out = sum(((sim_ltt - live_ltt) > quant[2,]) | ((sim_ltt - live_ltt) < quant[1,]))

        lines(ltt_times, quant[1,], col = "deepskyblue3", lwd = 2)
        lines(ltt_times, quant[2,], col = "deepskyblue3", lwd = 2)
        lines(ltt_times, sim_ltt - live_ltt, col = "orange1", lwd = 6, lty = 1)
        lines(ltt_times, sim_ltt - live_ltt, col = "orange3", lwd = 4, lty = 2)
        lines(ltt_times, maps[id_ltt] - live_ltt, col = "cadetblue1", lwd = 5, lty = 1)
        lines(ltt_times, maps[id_ltt] - live_ltt, col = "cadetblue4", lwd = 3, lty = 6)

        lines(ltt_times, means[id_ltt] - live_ltt, col = "darkseagreen1", lwd = 5, lty = 1)
        lines(ltt_times, means[id_ltt] - live_ltt, col = "darkseagreen4", lwd = 3, lty = 6)


        axis(1)
        axis(2)#, at = log(y_lab), lab = y_lab)
        title(main = c(cor((means[id_ltt] - live_ltt), (sim_ltt - live_ltt)),
            cor((maps[id_ltt] - live_ltt), (sim_ltt - live_ltt))))
        #print(log(maps[id_ltt]) - maps_log[id_ltt])
    """)
end

function plot_LTT_chain(chains, complete_tree::Tree, extant_tree::Tree, y_max::Int64, ltt_times::Array{Float64,1}, maps::Array{Float64,1} ; n_ltt = 100, alpha_col = 0.05, burn = 0.25, thin = 1)
    npar = extant_tree.n_nodes + 3 + (extant_tree.n_nodes+1)/2

    pruned_tree = root_tree(complete_tree)
    sim_ltt = LTT(pruned_tree, ltt_times)[2]
    new_tree = Tree(extant_tree.offsprings, 0., extant_tree.attributes)
    live_ltt = LTT(new_tree, ltt_times)[2]

    @rput chains
    @rput maps
    @rput npar
    @rput ltt_times
    @rput n_ltt
    @rput y_max
    @rput sim_ltt
    @rput alpha_col
    @rput burn
    @rput thin
    @rput live_ltt

    reval("""
        require(coda)
        n_row = length(chains[[1]][[1]])
        ini = floor(burn * n_row)
        if(ini == 0) ini = 1
        it = seq(ini, n_row, thin)
        npar2 = length(chains[[1]])
        map_chains = mcmc.list(lapply(1:3, function(i){
            mcmc(sapply(1:npar2, function(j){
                if(j<4 || (j>npar )){#}& j<(npar +3))) {
                    chains[[i]][[j]][it]
                }else{
                    log(chains[[i]][[j]][it])
                }
        }))}))
        unlist_chains = (sapply(1:npar2, function(k){
            c(map_chains[[1]][,k], map_chains[[2]][,k], map_chains[[3]][,k])
            }))

        #maps_log=sapply(1:npar2, function(i){D=density(unlist_chains[,i]);
        #    return(D[[1]][which.max(D[[2]])])})

        means=sapply(1:npar2, function(i){if(i>= npar +3){log(mean(unlist_chains[,i]))}else{mean(unlist_chains[,i])}})
    """)

    reval("""
        library(scales)
        plot(100000, axes = F, xlim = range(ltt_times), ylim = log(c(2,y_max)), xlab = "time", ylab = "n_lineage")
        y_lab = c(2,5,10,50,100,200,500,1000,2000,5000,10000,20000,50000)
        y_lab = y_lab[y_lab<=y_max]
        lines(ltt_times, log(live_ltt), col = "black", lwd = 6, lty = 1)

        id_plot = unique(floor(seq(ini, n_row, length.out = n_ltt)))#unique(floor(seq(2,length(chains[[1]][[1]]), length.out = n_ltt)))
        id_ltt = (npar+3):length(chains[[1]])
        Ys = c()
        for (i in 1:3){
            for (j in id_plot){
                y = sapply(id_ltt, function(k){log(chains[[i]][[k]][j])})
                Ys = rbind(Ys,y)
                lines(ltt_times, y, col = alpha("deepskyblue2", alpha = alpha_col), lwd = 2)
                }
            }
        quant = sapply(1:ncol(Ys), function(i){quantile(Ys[,i], probs = c(0.05,0.95))})

        lines(ltt_times, quant[1,], col = "deepskyblue3", lwd = 2)
        lines(ltt_times, quant[2,], col = "deepskyblue3", lwd = 2)
        lines(ltt_times, log(sim_ltt), col = "orange1", lwd = 6, lty = 1)
        lines(ltt_times, log(sim_ltt), col = "orange3", lwd = 4, lty = 2)
        lines(ltt_times, log(maps[id_ltt]), col = "cadetblue1", lwd = 5, lty = 1)
        lines(ltt_times, log(maps[id_ltt]), col = "cadetblue4", lwd = 3, lty = 6)

        lines(ltt_times, means[id_ltt], col = "darkseagreen1", lwd = 5, lty = 1)
        lines(ltt_times, means[id_ltt], col = "darkseagreen4", lwd = 3, lty = 6)

        axis(1)
        axis(2, at = log(y_lab), lab = y_lab)

        #print(log(maps[id_ltt]) - maps_log[id_ltt])
    """)

    #@rget means
    #@rget id_ltt
    #println(id_ltt)
    #return means[id_ltt]
end

function plot_LTT_chain_extinct(chains, complete_tree::Tree, extant_tree::Tree, y_max::T, ltt_times::Array{Float64,1}, maps::Array{Float64,1} ; n_ltt = 100, alpha_col = 0.05, burn = 0.25, thin = 1)where {T <: Number}

    npar = extant_tree.n_nodes + 3 + (extant_tree.n_nodes+1)/2

    pruned_tree = root_tree(complete_tree)
    sim_ltt = LTT(pruned_tree, ltt_times)[2]
    new_tree = Tree(extant_tree.offsprings, 0., extant_tree.attributes)
    live_ltt = LTT(new_tree, ltt_times)[2]

    @rput chains
    @rput maps
    @rput npar
    @rput ltt_times
    @rput n_ltt
    @rput y_max
    @rput sim_ltt
    @rput alpha_col
    @rput live_ltt
    @rput burn

    reval("""
        require(coda)
        n_row = length(chains[[1]][[1]])
        ini = floor(burn * n_row)
        if(ini == 0) ini = 1
        it = seq(ini, n_row, thin)
        npar2 = length(chains[[1]])
        map_chains = mcmc.list(lapply(1:3, function(i){
            mcmc(sapply(1:npar2, function(j){
                if(j<4 || (j>npar )){#}& j<(npar +3))) {
                    chains[[i]][[j]][it]
                }else{
                    log(chains[[i]][[j]][it])
                }
        }))}))
        unlist_chains = (sapply(1:npar2, function(k){
            c(map_chains[[1]][,k], map_chains[[2]][,k], map_chains[[3]][,k])
            }))

        #maps_log=sapply(1:npar2, function(i){D=density(unlist_chains[,i]);
        #    return(D[[1]][which.max(D[[2]])])})

        means=sapply(1:npar2, function(i){if(i>= npar){(mean(unlist_chains[,i]))}else{mean(unlist_chains[,i])}})
    """)

    reval("""
        library(scales)
        plot(100000, axes = F, xlim = range(ltt_times), ylim = range(c(-1,y_max, 1.5*(sim_ltt - live_ltt))), xlab = "time", ylab = "n_lineage")
        y_lab = c(2,5,10,50,100,200,500,1000,2000,5000,10000,20000,50000)
        y_lab = y_lab[y_lab<=y_max]
        id_plot = unique(floor(seq(ini, n_row, length.out = n_ltt)))#unique(floor(seq(2,length(chains[[1]][[1]]), length.out = n_ltt)))
        id_ltt = (npar+3):length(chains[[1]])
        Ys = c()
        for (i in 1:3){
            for (j in id_plot){
                y = sapply(id_ltt, function(k){(chains[[i]][[k]][j])- live_ltt[k -npar-2]})
                Ys = rbind(Ys,y)
                lines(ltt_times, y, col = alpha("deepskyblue2", alpha = alpha_col), lwd = 2)
                }
            }
        quant = sapply(1:ncol(Ys), function(i){quantile(Ys[,i], probs = c(0.05,0.95))})
        n_out = sum(((sim_ltt - live_ltt) > quant[2,]) | ((sim_ltt - live_ltt) < quant[1,]))

        lines(ltt_times, quant[1,], col = "deepskyblue3", lwd = 2)
        lines(ltt_times, quant[2,], col = "deepskyblue3", lwd = 2)
        lines(ltt_times, sim_ltt - live_ltt, col = "orange1", lwd = 6, lty = 1)
        lines(ltt_times, sim_ltt - live_ltt, col = "orange3", lwd = 4, lty = 2)
        lines(ltt_times, maps[id_ltt] - live_ltt, col = "cadetblue1", lwd = 5, lty = 1)
        lines(ltt_times, maps[id_ltt] - live_ltt, col = "cadetblue4", lwd = 3, lty = 6)

        lines(ltt_times, means[id_ltt] - live_ltt, col = "darkseagreen1", lwd = 5, lty = 1)
        lines(ltt_times, means[id_ltt] - live_ltt, col = "darkseagreen4", lwd = 3, lty = 6)


        axis(1)
        axis(2)#, at = log(y_lab), lab = y_lab)
        title(main = c(cor((means[id_ltt] - live_ltt), (sim_ltt - live_ltt)),
            cor((maps[id_ltt] - live_ltt), (sim_ltt - live_ltt))))
        #print(log(maps[id_ltt]) - maps_log[id_ltt])
    """)
end

function plot_LTT_chain(chains, sim_ltt::Array{Int64,1}, extant_tree::Tree, y_max::T, ltt_times::Array{Float64,1}, maps::Array{Float64,1} ; n_ltt = 100, alpha_col = 0.05, burn = 0.25, thin = 1)where {T <: Number}
    npar = extant_tree.n_nodes + 3 + (extant_tree.n_nodes+1)/2

    new_tree = Tree(extant_tree.offsprings, 0., extant_tree.attributes)
    live_ltt = LTT(new_tree, ltt_times)[2]

    @rput chains
    @rput maps
    @rput npar
    @rput ltt_times
    @rput n_ltt
    @rput y_max
    @rput sim_ltt
    @rput alpha_col
    @rput burn
    @rput thin
    @rput live_ltt

    reval("""
        require(coda)
        n_row = length(chains[[1]][[1]])
        ini = floor(burn * n_row)
        if(ini == 0) ini = 1
        it = seq(ini, n_row, thin)
        npar2 = length(chains[[1]])
        map_chains = mcmc.list(lapply(1:3, function(i){
            mcmc(sapply(1:npar2, function(j){
                if(j<4 || (j>npar )){#}& j<(npar +3))) {
                    chains[[i]][[j]][it]
                }else{
                    log(chains[[i]][[j]][it])
                }
        }))}))
        unlist_chains = (sapply(1:npar2, function(k){
            c(map_chains[[1]][,k], map_chains[[2]][,k], map_chains[[3]][,k])
            }))

        #maps_log=sapply(1:npar2, function(i){D=density(unlist_chains[,i]);
        #    return(D[[1]][which.max(D[[2]])])})

        #means=sapply(1:npar2, function(i){if(i>= npar +3){log(mean(unlist_chains[,i]))}else{mean(unlist_chains[,i])}})
        means=sapply(1:npar2, function(i){if(i>= npar +3){(mean(unlist_chains[,i]))}else{mean(unlist_chains[,i])}})
        means = log(means)
    """)

    reval("""
        library(scales)
        plot(100000, axes = F, xlim = range(ltt_times), ylim = log(c(2,y_max)), xlab = "time", ylab = "n_lineage")
        y_lab = c(2,5,10,50,100,200,500,1000,2000,5000,10000,20000,50000)
        y_lab = y_lab[y_lab<=y_max]
        lines(ltt_times, log(live_ltt), col = "black", lwd = 6, lty = 1)

        id_plot = unique(floor(seq(ini, n_row, length.out = n_ltt)))#unique(floor(seq(2,length(chains[[1]][[1]]), length.out = n_ltt)))
        id_ltt = (npar+3):length(chains[[1]])
        Ys = c()
        for (i in 1:3){
            for (j in id_plot){
                y = sapply(id_ltt, function(k){log(chains[[i]][[k]][j])})
                Ys = rbind(Ys,y)
                lines(ltt_times, y, col = alpha("deepskyblue2", alpha = alpha_col), lwd = 2)
                }
            }
        quant = sapply(1:ncol(Ys), function(i){quantile(Ys[,i], probs = c(0.05,0.95))})

        lines(ltt_times, quant[1,], col = "deepskyblue3", lwd = 2)
        lines(ltt_times, quant[2,], col = "deepskyblue3", lwd = 2)
        lines(ltt_times, log(sim_ltt), col = "orange1", lwd = 6, lty = 1)
        lines(ltt_times, log(sim_ltt), col = "orange3", lwd = 4, lty = 2)
        lines(ltt_times, log(maps[id_ltt]), col = "cadetblue1", lwd = 5, lty = 1)
        lines(ltt_times, log(maps[id_ltt]), col = "cadetblue4", lwd = 3, lty = 6)

        lines(ltt_times, means[id_ltt], col = "darkseagreen1", lwd = 5, lty = 1)
        lines(ltt_times, means[id_ltt], col = "darkseagreen4", lwd = 3, lty = 6)

        axis(1)
        axis(2, at = log(y_lab), lab = y_lab)

        #print(log(maps[id_ltt]) - maps_log[id_ltt])
    """)

end

function plot_LTT_chain_extinct(chains, sim_ltt::Array{Int64,1}, extant_tree::Tree, y_max::T, ltt_times::Array{Float64,1}, maps::Array{Float64,1} ; n_ltt = 100, alpha_col = 0.05, burn = 0.25, thin = 1)where {T <: Number}

    npar = extant_tree.n_nodes + 3 + (extant_tree.n_nodes+1)/2


    new_tree = Tree(extant_tree.offsprings, 0., extant_tree.attributes)
    live_ltt = LTT(new_tree, ltt_times)[2]

    @rput chains
    @rput maps
    @rput npar
    @rput ltt_times
    @rput n_ltt
    @rput y_max
    @rput sim_ltt
    @rput alpha_col
    @rput live_ltt
    @rput burn

    reval("""
        require(coda)
        n_row = length(chains[[1]][[1]])
        ini = floor(burn * n_row)
        if(ini == 0) ini = 1
        it = seq(ini, n_row, thin)
        npar2 = length(chains[[1]])
        map_chains = mcmc.list(lapply(1:3, function(i){
            mcmc(sapply(1:npar2, function(j){
                if(j<4 || (j>npar )){#}& j<(npar +3))) {
                    chains[[i]][[j]][it]
                }else{
                    log(chains[[i]][[j]][it])
                }
        }))}))
        unlist_chains = (sapply(1:npar2, function(k){
            c(map_chains[[1]][,k], map_chains[[2]][,k], map_chains[[3]][,k])
            }))

        #maps_log=sapply(1:npar2, function(i){D=density(unlist_chains[,i]);
        #    return(D[[1]][which.max(D[[2]])])})

        means=sapply(1:npar2, function(i){if(i>= npar){(mean(unlist_chains[,i]))}else{mean(unlist_chains[,i])}})
    """)

    reval("""
        library(scales)
        plot(100000, axes = F, xlim = range(ltt_times), ylim = range(c(-1,y_max, 1.5*(sim_ltt - live_ltt))), xlab = "time", ylab = "n_lineage")
        y_lab = c(2,5,10,50,100,200,500,1000,2000,5000,10000,20000,50000)
        y_lab = y_lab[y_lab<=y_max]
        id_plot = unique(floor(seq(ini, n_row, length.out = n_ltt)))#unique(floor(seq(2,length(chains[[1]][[1]]), length.out = n_ltt)))
        id_ltt = (npar+3):length(chains[[1]])
        Ys = c()
        for (i in 1:3){
            for (j in id_plot){
                y = sapply(id_ltt, function(k){(chains[[i]][[k]][j])- live_ltt[k -npar-2]})
                Ys = rbind(Ys,y)
                lines(ltt_times, y, col = alpha("deepskyblue2", alpha = alpha_col), lwd = 2)
                }
            }
        quant = sapply(1:ncol(Ys), function(i){quantile(Ys[,i], probs = c(0.05,0.95))})
        n_out = sum(((sim_ltt - live_ltt) > quant[2,]) | ((sim_ltt - live_ltt) < quant[1,]))

        lines(ltt_times, quant[1,], col = "deepskyblue3", lwd = 2)
        lines(ltt_times, quant[2,], col = "deepskyblue3", lwd = 2)
        lines(ltt_times, sim_ltt - live_ltt, col = "orange1", lwd = 6, lty = 1)
        lines(ltt_times, sim_ltt - live_ltt, col = "orange3", lwd = 4, lty = 2)
        lines(ltt_times, maps[id_ltt] - live_ltt, col = "cadetblue1", lwd = 5, lty = 1)
        lines(ltt_times, maps[id_ltt] - live_ltt, col = "cadetblue4", lwd = 3, lty = 6)

        lines(ltt_times, means[id_ltt] - live_ltt, col = "darkseagreen1", lwd = 5, lty = 1)
        lines(ltt_times, means[id_ltt] - live_ltt, col = "darkseagreen4", lwd = 3, lty = 6)


        axis(1)
        axis(2)#, at = log(y_lab), lab = y_lab)
        title(main = c(cor((means[id_ltt] - live_ltt), (sim_ltt - live_ltt)),
            cor((maps[id_ltt] - live_ltt), (sim_ltt - live_ltt))))
        #print(log(maps[id_ltt]) - maps_log[id_ltt])
    """)
end

function plot_LTT_chain(mcmc_sampler, extant_tree::Tree, y_max::T ; n_ltt = 100, alpha_col = 0.05, burn = 0.25, thin = 1) where {T <: Number}
    chains = mcmc_sampler[1][1]
    npar = mcmc_sampler[1][4][1].n_nodes + 3 + (mcmc_sampler[1][4][1].n_nodes+1)/2
    ltt_times = mcmc_sampler[1][6]
    maps = mcmc_sampler[2]

    new_tree = Tree(extant_tree.offsprings, 0., extant_tree.attributes)
    live_ltt = LTT(new_tree, ltt_times)[2]

    @rput chains
    @rput maps
    @rput npar
    @rput ltt_times
    @rput n_ltt
    @rput y_max
    @rput alpha_col
    @rput burn
    @rput thin
    @rput live_ltt

    reval("""
        require(coda)
        n_row = length(chains[[1]][[1]])
        ini = floor(burn * n_row)
        if(ini == 0) ini = 1
        it = seq(ini, n_row, thin)
        npar2 = length(chains[[1]])
        map_chains = mcmc.list(lapply(1:3, function(i){
            mcmc(sapply(1:npar2, function(j){
                if(j<4 || (j>npar )){#}& j<(npar +3))) {
                    chains[[i]][[j]][it]
                }else{
                    log(chains[[i]][[j]][it])
                }
        }))}))
        unlist_chains = (sapply(1:npar2, function(k){
            c(map_chains[[1]][,k], map_chains[[2]][,k], map_chains[[3]][,k])
            }))

        #maps_log=sapply(1:npar2, function(i){D=density(unlist_chains[,i]);
        #    return(D[[1]][which.max(D[[2]])])})

        means=sapply(1:npar2, function(i){if(i>= npar +3){log(mean(unlist_chains[,i]))}else{mean(unlist_chains[,i])}})
    """)

    reval("""
        library(scales)
        plot(100000, axes = F, xlim = range(ltt_times), ylim = log(c(2,y_max)), xlab = "time", ylab = "n_lineage")
        y_lab = c(2,5,10,50,100,200,500,1000,2000,5000,10000,20000,50000)
        y_lab = y_lab[y_lab<=y_max]
        lines(ltt_times, log(live_ltt), col = "black", lwd = 6, lty = 1)

        id_plot = unique(floor(seq(ini, n_row, length.out = n_ltt)))#unique(floor(seq(2,length(chains[[1]][[1]]), length.out = n_ltt)))
        id_ltt = (npar+3):length(chains[[1]])
        Ys = c()
        for (i in 1:3){
            for (j in id_plot){
                y = sapply(id_ltt, function(k){log(chains[[i]][[k]][j])})
                Ys = rbind(Ys,y)
                lines(ltt_times, y, col = alpha("deepskyblue2", alpha = alpha_col), lwd = 2)
                }
            }
        quant = sapply(1:ncol(Ys), function(i){quantile(Ys[,i], probs = c(0.05,0.95))})

        lines(ltt_times, quant[1,], col = "deepskyblue3", lwd = 2)
        lines(ltt_times, quant[2,], col = "deepskyblue3", lwd = 2)
        lines(ltt_times, log(maps[id_ltt]), col = "cadetblue1", lwd = 5, lty = 1)
        lines(ltt_times, log(maps[id_ltt]), col = "cadetblue4", lwd = 3, lty = 6)

        lines(ltt_times, means[id_ltt], col = "darkseagreen1", lwd = 5, lty = 1)
        lines(ltt_times, means[id_ltt], col = "darkseagreen4", lwd = 3, lty = 6)

        axis(1)
        axis(2, at = log(y_lab), lab = y_lab)

        #print(log(maps[id_ltt]) - maps_log[id_ltt])
    """)
end

function plot_LTT_chain(chains, extant_tree::Tree, y_max::T, maps; n_ltt = 100, alpha_col = 0.05, burn = 0.25, thin = 1) where {T <: Number}
    npar = extant_tree.n_nodes + 3 + (extant_tree.n_nodes+1)/2
    ltt_steps = Int64(length(chains[1]) - npar-3 )
    new_tree = Tree(extant_tree.offsprings, 0., extant_tree.attributes)
    root_depth = maximum(node_depths(new_tree)) * 1.01
    ltt_times = [0:ltt_steps...] * root_depth/ltt_steps
    println(ltt_steps)
    live_ltt = LTT(new_tree, ltt_times)[2]

    @rput chains
    @rput maps
    @rput npar
    @rput ltt_times
    @rput n_ltt
    @rput y_max
    @rput alpha_col
    @rput burn
    @rput thin
    @rput live_ltt

    reval("""
        require(coda)
        n_row = length(chains[[1]][[1]])
        ini = floor(burn * n_row)
        if(ini == 0) ini = 1
        it = seq(ini, n_row, thin)
        npar2 = length(chains[[1]])
        map_chains = mcmc.list(lapply(1:3, function(i){
            mcmc(sapply(1:npar2, function(j){
                if(j<4 || (j>npar )){#}& j<(npar +3))) {
                    chains[[i]][[j]][it]
                }else{
                    log(chains[[i]][[j]][it])
                }
        }))}))
        unlist_chains = (sapply(1:npar2, function(k){
            c(map_chains[[1]][,k], map_chains[[2]][,k], map_chains[[3]][,k])
            }))

        #maps_log=sapply(1:npar2, function(i){D=density(unlist_chains[,i]);
        #    return(D[[1]][which.max(D[[2]])])})

        means=sapply(1:npar2, function(i){if(i>= npar +3){log(mean(unlist_chains[,i]))}else{mean(unlist_chains[,i])}})
    """)

    reval("""
        library(scales)
        plot(100000, axes = F, xlim = range(ltt_times), ylim = log(c(2,y_max)), xlab = "time", ylab = "n_lineage")
        y_lab = c(2,5,10,50,100,200,500,1000,2000,5000,10000,20000,50000)
        y_lab = y_lab[y_lab<=y_max]
        lines(ltt_times, log(live_ltt), col = "black", lwd = 6, lty = 1)

        id_plot = unique(floor(seq(ini, n_row, length.out = n_ltt)))#unique(floor(seq(2,length(chains[[1]][[1]]), length.out = n_ltt)))
        id_ltt = (npar+3):length(chains[[1]])
        Ys = c()
        for (i in 1:3){
            for (j in id_plot){
                y = sapply(id_ltt, function(k){log(chains[[i]][[k]][j])})
                Ys = rbind(Ys,y)
                lines(ltt_times, y, col = alpha("deepskyblue2", alpha = alpha_col), lwd = 2)
                }
            }
        quant = sapply(1:ncol(Ys), function(i){quantile(Ys[,i], probs = c(0.05,0.95))})

        lines(ltt_times, quant[1,], col = "deepskyblue3", lwd = 2)
        lines(ltt_times, quant[2,], col = "deepskyblue3", lwd = 2)
        lines(ltt_times, log(maps[id_ltt]), col = "cadetblue1", lwd = 5, lty = 1)
        lines(ltt_times, log(maps[id_ltt]), col = "cadetblue4", lwd = 3, lty = 6)

        lines(ltt_times, means[id_ltt], col = "darkseagreen1", lwd = 5, lty = 1)
        lines(ltt_times, means[id_ltt], col = "darkseagreen4", lwd = 3, lty = 6)

        axis(1)
        axis(2, at = log(y_lab), lab = y_lab)

        #print(log(maps[id_ltt]) - maps_log[id_ltt])
    """)
end

function plot_LTT_chain_extinct(mcmc_sampler, extant_tree::Tree, y_max::T ; n_ltt = 100, alpha_col = 0.05, burn = 0.25, thin = 1) where {T <: Number}
    chains = mcmc_sampler[1][1]
    npar = mcmc_sampler[1][4][1].n_nodes + 3 + (mcmc_sampler[1][4][1].n_nodes+1)/2
    ltt_times = mcmc_sampler[1][6]
    maps = mcmc_sampler[2]

    new_tree = Tree(extant_tree.offsprings, 0., extant_tree.attributes)
    live_ltt = LTT(new_tree, ltt_times)[2]

    @rput chains
    @rput maps
    @rput npar
    @rput ltt_times
    @rput n_ltt
    @rput y_max
    @rput alpha_col
    @rput live_ltt
    @rput burn

    reval("""
        require(coda)
        n_row = length(chains[[1]][[1]])
        ini = floor(burn * n_row)
        if(ini == 0) ini = 1
        it = seq(ini, n_row, thin)
        npar2 = length(chains[[1]])
        map_chains = mcmc.list(lapply(1:3, function(i){
            mcmc(sapply(1:npar2, function(j){
                if(j<4 || (j>npar )){#}& j<(npar +3))) {
                    chains[[i]][[j]][it]
                }else{
                    log(chains[[i]][[j]][it])
                }
        }))}))
        unlist_chains = (sapply(1:npar2, function(k){
            c(map_chains[[1]][,k], map_chains[[2]][,k], map_chains[[3]][,k])
            }))

        #maps_log=sapply(1:npar2, function(i){D=density(unlist_chains[,i]);
        #    return(D[[1]][which.max(D[[2]])])})

        means=sapply(1:npar2, function(i){if(i>= npar){(mean(unlist_chains[,i]))}else{mean(unlist_chains[,i])}})
    """)

    reval("""
        library(scales)
        plot(100000, axes = F, xlim = range(ltt_times), ylim = c(-1,y_max), xlab = "time", ylab = "n_lineage")
        y_lab = c(2,5,10,50,100,200,500,1000,2000,5000,10000,20000,50000)
        y_lab = y_lab[y_lab<=y_max]
        id_plot = unique(floor(seq(ini, n_row, length.out = n_ltt)))#unique(floor(seq(2,length(chains[[1]][[1]]), length.out = n_ltt)))
        id_ltt = (npar+3):length(chains[[1]])
        Ys = c()
        for (i in 1:3){
            for (j in id_plot){
                y = sapply(id_ltt, function(k){(chains[[i]][[k]][j])- live_ltt[k -npar-2]})
                Ys = rbind(Ys,y)
                lines(ltt_times, y, col = alpha("deepskyblue2", alpha = alpha_col), lwd = 2)
                }
            }
        quant = sapply(1:ncol(Ys), function(i){quantile(Ys[,i], probs = c(0.05,0.95))})

        lines(ltt_times, quant[1,], col = "deepskyblue3", lwd = 2)
        lines(ltt_times, quant[2,], col = "deepskyblue3", lwd = 2)
        lines(ltt_times, maps[id_ltt] - live_ltt, col = "cadetblue1", lwd = 5, lty = 1)
        lines(ltt_times, maps[id_ltt] - live_ltt, col = "cadetblue4", lwd = 3, lty = 6)

        lines(ltt_times, means[id_ltt] - live_ltt, col = "darkseagreen1", lwd = 5, lty = 1)
        lines(ltt_times, means[id_ltt] - live_ltt, col = "darkseagreen4", lwd = 3, lty = 6)


        axis(1)
        axis(2)#, at = log(y_lab), lab = y_lab)
    """)
end
