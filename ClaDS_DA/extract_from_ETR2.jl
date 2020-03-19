function extract_S(tree::Tree ,edge_trees::Array{EdgeTreeRates2,1})

    S::Float64 = 0.
    for et in edge_trees
        S += et.stem_rate[1] * et.effective_bl
    end

    return S

end

function extract_nextinct(tree::Tree ,edge_trees::Array{EdgeTreeRates2,1})

    n::Int64 = 2 - n_tip(tree)
    i=0
    for et in edge_trees
        n += n_tip(et.tree) - et.tip_number
    end

    return n
end


function n_tip(tree::Tree ,edge_trees::Array{EdgeTreeRates2,1})
    n = n_tip(tree)
    for i in 1:length(edge_trees)
        n += n_tip(edge_trees[i].tree) - 1
    end
    return n
end

function n_extant_tips(tree::Tree ,edge_trees::Array{EdgeTreeRates2,1})

    n::Int64 = 0

    for et in edge_trees
        n += et.tip_number
    end

    return n
end

function extract_relative_rates(tree::Tree, edge_trees::Array{EdgeTreeRates2,1}, rates::Array{Float64,1}; id = 1)

    function aux(subtree, parent_rate, x)
        if subtree.n_nodes > 0
            if length(subtree.offsprings) == 0
                pushfirst!(x,log(subtree.attributes[id]) - parent_rate)
            else
                log_rate = deepcopy(log(subtree.attributes[id]))
                aux(subtree.offsprings[2], log_rate, x)
                aux(subtree.offsprings[1], log_rate, x)
                pushfirst!(x,log(subtree.attributes[id]) - parent_rate)
            end
        end
    end

    relative_rates = Array{Float64,1}(undef,0)

    i = 0
    for et in edge_trees
        i+= 1
        parent_edge = et.parent_edge

        parent_rate = rates[1]
        if parent_edge > 0
            parent_rate = edge_trees[parent_edge].tip_rate * rates[parent_edge + 1]
        end

        r = log(rates[i + 1]) - log(parent_rate)
        pushfirst!(relative_rates,r)
        if length(et.tree.offsprings) > 0
            aux(et.tree.offsprings[1], 0., relative_rates)
            aux(et.tree.offsprings[2], 0., relative_rates)
        end
    end
    return relative_rates
end

function extract_relative_rates(tree::Tree, edge_trees::Array{EdgeTreeRates2,1}; id = 1)
    function aux(subtree, parent_rate, x)
        if subtree.n_nodes > 0
            if length(subtree.offsprings) == 0
                pushfirst!(x,log(subtree.attributes[id]) - parent_rate)
            else
                log_rate = deepcopy(log(subtree.attributes[id]))
                aux(subtree.offsprings[2], log_rate, x)
                aux(subtree.offsprings[1], log_rate, x)
                pushfirst!(x,log(subtree.attributes[id]) - parent_rate)
            end
        end
    end

    relative_rates = Array{Float64,1}(undef,0)

    for et in edge_trees
        parent_edge = et.parent_edge

        parent_rate = tree.attributes[1]
        if parent_edge > 0
            parent_rate = edge_trees[parent_edge].tip_rate * edge_trees[parent_edge].stem_rate[1]
        end


        r = log(et.stem_rate[1]) - log(parent_rate)
        pushfirst!(relative_rates,r)

        if length(et.tree.offsprings) > 10000000
            aux(et.tree.offsprings[1], 0., relative_rates)
            aux(et.tree.offsprings[2], 0., relative_rates)
        end
    end

    return relative_rates
end

function extract_partial_relative_rates(tree::Tree, edge_trees::Array{EdgeTreeRates2,1}; id = 1)

    relative_rates = Array{Float64,1}(undef,0)

    for et in edge_trees
        parent_edge = et.parent_edge

        parent_rate = tree.attributes[1]
        if parent_edge > 0
            parent_rate = edge_trees[parent_edge].tip_rate * edge_trees[parent_edge].stem_rate[1]
        end


        r = log(et.stem_rate[1]) - log(parent_rate)
        pushfirst!(relative_rates,r)
    end

    return relative_rates
end

function extract_tip_rates(tree::Tree, edge_trees::Array{EdgeTreeRates2,1},
    tips_id::Array{Bool,1}, rates::Array{Float64,1} ; id = 1, return_extinct = true)

    tip_rates = Array{Float64,1}(undef,0)
    for i in 1:length(edge_trees)
        if tips_id[i]
            r = sample(edge_trees[i].tip_rates) * edge_trees[i].stem_rate[1]
            push!(tip_rates,r)
        end
    end
    return tip_rates
end

function node_depths_base(tree::Tree, edge_trees::Array{EdgeTreeRates2,1})
    function aux(subtree, node_times, current_depth, att)
        if subtree.n_nodes < 2
            pushfirst!(node_times, current_depth)
        else
            aux(subtree.offsprings[2],node_times, current_depth + subtree.branch_length/att, att)
            aux(subtree.offsprings[1],node_times, current_depth + subtree.branch_length/att, att)
            pushfirst!(node_times,current_depth)
        end
    end

    node_times = node_depths_base(tree)[2:end]
    node_times_full = Array{Float64,1}(undef,0)
    for i in length(edge_trees):-1:1
        aux(edge_trees[i].tree, node_times_full, node_times[i], edge_trees[i].scale)
    end
    #pushfirst!(node_times,0.)
    return node_times_full
end


function node_depths(tree::Tree, edge_trees::Array{EdgeTreeRates2,1})
    function aux(subtree, node_times, current_depth, sa)
        if subtree.n_nodes < 2
            pushfirst!(node_times, current_depth + subtree.branch_length/sa)
        else
            aux(subtree.offsprings[2],node_times, current_depth + subtree.branch_length/sa, sa)
            aux(subtree.offsprings[1],node_times, current_depth + subtree.branch_length/sa, sa)
            pushfirst!(node_times,current_depth+ subtree.branch_length/sa)
        end
    end

    node_times = node_depths_base(tree)[2:end]
    node_times_full = Array{Float64,1}(undef,0)
    for i in length(edge_trees):-1:1
        aux(edge_trees[i].tree, node_times_full, node_times[i], edge_trees[i].scale)
    end
    #pushfirst!(node_times,0.)
    return node_times_full
end

function extract_rates(tree::Tree ,edge_trees::Array{EdgeTreeRates2,1}; id = 1)
    function aux(subtree, rates,ri)
        if subtree.n_nodes == 0
            #return []
        elseif length(subtree.offsprings) == 0
            pushfirst!(rates,subtree.attributes[id]*ri)
        else
            aux(subtree.offsprings[2],rates,ri)
            aux(subtree.offsprings[1],rates,ri)
            pushfirst!(rates,subtree.attributes[id]*ri)
        end
    end

    rates =  Array{Float64,1}(undef,0)
    for i in length(edge_trees):-1:1
        aux(edge_trees[i].tree, rates, edge_trees[i].stem_rate[1])
    end
    pushfirst!(rates, tree.attributes[1])
    return rates
end
