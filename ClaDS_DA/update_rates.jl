function update_rates(tree::Tree, rates::Array{T,1} ; id = 1) where {T<:Number}
    function aux(sub_tree, sub_rates)
        if sub_tree.n_nodes < 2
            root_rate = popfirst!(sub_rates)
            current_rates = sub_tree.attributes
            current_rates[id] = root_rate
            return Tree(sub_tree.offsprings, sub_tree.branch_length, current_rates, sub_tree.n_nodes, sub_tree.label)
        else
            root_rate = popfirst!(sub_rates)
            t1 = aux(sub_tree.offsprings[1], sub_rates)
            t2 = aux(sub_tree.offsprings[2], sub_rates)
            current_rates = sub_tree.attributes
            current_rates[id] = root_rate
            return Tree([t1,t2], sub_tree.branch_length, current_rates, sub_tree.n_nodes, sub_tree.label)
        end
    end

    new_rates = deepcopy(rates)
    aux(tree, new_rates)
end

function update_rates!(tree::Tree, rates::Array{T,1} ; id = 1) where {T<:Number}
    function aux(sub_tree, sub_rates)
        if sub_tree.n_nodes < 2
            root_rate = popfirst!(sub_rates)
            sub_tree.attributes[id] = root_rate
        else
            root_rate = popfirst!(sub_rates)
            sub_tree.attributes[id] = root_rate
            t1 = sub_tree.offsprings[1]
            aux(t1, sub_rates)
            t2 = sub_tree.offsprings[2]
            aux(t2, sub_rates)
            #return Tree([t1,t2], sub_tree.branch_length, root_rate, sub_tree.n_nodes)
        end
    end

    new_rates = deepcopy(rates)
    aux(tree, new_rates)
end

function update_rates(tree::Tree, rates::Array{T,1}, edge_trees ; id = 1) where {T<:Number}
    tree = update_rates(tree,rates, id = id)

    for i in 1:length(edge_trees)
        edge_trees[i].tree.attributes[id] = rates[i+1]
    end

    return tree, edge_trees
end

function update_rates!(tree::Tree, rates::Array{T,1}, edge_trees ; id = 1) where {T<:Number}
    update_rates!(tree,rates, id = id)

    for i in 1:length(edge_trees)
        edge_trees[i].tree.attributes[id] = rates[i+1]
    end

    return tree, edge_trees
end

function update_rates(tree::Tree, rate::T ; id = id) where {T<:Number}
    function aux(sub_tree)
        current_rate = sub_tree.attributes
        current_rate[id] *= rate
        if sub_tree.n_nodes < 2
            return Tree(sub_tree.offsprings, sub_tree.branch_length, current_rate, sub_tree.n_nodes, sub_tree.extant)
        else
            return Tree([aux(sub_tree.offsprings[1]),aux(sub_tree.offsprings[2])], sub_tree.branch_length, current_rate, sub_tree.n_nodes, sub_tree.extant)
        end
    end

    aux(tree)
end

function update_rates!(tree::Tree, rate::T ; id = 1) where {T<:Number}
    function aux(sub_tree)
        if sub_tree.n_nodes < 2
            sub_tree.attributes[id] *= rate
        else
            sub_tree.attributes[id] *= rate
            t1 = sub_tree.offsprings[1]
            aux(t1)
            t2 = sub_tree.offsprings[2]
            aux(t2)
        end
    end

    aux(tree)
end

function update_rates(tree::Tree, rate::T, edge_trees::Union{Array{EdgeTree,1},Array{EdgeTreeRates,1}} ; id = 1) where {T<:Number}
    tree = update_rates(tree,rate, id = id)

    for i in 1:length(edge_trees)
        subtree = edge_trees[i].tree
        tip = edge_trees[i].tip_id
        n = edge_trees[i].n
        subtree = update_rates(subtree,rate, id = id)
        edge_trees[i] = EdgeTree(subtree,tip, n)
    end

    return tree, edge_trees
end
