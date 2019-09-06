struct Tree
    offsprings::Array{Tree,1}
    branch_length::Float64
    attributes::Array{T,1} where {T<:Number}
    n_nodes::Int64
    extant::Bool
    label::String
end

function Tree(offsprings::Array{Tree,1}, branch_length::Float64, attributes::T, n_nodes::Int64, extant::Bool, label::String) where {T<:Number}
    Tree(offsprings, branch_length, [attributes], n_nodes, extant, label)
end

function Tree(offsprings::Array{Tree,1}, branch_length::Float64, attributes::T, n_nodes::Int64, extant::Bool) where {T<:Number}
    Tree(offsprings, branch_length, attributes, n_nodes, extant, "")
end

function Tree(offsprings::Array{Tree,1}, branch_length::Float64, attributes::Array{T,1}, n_nodes::Int64, extant::Bool) where {T<:Number}
    Tree(offsprings, branch_length, attributes, n_nodes, extant, "")
end

function Tree(offsprings::Array{Tree,1}, branch_length::Float64, attributes::Array{T,1}, n_nodes::Int64) where {T<:Number}
    Tree(offsprings, branch_length, attributes, n_nodes, true)
end

function Tree(offsprings::Array{Tree,1}, branch_length::Float64, attributes::T, n_nodes::Int64) where {T<:Number}
    Tree(offsprings, branch_length, [attributes], n_nodes, true)
end

function Tree(offsprings::Array{Tree,1}, branch_length::Float64, attributes::Array{T,1}, n_nodes::Int64, label::String) where {T<:Number}
    Tree(offsprings, branch_length, attributes, n_nodes, true, label)
end

function Tree(offsprings::Array{Tree,1}, branch_length::Float64, attributes::T, n_nodes::Int64, label::String) where {T<:Number}
    Tree(offsprings, branch_length, [attributes], n_nodes, true, label)
end

function Tree(offsprings::Array{Tree,1}, branch_length::Float64, attributes)
    n_nodes = length(branch_length)
    extant = false
    for i in offsprings
        n_nodes += i.n_nodes
        extant = extant || i.extant
    end
    extant = extant || n_nodes == 1
    Tree(offsprings, branch_length, attributes, n_nodes, extant)
end

function Tree(offsprings::Array{Tree,1}, branch_length::Float64, attributes, label::String)
    #print("here")
    n_nodes = length(branch_length)
    extant = false
    for i in offsprings
        n_nodes += i.n_nodes
        extant = extant || i.extant
    end
    extant = extant || n_nodes == 1
    Tree(offsprings, branch_length, attributes, n_nodes, extant, label)
end

function Tree(offsprings::Array{Tree,1}, branch_length::Float64, attributes, extant::Bool)
    n_nodes = length(branch_length)
    for i in offsprings
        n_nodes += i.n_nodes
    end
    Tree(offsprings, branch_length, attributes, n_nodes, extant)
end

function Tree(offsprings::Array{Tree,1}, branch_length::Float64, attributes, extant::Bool, label::String)
    n_nodes = length(branch_length)
    for i in offsprings
        n_nodes += i.n_nodes
    end
    Tree(offsprings, branch_length, attributes, n_nodes, extant, label)
end

function Tree()
    Tree(Array{Tree,1}(undef,0), 0., Array{Float64,1}(undef,0), 0, false)
end

struct EdgeTree
    tree::Tree
    tip_id::Int64
    n::Int64
end

struct EdgeTreeRates
    tree::Tree
    tip_id::Int64 #for tip edges it is the number of alive species
    rate::Float64
    rates::Array{Float64,1}
end

function rasterize(tree::Tree)
    function aux(subtree)
        if length(subtree.offsprings) == 0
            return subtree
        else
            left = aux(subtree.offsprings[1])
            right = aux(subtree.offsprings[2])

            if left.n_nodes < right.n_nodes
                return Tree([left,right], subtree.branch_length,subtree.attributes,subtree.n_nodes,subtree.extant,subtree.label)
            else
                return Tree([right,left], subtree.branch_length,subtree.attributes,subtree.n_nodes,subtree.extant,subtree.label)
            end
        end
    end

    return aux(tree)
end

function n_tip(tree::Tree)
    return Int64((tree.n_nodes+1)/2)
end

function n_tip(tree::Tree ,edge_trees::Union{Array{EdgeTree,1},Array{EdgeTreeRates,1}})
    n = n_tip(tree)
    for i in 1:length(edge_trees)
        n += n_tip(edge_trees[i].tree) - 1
    end
    return n
end

function n_extant_tips(tree::Tree)
    function aux(subtree)
        if subtree.n_nodes < 2
            if subtree.extant
                return 1
            else
                return  0
            end
        else
            return aux(subtree.offsprings[1]) + aux(subtree.offsprings[2])
        end
    end

    aux(tree)
end

function n_extant_tips(tree::Tree ,edge_trees::Union{Array{EdgeTree,1},Array{EdgeTreeRates,1}})
    n = n_extant_tips(tree)
    for i in 1:length(edge_trees)
        n += n_extant_tips(edge_trees[i].tree) - 1
    end
    return n
end

function n_extinct(tree::Tree)
    n_tip(tree) - n_extant_tips(tree)
end

function n_extinct(tree::Tree ,edge_trees::Union{Array{EdgeTree,1},Array{EdgeTreeRates,1}})
    n = 0
    for i in 1:length(edge_trees)
        n += n_extinct(edge_trees[i].tree)
    end
    return n
end

function is_tip(tree::Tree)
    length(tree.offsprings) == 0
end

function tips(tree::Tree)
    function aux(subtree, x)
        if length(subtree.offsprings) == 0
            pushfirst!(x,true)
        else
            aux(subtree.offsprings[2], x)
            aux(subtree.offsprings[1], x)
            pushfirst!(x,false)
        end
    end

    x = Array{Bool,1}(undef,0)
    aux(tree, x)
    return x
end

function extant_tips(tree::Tree)
    function aux(subtree, x)
        if length(subtree.offsprings) == 0 && subtree.extant
            pushfirst!(x,true)
        elseif length(subtree.offsprings) == 0
            pushfirst!(x,false)
        else
            aux(subtree.offsprings[2], x)
            aux(subtree.offsprings[1], x)
            pushfirst!(x,false)
        end
    end

    x = Array{Bool,1}(undef,0)
    aux(tree, x)
    return x
end

function n_left(tree::Tree)
    function aux(subtree, x)
        if length(subtree.offsprings) == 0
            pushfirst!(x,0)
        else
            aux(subtree.offsprings[2], x)
            aux(subtree.offsprings[1], x)
            pushfirst!(x,subtree.offsprings[1].n_nodes)
        end
    end

    x = Array{Int64,1}(undef,0)
    aux(tree, x)
    return x
end

function extract_rates(tree::Tree; id=1)
    function aux(subtree, rates)
        if subtree.n_nodes == 0
            #return []
        elseif length(subtree.offsprings) == 0
            pushfirst!(rates,subtree.attributes[id])
        else
            aux(subtree.offsprings[2],rates)
            aux(subtree.offsprings[1],rates)
            pushfirst!(rates,subtree.attributes[id])
        end
    end

    rates = Array{Float64,1}(undef,0)
    aux(tree, rates)
    return(rates)
end

function extract_rates(tree::Tree ,edge_trees::Union{Array{EdgeTree,1},Array{EdgeTreeRates,1}}; id = 1)
    function aux(subtree, rates)
        if subtree.n_nodes == 0
            #return []
        elseif length(subtree.offsprings) == 0
            pushfirst!(rates,subtree.attributes[id])
        else
            aux(subtree.offsprings[2],rates)
            aux(subtree.offsprings[1],rates)
            pushfirst!(rates,subtree.attributes[id])
        end
    end

    rates =  Array{Float64,1}(undef,0)
    for i in length(edge_trees):-1:1
        aux(edge_trees[i].tree, rates)
    end
    pushfirst!(rates, tree.attributes[1])
    return rates
end

function extract_branch_lengths(tree::Tree)
    function aux(subtree, bl)
        if subtree.n_nodes == 0
            #return []
        elseif length(subtree.offsprings) == 0
            pushfirst!(bl,subtree.branch_length)
        else
            aux(subtree.offsprings[2],bl)
            aux(subtree.offsprings[1],bl)
            pushfirst!(bl,subtree.branch_length)
        end
    end

    bl =  Array{Float64,1}(undef,0)
    aux(tree, bl)
    return bl
end

function extract_branch_lengths(tree::Tree ,edge_trees::Union{Array{EdgeTree,1},Array{EdgeTreeRates,1}})

    function aux(subtree, bl)
        if subtree.n_nodes == 0
            #return []
        elseif length(subtree.offsprings) == 0
            pushfirst!(bl,subtree.branch_length)
        else
            aux(subtree.offsprings[2],bl)
            aux(subtree.offsprings[1],bl)
            pushfirst!(bl,subtree.branch_length)
        end
    end

    bl =  Array{Float64,1}(undef,0)
    for i in length(edge_trees):-1:1
        aux(edge_trees[i].tree,bl)
    end
    pushfirst!(bl, tree.branch_length)
    return bl
end

function which_extant(tree::Tree)
    function aux(subtree, x)
        if subtree.n_nodes == 0
        elseif length(subtree.offsprings) == 0
            if subtree.extant
                pushfirst!(x, 1)
            else
                pushfirst!(x, 0)
            end
        else
            aux(subtree.offsprings[2],x)
            aux(subtree.offsprings[1],x)
            if subtree.extant
                pushfirst!(extant, 1)
            else
                pushfirst!(extant, 0)
            end
        end
    end

    x = Array{Int64,1}(undef,0)
    aux(tree, x)
    return x
end

function get_parent_edge(tree::Tree, edge_id::Int64)
    if edge_id == 0
        return NaN
    end

    function aux(subtree, edge)
        n_edges_left = subtree.offsprings[1].n_nodes + 1
        if edge == 2 || edge == (n_edges_left + 1)
            return 1
        elseif edge > (n_edges_left + 1)
            return n_edges_left + aux(subtree.offsprings[2], edge - n_edges_left)
        else
            return 1 + aux(subtree.offsprings[1], edge - 1)
        end
    end

    aux(tree, edge_id + 1) - 1
end

function get_parent_edges(tree::Tree)
    parents = Array{Int64,1}(undef,0)
    for i in 1:(tree.n_nodes-1)
        push!(parents,get_parent_edge(tree, i))
    end
    return parents
end

function get_parent_rate(tree::Tree, edge_id::Int64, edge_trees::Union{Array{EdgeTree,1},Array{EdgeTreeRates,1}})
    parent_edge = get_parent_edge(tree, edge_id)
    if parent_edge == 0
        return tree.attributes[1]
    else
        if edge_trees[parent_edge].tree.n_nodes < 2
            return edge_trees[parent_edge].tree.attributes[1]
        else
            return extract_tip_rates(edge_trees[parent_edge].tree, edge_trees[parent_edge].tip_id)
            #tip_rates = extract_tip_rates(edge_trees[parent_edge].tree)
            #return tip_rates[edge_trees[parent_edge].tip_id]
        end
    end
end

function get_parent_rate(tree::Tree, edge_id::Int64, edge_trees::Union{Array{EdgeTree,1},Array{EdgeTreeRates,1}}, rates::Array{Float64,1})
    parent_edge = get_parent_edge(tree, edge_id)
    if parent_edge == 0
        return rates[1]
    else
        if edge_trees[parent_edge].tree.n_nodes < 2
            return rates[parent_edge+1]#
        else
            return extract_tip_rates(edge_trees[parent_edge].tree, edge_trees[parent_edge].tip_id)
            #tip_rates = extract_tip_rates(edge_trees[parent_edge].tree)
            #return tip_rates[edge_trees[parent_edge].tip_id]
        end
    end
end

function get_parent_rate(tree::Tree, edge_trees::Union{Array{EdgeTree,1},Array{EdgeTreeRates,1}}, rates::Array{Float64,1}, parent_edge::Int64)
    if parent_edge == 0
        return rates[1]
    else
        if edge_trees[parent_edge].tree.n_nodes < 2
            return rates[parent_edge+1]#
        else
            return extract_tip_rates(edge_trees[parent_edge].tree, edge_trees[parent_edge].tip_id)
            #tip_rates = extract_tip_rates(edge_trees[parent_edge].tree)
            #return tip_rates[edge_trees[parent_edge].tip_id]
        end
    end
end

function get_parent_rates(tree::Tree; id = 1)
    function aux(subtree, root_rate, rates)
        if subtree.n_nodes == 0
            #return []
        elseif length(subtree.offsprings) == 0
            pushfirst!(rates,root_rate)
        else
            new_root = deepcopy(subtree.attributes[id])
            aux(subtree.offsprings[2],new_root,rates)
            aux(subtree.offsprings[1],new_root,rates)
            pushfirst!(rates,root_rate)
        end
    end

    root_rate = tree.attributes[id]
    x = Array{Float64, 1}(undef,0)
    aux(tree,root_rate,x)
    return(x)
end

function extract_relative_rates(tree::Tree; id = 1)
    function aux(subtree, parent_rate, x)
        if subtree.n_nodes == 0
        elseif length(subtree.offsprings) == 0
            pushfirst!(x,log(subtree.attributes[id]) - parent_rate)
        else
            log_rate = deepcopy(log(subtree.attributes[id]))
            aux(subtree.offsprings[2], log_rate, x)
            aux(subtree.offsprings[1], log_rate, x)
            pushfirst!(x,log(subtree.attributes[id]) - parent_rate)
        end
    end

    x = Array{Float64,1}(undef,0)
    aux(tree, 0., x)
    return x
end

function extract_relative_rates(tree::Tree, edge_trees, rates; id = 1)
    function aux(subtree, parent_rate, x)
        if subtree.n_nodes == 0
        elseif length(subtree.offsprings) == 0
            #print("  $(log(subtree.attributes[id]) - parent_rate)  ")
            pushfirst!(x,log(subtree.attributes[id]) - parent_rate)
        else
            log_rate = deepcopy(log(subtree.attributes[id]))
            #print("  $log_rate  ")
            aux(subtree.offsprings[2], log_rate, x)
            aux(subtree.offsprings[1], log_rate, x)
            pushfirst!(x,log(subtree.attributes[id]) - parent_rate)
        end
    end

    relative_rates = []
    for edge_id in 1:length(edge_trees)
        parent_rate = log(get_parent_rate(tree, edge_id, edge_trees))
        #print("  ; $edge_id  $parent_rate  $(length(relative_rates))")
        r = log(rates[edge_id+1]) - parent_rate
        #print(" $r ")
        if isnan(r) || r == -Inf
            println("$edge_id, $r, $(rates[edge_id+1]), $(get_parent_rate(tree, edge_id, edge_trees))")
        end
        if edge_trees[edge_id].tree.n_nodes == 1
            #print("go to 1")
            pushfirst!(relative_rates,r)
        else
            #print("go to 2")
            aux(edge_trees[edge_id].tree,parent_rate, relative_rates)
        end
    end
    return relative_rates
end

function extract_relative_rates(tree::Tree, edge_trees; id = 1)
    rates = extract_rates(tree)
    extract_relative_rates(tree, edge_trees, rates, id = id)
end

function make_ape(tree::Tree ; id = 1)
    n_nodes = tree.n_nodes
    node_id = (n_nodes + 3)/2

    function aux(subtree, root, next_node, tip)
        if subtree.n_nodes == 0
            return Array{Int64,2}(undef,0,2), Array{Float64,1}(undef,0), Array{Float64,1}(undef,0), next_node, tip, Array{String,1}(undef,0)
        elseif length(subtree.offsprings) == 0
            return [root tip], subtree.branch_length, subtree.attributes[id], next_node, tip+1, subtree.label
        else
            tupple1 = aux(subtree.offsprings[1], next_node, next_node+1, tip)
            tupple2 = aux(subtree.offsprings[2], next_node, tupple1[4], tupple1[5])
            if subtree.branch_length == 0.
                return [tupple1[1]; tupple2[1]],
                    [tupple1[2]; tupple2[2]],
                    [tupple1[3]; tupple2[3]],
                    tupple2[4], tupple2[5],
                    [tupple1[6]; tupple2[6]]
            else
                return [[root next_node]; tupple1[1]; tupple2[1]],      # edges
                    [subtree.branch_length; tupple1[2]; tupple2[2]],    # branch lengths
                    [subtree.attributes[id]; tupple1[3]; tupple2[3]],   # rates
                    tupple2[4], tupple2[5],                             # auxiliary variables
                    [tupple1[6]; tupple2[6]]                            # tip labels
            end
        end
    end

    tupple = aux(tree, -1, node_id, 1)
    return tupple[1], tupple[2], tupple[3], tupple[6]
end

function Rsave_tree(tree::Tree ; file = "tree.Rdata")
    edges, branch_lengths, rates, tip_labels = make_ape(tree)
    ntip = (tree.n_nodes + 1)/2

    @rput edges
    @rput branch_lengths
    @rput rates
    @rput ntip
    @rput file
    @rput tip_labels

    reval("""
        tree = list(edge = edges, Nnode = ntip - 1, edge.lengths = branch_lengths, tip.labels = tip_labels)
        class(tree) = "phylo"
        save(tree, rates, file = file)
    """)
end

function get_daughter_edges(tree, edge_id, lefts)

    function aux(subtree, edge, add, x)
        if subtree.n_nodes <= 1
            pushfirst!(x, edge + add)
        elseif edge == 0
            aux(subtree.offsprings[2], 0, add + 1 + lefts[add + 1], x)
            aux(subtree.offsprings[1], 0, add + 1, x)
            pushfirst!(x, add)
        elseif edge < lefts[add + 1]
            aux(subtree.offsprings[1], edge - 1, add + 1, x)
        else
            aux(subtree.offsprings[2], edge - 1 - lefts[add + 1], add + lefts[add + 1] + 1, x)
        end
    end

    x = Array{Int64,1}(undef,0)
    aux(tree, edge_id, 0,x)
    return x
end

function get_daughter_edges(tree, edge_id)
    lefts = n_left(tree)

    get_daughter_edges(tree, edge_id, lefts)
end

function daughters(tree, lefts)
    daughter_edges = Array{Array{Int64,1},1}(undef,tree.n_nodes-1)
    for i in 2:tree.n_nodes
        daughter_edges[i-1] = get_daughter_edges(tree, i-1, lefts)
    end

    return daughter_edges
end

function daughters(tree)
    lefts = n_left(tree)
    daughter_edges = Array{Array{Int64,1},1}(undef,tree.n_nodes)
    for i in 2:tree.n_nodes
        daughter_edges[i-1] = get_daughter_edges(tree, i-1, lefts)
    end

    return daughter_edges
end

function get_daughter_rates(tree, edge_id ; rate_id = 1)

    function aux(subtree, edge)
        if edge == 1
            if subtree.n_nodes < 2
                return [NaN,NaN]
            else
                return [subtree.offsprings[1].attributes[rate_id] ; subtree.offsprings[2].attributes[rate_id]]
            end
        else
            t11 = subtree.offsprings[1]
            t12 = subtree.offsprings[2]
            if t11.n_nodes < (edge-1)
                return aux(t12, edge-1-t11.n_nodes)
            else
                return aux(t11, edge-1)
            end
        end
    end

    aux(tree, edge_id + 1)
end

function get_daughter_rates(tree, edge_id, edge_trees ; rate_id = 1)

    if edge_id == 0
        return [tree.offsprings[1].attributes[rate_id] ; tree.offsprings[2].attributes[rate_id]]
    elseif edge_trees[edge_id].tree.n_nodes > 1
        return [edge_trees[edge_id].tree.offsprings[1].attributes[rate_id] ; edge_trees[edge_id].tree.offsprings[2].attributes[rate_id]]
    end

    function aux(subtree, edge)
        if edge == 1
            if subtree.n_nodes < 2
                return [NaN,NaN]
            else
                return [subtree.offsprings[1].attributes[rate_id] ; subtree.offsprings[2].attributes[rate_id]]
            end
        else
            t11 = subtree.offsprings[1]
            t12 = subtree.offsprings[2]
            if t11.n_nodes < (edge-1)
                return aux(t12, edge-1-t11.n_nodes)
            else
                return aux(t11, edge-1)
            end
        end
    end

    aux(tree, edge_id + 1)
end

function get_daughter_rates(tree, edge_id, edge_trees, rates, lefts ; rate_id = 1, with_edge_tree = true)

    edge = edge_id +1
    if edge_id == 0
        return [rates[edge+1]; rates[edge+1+lefts[edge]]]
    elseif edge_trees[edge_id].tree.n_nodes > 1 && with_edge_tree
        return [edge_trees[edge_id].tree.offsprings[1].attributes[rate_id] ; edge_trees[edge_id].tree.offsprings[2].attributes[rate_id]]
    elseif lefts[edge]==0
        return [NaN,NaN]
    else
        return [rates[edge+1]; rates[edge+1+lefts[edge]]]
    end
end

function get_daughter_rates(tree, edge_id, edge_trees, rates ; rate_id = 1)

    lefts = n_left(tree)
    get_daughter_rates(tree, edge_id, edge_trees, rates, lefts, rate_id = rate_id)
end

function get_node_depth(tree, edge_id)
    function aux(subtree, edge)
        if subtree.n_nodes < 2
            return 0.
        else
            if edge == 1
                return aux(subtree.offsprings[1], edge) + subtree.offsprings[1].branch_length
            elseif (edge - 1) > subtree.offsprings[1].n_nodes
                return aux(subtree.offsprings[2], edge - 1 - subtree.offsprings[1].n_nodes)
            else
                return aux(subtree.offsprings[1], edge - 1)
            end
        end
    end

    aux(tree, edge_id + 1)
end

function get_base_node_depth(tree, edge_id)
    function aux(subtree, edge)
        if subtree.n_nodes < 2
            return subtree.branch_length
        else
            if edge == 1
                return aux(subtree.offsprings[1], edge) + subtree.branch_length
            elseif (edge - 1) > subtree.offsprings[1].n_nodes
                return aux(subtree.offsprings[2], edge - 1 - subtree.offsprings[1].n_nodes)
            else
                return aux(subtree.offsprings[1], edge - 1)
            end
        end
    end

    aux(tree, edge_id + 1)
end

function extract_tip_rates(tree::Tree ; id = 1, return_extinct = true)

    function aux(subtree, x)
        if subtree.n_nodes == 0
        elseif length(subtree.offsprings) == 0
            if subtree.extant || return_extinct
                pushfirst!(x,subtree.attributes[id])
            else
                pushfirst!(x,NaN)
            end
        else
            aux(subtree.offsprings[2], x)
            aux(subtree.offsprings[1], x)
        end
    end

    x = Array{Float64,1}(undef,0)
    aux(tree,x)
    return x
end

function extract_tip_rates_light(tree::Tree ; id = 1, return_extinct = true)

    function aux(subtree, x)
        if subtree.n_nodes == 0
        elseif length(subtree.offsprings) == 0
            if subtree.extant || return_extinct
                pushfirst!(x,subtree.attributes[id])
            end
        else
            aux(subtree.offsprings[2], x)
            aux(subtree.offsprings[1], x)
        end
    end

    x = Array{Float64,1}(undef,0)
    aux(tree,x)
    return x
end

function extract_live_tip_rates(tree::Tree ; id = 1)

    function aux(subtree, x)
        if length(subtree.offsprings) == 0
            if subtree.extant
                pushfirst!(x,subtree.attributes[id])
            end
        else
            aux(subtree.offsprings[2], x)
            aux(subtree.offsprings[1], x)
        end
    end

    x = Array{Float64,1}(undef,0)
    aux(tree,x)
    return x
end

function extract_tip_rates_bu1(tree::Tree, edge_trees::Union{Array{EdgeTree,1},Array{EdgeTreeRates,1}}, tips_id::Array{Bool,1} ; id = 1, return_extinct = true)
    function aux(subtree,x)
        if length(subtree.offsprings)==0
            push!(x,subtree.attributes[id])
        else
            if subtree.offsprings[1].extant
                aux(subtree.offsprings[1],x)
            else
                aux(subtree.offsprings[2],x)
            end
        end
    end
    tip_rates = Array{Float64,1}(undef,0)
    for i in 1:length(edge_trees)
        if tips_id[i]
            aux(edge_trees[i].tree, tip_rates)
        end
    end
    return tip_rates
end



function extract_tip_rates(tree::Tree, edge_trees::Array{EdgeTree,1}, tips_id::Array{Bool,1}, rates::Array{Float64,1} ; id = 1, return_extinct = true)

    tip_rates = Array{Float64,1}(undef,0)
    for i in 1:length(edge_trees)
        if tips_id[i]
            push!(tip_rates,sample(extract_live_tip_rates(edge_trees[i].tree)))
        end
    end
    return tip_rates
end

function extract_tip_rates(tree::Tree, edge_trees::Array{EdgeTreeRates,1},
    tips_id::Array{Bool,1}, rates::Array{Float64,1} ; id = 1, return_extinct = true)
    tip_rates = Array{Float64,1}(undef,0)
    for i in 1:length(edge_trees)
        if tips_id[i]
            r = sample(edge_trees[i].rates)*rates[i+1]/edge_trees[i].rate
            push!(tip_rates,r)
        end
    end
    return tip_rates
end

function extract_tip_rates(tree::Tree, tip_id::Int64 ; id = 1, return_extinct = true)

    function aux(subtree, tip)
        #println(tip)
        if subtree.n_nodes == 0
        elseif length(subtree.offsprings) == 0
            if subtree.extant || return_extinct
                return subtree.attributes[id]
            else
                return NaN
            end
        elseif ((subtree.offsprings[1].n_nodes+1)/1) >= tip
            #println(" $((subtree.offsprings[1].n_nodes+1)/2) $tip $(tip - (subtree.offsprings[1].n_nodes+1)/2) go1")
            aux(subtree.offsprings[1], tip)
        else
            #println(" $((subtree.offsprings[1].n_nodes+1)/2) $tip $(tip - (subtree.offsprings[1].n_nodes+1)/2) go2")
            aux(subtree.offsprings[2], tip - (subtree.offsprings[1].n_nodes+1)/1)
        end
    end

    return aux(tree, tip_id * 2)
end

function cut_tree(tree::Tree, time; prune_extinct = false) #time from the root

    function aux(subtree, sub_time)
        if subtree.n_nodes < 2
            if subtree.branch_length < sub_time
                if prune_extinct && ! subtree.extant
                    return Tree()
                else
                    return subtree
                end
            else
                return Tree(subtree.offsprings, sub_time, subtree.attributes, subtree.n_nodes, true)
            end
        else
            if subtree.branch_length > sub_time
                return Tree(Array{Tree,1}(undef,0), sub_time, subtree.attributes)
            else
                left_tree = aux(subtree.offsprings[1], sub_time - subtree.branch_length)
                right_tree = aux(subtree.offsprings[2], sub_time - subtree.branch_length)
                if ! left_tree.extant && prune_extinct
                    if ! right_tree.extant
                        return Tree()
                    else
                        return Tree(right_tree.offsprings, right_tree.branch_length + subtree.branch_length,
                            subtree.attributes)
                    end
                elseif ! right_tree.extant && prune_extinct
                    return Tree(left_tree.offsprings, left_tree.branch_length + subtree.branch_length,
                        subtree.attributes)
                else
                    return Tree([left_tree, right_tree], subtree.branch_length, subtree.attributes)
                end
            end

        end
    end

    aux(tree, time)
end
