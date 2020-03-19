#=
New DA structure that allows minimal manipulation
=#

struct EdgeTreeRates2
    tree::Tree                          # the grafted tree
    stem_rate::Array{Float64,1}         # λ_i : every other rates are given divided by λ_i
    tip_rate::Float64                   # the tip rate
    tip_id::Int64                       # tip_id
    tip_rates::Array{Float64,1}         # rates at t_i (for internal branches)
    effective_bl::Float64               # t_j * λ_j / λ_i
    tip_number::Int64                   # extant tips
    scale::Float64
    #relative_rates::Array{Float64,1}
    depth::Float64
    stem_depth::Float64
    parent_edge::Int64
end

function EdgeTreeRates2(tree::Tree, stem_rate::Float64, tip_rate::Float64,
    tip_id::Int64, tip_rates::Array{Float64,1}, tip_number::Int64, depth::Float64,stem_depth::Float64, parent_edge::Int64)
    effective_bl = sum(extract_branch_lengths(tree) .* extract_rates(tree)) / stem_rate
    return EdgeTreeRates2(tree, [stem_rate], tip_rate, tip_id, tip_rates,
        effective_bl, tip_number, depth, stem_depth,parent_edge)
end

function EdgeTreeRates2(tree::Tree, stem_rate::Float64, tip_rate::Float64,
    tip_id::Int64, tip_rates::Array{Float64,1},
    effective_bl::Float64, tip_number::Int64,depth::Float64,stem_depth::Float64,
    parent_edge::Int64)
    return EdgeTreeRates2(tree, [stem_rate], tip_rate, tip_id, tip_rates, effective_bl, tip_number, depth,stem_depth, parent_edge)
end

function EdgeTreeRates2(tree::Tree, stem_rate::Array{Float64,1}, tip_rate::Float64,
    tip_id::Int64, tip_rates::Array{Float64,1},
    effective_bl::Float64, tip_number::Int64,
    depth::Float64,stem_depth::Float64,parent_edge::Int64)

    return EdgeTreeRates2(tree, stem_rate, tip_rate, tip_id, tip_rates, effective_bl,
        tip_number, stem_rate[1], depth,stem_depth, parent_edge)
end

#=,
function EdgeTreeRates2(tree::Tree, stem_rate::Array{Float64,1}, tip_rate::Float64,
    tip_id::Int64, tip_rates::Array{Float64,1},
    effective_bl::Array{Float64,1}, tip_number::Int64,depth::Float64,
    parent_edge::Int64)

    relative_rates = extract_relative_rates(tree)[2:end]
    return EdgeTreeRates2(tree, [stem_rate], tip_rate, tip_id, tip_rates, effective_bl,tip_number, relative_rates,depth, parent_edge)
end
=#

function init_edge_tree_rates2(tree, rates)
    bl = extract_branch_lengths(tree)
    edge_trees = Array{EdgeTreeRates2,1}(undef,0)
    lefts = n_left(tree)
    for i in 2:tree.n_nodes
        nd = get_node_depth(tree, i-1)
        parent_edge = get_parent_edge(tree, i-1)

        ne = 0
        if lefts[i]==0
            ne = 1
        end

        push!(edge_trees,
            EdgeTreeRates2(Tree(Array{Tree,1}(undef,0), bl[i] * rates[i], 1.),       # tree
                rates[i],                                                           # stem_rate
                1.,                                                                 # tip_rate
                1,                                                                  # tip_id
                [1.],                                                               # tip_rates
                bl[i],                                                              # effective_bl
                ne,                                                                  # tip_number
                nd,
                nd + bl[i],
                parent_edge))
    end
    return edge_trees
end
