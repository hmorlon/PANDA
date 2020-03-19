function update_edges_ETR2!(tree::Tree, edge_trees::Array{EdgeTreeRates2,1}, σ::Float64, α::Float64, ε::Float64, rates::Array{Float64,1},
    lefts::Array{Int64,1}; it_rates = 1)

    ε = draw_ε_crown(tree, edge_trees, lefts)

    for j in 1:it_rates

        draw_λi_quad!(rates, edge_trees, σ, α, ε, tree, lefts)
        relative_rates = extract_relative_rates(tree, edge_trees, rates)

        σ = draw_σ(relative_rates, α, β0 = 0.05, α0 = 0.5)
        α = draw_α(relative_rates, σ,  α_0 = -0.05, σ_0 = 0.1)

        ε = draw_ε_crown(tree, edge_trees, lefts)
    end

    return ε, σ, α
end

function change_bl(tree::Tree, scale::Float64)
    function aux(sub_tree)
        if sub_tree.n_nodes < 2
            return Tree(sub_tree.offsprings,
                sub_tree.branch_length / scale,
                sub_tree.attributes.*scale,
                sub_tree.n_nodes,
                sub_tree.extant)
        else
            left = aux(sub_tree.offsprings[1])
            right = aux(sub_tree.offsprings[2])
            return Tree([left,right],
                sub_tree.branch_length / scale,
                sub_tree.attributes.*scale,
                sub_tree.n_nodes,
                sub_tree.extant)
        end
    end

    aux(tree)
end

function change_bl(tree::Tree, scale::Float64, stem_rate::Float64)
    function aux(sub_tree)
        if sub_tree.n_nodes < 2
            return Tree(sub_tree.offsprings,
                sub_tree.branch_length / scale,
                sub_tree.attributes.*stem_rate,
                sub_tree.n_nodes,
                sub_tree.extant)
        else
            left = aux(sub_tree.offsprings[1])
            right = aux(sub_tree.offsprings[2])
            return Tree([left,right],
                sub_tree.branch_length / scale,
                sub_tree.attributes.*stem_rate,
                sub_tree.n_nodes,
                sub_tree.extant)
        end
    end

    aux(tree)
end

function graft_to_edge_ETR2(tree::Tree, sub_tree::Tree, edge::Int64, tip::Int64, scale::Float64)

    if sub_tree.n_nodes < 2
        return Tree(tree.offsprings, tree.branch_length + sub_tree.branch_length, tree.attributes, tree.n_nodes, sub_tree.extant)
    end

    function aux(t1, t2, edge_id)
        if t1.n_nodes == 1
            return t2
        elseif edge_id == 1
            new_t1 = Tree(t1.offsprings, 0., t1.attributes, t1.n_nodes, t1.extant)
            return graft_to_tip(t2, new_t1, tip)
        else
            t11 = t1.offsprings[1]
            t12 = t1.offsprings[2]
            if t11.n_nodes < (edge_id-1)
                t12 = aux(t12, t2, edge_id-1-t11.n_nodes)
                return Tree([t11, t12], t1.branch_length, t1.attributes)
            else
                t11 = aux(t11, t2, edge_id-1)
                return Tree([t11, t12], t1.branch_length, t1.attributes)
            end
        end
    end

    return aux(tree, sub_tree, edge+1)
end

function graft_edge_trees(tree, edge_trees::Array{EdgeTreeRates2,1})

    function aux(subtree, sub_edge_trees)
        if subtree.n_nodes < 2
            return change_bl(sub_edge_trees[1].tree, sub_edge_trees[1].scale, sub_edge_trees[1].stem_rate[1])
        else
            n_edges_left = subtree.offsprings[1].n_nodes + 1
            sub_left = deepcopy(sub_edge_trees[2:n_edges_left])
            sub_right = deepcopy(sub_edge_trees[(n_edges_left+1):subtree.n_nodes])
            subtree_left = aux(subtree.offsprings[1], sub_left)
            subtree_right = aux(subtree.offsprings[2], sub_right)
            return graft_to_edge(Tree([subtree_left, subtree_right], 0., subtree.attributes),
                change_bl(sub_edge_trees[1].tree, sub_edge_trees[1].scale, sub_edge_trees[1].stem_rate[1]), 0, sub_edge_trees[1].tip_id)
        end
    end

    n_edges_left = tree.offsprings[1].n_nodes
    left = deepcopy(edge_trees[1:n_edges_left])
    right = deepcopy(edge_trees[(n_edges_left+1):(tree.n_nodes - 1)])
    tree_left = aux(tree.offsprings[1], left)
    tree_right = aux(tree.offsprings[2], right)
    return Tree([tree_left, tree_right], tree.branch_length, tree.attributes)
end
