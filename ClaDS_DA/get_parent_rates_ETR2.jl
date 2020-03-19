
function get_parent_rate(tree::Tree, edge_id::Int64, edge_trees::Array{EdgeTreeRates2,1})
    parent_edge = edge_trees[edge_id].parent_edge

    if parent_edge == 0
        return tree.attributes[1]
    else
        return edge_trees[parent_edge].tip_rate * edge_trees[parent_edge].stem_rate[1]
    end
end

function get_parent_rate(tree::Tree, edge_id::Int64, edge_trees::Array{EdgeTreeRates2,1}, rates::Array{Float64,1})
    parent_edge = edge_trees[edge_id].parent_edge

    if parent_edge == 0
        return rates[1]
    else
        return edge_trees[parent_edge].tip_rate * rates[parent_edge+1]
    end
end
