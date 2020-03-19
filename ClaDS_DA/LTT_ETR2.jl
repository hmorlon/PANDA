function LTT(tree::Tree, edge_trees::Array{EdgeTreeRates2,1}, times::Array{Float64,1})
    n = length(times)

    function aux(subtree::Tree, from_root::Float64, id_time::Int64, ltt::Array{Int64,1}, scale::Float64)
        if id_time <= n
            if (subtree.branch_length/scale + from_root) <= times[id_time]
                if subtree.n_nodes < 2
                    if !subtree.extant
                        ltt[id_time] += -1
                    end
                else
                    ltt[id_time] += 1
                    aux(subtree.offsprings[2],from_root+ subtree.branch_length/scale, id_time, ltt, scale)
                    aux(subtree.offsprings[1],from_root+ subtree.branch_length/scale, id_time, ltt, scale)
                end
            else
                aux(subtree,from_root, id_time+1, ltt, scale)
            end
        end

    end

    ltt = fill(0, length(times))
    ltt[1] = 1
    aux(tree, 0., 1, ltt, 1.)
    #println("1 $ltt ")
    i=0
    md = edge_trees[1].stem_depth
    for et in edge_trees# in 1:length(edge_trees)
        if et.tree.n_nodes > 1
            #print("$i $(live_nd[i+1] + maximum(node_depths(edge_trees[i][1]))) ; ")
            aux(et.tree, md - et.stem_depth, 1, ltt,et.scale)
        end
    end
    #println("2 $ltt ")
    ltt = cumsum(ltt, dims = 1)
    #println("3 $ltt ")

    return times, ltt
end
