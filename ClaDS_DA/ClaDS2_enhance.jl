
function graft_to_tip(tree,subtree, tip_id)

    function aux(t1, t2, tip)
        if t1.n_nodes == 1
            return Tree(t2.offsprings, t2.branch_length + t1.branch_length, t1.attributes, t2.n_nodes, t2.extant)
        else
            t11 = t1.offsprings[1]
            t12 = t1.offsprings[2]
            t11_tips = (t11.n_nodes + 1) / 2
            if t11_tips < tip
                t12 = aux(t12, t2, tip - t11_tips)
                return Tree([t11, t12], t1.branch_length, t1.attributes)
            else
                t11 = aux(t11, t2, tip)
                return Tree([t11, t12], t1.branch_length, t1.attributes)
            end
        end
    end

    aux(tree,subtree, tip_id)
end

function graft_to_edge(tree, sub_tree, edge, tip)

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

function graft_to_tips(tree, subtree_list, tips)
    if length(tips) == 0
        return tree
    end

    ntips = Int64((tree.n_nodes + 1) / 2)
    subtrees = fill(Tree(),ntips)
    #println("$tips, $ntips")

    for i in 1:ntips
        for j in 1:length(tips)
            if i==tips[j]
                subtrees[i] = subtree_list[j]
                break
            end
        end
    end


    function aux(t1, t2)
        if t1.n_nodes <= 1
            #=if ! t2[1].extant
                println("is it ok ? $(t2[1]) ; $(t1.extant)")
            end=#
            return Tree(t2[1].offsprings, t2[1].branch_length + t1.branch_length, t1.attributes, t2[1].extant)
        else
            t11 = t1.offsprings[1]
            t12 = t1.offsprings[2]
            t11_tips = Int64((t11.n_nodes + 1) / 2)

            t11 = aux(t11, t2[1:t11_tips])
            t12 = aux(t12, t2[(t11_tips+1):end])

            return Tree([t11, t12], t1.branch_length, t1.attributes)
        end
    end

    new_tree = aux(tree,subtrees)
    return new_tree
end

function init_edge_tree(tree, rates)
    bl = extract_branch_lengths(tree)
    edge_trees = Array{EdgeTree,1}(undef,0)

    for i in 2:tree.n_nodes
        push!(edge_trees, EdgeTree(Tree(Array{Tree,1}(undef,0), bl[i], rates[i]), 1, 1))
    end
    return edge_trees
end

function change_edge_trees_rejection!(edge_trees, tree, edge_id, σ, α, ε, f, rates, branch_lengths ; do_tips = true, max_try = 1_000, max_node_number = 50, keep_if_any = true)
    trial = 0

    daughter_rates = get_daughter_rates(tree, edge_id)

    if isnan(daughter_rates[1])
        if do_tips
            n_former = n_extant_tips(edge_trees[edge_id].tree)
            int_tree = Tree()
            while trial < max_try
                trial += 1

                int_tree, unsampled, n_alive = sim_ClaDS2_time_unsampled(branch_lengths[edge_id], σ, α, ε, rates[edge_id],return_if_sampled = true,
                max_node_number = max_node_number, not_sampled = false, sampling_proba = f, return_if_extinct = false, make_tree = false)
                if unsampled
                    edge_trees[edge_id] = EdgeTree(int_tree, 1, 1)
                    return edge_trees[edge_id]
                end
            end
        end
        return edge_trees[edge_id]
    else
        log_α = log(α)
        node_depth = get_node_depth(tree, edge_id)
        if keep_if_any
            if edge_trees[edge_id].tree.n_nodes < 2
                former_parent_rate = rates[edge_id]
            else
                edge_tree_tip_rates = extract_tip_rates(edge_trees[edge_id].tree, return_extinct = false)
                former_parent_rate = edge_tree_tip_rates[edge_trees[edge_id].tip_id]
            end

            former_lik = 0. -
                ((log(former_parent_rate) + log_α - log(daughter_rates[1]))^2 + (log(former_parent_rate) + log_α - log(daughter_rates[2]))^2)/(2 * σ^2) + log(former_parent_rate)
        else
            former_lik = (log(daughter_rates[1]) + log(daughter_rates[1]))/2 -log_α+(σ^2)/2
        end

        former_lik += log(edge_trees[edge_id].n)

        while trial < max_try
            trial += 1

            int_tree =
                sim_ClaDS2_time(branch_lengths[edge_id], σ, α, ε, rates[edge_id],prune_extinct = false,
                max_node_number = max_node_number, return_if_extinct = false, make_tree = true, return_if_max = false)
            nalive = n_extant_tips(int_tree)
            sampled = []
            tip_trees = []
            graft_tips = []
            keep = false
            keep_tip = false
            add_to_tip_id = 0
            new_tip_id = NaN
            new_lik = 0.
            new_parent_rate = 0.
             any_kept = false
             w = Array{Float64,1}(undef,0)
            if (0 < nalive) & (int_tree.branch_length >= 0)
                int_tip_rates = extract_tip_rates(int_tree, return_extinct = false)
                for new_tip in sample(1:nalive)
                    k = 0
                    sampled = []
                    tip_trees = []
                    graft_tips = []
                    keep = false
                    keep_tip = false
                    add_to_tip_id = 0
                    new_tip_id = NaN
                    new_lik = 0.
                    new_parent_rate = 0.

                    for i in 1:length(int_tip_rates)
                        if ! isnan(int_tip_rates[i])
                            k += 1
                            if k == new_tip
                                new_parent_rate = int_tip_rates[i]
                                new_lik = 0. - ((log(new_parent_rate) +
                                    log_α - log(daughter_rates[1]))^2 + (log(new_parent_rate) +
                                    log_α - log(daughter_rates[2]))^2)/(2 * σ^2)  +log(new_parent_rate) + log(nalive)
                                u = rand()
                                keep = log(u) < new_lik - former_lik
                                keep_tip = deepcopy(keep)
                                new_tip_id = i
                            end
                        end
                    end

                    if keep | keep_if_any
                        keep = true
                        k = 0
                        i = 0
                        while (keep) & (i < length(int_tip_rates))
                            i += 1
                            if ! isnan(int_tip_rates[i])
                                k += 1
                                if k != new_tip
                                    tip_tree, unsampled, n_alive = sim_ClaDS2_time_unsampled(node_depth, σ, α, ε, int_tip_rates[i],
                                        max_node_number = max_node_number, not_sampled = true, sampling_proba = f, return_if_sampled = true, make_tree = true)
                                    if unsampled
                                        push!(tip_trees, tip_tree)
                                        push!(graft_tips, i)
                                        if i < new_tip_id
                                            add_to_tip_id += n_tip(tip_tree) - 1
                                        end
                                    else
                                        keep = false
                                    end
                                end
                            end
                        end
                    end

                    if keep & keep_tip
                        edge_trees[edge_id] = EdgeTree(graft_to_tips(int_tree, [Tree(Array{Tree,1}(undef,0),0., [0.], 0, true); tip_trees],
                            [new_tip_id; graft_tips]), new_tip_id + add_to_tip_id, nalive)

                        return  edge_trees[edge_id]
                    elseif keep & keep_if_any

                        return  edge_trees[edge_id]
                    end
                end
            end
        end
        println(trial)
    end
end

function change_edge_trees_MH!(edge_trees, tree, edge_id, σ, α, ε, f, rates, branch_lengths ; do_tips = true, max_try = 1_000, max_node_number = 50, keep_if_any = true)
    trial = 0

    daughter_rates = get_daughter_rates(tree, edge_id)

    if isnan(daughter_rates[1])
        if do_tips
            int_tree = Tree()
            n_alive = edge_trees[edge_id].n
            if n_alive == 0
                former_lik = -Inf
            elseif n_alive == 1
                former_lik = log(f)
            else
                former_lik = log(n_alive) + log(f) + (n_alive - 1) * log(1-f)
            end
            while trial < 1#max_try
                trial += 1

                int_tree, unsampled, n_alive = sim_ClaDS2_time_unsampled(branch_lengths[edge_id], σ, α, ε, rates[edge_id],return_if_sampled = true,
                max_node_number = max_node_number, not_sampled = false, sampling_proba = f, return_if_extinct = false, make_tree = false)
                if n_alive == 0
                    log_sampling_proba = -Inf
                elseif n_alive == 1
                    log_sampling_proba = log(f)
                else
                    log_sampling_proba = log(n_alive) + log(f) + (n_alive - 1) * log(1-f)
                end
                u = log(rand())
                unsampled = u < (log_sampling_proba-former_lik)
                if unsampled
                    edge_trees[edge_id] = EdgeTree(int_tree, n_alive, n_alive)
                    return edge_trees[edge_id]
                end
            end
        end
        return edge_trees[edge_id]
    else
        log_α = log(α)
        node_depth = get_node_depth(tree, edge_id)
        if keep_if_any
            if edge_trees[edge_id].tree.n_nodes < 2
                former_parent_rate = rates[edge_id]
            else
                edge_tree_tip_rates = extract_tip_rates(edge_trees[edge_id].tree, return_extinct = false)
                former_parent_rate = edge_tree_tip_rates[edge_trees[edge_id].tip_id]
            end

            former_lik = 0. -
                ((log(former_parent_rate) + log_α - log(daughter_rates[1]))^2 + (log(former_parent_rate) + log_α - log(daughter_rates[2]))^2)/(2 * σ^2) + log(former_parent_rate)
        else
            former_lik = (log(daughter_rates[1]) + log(daughter_rates[1]))/2 -log_α+(σ^2)/2
        end

        former_lik += log(edge_trees[edge_id].n)

        while trial < max_try
            trial += 1

            int_tree =
                sim_ClaDS2_time(branch_lengths[edge_id], σ, α, ε, rates[edge_id],prune_extinct = false,
                max_node_number = max_node_number, return_if_extinct = false, make_tree = true, return_if_max = false)
            nalive = n_extant_tips(int_tree)
            sampled = []
            tip_trees = []
            graft_tips = []
            keep = false
            keep_tip = false
            add_to_tip_id = 0
            new_tip_id = NaN
            new_lik = 0.
            new_parent_rate = 0.
             any_kept = false
             w = Array{Float64,1}(undef,0)
            if (0 < nalive) & (int_tree.branch_length >= 0)
                new_lik_rate = 0.

                int_tip_rates = extract_tip_rates(int_tree, return_extinct = false)
                for new_tip in sample(1:nalive)
                    k = 0
                    sampled = []
                    tip_trees = []
                    graft_tips = []
                    keep = false
                    keep_tip = false
                    add_to_tip_id = 0
                    new_tip_id = NaN
                    new_lik = 0.
                    new_parent_rate = 0.

                    for i in 1:length(int_tip_rates)
                        if ! isnan(int_tip_rates[i])
                            k += 1
                            if k == new_tip
                                new_parent_rate = int_tip_rates[i]
                                new_lik_rate = 0. - ((log(new_parent_rate) +
                                    log_α - log(daughter_rates[1]))^2 + (log(new_parent_rate) +
                                    log_α - log(daughter_rates[2]))^2)/(2 * σ^2)  +log(new_parent_rate) + log(nalive)
                                #u = rand()
                                #keep = log(u) < new_lik - former_lik
                                #keep_tip = deepcopy(keep)
                                new_tip_id = i
                            end
                        end
                    end

                    if true#keep | keep_if_any
                        former_tip_number = n_extant_tips(edge_trees[edge_id].tree) - 1
                        if former_tip_number == 0
                            former_lik += 0
                        else
                            former_lik += former_tip_number * log(1-f)
                        end
                        u = log(rand())

                        keep = true
                        new_lik = new_lik_rate
                        keep = u < (new_lik-former_lik)

                        k = 0
                        i = 0
                        ntips = 0
                        while (keep) & (i < length(int_tip_rates))
                            i += 1
                            if ! isnan(int_tip_rates[i])
                                k += 1
                                if k != new_tip
                                    tip_tree, unsampled, n_alive = sim_ClaDS2_time_unsampled(node_depth, σ, α, ε, int_tip_rates[i],
                                        max_node_number = max_node_number, not_sampled = true, sampling_proba = f, return_if_sampled = true, make_tree = true)
                                    if true#unsampled
                                        ntips += n_alive
                                        if ntips == 0
                                            new_lik = 0
                                        elseif ntips == Inf
                                            new_lik = -Inf
                                        else
                                            new_lik = new_lik_rate + ntips * log(1-f)
                                        end
                                        keep = u < (new_lik-former_lik)
                                        push!(tip_trees, tip_tree)
                                        push!(graft_tips, i)
                                        if i < new_tip_id
                                            add_to_tip_id += n_tip(tip_tree) - 1
                                        end
                                    else
                                        keep = false
                                    end
                                end
                            end
                        end

                    end

                    if keep #& keep_tip
                        edge_trees[edge_id] = EdgeTree(graft_to_tips(int_tree, [Tree(Array{Tree,1}(undef,0),0., [0.], 0, true); tip_trees],
                            [new_tip_id; graft_tips]), new_tip_id + add_to_tip_id, nalive)

                        return  edge_trees[edge_id]
                    elseif true #& keep_if_any

                        return  edge_trees[edge_id]
                    end
                end
            end
        end
        println(trial)
    end
end

function update_edge_trees!(edge_trees, tree, σ, α, ε, fs, rates_all, branch_lengths_all ; do_tips = true, max_try = 1_000, max_node_number = 50, keep_if_any = true, enhance_method = "reject")
    n_edges = length(edge_trees)
    rates = rates_all[2:end]
    branch_lengths = branch_lengths_all[2:end]

    for edge_id in n_edges:-1:1
        f = fs[edge_id]
        if enhance_method == "reject"
            change_edge_trees_rejection!(edge_trees, tree, edge_id, σ, α, ε, f, rates, branch_lengths,
                max_try = max_try, max_node_number = max_node_number, keep_if_any = keep_if_any, do_tips = do_tips);
        elseif enhance_method == "MH"
            change_edge_trees_MH!(edge_trees, tree, edge_id, σ, α, ε, f, rates, branch_lengths,
                max_try = max_try, max_node_number = max_node_number, keep_if_any = keep_if_any, do_tips = do_tips);
        elseif enhance_method == "rates"
            change_edge_trees_rates!(edge_trees, tree, edge_id, σ, α, ε, f, rates, branch_lengths,
                max_try = max_try, max_node_number = max_node_number, keep_if_any = keep_if_any, do_tips = do_tips);
        elseif enhance_method == "MHrates"
            change_edge_trees_rateMH!(edge_trees, tree, edge_id, σ, α, ε, f, rates, branch_lengths,
                max_try = max_try, max_node_number = max_node_number, keep_if_any = keep_if_any, do_tips = do_tips);
        end

    end
end

function update_edge_trees!(edge_trees, tree, σ, α, ε, fs; max_try = 1_000, max_node_number = 50, keep_if_any = true, do_tips = true, enhance_method = "reject")
    rates_all = extract_rates(tree)
    branch_lengths_all = extract_branch_lengths(tree)

    update_edge_trees!(edge_trees, tree, σ, α, ε, fs, rates_all, branch_lengths_all, max_try = max_try,
        max_node_number = max_node_number, keep_if_any = keep_if_any, do_tips = do_tips, enhance_method = enhance_method)
end

function graft_edge_trees(tree, edge_trees)

    function aux(subtree, sub_edge_trees)
        if subtree.n_nodes < 2
            return sub_edge_trees[1].tree
        else
            n_edges_left = subtree.offsprings[1].n_nodes + 1
            sub_left = deepcopy(sub_edge_trees[2:n_edges_left])
            sub_right = deepcopy(sub_edge_trees[(n_edges_left+1):subtree.n_nodes])
            subtree_left = aux(subtree.offsprings[1], sub_left)
            subtree_right = aux(subtree.offsprings[2], sub_right)
            return graft_to_edge(Tree([subtree_left, subtree_right], 0., subtree.attributes), sub_edge_trees[1].tree, 0, sub_edge_trees[1].tip_id)
        end
    end

    n_edges_left = tree.offsprings[1].n_nodes
    left = deepcopy(edge_trees[1:n_edges_left])
    right = deepcopy(edge_trees[(n_edges_left+1):(tree.n_nodes - 1)])
    tree_left = aux(tree.offsprings[1], left)
    tree_right = aux(tree.offsprings[2], right)
    return Tree([tree_left, tree_right], tree.branch_length, tree.attributes)
end
