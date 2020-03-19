function init_edge_tree_rates(tree, rates)
    bl = extract_branch_lengths(tree)
    edge_trees = Array{EdgeTreeRates,1}(undef,0)

    for i in 2:tree.n_nodes
        push!(edge_trees, EdgeTreeRates(Tree(Array{Tree,1}(undef,0), bl[i], rates[i]), 1, rates[i], [rates[i]]))
    end
    return edge_trees
end

function change_edge_trees_rates_bu!(edge_trees, tree, edge_id, σ, α, ε, f, rates, branch_lengths ; do_tips = true, max_try = 1_000, max_node_number = 50, keep_if_any = true)
    trial = 0

    daughter_rates = get_daughter_rates(tree, edge_id)

    if isnan(daughter_rates[1])
        if do_tips
            int_tree = Tree()
            while trial < max_try
                trial += 1

                int_tree, unsampled, n_alive = sim_ClaDS2_time_unsampled(branch_lengths[edge_id], σ, α, ε, rates[edge_id],return_if_sampled = true,
                max_node_number = max_node_number, not_sampled = false, sampling_proba = f, return_if_extinct = false, make_tree = false)
                if unsampled
                    edge_trees[edge_id] = EdgeTreeRates(int_tree, 1, rates[edge_id], [rates[edge_id]])
                    return edge_trees[edge_id]
                end
            end
        end
        return edge_trees[edge_id]
    else
        log_α = log(α)
        node_depth = get_node_depth(tree, edge_id)
        #=if keep_if_any
            if edge_trees[edge_id].tree.n_nodes < 2
                former_parent_rate = rates[edge_id]
            else
                edge_tree_tip_rates = extract_tip_rates(edge_trees[edge_id].tree, return_extinct = false)
                former_parent_rate = edge_tree_tip_rates[edge_trees[edge_id].tip_id]
            end

            former_lik = 0. #-
                #((log(former_parent_rate) + log_α - log(daughter_rates[1]))^2 + (log(former_parent_rate) + log_α - log(daughter_rates[2]))^2)/(2 * σ^2) + log(former_parent_rate)
        else
            former_lik = (log(daughter_rates[1]) + log(daughter_rates[1]))/2 -log_α+(σ^2)/2
        end=#

        sum_rates = 0.
        #ratio = former_parent_rate / edge_trees[edge_id].rate
        function nu(former_rate)#, former_ratio)
            if isnan(former_rate)
                return 0.
            end
            fr = former_rate#*former_ratio
            log_fr = log(fr)
            return fr * exp((-(log_fr + log_α - log(daughter_rates[1]))^2 - (log_fr + log_α - log(daughter_rates[2]))^2)/(2 * σ^2))
        end



        while trial < max_try
            trial += 1

            int_tree, int_tip_rates =
                sim_ClaDS2_time(branch_lengths[edge_id], σ, α, ε, rates[edge_id],prune_extinct = false,
                max_node_number = max_node_number, return_if_extinct = false, make_tree = true, return_if_max = false, return_na = true)
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
            if (0 < nalive) & (int_tree.branch_length >= 0)
                #int_tip_rates = extract_tip_rates(int_tree, return_extinct = false)
                #println(int_tip_rates)
                if sum_rates==0
                    tr =
                    for r in tr#edge_trees[edge_id].rates
                        sum_rates += nu(r)#, ratio)
                    end
                end
                former_lik += + log(sum_rates)
                w = Array{Float64,1}(undef,length(int_tip_rates))
                Ws = 0.
                for r in 1:length(int_tip_rates)
                    w[r] = nu(int_tip_rates[r],1.)
                    Ws += w[r]
                end
                for new_tip in sample(1:(length(int_tip_rates)), Weights(w))
                    #println(w)
                    #println("$Ws $sum_rates")
                    #println(new_tip)
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
                            if i == new_tip
                                new_parent_rate = int_tip_rates[i]
                                new_lik = 0. + log(Ws)
                                #- ((log(new_parent_rate) +
                                    #log_α - log(daughter_rates[1]))^2 + (log(new_parent_rate) +
                                    #log_α - log(daughter_rates[2]))^2)/(2 * σ^2)  +log(new_parent_rate) + log(nalive)
                                u = rand()
                                keep = log(u) < new_lik - former_lik
                                #println("$σ $((Ws)) $sum_rates $keep $ratio")
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
                                if i != new_tip
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
                        edge_trees[edge_id] = EdgeTreeRates(graft_to_tips(int_tree, [Tree(Array{Tree,1}(undef,0),0., [0.], 0, true); tip_trees],
                            [new_tip_id; graft_tips]), new_tip_id + add_to_tip_id, int_tip_rates[new_tip], int_tip_rates)

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

function change_edge_trees_rateMH_bu!(edge_trees, tree, edge_id, σ, α, ε, f, rates, branch_lengths ; do_tips = true, max_try = 1_000, max_node_number = 50, keep_if_any = true)
    trial = 0

    daughter_rates = get_daughter_rates(tree, edge_id)

    if isnan(daughter_rates[1])
        if do_tips
            int_tree = Tree()
            n_alive = edge_trees[edge_id].tip_id#n_extant_tips(edge_trees[edge_id].tree)
            if n_alive == 0
                former_lik = -Inf
            elseif n_alive == 1
                former_lik = log(f)
            else
                former_lik = log(n_alive) + log(f) + (n_alive - 1) * log(1-f)
            end
            while trial < 1#max_try
                trial += 1

                int_tree, unsampled, n_alive, tip_rates = sim_ClaDS2_time_unsampled_rates(branch_lengths[edge_id], σ, α, ε,
                    rates[edge_id],return_if_sampled = true,
                    max_node_number = max_node_number, not_sampled = false,
                    sampling_proba = f, return_if_extinct = false, make_tree = false)
                if n_alive == 0 #| (int_tree.branch_length <=0)
                    #log_sampling_proba = -Inf
                    return edge_trees[edge_id]
                elseif n_alive == 1
                    log_sampling_proba = log(f)
                else
                    log_sampling_proba = log(n_alive) + log(f) + (n_alive - 1) * log(1-f)
                end
                u = log(rand())
                unsampled = u < (log_sampling_proba-former_lik)
                if unsampled
                    edge_trees[edge_id] = EdgeTreeRates(int_tree, n_alive, rates[edge_id], tip_rates)
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

            former_lik = 0. #-
                #((log(former_parent_rate) + log_α - log(daughter_rates[1]))^2 + (log(former_parent_rate) + log_α - log(daughter_rates[2]))^2)/(2 * σ^2) + log(former_parent_rate)
        else
            former_lik = (log(daughter_rates[1]) + log(daughter_rates[1]))/2 -log_α+(σ^2)/2
        end

        sum_rates = 0.
        ratio = former_parent_rate / edge_trees[edge_id].rate
        function nu(former_rate, former_ratio)
            if isnan(former_rate)
                return 0.
            end
            fr = former_rate*former_ratio
            log_fr = log(fr)
            return fr * exp((-(log_fr + log_α - log(daughter_rates[1]))^2 - (log_fr + log_α - log(daughter_rates[2]))^2)/(2 * σ^2))
        end

        for r in edge_trees[edge_id].rates
            sum_rates += nu(r, ratio)
        end
        former_lik += + log(sum_rates)

        while trial < max_try
            trial += 1
            int_tree, int_tip_rates, which_alive =
                sim_ClaDS2_time_rates(branch_lengths[edge_id], σ, α, ε, rates[edge_id],prune_extinct = false,
                max_node_number = max_node_number, make_tree = true, return_if_max = true, return_na = true)

            nalive = length(int_tip_rates)
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
            if int_tree.branch_length<=0
                return  edge_trees[edge_id]

           elseif (0 < nalive) #& (int_tree.branch_length >= 0)
               #int_tip_rates = extract_tip_rates(int_tree, return_extinct = false)
               #println(int_tip_rates)
               w = Array{Float64,1}(undef,0)
               Ws = 0.
               for r in 1:length(int_tip_rates)
                   #if !isnan(int_tip_rates[r])
                   push!(w,nu(int_tip_rates[r],1.))
                   Ws += w[end]
                   #end
               end
               for new_tip in sample(1:nalive, Weights(w))
                   current_n = nalive
                   k = 0
                   sampled = []
                   tip_trees = []
                   graft_tips = []
                   keep = false
                   keep_tip = false
                   add_to_tip_id = 0
                   new_tip_id = NaN
                   new_tip_id_rate = NaN
                   new_lik = 0.
                   new_parent_rate = 0.
                   new_lik_rate = 0.
                  # println(new_tip)

                   for i in 1:length(which_alive)
                       if which_alive[i]
                           k += 1
                           if k == new_tip
                               new_parent_rate = int_tip_rates[k]
                               new_lik_rate = 0. + log(Ws)
                               u = rand()
                               keep = log(u) < new_lik - former_lik
                               keep_tip = deepcopy(keep)
                               new_tip_id = i
                               new_tip_id_rate = k
                           end
                       end
                   end

                    if true
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
                        n_alive = 0
                        current_lik = new_lik_rate
                        while (keep) & (i < length(which_alive))
                            i += 1
                            if which_alive[i]
                                k += 1
                                if k != new_tip
                                    tip_tree = sim_ClaDS2_time(node_depth, σ, α, ε, int_tip_rates[k],#u,current_lik,log(1-f),
                                        max_node_number = max_node_number)
                                    if tip_tree.branch_length < 0
                                        n_alive = Inf
                                    else
                                        n_alive = n_extant_tips(tip_tree)
                                    end
                                    if true#unsampled
                                        ntips += n_alive
                                        current_n += -1#n_alive - 1
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
                        edge_trees[edge_id] = EdgeTreeRates(graft_to_tips(int_tree, [Tree(Array{Tree,1}(undef,0),0., [0.], 0, true); tip_trees],
                            [new_tip_id; graft_tips]), new_tip_id + add_to_tip_id, int_tip_rates[new_tip_id_rate], int_tip_rates)

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


function change_edge_trees_rateMH!(edge_trees, tree, edge_id, σ, α, ε, f, rates, branch_lengths ;
    do_tips = true, max_try = 1_000, max_node_number = 50, keep_if_any = true, rep = 1)
    trial = 0

    daughter_rates = get_daughter_rates(tree, edge_id)

    if isnan(daughter_rates[1])
        if do_tips
            int_tree = Tree()
            n_alive = edge_trees[edge_id].tip_id#n_extant_tips(edge_trees[edge_id].tree)
            if n_alive == 0
                former_lik = -Inf
            elseif n_alive == 1
                former_lik = log(f)
            else
                former_lik = log(n_alive) + log(f) + (n_alive - 1) * log(1-f)
            end
            while trial < rep
                trial += 1
                u = log(rand())

                int_tree, unsampled, n_alive, tip_rates = sim_ClaDS2_time_unsampled_rates(branch_lengths[edge_id], σ, α, ε,
                    rates[edge_id],u,-1*former_lik,return_if_sampled = true,
                    max_node_number = max_node_number, not_sampled = false,
                    sampling_proba = f, return_if_extinct = false, make_tree = false)
                if unsampled
                    edge_trees[edge_id] = EdgeTreeRates(int_tree, n_alive, rates[edge_id], tip_rates)
                    former_lik = log(n_alive) + log(f) + (n_alive - 1) * log(1-f)
                    #return edge_trees[edge_id]
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

            former_lik = 0. #-
                #((log(former_parent_rate) + log_α - log(daughter_rates[1]))^2 + (log(former_parent_rate) + log_α - log(daughter_rates[2]))^2)/(2 * σ^2) + log(former_parent_rate)
        else
            former_lik = (log(daughter_rates[1]) + log(daughter_rates[1]))/2 -log_α+(σ^2)/2
        end

        sum_rates = 0.
        ratio = former_parent_rate / edge_trees[edge_id].rate
        function nu(former_rate, former_ratio)
            if isnan(former_rate)
                return 0.
            end
            fr = former_rate*former_ratio
            log_fr = log(fr)
            return fr * exp((-(log_fr + log_α - log(daughter_rates[1]))^2 - (log_fr + log_α - log(daughter_rates[2]))^2)/(2 * σ^2))
        end

        for r in edge_trees[edge_id].rates
            sum_rates += nu(r, ratio)
        end
        former_lik += + log(sum_rates)
        it = 0

        while it < rep
            it += 1
            #println("$edge_id ; $it")
            trial = 0
            while trial < max_try
                nalive = 0
                int_tree = Tree()
                int_tip_rates = Array{Float64,1}(undef,0)
                which_alive = Array{Bool,1}(undef,0)
                while (nalive == 0) && (trial < max_try)
                    trial += 1
                    int_tree, int_tip_rates, which_alive =
                        sim_ClaDS2_time_rates(branch_lengths[edge_id], σ, α, ε, rates[edge_id],prune_extinct = false,
                        max_node_number = max_node_number, make_tree = true, return_if_max = true, return_na = true)

                    if int_tree.branch_length<=0
                        nalive = Inf
                    else
                        nalive = length(int_tip_rates)
                    end
                end
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
                if int_tree.branch_length<=0
                    break
                    #return  edge_trees[edge_id]

               elseif (0 < nalive) #& (int_tree.branch_length >= 0)
                   #int_tip_rates = extract_tip_rates(int_tree, return_extinct = false)
                  # print("$edge_id : $nalive")
                   w = Array{Float64,1}(undef,0)
                   Ws = 0.
                   for r in 1:length(int_tip_rates)
                       #if !isnan(int_tip_rates[r])
                       push!(w,nu(int_tip_rates[r],1.))
                       Ws += w[end]
                       #end
                   end
                   for new_tip in sample(1:nalive, Weights(w))
                       current_n = nalive
                       k = 0
                       sampled = []
                       tip_trees = []
                       graft_tips = []
                       keep = false
                       keep_tip = false
                       add_to_tip_id = 0
                       new_tip_id = NaN
                       new_tip_id_rate = NaN
                       new_lik = 0.
                       new_parent_rate = 0.
                       new_lik_rate = 0.
                      # println(new_tip)

                       for i in 1:length(which_alive)
                           if which_alive[i]
                               k += 1
                               if k == new_tip
                                   new_parent_rate = int_tip_rates[k]
                                   new_lik_rate = 0. + log(Ws)
                                   new_tip_id = i
                                   new_tip_id_rate = k
                               end
                           end
                       end

                        if true
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
                            n_alive = 0
                            current_lik = new_lik_rate - former_lik
                            while (keep) & (i < length(which_alive))
                                i += 1
                                if which_alive[i]
                                    k += 1
                                    if k != new_tip
                                        tip_tree, current_lik = sim_ClaDS2_time(node_depth, σ, α, ε, int_tip_rates[k],u,current_lik,log(1-f),
                                            max_node_number = max_node_number)
                                        #p
                                        if tip_tree.branch_length < 0
                                            keep = false
                                        else
                                            push!(tip_trees, tip_tree)
                                            push!(graft_tips, i)
                                            if i < new_tip_id
                                                add_to_tip_id += n_tip(tip_tree) - 1
                                            end
                                        end
                                    end
                                end
                            end
                        end

                        if keep #& keep_tip
                            edge_trees[edge_id] = EdgeTreeRates(graft_to_tips(int_tree, [Tree(Array{Tree,1}(undef,0),0., [0.], 0, true); tip_trees],
                                [new_tip_id; graft_tips]), new_tip_id + add_to_tip_id, int_tip_rates[new_tip_id_rate], int_tip_rates)
                            #int_tip_rates[new_tip_id_rate]
                            former_lik = log(Ws)
                            trial = max_try
                            #println(" kept  ; $Ws")
                            #println("$edge_id ; $it ; $former_lik")
                            #break
                            #return  edge_trees[edge_id]
                        elseif true #& keep_if_any
                            trial = max_try
                            #println(" reject  ; $Ws")
                            #break
                            #return  edge_trees[edge_id]
                        end
                    end
                end
            end
        end
        return edge_trees[edge_id]
        println(trial)
    end
end


function change_edge_trees_rates!(edge_trees, tree, edge_id, σ, α, ε, f, rates, branch_lengths ; do_tips = true, max_try = 1_000, max_node_number = 50, keep_if_any = true)
    trial = 0

    daughter_rates = get_daughter_rates(tree, edge_id)

    if isnan(daughter_rates[1])
        if do_tips
            int_tree = Tree()
            n_alive = edge_trees[edge_id].tip_id#n_extant_tips(edge_trees[edge_id].tree)
            former_lik = 0
            while trial < max_try
                trial += 1
                u = log(rand())

                int_tree, unsampled, n_alive, tip_rates = sim_ClaDS2_time_unsampled_rates(branch_lengths[edge_id], σ, α, ε,
                    rates[edge_id],u,-1*former_lik,return_if_sampled = true,
                    max_node_number = max_node_number, not_sampled = false,
                    sampling_proba = f, return_if_extinct = false, make_tree = false)
                if unsampled
                    edge_trees[edge_id] = EdgeTreeRates(int_tree, n_alive, rates[edge_id], tip_rates)
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

            former_lik = 0. #-
                #((log(former_parent_rate) + log_α - log(daughter_rates[1]))^2 + (log(former_parent_rate) + log_α - log(daughter_rates[2]))^2)/(2 * σ^2) + log(former_parent_rate)
        else
            former_lik = (log(daughter_rates[1]) + log(daughter_rates[1]))/2 -log_α+(σ^2)/2
        end

        sum_rates = 0.
        ratio = former_parent_rate / edge_trees[edge_id].rate
        function nu(former_rate, former_ratio)
            if isnan(former_rate)
                return 0.
            end
            fr = former_rate*former_ratio
            log_fr = log(fr)
            return fr * exp((-(log_fr + log_α - log(daughter_rates[1]))^2 - (log_fr + log_α - log(daughter_rates[2]))^2)/(2 * σ^2))
        end

        for r in edge_trees[edge_id].rates
            sum_rates += nu(r, ratio)
        end
        former_lik += + log(sum_rates)

        while trial < max_try
            trial += 1
            int_tree, int_tip_rates, which_alive =
                sim_ClaDS2_time_rates(branch_lengths[edge_id], σ, α, ε, rates[edge_id],prune_extinct = false,
                max_node_number = max_node_number, make_tree = true, return_if_max = true, return_na = true)
            nalive = length(int_tip_rates)
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
            if int_tree.branch_length<=0
                return  edge_trees[edge_id]

           elseif (0 < nalive) #& (int_tree.branch_length >= 0)
               #int_tip_rates = extract_tip_rates(int_tree, return_extinct = false)
               #println(int_tip_rates)
               w = Array{Float64,1}(undef,0)
               Ws = 0.
               for r in 1:length(int_tip_rates)
                   #if !isnan(int_tip_rates[r])
                   push!(w,nu(int_tip_rates[r],1.))
                   Ws += w[end]
                   #end
               end
               for new_tip in sample(1:nalive, Weights(w))
                   current_n = nalive
                   k = 0
                   sampled = []
                   tip_trees = []
                   graft_tips = []
                   keep = false
                   keep_tip = false
                   add_to_tip_id = 0
                   new_tip_id = NaN
                   new_tip_id_rate = NaN
                   new_lik = 0.
                   new_parent_rate = 0.
                   new_lik_rate = 0.
                  # println(new_tip)

                  for i in 1:length(which_alive)
                      if which_alive[i]
                          k += 1
                          if k == new_tip
                              new_parent_rate = int_tip_rates[k]
                              new_lik_rate = 0. + log(Ws)
                              u = rand()
                              keep = log(u) < new_lik_rate - former_lik
                              keep_tip = deepcopy(keep)
                              new_tip_id = i
                              new_tip_id_rate = k
                          end
                      end
                  end

                    if true
                        former_tip_number = n_extant_tips(edge_trees[edge_id].tree) - 1

                        u = log(rand())

                        keep = true
                        new_lik = 0#new_lik_rate
                        keep = u < new_lik#-former_lik)

                        k = 0
                        i = 0
                        ntips = 0
                        n_alive = 0
                        current_lik = 0#new_lik_rate# - former_lik
                        while (keep) & (i < length(which_alive))
                            i += 1
                            if which_alive[i]
                                k += 1
                                if k != new_tip
                                    tip_tree, current_lik = sim_ClaDS2_time(node_depth, σ, α, ε, int_tip_rates[k],u,current_lik,log(1-f),
                                        max_node_number = max_node_number)
                                    #p
                                    if tip_tree.branch_length < 0
                                        keep = false
                                    else
                                        push!(tip_trees, tip_tree)
                                        push!(graft_tips, i)
                                        if i < new_tip_id
                                            add_to_tip_id += n_tip(tip_tree) - 1
                                        end
                                    end
                                end
                            end
                        end
                    end

                    if keep

                        if keep_tip
                            edge_trees[edge_id] = EdgeTreeRates(graft_to_tips(int_tree, [Tree(Array{Tree,1}(undef,0),0., [0.], 0, true); tip_trees],
                                [new_tip_id; graft_tips]), new_tip_id + add_to_tip_id, int_tip_rates[new_tip_id_rate], int_tip_rates)

                                return  edge_trees[edge_id]
                        else
                            return  edge_trees[edge_id]

                        end
                    end
                end
            end
        end

        return  edge_trees[edge_id]

        println(trial)
    end
end

function change_edge_trees_ratesR!(edge_trees, tree, edge_id, σ, α, ε, f, rates, branch_lengths ; do_tips = true, max_try = 1_000, max_node_number = 50, keep_if_any = true)
    trial = 0

    daughter_rates = get_daughter_rates(tree, edge_id)

    if isnan(daughter_rates[1])
        if do_tips
            int_tree = Tree()
            n_alive = edge_trees[edge_id].tip_id#n_extant_tips(edge_trees[edge_id].tree)
            former_lik = 0
            while trial < max_try
                trial += 1
                u = log(rand())

                int_tree, unsampled, n_alive, tip_rates = sim_ClaDS2_time_unsampled_rates(branch_lengths[edge_id], σ, α, ε,
                    rates[edge_id],u,-1*former_lik,return_if_sampled = true,
                    max_node_number = max_node_number, not_sampled = false,
                    sampling_proba = f, return_if_extinct = false, make_tree = false)
                if unsampled
                    edge_trees[edge_id] = EdgeTreeRates(int_tree, n_alive, rates[edge_id], tip_rates)
                    return edge_trees[edge_id]
                end
            end
        end
        return edge_trees[edge_id]
    else
        log_α = log(α)
        node_depth = get_node_depth(tree, edge_id)
        former_lik = (log(daughter_rates[1]) + log(daughter_rates[1]))/2 -log_α+(σ^2)/2

        sum_rates = 0.
        #ratio = former_parent_rate / edge_trees[edge_id].rate
        function nu(former_rate, former_ratio)
            if isnan(former_rate)
                return 0.
            end
            fr = former_rate*former_ratio
            log_fr = log(fr)
            #return fr * exp((-(log_fr + log_α - log(daughter_rates[1]))^2 - (log_fr + log_α - log(daughter_rates[2]))^2)/(2 * σ^2))
            return log_fr + ((-(log_fr + log_α - log(daughter_rates[1]))^2 - (log_fr + log_α - log(daughter_rates[2]))^2)/(2 * σ^2))
        end

        #for r in edge_trees[edge_id].rates
        #    sum_rates += nu(r, ratio)
        #end
        #former_lik += + log(sum_rates)

        while trial < max_try
            trial += 1
            int_tree, int_tip_rates, which_alive =
                sim_ClaDS2_time_rates(branch_lengths[edge_id], σ, α, ε, rates[edge_id],prune_extinct = false,
                max_node_number = max_node_number, make_tree = true, return_if_max = true, return_na = true)
            nalive = length(int_tip_rates)
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
            if int_tree.branch_length<=0
                return  edge_trees[edge_id]

           elseif (0 < nalive) #& (int_tree.branch_length >= 0)
               #int_tip_rates = extract_tip_rates(int_tree, return_extinct = false)
               #println(int_tip_rates)
               #w = Array{Float64,1}(undef,0)
               #Ws = 0.
               #for r in 1:length(int_tip_rates)
                   #if !isnan(int_tip_rates[r])
                   #push!(w,nu(int_tip_rates[r],1.))
                   #Ws += w[end]
                   #end
               #end
               for new_tip in sample(1:nalive)
                   current_n = nalive
                   k = 0
                   sampled = []
                   tip_trees = []
                   graft_tips = []
                   keep = false
                   keep_tip = false
                   add_to_tip_id = 0
                   new_tip_id = NaN
                   new_tip_id_rate = NaN
                   new_lik = 0.
                   new_parent_rate = 0.
                   new_lik_rate = 0.
                  # println(new_tip)

                  for i in 1:length(which_alive)
                      if which_alive[i]
                          k += 1
                          if k == new_tip
                              new_parent_rate = int_tip_rates[k]
                              new_lik_rate = nu(new_parent_rate, 1.)
                              u = rand()
                              keep = log(u) < new_lik_rate - former_lik
                              keep_tip = deepcopy(keep)
                              new_tip_id = i
                              new_tip_id_rate = k
                          end
                      end
                  end

                    if true
                        former_tip_number = n_extant_tips(edge_trees[edge_id].tree) - 1

                        u = log(rand())

                        keep = true
                        new_lik = new_lik_rate
                        keep = u < (new_lik-former_lik)

                        k = 0
                        i = 0
                        ntips = 0
                        n_alive = 0
                        current_lik = (new_lik-former_lik)
                        while (keep) & (i < length(which_alive))
                            i += 1
                            if which_alive[i]
                                k += 1
                                if k != new_tip
                                    tip_tree, current_lik = sim_ClaDS2_time(node_depth, σ, α, ε, int_tip_rates[k],u,current_lik,log(1-f),
                                        max_node_number = max_node_number)
                                    #p
                                    if tip_tree.branch_length < 0
                                        keep = false
                                    else
                                        push!(tip_trees, tip_tree)
                                        push!(graft_tips, i)
                                        if i < new_tip_id
                                            add_to_tip_id += n_tip(tip_tree) - 1
                                        end
                                    end
                                end
                            end
                        end
                    end

                    if keep
                        edge_trees[edge_id] = EdgeTreeRates(graft_to_tips(int_tree, [Tree(Array{Tree,1}(undef,0),0., [0.], 0, true); tip_trees],
                            [new_tip_id; graft_tips]), new_tip_id + add_to_tip_id, int_tip_rates[new_tip_id_rate], int_tip_rates)

                        return  edge_trees[edge_id]
                    end
                end
            end
        end

        return  edge_trees[edge_id]

        println(trial)
    end
end

function change_edge_trees_rrMH_old!(edge_trees, tree, edge_id, σ, α, ε, f, rates, branch_lengths, lefts;
    do_tips = true, max_try = 1_000, max_node_number = 50, keep_if_any = true, rep = 1)
    trial = 0
    if lefts[edge_id+1]>0
        daughter_rates = [rates[edge_id+2]; rates[edge_id+2+lefts[edge_id+1]]]
    else
        daughter_rates =[NaN,NaN]
    end
    parent_rate = get_parent_rate(tree, edge_id, edge_trees, rates)
    lambda_law = LogNormal(log(α),σ)

    if isnan(daughter_rates[1])
        function new_rates2(l)
            l * rand(lambda_law, 1)[1]
        end
        if do_tips
            int_tree = Tree()
            n_alive = edge_trees[edge_id].tip_id
            if n_alive == 0
                former_lik = -Inf
            elseif n_alive == 1
                former_lik = log(f)
            else
                former_lik = log(n_alive) + log(f) + (n_alive - 1) * log(1-f)
            end
            while trial < rep
                trial += 1
                ext = 0
                while ext < max_try
                    ext += 1
                    u = log(rand())
                    λ = new_rates2(parent_rate)
                    int_tree, unsampled, n_alive, tip_rates = sim_ClaDS2_time_unsampled_rates(branch_lengths[edge_id], σ, α, ε,
                        λ,u,-1*former_lik,return_if_sampled = true,
                        max_node_number = max_node_number, not_sampled = false,
                        sampling_proba = f, return_if_extinct = true, make_tree = false)
                    if (unsampled && int_tree.extant)
                        rates[edge_id+1] = λ
                        edge_trees[edge_id] = EdgeTreeRates(int_tree, n_alive, λ, tip_rates)
                        former_lik = log(n_alive) + log(f) + (n_alive - 1) * log(1-f)
                    elseif n_alive > 0
                        break
                    end
                end
            end
        end
        return edge_trees[edge_id]
    else
        function new_rates3(l)
            l * rand(lambda_law, 1)[1]
        end
        log_α = log(α)
        node_depth = get_node_depth(tree, edge_id)
        ratio = 1.
        #tr = Array{Float64,1}(undef,0)
        if edge_trees[edge_id].tree.n_nodes < 2
            #former_parent_rate = rates[edge_id + 1]
            #tr = [rates[edge_id + 1]]
            ratio = 1.
        else
            #edge_tree_tip_rates = extract_tip_rates(edge_trees[edge_id].tree, return_extinct = false)
            #former_parent_rate = edge_tree_tip_rates[edge_trees[edge_id].tip_id]
            ratio = rates[edge_id + 1] / edge_trees[edge_id].rate
            #println("$former_parent_rate ; $(edge_trees[edge_id].rates*ratio)")
        end
        former_lik = 0.
        sum_rates = 0.
        function nu(former_rate, former_ratio)
            if isnan(former_rate)
                return 0.
            end
            fr = former_rate*former_ratio
            log_fr = log(fr)
            #println("$edge_id, $fr, $former_ratio, $(daughter_rates./α)")
            return fr * exp((-(log_fr + log_α - log(daughter_rates[1]))^2 - (log_fr + log_α - log(daughter_rates[2]))^2)/(2 * σ^2))
        end

        for r in edge_trees[edge_id].rates
            sum_rates += nu(r, ratio)
        end
        #println(sum_rates )
        former_lik += + log(sum_rates)
        it = 0

        while it < rep
            it += 1
            ##println("$edge_id ; $it")
            trial = 0
            while trial < max_try
                nalive = 0
                int_tree = Tree()
                int_tip_rates = Array{Float64,1}(undef,0)
                which_alive = Array{Bool,1}(undef,0)
                λ = 1.

                while (nalive == 0) && (trial < max_try)
                    λ = new_rates3(parent_rate)
                    #println("$λ  $(rates[edge_id+1]) $parent_rate $α")
                    trial += 1
                    int_tree, int_tip_rates, which_alive =
                        sim_ClaDS2_time_rates(branch_lengths[edge_id], σ, α, ε, λ,prune_extinct = false,
                        max_node_number = max_node_number, make_tree = true, return_if_max = true, return_na = true)

                    if int_tree.branch_length<=0
                        nalive = Inf
                    else
                        nalive = length(int_tip_rates)
                    end
                end
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
                if int_tree.branch_length<=0
                    ##println("$edge_id $trial $it $rep")
                    break
                    #return  edge_trees[edge_id]

               elseif (0 < nalive) #& (int_tree.branch_length >= 0)
                   #int_tip_rates = extract_tip_rates(int_tree, return_extinct = false)
                  # #print("$edge_id : $nalive")
                   w = Array{Float64,1}(undef,0)
                   Ws = 0.
                   for r in 1:length(int_tip_rates)
                       #if !isnan(int_tip_rates[r])
                       push!(w,nu(int_tip_rates[r],1.))
                       Ws += w[end]
                       #end
                   end
                   #println("$w ; $Ws")
                   for new_tip in sample(1:nalive, Weights(w))
                       current_n = nalive
                       k = 0
                       sampled = []
                       tip_trees = []
                       graft_tips = []
                       keep = false
                       keep_tip = false
                       add_to_tip_id = 0
                       new_tip_id = NaN
                       new_tip_id_rate = NaN
                       new_lik = 0.
                       new_parent_rate = 0.
                       new_lik_rate = 0.
                      # #println(new_tip)

                       for i in 1:length(which_alive)
                           if which_alive[i]
                               k += 1
                               if k == new_tip
                                   new_parent_rate = int_tip_rates[k]
                                   new_lik_rate = 0. + log(Ws)
                                   new_tip_id = i
                                   new_tip_id_rate = k
                               end
                           end
                       end

                        if true
                            former_tip_number = n_extant_tips(edge_trees[edge_id].tree) - 1
                            if former_tip_number == 0
                                former_lik += 0
                            else
                                former_lik += former_tip_number * log(1-f)
                            end
                            u = log(rand())

                            keep = true
                            new_lik = new_lik_rate
                            #println(new_lik-former_lik)
                            keep = u < (new_lik-former_lik)

                            k = 0
                            i = 0
                            ntips = 0
                            n_alive = 0
                            current_lik = new_lik_rate - former_lik
                            ##println("$edge_id ; $current_lik ; $former_lik")
                            #println(keep)
                            while (keep) & (i < length(which_alive))
                                i += 1
                                if which_alive[i]
                                    k += 1
                                    #println("$i, $k, $new_tip")
                                    if k != new_tip
                                        #println(" rate $(int_tip_rates[k])")
                                        tip_tree, current_lik = sim_ClaDS2_time(node_depth, σ, α, ε, int_tip_rates[k],u,current_lik,log(1-f),
                                            max_node_number = max_node_number)
                                        #p
                                        if tip_tree.branch_length < 0
                                            keep = false
                                        else
                                            push!(tip_trees, tip_tree)
                                            push!(graft_tips, i)
                                            if i < new_tip_id
                                                add_to_tip_id += n_tip(tip_tree) - 1
                                            end
                                        end
                                    end
                                end
                            end
                        end
                        if keep #& keep_tip
                            rates[edge_id+1] = λ
                            edge_trees[edge_id] = EdgeTreeRates(graft_to_tips(int_tree, [Tree(Array{Tree,1}(undef,0),0., [0.], 0, true); tip_trees],
                                [new_tip_id; graft_tips]), new_tip_id + add_to_tip_id, λ, int_tip_rates)
                            former_lik = log(Ws)
                            trial = max_try
                        elseif true
                            trial = max_try
                        end
                    end
                end
            end
        end
        return edge_trees[edge_id]
        println(trial)
    end
end
