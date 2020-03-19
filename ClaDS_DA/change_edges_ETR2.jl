function change_edge_trees_rrMH!(edge_trees, tree, edge_id, σ, α, ε, f, rates, branch_lengths, lefts;
    max_try = 1_000, max_node_number = 50, rep = 1)

    log_α = log(α)
    lambda_law = LogNormal(log_α,σ)                                                        # relative new rate density

    trial = 0
    parent_rate = get_parent_rate(tree, edge_id, edge_trees, rates)
    #pr2 = rates[1]
    #if edge_trees[edge_id].parent_edge > 0
    #    pr2 = edge_trees[edge_trees[edge_id].parent_edge].stem_rate[1] *edge_trees[edge_trees[edge_id].parent_edge].tip_rate             # the parent rate in the augmented tree
    #end
    #println("$pr2, $parent_rate")
    depth = edge_trees[edge_id].depth
    sd = edge_trees[edge_id].stem_depth

    parent_edge = edge_trees[edge_id].parent_edge

    if lefts[edge_id+1] == 0                                                                   # meaning it's a tip branch

        function new_rates_tip(l)                                                                   # new rate density
            l * rand(lambda_law, 1)[1]
        end

        int_tree = Tree()                                                                       # initialize the tree
        n_alive = edge_trees[edge_id].tip_number                                                    # tip number in the current MCMC state

        if n_alive == 0                                                                         # compute the former likelihood of having exactly one tip in the reconstructed tree
            former_lik = -Inf
        elseif n_alive == 1
            former_lik = log(f)
        else
            former_lik = log(n_alive) + log(f) + (n_alive - 1) * log(1-f)
        end

        while trial < rep

            trial += 1
            ext = 0

            while ext < max_try                                                                 # conditionned to the non-extinction (if you don't want it put max_try to 1)
                ext += 1

                u = log(rand())
                λ = new_rates_tip(parent_rate)                                                      # draw the new stem_rate
                int_tree, unsampled, nalive, tip_rates =
                    sim_ClaDS2_time_unsampled_rates(branch_lengths[edge_id]*λ, σ, α, ε,
                        1.,u,-1*former_lik,
                        max_node_number = max_node_number,
                        sampling_proba = f, return_if_extinct = true)
                #println("$edge_id ; $trial , $ext , $nalive, $unsampled, $(unsampled && nalive > 0 )")

                if (unsampled && nalive > 0 )                                                       # one tip sampled exactly
                    rates[edge_id+1] = λ
                    edge_trees[edge_id] =
                        EdgeTreeRates2(int_tree, λ, 1., 1,
                            tip_rates, nalive, depth,sd, parent_edge)

                    #former_lik = log(nalive) + log(f) + (nalive - 1) * log(1-f)               # new likelihood
                    break
                elseif nalive > 0
                    break
                end
                #println(ext)
            end
        end

        return edge_trees[edge_id]


    else                                                                                    # meaning it's an internal branch
        daughter_rates = [rates[edge_id+2]; rates[edge_id+2+lefts[edge_id+1]]]                  # the daugther rates in that case
        stem_rate = rates[edge_id+1]
        function new_rates_int(l)                                                                   # new rate density
            l * rand(lambda_law, 1)[1]
        end

        function nu(former_rate, former_ratio)                                                  # speciating at t_i and giving birth to the new  daugther rates
            if isnan(former_rate)
                return 0.
            end
            fr = former_rate*former_ratio
            log_fr = log(fr)
            return fr * exp((-(log_fr + log_α - log(daughter_rates[1]))^2 -
                (log_fr + log_α - log(daughter_rates[2]))^2)/(2 * σ^2))
        end

        sum_rates = 0.
        for r in edge_trees[edge_id].tip_rates
            sum_rates += nu(r, stem_rate)
        end
        former_lik_rate = log(sum_rates)                                                        # Σ λνν
        former_tip_number = edge_trees[edge_id].tip_number
        former_lik = former_lik_rate

        if former_tip_number == 0
            former_lik += 0
        else
            former_lik += former_tip_number * log(1-f)                                          # no sampled tip
        end

        it = 0
        λ = 1.

        while it < rep
            it += 1
            trial = 0

            while trial < max_try
                nalive = 0
                int_tree = Tree()
                int_tip_rates = Array{Float64,1}(undef,0)
                which_alive = Array{Bool,1}(undef,0)

                while (nalive == 0) && (trial < max_try)

                    trial += 1
                    λ = new_rates_int(parent_rate)

                    int_tree, int_tip_rates, which_alive =
                        sim_ClaDS2_time_rates(branch_lengths[edge_id]*λ, σ, α, ε, 1.,
                            prune_extinct = false, max_node_number = max_node_number,
                            return_if_max = true, return_na = true)

                    if int_tree.branch_length<=0
                        nalive = Inf
                    else
                        nalive = length(int_tip_rates)
                    end
                    #=if nalive>1
                        println("$edge_id ; $trial $nalive $ε $(branch_lengths[edge_id]*λ) $σ $α")
                        println(int_tip_rates)
                        plot_ClaDS(int_tree)
                    end=#
                end


                if nalive == Inf
                    break
                end

               if (0 < nalive)

                   w = Array{Float64,1}(undef,0)
                   Ws = 0.
                   for r in 1:length(int_tip_rates)
                       push!(w,nu(int_tip_rates[r],λ))
                       Ws += w[end]
                   end
                   new_lik = log(Ws)                                                        # Σ λνν

                   new_tip = sample(1:nalive, Weights(w))


                   new_tip_id = NaN
                   new_parent_rate = 0.
                   former_tip_number = 0

                   k = 0
                   for i in 1:length(which_alive)
                       if which_alive[i]
                           k += 1
                           if k == new_tip
                               new_parent_rate = int_tip_rates[k]
                               new_tip_id = i
                               #which_alive[i] = false
                           end
                       end
                   end

                   u = log(rand())

                   keep = u < (new_lik-former_lik)


                   ntips = 0
                   n_alive = 0
                   current_lik = new_lik - former_lik

                   k = 0
                   i = 0
                   tip_trees = []
                   graft_tips = []
                   add_to_tip_id = 0
                   new_tip_number = 0

                   while keep & (i < length(which_alive))
                       i += 1

                       if which_alive[i]
                           k += 1
                           if k != new_tip
                               tip_tree, current_lik =
                                    sim_ClaDS2_time(depth*λ, σ, α, ε,
                                        int_tip_rates[k],u,
                                        current_lik,log(1-f),
                                        max_node_number = max_node_number)

                                   if tip_tree.branch_length < 0                            # more than max_node_number nodes, or sampled
                                       keep = false
                                   else
                                       push!(tip_trees, tip_tree)
                                       push!(graft_tips, i)
                                       nt = n_extant_tips(tip_tree)
                                       new_tip_number += nt

                                       if i < new_tip_id
                                           add_to_tip_id += n_tip(tip_tree) - 1
                                       end
                                   end

                                   #=if (keep && (tip_tree.n_nodes > 30))
                                       ep = ε
                                       nn = tip_tree.n_nodes
                                       @rput ep
                                       @rput nn
                                       plot_ClaDS(tip_tree)
                                       R"title(main = c(ep,nn))"
                                       println("$nalive, $(var(extract_relative_rates(int_tree))) : $(σ^2) ; $(mean(extract_relative_rates(int_tree))) : $(log(α))
                                        $(new_parent_rate * λ) , $daughter_rates, $(log(Ws) -former_lik)")
                                  end=#
                              end
                         end

                   end

                   if keep
                       rates[edge_id+1] = λ
                       to_graft = graft_to_tips(int_tree,
                            [Tree(Array{Tree,1}(undef,0),0., [0.], 0, true); tip_trees],    # the tree
                            [new_tip_id; graft_tips])
                       edge_trees[edge_id] =
                            EdgeTreeRates2(to_graft,
                                λ,                                                              # stem_rate
                                new_parent_rate,                                                # tip_rate
                                new_tip_id + add_to_tip_id,                                     # tip_id
                                int_tip_rates,                                                  # tip_rates
                                new_tip_number,
                                depth, sd, parent_edge)

                        if it < rep
                            former_lik = log(Ws)

                            if new_tip_number == 0
                                former_lik += 0
                            else
                                former_lik += new_tip_number * log(1-f)                                          # no sampled tip
                            end
                        end
                   end

                   trial = max_try
               end
           end
        end
        return edge_trees[edge_id]
    end
end

function change_edge_trees_rr!(edge_trees, tree, edge_id, σ, α, ε, f, rates, branch_lengths, lefts;
    max_try = 1_000, max_node_number = 50, rep = 1)

    log_α = log(α)
    lambda_law = LogNormal(log_α,σ)                                                        # relative new rate density

    trial = 0
    #parent_rate = get_parent_rate(tree, edge_id, edge_trees, rates)
    rate_i = rates[edge_id+1]
    #pr2 = rates[1]
    #if edge_trees[edge_id].parent_edge > 0
    #    pr2 = edge_trees[edge_trees[edge_id].parent_edge].stem_rate[1] *edge_trees[edge_trees[edge_id].parent_edge].tip_rate             # the parent rate in the augmented tree
    #end
    #println("$pr2, $parent_rate")
    depth = edge_trees[edge_id].depth
    sd = edge_trees[edge_id].stem_depth

    parent_edge = edge_trees[edge_id].parent_edge

    if lefts[edge_id+1] == 0                                                                   # meaning it's a tip branch

        int_tree = Tree()                                                                       # initialize the tree
        n_alive = edge_trees[edge_id].tip_number                                                    # tip number in the current MCMC state

        if n_alive == 0                                                                         # compute the former likelihood of having exactly one tip in the reconstructed tree
            former_lik = -Inf
        elseif n_alive == 1
            former_lik = log(f)
        else
            former_lik = log(n_alive) + log(f) + (n_alive - 1) * log(1-f)
        end

        while trial < rep

            trial += 1
            ext = 0

            while ext < max_try                                                                 # conditionned to the non-extinction (if you don't want it put max_try to 1)
                ext += 1

                u = log(rand())
                λ = rate_i#new_rates_tip(parent_rate)                                                      # draw the new stem_rate
                int_tree, unsampled, nalive, tip_rates =
                    sim_ClaDS2_time_unsampled_rates(branch_lengths[edge_id]*λ, σ, α, ε,
                        1.,u,-1*former_lik,
                        max_node_number = max_node_number,
                        sampling_proba = f, return_if_extinct = true)
                #println("$edge_id ; $trial , $ext , $nalive, $unsampled, $(unsampled && nalive > 0 )")

                if (unsampled && nalive > 0 )                                                       # one tip sampled exactly
                    rates[edge_id+1] = λ
                    edge_trees[edge_id] =
                        EdgeTreeRates2(int_tree, λ, 1., 1,
                            tip_rates, nalive, depth,sd, parent_edge)

                    former_lik = log(nalive) + log(f) + (nalive - 1) * log(1-f)               # new likelihood
                    break
                elseif nalive > 0
                    break
                end
                #println(ext)
            end
        end

        return edge_trees[edge_id]


    else                                                                                    # meaning it's an internal branch
        daughter_rates = [rates[edge_id+2]; rates[edge_id+2+lefts[edge_id+1]]]                  # the daugther rates in that case
        stem_rate = rates[edge_id+1]

        function nu(former_rate, former_ratio)                                                  # speciating at t_i and giving birth to the new  daugther rates
            if isnan(former_rate)
                return 0.
            end
            fr = former_rate*former_ratio
            log_fr = log(fr)
            return fr * exp((-(log_fr + log_α - log(daughter_rates[1]))^2 -
                (log_fr + log_α - log(daughter_rates[2]))^2)/(2 * σ^2))
        end

        sum_rates = 0.
        for r in edge_trees[edge_id].tip_rates
            sum_rates += nu(r, stem_rate)
        end
        former_lik_rate = log(sum_rates)                                                        # Σ λνν
        former_tip_number = edge_trees[edge_id].tip_number
        former_lik = former_lik_rate

        if former_tip_number == 0
            former_lik += 0
        else
            former_lik += former_tip_number * log(1-f)                                          # no sampled tip
        end

        it = 0
        λ = 1.

        while it < rep
            it += 1
            trial = 0

            while trial < max_try
                nalive = 0
                int_tree = Tree()
                int_tip_rates = Array{Float64,1}(undef,0)
                which_alive = Array{Bool,1}(undef,0)

                while (nalive == 0) && (trial < max_try)

                    trial += 1
                    λ = rate_i#new_rates_int(parent_rate)

                    int_tree, int_tip_rates, which_alive =
                        sim_ClaDS2_time_rates(branch_lengths[edge_id]*λ, σ, α, ε, 1.,
                            prune_extinct = false, max_node_number = max_node_number,
                            return_if_max = true, return_na = true)

                    if int_tree.branch_length<=0
                        nalive = Inf
                    else
                        nalive = length(int_tip_rates)
                    end
                    #=if nalive>1
                        println("$edge_id ; $trial $nalive $ε $(branch_lengths[edge_id]*λ) $σ $α")
                        println(int_tip_rates)
                        plot_ClaDS(int_tree)
                    end=#
                end


                if nalive == Inf
                    break
                end

               if (0 < nalive)

                   w = Array{Float64,1}(undef,0)
                   Ws = 0.
                   for r in 1:length(int_tip_rates)
                       push!(w,nu(int_tip_rates[r],λ))
                       Ws += w[end]
                   end
                   new_lik = log(Ws)                                                        # Σ λνν

                   new_tip = sample(1:nalive, Weights(w))


                   new_tip_id = NaN
                   new_parent_rate = 0.
                   former_tip_number = 0

                   k = 0
                   for i in 1:length(which_alive)
                       if which_alive[i]
                           k += 1
                           if k == new_tip
                               new_parent_rate = int_tip_rates[k]
                               new_tip_id = i
                               #which_alive[i] = false
                           end
                       end
                   end

                   u = log(rand())

                   keep = u < (new_lik-former_lik)


                   ntips = 0
                   n_alive = 0
                   current_lik = new_lik - former_lik

                   k = 0
                   i = 0
                   tip_trees = []
                   graft_tips = []
                   add_to_tip_id = 0
                   new_tip_number = 0

                   while keep & (i < length(which_alive))
                       i += 1

                       if which_alive[i]
                           k += 1
                           if k != new_tip
                               tip_tree, current_lik =
                                    sim_ClaDS2_time(depth*λ, σ, α, ε,
                                        int_tip_rates[k],u,
                                        current_lik,log(1-f),
                                        max_node_number = max_node_number)

                                   if tip_tree.branch_length < 0                            # more than max_node_number nodes, or sampled
                                       keep = false
                                   else
                                       push!(tip_trees, tip_tree)
                                       push!(graft_tips, i)
                                       nt = n_extant_tips(tip_tree)
                                       new_tip_number += nt

                                       if i < new_tip_id
                                           add_to_tip_id += n_tip(tip_tree) - 1
                                       end
                                   end

                                   #=if (keep && (tip_tree.n_nodes > 30))
                                       ep = ε
                                       nn = tip_tree.n_nodes
                                       @rput ep
                                       @rput nn
                                       plot_ClaDS(tip_tree)
                                       R"title(main = c(ep,nn))"
                                       println("$nalive, $(var(extract_relative_rates(int_tree))) : $(σ^2) ; $(mean(extract_relative_rates(int_tree))) : $(log(α))
                                        $(new_parent_rate * λ) , $daughter_rates, $(log(Ws) -former_lik)")
                                  end=#
                              end
                         end

                   end

                   if keep
                       rates[edge_id+1] = λ
                       to_graft = graft_to_tips(int_tree,
                            [Tree(Array{Tree,1}(undef,0),0., [0.], 0, true); tip_trees],    # the tree
                            [new_tip_id; graft_tips])
                       edge_trees[edge_id] =
                            EdgeTreeRates2(to_graft,
                                λ,                                                              # stem_rate
                                new_parent_rate,                                                # tip_rate
                                new_tip_id + add_to_tip_id,                                     # tip_id
                                int_tip_rates,                                                  # tip_rates
                                new_tip_number,
                                depth, sd, parent_edge)

                        if it < rep
                            former_lik = log(Ws)

                            if new_tip_number == 0
                                former_lik += 0
                            else
                                former_lik += new_tip_number * log(1-f)                                          # no sampled tip
                            end
                        end
                   end

                   trial = max_try
               end
           end
        end
        return edge_trees[edge_id]
    end
end

#=
function change_edge_trees_rr!(edge_trees, tree, edge_id, σ, α, ε, f, rates, branch_lengths, lefts;
    max_try = 1_000, max_node_number = 50, rep = 1)

    log_α = log(α)
    lambda_law = LogNormal(log_α,σ)                                                        # relative new rate density

    trial = 0
    parent_rate = get_parent_rate(tree, edge_id, edge_trees, rates)
    #pr2 = rates[1]
    #if edge_trees[edge_id].parent_edge > 0
    #    pr2 = edge_trees[edge_trees[edge_id].parent_edge].stem_rate[1] *edge_trees[edge_trees[edge_id].parent_edge].tip_rate             # the parent rate in the augmented tree
    #end
    #println("$pr2, $parent_rate")
    depth = edge_trees[edge_id].depth
    sd = edge_trees[edge_id].stem_depth

    parent_edge = edge_trees[edge_id].parent_edge

    if lefts[edge_id+1] == 0                                                                   # meaning it's a tip branch

        function new_rates_tip(l)                                                                   # new rate density
            l * rand(lambda_law, 1)[1]
        end

        int_tree = Tree()                                                                       # initialize the tree
        n_alive = edge_trees[edge_id].tip_number                                                    # tip number in the current MCMC state

        former_lik = 0

        while trial < rep

            trial += 1
            ext = 0

            while ext < max_try                                                                 # conditionned to the non-extinction (if you don't want it put max_try to 1)
                ext += 1

                u = log(rand())
                λ = new_rates_tip(parent_rate)                                                      # draw the new stem_rate
                int_tree, unsampled, nalive, tip_rates =
                    sim_ClaDS2_time_unsampled_rates(branch_lengths[edge_id]*λ, σ, α, ε,
                        1.,u,-1*former_lik,
                        max_node_number = max_node_number,
                        sampling_proba = f, return_if_extinct = true)
                #println("$edge_id ; $trial , $ext , $nalive, $unsampled, $(unsampled && nalive > 0 )")

                if (unsampled && nalive > 0 )                                                       # one tip sampled exactly
                    rates[edge_id+1] = λ
                    edge_trees[edge_id] =
                        EdgeTreeRates2(int_tree, λ, 1., 1,
                            tip_rates, nalive, depth,sd, parent_edge)

                    former_lik = 0#log(nalive) + log(f) + (nalive - 1) * log(1-f)               # new likelihood
                    break
                elseif nalive > 0
                    #break
                end
                #println(ext)
            end
        end

        return edge_trees[edge_id]


    else                                                                                    # meaning it's an internal branch
        daughter_rates = [rates[edge_id+2]; rates[edge_id+2+lefts[edge_id+1]]]                  # the daugther rates in that case
        stem_rate = rates[edge_id+1]
        function new_rates_int(l)                                                                   # new rate density
            l * rand(lambda_law, 1)[1]
        end

        function nu(former_rate, former_ratio)                                                  # speciating at t_i and giving birth to the new  daugther rates
            if isnan(former_rate)
                return 0.
            end
            fr = former_rate*former_ratio
            log_fr = log(fr)
            return fr * exp((-(log_fr + log_α - log(daughter_rates[1]))^2 -
                (log_fr + log_α - log(daughter_rates[2]))^2)/(2 * σ^2))
        end

        sum_rates = 0.
        for r in edge_trees[edge_id].tip_rates
            sum_rates += nu(r, stem_rate)
        end
        former_lik_rate = log(sum_rates)                                                        # Σ λνν
        former_tip_number = edge_trees[edge_id].tip_number
        former_lik = former_lik_rate

        #=if former_tip_number == 0
            former_lik += 0
        else
            former_lik += former_tip_number * log(1-f)                                          # no sampled tip
        end=#

        it = 0
        λ = 1.

        while it < rep
            it += 1
            trial = 0

            while trial < max_try
                nalive = 0
                int_tree = Tree()
                int_tip_rates = Array{Float64,1}(undef,0)
                which_alive = Array{Bool,1}(undef,0)

                while (nalive == 0) && (trial < max_try)

                    trial += 1
                    λ = new_rates_int(parent_rate)

                    int_tree, int_tip_rates, which_alive =
                        sim_ClaDS2_time_rates(branch_lengths[edge_id]*λ, σ, α, ε, 1.,
                            prune_extinct = false, max_node_number = max_node_number,
                            return_if_max = true, return_na = true)

                    if int_tree.branch_length<=0
                        nalive = Inf
                    else
                        nalive = length(int_tip_rates)
                    end
                    #=if nalive>1
                        println("$edge_id ; $trial $nalive $ε $(branch_lengths[edge_id]*λ) $σ $α")
                        println(int_tip_rates)
                        plot_ClaDS(int_tree)
                    end=#
                end


                if nalive == Inf
                    #break
                elseif (0 < nalive)

                   w = Array{Float64,1}(undef,0)
                   Ws = 0.
                   for r in 1:length(int_tip_rates)
                       push!(w,nu(int_tip_rates[r],λ))
                       Ws += w[end]
                   end
                   new_lik = log(Ws)                                                        # Σ λνν

                   new_tip = sample(1:nalive, Weights(w))


                   new_tip_id = NaN
                   new_parent_rate = 0.
                   former_tip_number = 0

                   k = 0
                   for i in 1:length(which_alive)
                       if which_alive[i]
                           k += 1
                           if k == new_tip
                               new_parent_rate = int_tip_rates[k]
                               new_tip_id = i
                               #which_alive[i] = false
                           end
                       end
                   end

                   u = log(rand())

                   keep = u < (new_lik-former_lik)


                   ntips = 0
                   n_alive = 0
                   current_lik = new_lik - former_lik

                   k = 0
                   i = 0
                   tip_trees = []
                   graft_tips = []
                   add_to_tip_id = 0
                   new_tip_number = 0

                   while keep & (i < length(which_alive))
                       i += 1

                       if which_alive[i]
                           k += 1
                           if k != new_tip
                               tip_tree, current_lik =
                                    sim_ClaDS2_time(depth*λ, σ, α, ε,
                                        int_tip_rates[k],u,
                                        current_lik,log(1-f),
                                        max_node_number = max_node_number)

                                   if tip_tree.branch_length < 0                            # more than max_node_number nodes, or sampled
                                       keep = false
                                   else
                                       push!(tip_trees, tip_tree)
                                       push!(graft_tips, i)
                                       nt = n_extant_tips(tip_tree)
                                       new_tip_number += nt

                                       if i < new_tip_id
                                           add_to_tip_id += n_tip(tip_tree) - 1
                                       end
                                   end

                                   #=if (keep && (tip_tree.n_nodes > 30))
                                       ep = ε
                                       nn = tip_tree.n_nodes
                                       @rput ep
                                       @rput nn
                                       plot_ClaDS(tip_tree)
                                       R"title(main = c(ep,nn))"
                                       println("$nalive, $(var(extract_relative_rates(int_tree))) : $(σ^2) ; $(mean(extract_relative_rates(int_tree))) : $(log(α))
                                        $(new_parent_rate * λ) , $daughter_rates, $(log(Ws) -former_lik)")
                                  end=#
                              end
                         end

                   end

                   if keep
                       rates[edge_id+1] = λ
                       to_graft = graft_to_tips(int_tree,
                            [Tree(Array{Tree,1}(undef,0),0., [0.], 0, true); tip_trees],    # the tree
                            [new_tip_id; graft_tips])
                       edge_trees[edge_id] =
                            EdgeTreeRates2(to_graft,
                                λ,                                                              # stem_rate
                                new_parent_rate,                                                # tip_rate
                                new_tip_id + add_to_tip_id,                                     # tip_id
                                int_tip_rates,                                                  # tip_rates
                                new_tip_number,
                                depth, sd, parent_edge)

                        trial = max_try
                        if it < rep
                            former_lik = log(Ws)

                            #=if new_tip_number == 0
                                former_lik += 0
                            else
                                former_lik += new_tip_number * log(1-f)                                          # no sampled tip
                            end=#
                        end
                   end

                   #
               end
           end
        end
        return edge_trees[edge_id]
    end
end
=#
