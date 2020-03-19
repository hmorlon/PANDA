function add_iter_ClaDS2_ETR2(sampler, n_reccord::Int64; thin = 1, fs = 1., plot_tree = 0, print_state = 0, quad = 1,
    max_node_number = 1_000, max_try = 100_000, it_edge_tree = 1, print_all = false, it_rates = 1, enhance_method = "reject")
    chain_s, param_s, edge_trees_s, tree_s, extant_branch_lengths, ltt_times, live_nd, mean_rates_chains = sampler

    n_chains = length(chain_s)
    tips_id = tips(tree_s[1])[2:end]
    lefts = n_left(tree_s[1])
    ntips = Int64((tree_s[1].n_nodes+1)/2)

    if plot_tree>0
        R"par(mfrow=c(3,4), mar=c(5,2,2,2))"
    end

    n_par = tree_s[1].n_nodes + 3

    for k in 1:n_chains
        chain = chain_s[k]
        tree = tree_s[k]

        σ = param_s[k][1]
        α = param_s[k][2]
        ε = param_s[k][3]
        rates = param_s[k][4:(tree.n_nodes + 3)]

        edge_trees = edge_trees_s[k]

        relative_rates = Array{Float64,1}(undef,0)

        for i in 1:n_reccord
            for j in 1:thin
                if quad>0
                    for l in 1:it_edge_tree
                        #println("$(edge_trees[5].stem_rate), $σ, $α, $ε")
                        update_edge_trees!(edge_trees, tree, σ, α, ε, fs, rates, extant_branch_lengths, enhance_method= enhance_method,
                            max_node_number = max_node_number, max_try = max_try, keep_if_any = true,
                            do_tips =  true)


                        if 1 < quad <4
                            relative_rates = extract_relative_rates(tree, edge_trees, rates)
                            #println("$(mean(relative_rates)) ; $(mean(extract_relative_rates(tree))) / $(extract_nextinct(tree,edge_trees))")
                            σ = draw_σ(relative_rates, α, β0 = 0.05, α0 = 0.5)
                            α = draw_α(relative_rates, σ,  α_0 = -0.05, σ_0 = 0.1)
                            ε = draw_ε_crown(tree, edge_trees, lefts)
                            draw_λ0_slicing!(rates, edge_trees, σ, α, ε, lefts)
                            tree.attributes[1] = rates[1]
                        elseif quad == 1
                            relative_rates = extract_relative_rates(tree, edge_trees, rates)
                            #println("$(mean(relative_rates)) ; $(mean(extract_relative_rates(tree))) / $(extract_nextinct(tree,edge_trees))")
                            draw_λ0_slicing!(rates, edge_trees, σ, α, ε, lefts)
                            σ = draw_σ2(relative_rates, α, β0 = 0.05, α0 = 0.5)
                            α = draw_α(relative_rates, σ,  α_0 = -0.05, σ_0 = 0.1)
                            ε = draw_ε_crown(tree, edge_trees, lefts)
                            draw_λ0!(tree::Tree, ε::Float64, rates::Array{Float64,1},
                                edge_trees::Array{EdgeTreeRates2,1})
                            tree.attributes[1] = rates[1]
                        elseif quad == 4
                            relative_rates = extract_relative_rates(tree, edge_trees, rates)
                            #println("$(mean(relative_rates)) ; $(mean(extract_relative_rates(tree))) / $(extract_nextinct(tree,edge_trees))")
                            α = draw_α(relative_rates, σ,  α_0 = -0.05, σ_0 = 0.1)
                            σ = draw_σ(relative_rates, α, β0 = 0.05, α0 = 0.5)
                            ε = draw_ε_crown(tree, edge_trees, lefts)
                            draw_λ0_slicing!(rates, edge_trees, σ, α, ε, lefts)
                            tree.attributes[1] = rates[1]
                        elseif quad == 5
                            relative_rates = extract_relative_rates(tree, edge_trees, rates)
                            #println("$(mean(relative_rates)) ; $(mean(extract_relative_rates(tree))) / $(extract_nextinct(tree,edge_trees))")
                            σ = draw_σ2(relative_rates, α, β0 = 0.05, α0 = 0.5)
                            α = draw_α(relative_rates, σ,  α_0 = -0.05, σ_0 = 0.1)
                            ε = draw_ε_crown(tree, edge_trees, lefts)
                            draw_λ0_slicing!(rates, edge_trees, σ, α, ε, lefts)
                            tree.attributes[1] = rates[1]
                        elseif quad == 6
                            relative_rates = extract_relative_rates(tree, edge_trees, rates)
                            #println("$(mean(relative_rates)) ; $(mean(extract_relative_rates(tree))) / $(extract_nextinct(tree,edge_trees))")
                            σ = draw_σ2(relative_rates, α, β0 = 0.05, α0 = 0.5)
                            α = draw_α(relative_rates, σ,  α_0 = -0.05, σ_0 = 0.1)
                            ε = draw_ε_crown(tree, edge_trees, lefts)
                        end
                    end
                else
                    update_edge_trees!(edge_trees, tree, σ, α, ε, fs, rates, extant_branch_lengths, enhance_method= enhance_method,
                        max_node_number = max_node_number, max_try = max_try, keep_if_any = true,rep = it_rates,
                        do_tips =  true)

                    if quad == 0
                        relative_rates = extract_relative_rates(tree, edge_trees, rates)
                        #println("$(mean(relative_rates)) ; $(mean(extract_relative_rates(tree))) / $(extract_nextinct(tree,edge_trees))")
                        σ = draw_σ2(relative_rates, α, β0 = 0.05, α0 = 0.5)
                        α = draw_α(relative_rates, σ,  α_0 = -0.05, σ_0 = 0.1)
                        ε = draw_ε_crown(tree, edge_trees, lefts)
                        draw_λ0_slicing!(rates, edge_trees, σ, α, ε, lefts)
                        tree.attributes[1] = rates[1]
                    end
                end

                if it_rates > 0
                    if quad == 1
                        ε, σ, α = update_edges_ETR2!(tree, edge_trees, σ, α, ε, rates,lefts,it_rates = it_rates )
                    elseif quad == 2
                        draw_λ0!(tree::Tree, ε::Float64, rates::Array{Float64,1},
                            edge_trees::Array{EdgeTreeRates2,1})
                    elseif quad == 3
                        for l in 1:it_rates
                            draw_λ0!(tree::Tree, ε::Float64, rates::Array{Float64,1},
                                edge_trees::Array{EdgeTreeRates2,1})
                            ε = draw_ε_crown(tree, edge_trees, lefts)
                        end
                    else
                        for l in 1:it_rates
                            draw_λ0!(tree::Tree, ε::Float64, rates::Array{Float64,1},
                                edge_trees::Array{EdgeTreeRates2,1})
                            ε = draw_ε_crown(tree, edge_trees, lefts)
                        end
                    end
                end

                if plot_tree>0
                    if i%plot_tree == 0 && j == thin
                        enhanced_tree = graft_edge_trees(tree, edge_trees)
                        plot_ClaDS(enhanced_tree)#, extant_edges[2:end], ln = false)
                    end
                end
                if print_state>0
                    if i%print_state == 0 && (j == thin || print_all)
                        relative_rates = extract_relative_rates(tree, edge_trees)
                        println("$i : sd $(sqrt(var(relative_rates))), mean $(mean(relative_rates)), σ $σ, α $α, ε $ε, λ0 $(rates[1]), n_ext $(extract_nextinct(tree,edge_trees)) $(n_tip(tree,edge_trees)) $(n_extant_tips(tree,edge_trees))")
                    end
                end

                #println(" 3 nnn")

            end
            #println("$σ ")

            tip_rates = extract_tip_rates(tree, edge_trees, tips_id, rates)
            param_s[k][1] = σ
            param_s[k][2] = α
            param_s[k][3] = ε
            param_s[k][4:n_par] = rates
            param_s[k][(n_par+1):(n_par+ntips)] = tip_rates
            param_s[k][n_par + ntips + 1] = extract_nextinct(tree,edge_trees)
            param_s[k][n_par + ntips + 2] = n_tip(tree,edge_trees) - extract_nextinct(tree,edge_trees)
            #enhanced_tree = graft_edge_trees(tree, edge_trees)
            #ltt = LTT(enhanced_tree, ltt_times)[2]
            ltt = LTT(tree, edge_trees, ltt_times)[2]
            param_s[k][(n_par + ntips + 3):end] = ltt
            edge_trees_s[k] = edge_trees
            add_to_chain!(chain_s[k], param_s[k])
            #println(time_rates(tree,edge_trees,ltt_times))
            #println(mean_rates_chains)
            for imr in 1:3
                #add_to_chain!(mean_rates_chains[k], time_rates(tree,edge_trees,ltt_times))
                push!(mean_rates_chains[k], time_rates(tree,edge_trees,ltt_times))
            end
            tree_s[k] = tree
        end
    end

    println(" ")
    return chain_s, param_s, edge_trees_s, tree_s, extant_branch_lengths, ltt_times, live_nd, mean_rates_chains
end

function run_ClaDS2(tree::Tree, n_reccord::Int64; ini_par = [], initialize_rates = 0, goal_gelman = 1.05,
    thin = 10, burn = 1/4, f = 1., plot_tree = 0, print_state = 0, max_node_number = 100, plot_chain = false,quad = 1,
    max_try = 10_000, it_edge_tree = 10, print_all = false, it_rates = 1, sampler = [], plot_burn = NaN, ltt_steps = 50,
    save_as_R = false, Rfile = "coda.Rdata", max_it_number = Inf, enhance_method = "MHrr", end_it = Inf, n_chains = 3)

    ntips = Int64((tree.n_nodes + 1)/2)
    if isnan(plot_burn)
        plot_burn = burn
    end

    if length(sampler) == 0
        sampler = initialize_ClaDS2_LTT(tree , ini_par = ini_par, initialize_rates = initialize_rates, ltt_steps = ltt_steps, enhance_method = enhance_method, n_chains = n_chains)
    end

    if length(f) == 1
        fs = fill(f,tree.n_nodes-1)
    elseif length(f) == ntips
        fs = sample_fractions(tree,f)[2:end]
    else
        fs = f
    end

    gelman = 10.
    MAPS=[]
    nit = 0
    while (gelman > goal_gelman) & (nit < end_it)
        sampler = add_iter_ClaDS2_ETR2(sampler, n_reccord, thin = thin, fs = fs, plot_tree = plot_tree, print_state = print_state,
            max_node_number = max_node_number, max_try = max_try, it_edge_tree = it_edge_tree,
            print_all = print_all, it_rates = it_rates, enhance_method = enhance_method, quad = quad)

        nit += n_reccord



        chains = sampler[1]
        npar = sampler[4][1].n_nodes + 3
        @rput chains
        @rput npar

        if goal_gelman == 0.
            id_gelman = 0
            gelman = 10
        else
            gelman = gelman_est(sampler[1], npar, burn = burn)
            id_gelman = gelman[1]
            gelman = gelman[2]
        end

        if plot_chain
            plot_coda(sampler, burn = plot_burn, id_par=[1,2,3,max(4,id_gelman)])
        else
            chains = sampler[1]
            @rput chains
            @rput thin
            @rput burn
            id_par = 1
            @rput id_par
            reval("""
                require(coda)
                n_row = length(chains[[1]][[1]])
                ini = floor(burn * n_row)
                if(ini == 0) ini = 1
                    it = seq(ini, n_row, 1)
            """)
        end

        @rput ntips
        @rput n_chains

        reval("""
            npar2 = length(chains[[1]])
            n_row = length(chains[[1]][[1]])
            ini = floor(burn * n_row)
            if(ini == 0) ini = 1
            it = c(seq(ini, n_row-1, 100),n_row)

            map_chains = mcmc.list(lapply(1:n_chains, function(i){
                mcmc(sapply(1:npar2, function(j){
                    if(j<=3 || j>(npar+ntips)) {
                        chains[[i]][[j]][it]
                    }else{
                        log(chains[[i]][[j]][it])
                    }
            }))}))

            unlist_chains = (sapply(1:npar2, function(k){
                x = c()
                for (n in 1:n_chains){
                    x = c(x,map_chains[[n]][,k] )
                    }
                return(x)
                }))

            MAPS=sapply(1:npar2, function(i){D=density(unlist_chains[,i]);
                return(D[[1]][which.max(D[[2]])])})
            #MAPS[2] = exp(MAPS[2])
        """)
        @rget MAPS
        if save_as_R
            chains = sampler[1]
            sampler_to_Rdata(tree, (sampler, MAPS), Rfile , sample_fraction = f, max_it_number = max_it_number)
        end
        println("iteration $(length(sampler[1][1][1])) ; gelman = $gelman for param $id_gelman")
        println("   σ = $(MAPS[1]), α = $(MAPS[2]), ε = $(MAPS[3]), λ0 = $(MAPS[4])")
        println("")
    end

    gelm = (0,0)
    npar = sampler[4][1].n_nodes + 3

    if goal_gelman > 0.
        gelm = gelman_est(sampler[1], npar, burn = burn)
    end

    @rput ntips
    @rput n_chains

    reval("""
        npar2 = length(chains[[1]])
        n_row = length(chains[[1]][[1]])
        ini = floor(burn * n_row)
        if(ini == 0) ini = 1
        it = seq(ini, n_row, 1)

        map_chains = mcmc.list(lapply(1:n_chains, function(i){
            mcmc(sapply(1:npar2, function(j){
                if(j==1 || j==3 || j>(npar+ntips)) {
                    chains[[i]][[j]][it]
                }else{
                    log(chains[[i]][[j]][it])
                }
        }))}))

        unlist_chains = (sapply(1:npar2, function(k){
            x = c()
            for (n in 1:n_chains){
                x = c(x,map_chains[[n]][,k] )
                }
            return(x)
            }))

        MAPS=sapply(1:npar2, function(i){D=density(unlist_chains[,i]);
            return(D[[1]][which.max(D[[2]])])})
        MAPS[2] = exp(MAPS[2])
    """)
    @rget MAPS
    return sampler, MAPS, gelm
end
