struct MCMClist
    chain1::Array{Array{Float64,1},1}
    chain2::Array{Array{Float64,1},1}
    chain3::Array{Array{Float64,1},1}
end

function MCMClist(chains::Array{Array{Array{Float64,1},1},1})
    chain1 = chains[1]
    chain2 = chains[2]
    chain3 = chains[3]
    MCMClist(chain1, chain2, chain3)
end

function chains_to_R_coda(chains ; save_chain = false, file = "coda_chain.Rdata", max_it_number = Inf)
    @rput chains
    @rput save_chain
    @rput file
    @rput max_it_number

    reval("""
        library(coda)
        it_number = length(chains[[1]][[1]])
        if(it_number > max_it_number){
            rows = unique(floor(seq(1,it_number,length.out = max_it_number)))
        }else{
            rows = 1:it_number
        }
        coda_chain = mcmc.list(lapply(1:3, function(i){
            mcmc(sapply(1:length(chains[[1]]), function(j){
                chains[[i]][[j]][rows]
            }))}))

        if(save_chain){
            save(coda_chain, it_number, file = file)
            }
    """)
end

function R_coda_to_chains(file)
    @rput file

    reval("""
        library(coda)
        load(file)
        list = lapply(coda_chain, function(a){
            lapply(1:ncol(a), function(i){
                a[,i]
            })
        })
    """)

    @rget list
    chains = Array{Array{Array{Float64,1},1},1}(undef,3)
    for i in 1:3
        chains[i] = [list[i][j] for j in 1:length(list[i])]
    end
    return chains
end

function sampler_to_Rdata(tree, sampler, file ; sample_fraction = 1., max_it_number = Inf)

    plot_tree = Tree(tree.offsprings, 0., tree.attributes, tree.n_nodes)
    edges, branch_lengths, rates, tip_labels = make_ape(tree)
    ntip = (tree.n_nodes + 1)/2
    @rput edges
    @rput branch_lengths
    @rput rates
    @rput ntip
    @rput tip_labels
    reval("""
        tree = list(edge = edges, Nnode = ntip - 1, edge.lengths = branch_lengths, tip.labels = tip_labels)
        class(tree) = "phylo"
    """)

    chains = sampler[1][1]
    MAPS = sampler[2]
    npar = tree.n_nodes + 3

    @rput chains
    @rput file
    @rput MAPS
    @rput npar
    @rput sample_fraction
    @rput max_it_number

    reval("""
        library(coda)
        it_number = length(chains[[1]][[1]])
        if(it_number > max_it_number){
            rows = unique(floor(seq(1,it_number,length.out = max_it_number)))
        }else{
            rows = 1:it_number
        }
        coda_chain = mcmc.list(lapply(1:3, function(i){
            mcmc(sapply(1:length(chains[[1]]), function(j){
                chains[[i]][[j]][rows]
            }))}))
        MAPS[4:(npar+ntip)] = exp(MAPS[4:(npar+ntip)])
        save(tree, coda_chain, MAPS, npar, sample_fraction, it_number, file = file)
    """)
end

function add_to_chain!(chain, param)
    for i in 1:length(param)
        push!(chain[i], param[i])
    end
end

function initialize_ClaDS2_LTT(tree::Tree ; ini_par = [], initialize_rates = 0, ltt_steps = 10, enhance_method = "reject") where {N}

    new_tree = Tree(deepcopy(tree.offsprings), 0., deepcopy(tree.attributes))

    root_depth = maximum(node_depths(new_tree)) * 1.01
    live_nd = node_depths_base(new_tree)
    println(root_depth)
    ltt_times = [0:ltt_steps...] * root_depth/ltt_steps
    ltt_extant = LTT(tree,ltt_times)
    if length(ini_par) == 0
        ini_par = [[[10. ^(1-i),1,0.5];fill(0.000001, tree.n_nodes); fill(0.000001, Int64((tree.n_nodes+1)/2));[0, n_extant_tips(tree)] ; ltt_extant[2]] for i in 1:3]
    end

    if (enhance_method == "rates") || (enhance_method == "MHrates")
        edge_trees_s = [init_edge_tree_rates(update_rates(new_tree, ini_par[i][4:(tree.n_nodes + 3)]), ini_par[i][4:(end-2)]) for i in 1:3]
    else
        edge_trees_s = [init_edge_tree(update_rates(new_tree, ini_par[i][4:(tree.n_nodes + 3)]), ini_par[i][4:(end-2)]) for i in 1:3]
    end

    for i in 1:3
        new_tree = Tree(deepcopy(tree.offsprings), 0., deepcopy(tree.attributes))
        update_rates!(new_tree, ini_par[i][4:(tree.n_nodes + 3)])

        rates = ini_par[i][4:(tree.n_nodes + 3)]
        ε = ini_par[i][3]
        σ = ini_par[i][1]
        α = ini_par[i][2]
        edge_trees = edge_trees_s[i]
        for j in 1:initialize_rates
            rates, ε, σ, α = update_edges_ini!(new_tree, edge_trees, σ, α, ε, rates, it_rates = 1, with_ε = false)
            update_rates!(new_tree, rates)
        end
        ini_par[i][4:(tree.n_nodes + 3)] = rates
        ini_par[i][3] = ε
        edge_trees_s[i] = edge_trees
        relative_rates = extract_relative_rates(tree, edge_trees, rates)[2:end]
    end

    extant_branch_lengths = extract_branch_lengths(tree)
    trees = [deepcopy(update_rates(new_tree, ini_par[i][4:end])) for i in 1:3]
    chains = [[[deepcopy(ini_par[i][j])] for j in 1:length(ini_par[1])] for i in 1:3]
    param = deepcopy(ini_par)

    println(" ")
    return chains, param, edge_trees_s, trees, extant_branch_lengths, ltt_times, live_nd
end

function add_iter_ClaDS2_LTT(sampler, n_reccord::Int64; thin = 1, fs = 1., plot_tree = 0, print_state = 0,
    max_node_number = 1_000, max_try = 100_000, it_edge_tree = 1, print_all = false, it_rates = 1, enhance_method = "reject")
    chain_s, param_s, edge_trees_s, tree_s, extant_branch_lengths, ltt_times, live_nd = sampler

    tips_id = tips(tree_s[1])[2:end]
    lefts = n_left(tree_s[1])
    ntips = Int64((tree_s[1].n_nodes+1)/2)
    if plot_tree>0
        R"par(mfrow=c(3,4), mar=c(5,2,2,2))"
    end
    n_par = tree_s[1].n_nodes + 3
    for k in 1:3
        chain = chain_s[k]
        tree = tree_s[k]

        σ = param_s[k][1]
        α = param_s[k][2]
        ε = param_s[k][3]
        rates = param_s[k][4:(tree.n_nodes + 3)]
        edge_trees = edge_trees_s[k]
        lefts = n_left(tree)
        daughter_edges = daughters(tree,lefts)
        relative_rates = []
        parents = get_parent_edges(tree)
        for i in 1:n_reccord
            for j in 1:thin
                for l in 1:it_edge_tree
                    update_edge_trees!(edge_trees, tree, σ, α, ε, fs, rates, extant_branch_lengths, enhance_method= enhance_method,
                        max_node_number = max_node_number, max_try = max_try, keep_if_any = true,
                        do_tips = (l==1))
                end

                if it_rates > 0
                    ε, σ, α = update_edges!(tree, edge_trees, σ, α, ε, rates, parents, with_ε = true, it_rates = it_rates )
                end

                if plot_tree>0
                    if i%plot_tree == 0 && j == thin
                        enhanced_tree = graft_edge_trees(tree, edge_trees)
                        plot_ClaDS(enhanced_tree)#, extant_edges[2:end], ln = false)
                    end
                end
                if print_state>0
                    if i%print_state == 0 && (j == thin || print_all)
                        relative_rates = extract_relative_rates(tree, edge_trees, rates)
                        println("$i : sd $(sqrt(var(relative_rates))), mean $(mean(relative_rates)), σ $σ, α $α, ε $ε, λ0 $(rates[1]), n_ext $(n_extinct(tree,edge_trees)) $(n_tip(tree,edge_trees))  $(n_extant_tips(tree,edge_trees))")
                    end
                end
            end

            tip_rates = extract_tip_rates(tree, edge_trees, tips_id, rates)
            param_s[k][1] = σ
            param_s[k][2] = α
            param_s[k][3] = ε
            param_s[k][4:n_par] = rates
            param_s[k][(n_par+1):(n_par+ntips)] = tip_rates
            param_s[k][n_par + ntips + 1] = n_extinct(tree,edge_trees)
            param_s[k][n_par + ntips + 2] = n_extant_tips(tree,edge_trees)
            ltt = LTT(tree, edge_trees, ltt_times, live_nd)[2]
            param_s[k][(n_par + ntips + 3):end] = ltt
            edge_trees_s[k] = edge_trees
            add_to_chain!(chain_s[k], param_s[k])
            tree_s[k] = tree
        end
    end

    println(" ")
    return chain_s, param_s, edge_trees_s, tree_s, extant_branch_lengths, ltt_times, live_nd
end


function plot_coda(sampler ; burn = 0, thin = 1, id_par = [1:4...])
    chains = sampler[1]
    @rput chains
    @rput thin
    @rput burn
    @rput id_par

    reval("""
        require(coda)
        n_row = length(chains[[1]][[1]])
        ini = floor(burn * n_row)
        if(ini == 0) ini = 1
            it = seq(ini, n_row, thin)

        plot_chains = mcmc.list(lapply(1:3, function(i){
        mcmc(sapply((id_par), function(j){
        if(j<4) {
            chains[[i]][[j]][it]
        }else{
            log(chains[[i]][[j]][it])
        }
        }))}))

        plot(plot_chains)
    """)
end

function R_gelman(sampler ; burn = 0, thin = 1)
    chains = sampler[1]
    npar = sampler[4][1].n_nodes + 3
    @rput chains
    @rput thin
    @rput burn
    @rput npar
    reval("""
        require(coda)
        n_row = length(chains[[1]][[1]])
        ini = floor(burn * n_row)
        if(ini == 0) ini = 1
        it = seq(ini, n_row, thin)

        gelman_chains = mcmc.list(lapply(1:3, function(i){
          mcmc(sapply(1:npar, function(j){
            if(j<4) {
              chains[[i]][[j]][it]
            }else{
              log(chains[[i]][[j]][it])
            }
          }))}))

        gelman = try(gelman.diag(gelman_chains))
        if (inherits(gelman, "try-error")){
            id_gelman = 1
            gelman = 10
        }else{
            id_gelman = which.max(gelman[[1]][,1])
            gelman = gelman[[1]][id_gelman,1]
        }

    """)

    @rget id_gelman
    @rget gelman

    return id_gelman, gelman
end

function sample_fractions(tree, sampling_at_tips)
    function aux(subtree, sub_sampling, x)
        if subtree.n_nodes < 2
            s = pop!(sub_sampling)
            pushfirst!(x,s)
            return s
        else
            s_left = aux(subtree.offsprings[2], sub_sampling,x)
            n_left = (subtree.offsprings[2].n_nodes+1)/2
            s_right = aux(subtree.offsprings[1], sub_sampling,x)
            n_right = (subtree.offsprings[1].n_nodes+1)/2
            s = (n_left + n_right)/(n_left/s_left + n_right/s_right)
            pushfirst!(x,s)

            return s
        end
    end

    fs = Array{Float64,1}(undef,0)
    aux(tree, deepcopy(sampling_at_tips), fs)
    return fs
end

function run_ClaDS_LTT(tree::Tree, n_reccord::Int64; ini_par = [], initialize_rates = 1_000, goal_gelman = 1.05,
    thin = 10, burn = 1/4, f = 1., plot_tree = 0, print_state = 0, max_node_number = 2_000, plot_chain = false,
    max_try = 100_000, it_edge_tree = 1, print_all = false, it_rates = 1, sampler = [], plot_burn = NaN, ltt_steps = 50,
    save_as_R = false, Rfile = "coda.Rdata", max_it_number = Inf, enhance_method = "MHrates", end_it = Inf)

    ntips = Int64((tree.n_nodes + 1)/2)
    if isnan(plot_burn)
        plot_burn = burn
    end

    if length(sampler) == 0
        sampler = initialize_ClaDS2_LTT(tree , ini_par = ini_par, initialize_rates = initialize_rates, ltt_steps = ltt_steps, enhance_method = enhance_method)
    end

    if length(f) == 1
        fs = fill(f,tree.n_nodes-1)
    else
        fs = sample_fractions(tree,f)[2:end]
    end

    gelman = 10.
    MAPS=[]
    nit = 0
    while (gelman > goal_gelman) & (nit < end_it)
        sampler = add_iter_ClaDS2_LTT(sampler, n_reccord, thin = thin, fs = fs, plot_tree = plot_tree, print_state = print_state,
            max_node_number = max_node_number, max_try = max_try, it_edge_tree = it_edge_tree,
            print_all = print_all, it_rates = it_rates, enhance_method = enhance_method)

        nit += n_reccord
        if plot_chain
            plot_coda(sampler, burn = plot_burn, id_par=[1,2,3,4])
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
                    it = seq(ini, n_row, thin)
            """)
        end


        if goal_gelman == 0.
            id_gelman = 0
            gelman = 10
            chains = sampler[1]
            npar = sampler[4][1].n_nodes + 3
            @rput chains
            @rput npar
        else
            gelman = R_gelman(sampler, burn = burn)
            id_gelman = gelman[1]
            gelman = gelman[2]
        end

        @rput ntips

        reval("""
            npar2 = length(chains[[1]])
            map_chains = mcmc.list(lapply(1:3, function(i){
                mcmc(sapply(1:npar2, function(j){
                    if(j<4 || j>(npar+ntips)) {
                        chains[[i]][[j]][it]
                    }else{
                        log(chains[[i]][[j]][it])
                    }
            }))}))

            unlist_chains = (sapply(1:npar2, function(k){
                c(map_chains[[1]][,k], map_chains[[2]][,k], map_chains[[3]][,k])
                }))

            MAPS=sapply(1:npar2, function(i){D=density(unlist_chains[,i]);
                return(D[[1]][which.max(D[[2]])])})
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

    if goal_gelman > 0.
        gelm = R_gelman(sampler, burn = burn)
    end

    return sampler, MAPS, gelm
end
