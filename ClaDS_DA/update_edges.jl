function update_edges_ini!(tree::Tree, edge_trees::Union{Array{EdgeTree,1},Array{EdgeTreeRates,1}}, σ::Float64, α::Float64, ε::Float64, rates::Array{Float64,1} ;
    it_rates=1, scales = 1_000, with_ε = true, ε_version = 1, update_hp = true)
    lefts = n_left(tree)

    enhanced_branch_lengths = extract_branch_lengths(tree, edge_trees)
    enhanced_rates = extract_rates(tree, edge_trees)
    n_extinct = n_tip(tree, edge_trees) - n_extant_tips(tree, edge_trees)
    update_edges_ini!(tree, edge_trees, σ, α, ε, rates, lefts, enhanced_branch_lengths, enhanced_rates, n_extinct , scales = scales , with_ε = with_ε, it_rates = it_rates, ε_version = ε_version, update_hp = update_hp)
end

function update_edges_ini!(tree::Tree, edge_trees::Union{Array{EdgeTree,1},Array{EdgeTreeRates,1}}, σ::Float64, α::Float64, ε::Float64, rates::Array{Float64,1},
    lefts::Array{Int64,1}, enhanced_branch_lengths::Array{Float64,1}, enhanced_rates::Array{Float64,1}, n_extinct::Int64;
    scales = 1_000, with_ε = true, it_rates = 1, ε_version = 1, update_hp =true) where{T<:Number}
    tree, rates = draw_λ0_slicing!(tree, ε, edge_trees, enhanced_branch_lengths, enhanced_rates)

    new_ε = ε

    for j in 1:it_rates

        if with_ε
            if ε_version == 1
                new_ε = draw_ε(tree, edge_trees)
            elseif ε_version == 2
                new_ε = draw_ε_v2(tree, edge_trees)
                rates = rates .* (1+ε)/(1+new_ε)
                enhanced_rates = enhanced_rates.* (1+ε)/(1+new_ε)
                tree, edge_trees = update_rates(tree, (1+ε)/(1+new_ε), edge_trees)
            end
        else
            new_ε = ε
        end

         ε = new_ε

        for i in 0:(length(rates)-1)
            draw_λi_slicing!(rates, edge_trees, i, σ, α, new_ε, tree, lefts)
        end
        for i in (length(rates)-1):-1:0
            draw_λi_slicing!(rates, edge_trees, i, σ, α, new_ε, tree, lefts)
        end

        if update_hp
            relative_rates = extract_relative_rates(tree, edge_trees, rates)[2:end]
            #tree = update_rates(tree, rates)
            #enhanced_tree = graft_edge_trees(tree, edge_trees)
            #relative_rates = extract_relative_rates(enhanced_tree, edge_trees)[2:end]
            σ = draw_σ(relative_rates, α)
            α = draw_α(relative_rates, σ)
        end
    end

    return rates, new_ε, σ, α
end

function update_edges!(tree::Tree, edge_trees::Union{Array{EdgeTree,1},Array{EdgeTreeRates,1}}, σ::Float64, α::Float64, ε::Float64, rates::Array{Float64,1}, parents::Array{Int64,1}; it_rates=1, scales = 1_000, with_ε = true, ε_version = 2, update_hp = true)
    lefts = n_left(tree)

    enhanced_branch_lengths = extract_branch_lengths(tree, edge_trees)
    enhanced_rates = extract_rates(tree, edge_trees)
    n_extinct = n_tip(tree, edge_trees) - n_extant_tips(tree, edge_trees)
    update_edges!(tree, edge_trees, σ, α, ε, rates, lefts, enhanced_branch_lengths, enhanced_rates, n_extinct, parents , it_rates = it_rates)
end

function update_edges!(tree::Tree, edge_trees::Union{Array{EdgeTree,1},Array{EdgeTreeRates,1}}, σ::Float64, α::Float64, ε::Float64, rates::Array{Float64,1},
    lefts::Array{Int64,1}, enhanced_branch_lengths::Array{Float64,1}, enhanced_rates::Array{Float64,1}, n_extinct::Int64, parents::Array{Int64,1};
    it_rates = 1)
    new_ε = ε
    tree, ratio = draw_λ0_slicing_new!(tree, ε, edge_trees, enhanced_branch_lengths, enhanced_rates)
    for j in 1:length(rates)
        rates[j] *= ratio
    end

    new_ε = draw_ε_crown(tree, edge_trees, lefts)
    ε = new_ε


    for j in 1:it_rates

        for i in 1:(length(rates)-1)
            draw_λi_rel!(rates, edge_trees, i, σ, α, new_ε, tree, lefts, parents)
        end
        draw_λi_slicing!(rates, edge_trees, 0, σ, α, new_ε, tree, lefts)


    end

    update_rates!(tree, rates)

    relative_rates = extract_relative_rates(tree, edge_trees, rates)[2:end]
    σ = draw_σ(relative_rates, α, β0 = 0.04, α0 = 0.5)
    α = draw_α(relative_rates, σ, σ_0 = 1)
    return new_ε, σ, α
end
