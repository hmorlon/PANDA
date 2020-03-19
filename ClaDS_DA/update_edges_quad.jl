function update_edges_quad!(tree::Tree, edge_trees::Union{Array{EdgeTree,1},Array{EdgeTreeRates,1}}, σ::Float64, α::Float64, ε::Float64, rates::Array{Float64,1},
    lefts::Array{Int64,1}, parents::Array{Int64,1};
    it_rates = 1)
    new_ε = ε
    #relative_rates = extract_relative_rates(tree, edge_trees, rates)[2:end]
    #σ = draw_σ(relative_rates, α, β0 = 0.01, α0 = 0.12)
    #α = draw_α(relative_rates, σ, σ_0 = 1)


    for j in 1:it_rates

        new_ε = draw_ε_crown(tree, edge_trees, lefts)
        ε = new_ε
        draw_λi_quad!(rates, edge_trees, σ, α, new_ε, tree, lefts, parents)

        relative_rates = extract_relative_rates(tree, edge_trees, rates)[2:end]
        σ = draw_σ(relative_rates, α, β0 = 0.01, α0 = 0.12)
        α = draw_α(relative_rates, σ, σ_0 = 1)

    end
    return new_ε, σ, α


end

function update_edges_quad!(tree::Tree, edge_trees::Array{EdgeTreeRates,1}, σ::Float64, α::Float64, ε::Float64, rates::Array{Float64,1}, parents::Array{Int64,1}; it_rates=1, scales = 1_000, with_ε = true, ε_version = 2, update_hp = true)
    lefts = n_left(tree)
    update_edges_quad!(tree, edge_trees, σ, α, ε, rates, lefts, parents , it_rates = it_rates)
end

function update_edges_quad_hp!(tree::Tree, edge_trees::Union{Array{EdgeTree,1},Array{EdgeTreeRates,1}}, σ::Float64, α::Float64, ε::Float64, rates::Array{Float64,1},
    lefts::Array{Int64,1}, parents::Array{Int64,1};
    it_rates = 1)
    relative_rates = extract_relative_rates(tree, edge_trees, rates)[2:end]
    σ = draw_σ(relative_rates, α, β0 = 0.05, α0 = 0.5)
    #α = draw_α(relative_rates, σ,  α_0 = -(σ^2/2), σ_0 = 0.2)
    α = draw_α(relative_rates, σ,  α_0 = -0.05, σ_0 = 0.1)

    new_ε = draw_ε_crown(tree, edge_trees, lefts)
    #new_ε = ε
    for j in 1:it_rates

        draw_λi_quad!(rates, edge_trees, σ, α, new_ε, tree, lefts, parents)
        relative_rates = extract_relative_rates(tree, edge_trees, rates)[2:end]
        σ = draw_σ(relative_rates, α, β0 = 0.05, α0 = 0.5)
        #α = draw_α(relative_rates, σ,  α_0 = -(σ^2/2), σ_0 = 0.2)
        α = draw_α(relative_rates, σ,  α_0 = -0.05, σ_0 = 0.1)

        new_ε = draw_ε_crown(tree, edge_trees, lefts)
        ε = new_ε
    end


    return new_ε, σ, α
end

function update_edges_quad_hpbu!(tree::Tree, edge_trees::Union{Array{EdgeTree,1},Array{EdgeTreeRates,1}}, σ::Float64, α::Float64, ε::Float64, rates::Array{Float64,1},
    lefts::Array{Int64,1}, parents::Array{Int64,1};
    it_rates = 1)
    relative_rates = extract_relative_rates(tree, edge_trees, rates)[2:end]
    σ = draw_σ(relative_rates, α, β0 = 0.05, α0 = 0.5)
    #α = draw_α(relative_rates, σ,  α_0 = -(σ^2/2), σ_0 = 0.2)
    α = draw_α(relative_rates, σ,  α_0 = -0.05, σ_0 = 0.1)

    new_ε = draw_ε_crown(tree, edge_trees, lefts)
    #new_ε = ε
    for j in 1:it_rates

        for i in 1:(length(rates)-1)
        #    draw_λi_rel!(rates, edge_trees, i, σ, α, new_ε, tree, lefts, parents)
            draw_λi_slicing!(rates, edge_trees, i, σ, α, new_ε, tree, lefts)
        end
        draw_λi_slicing!(rates, edge_trees, 0, σ, α, new_ε, tree, lefts)
        update_rates!(tree, rates)
        relative_rates = extract_relative_rates(tree, edge_trees, rates)[2:end]
        σ = draw_σ(relative_rates, α, β0 = 0.05, α0 = 0.5)
        #α = draw_α(relative_rates, σ,  α_0 = -(σ^2/2), σ_0 = 0.2)
        α = draw_α(relative_rates, σ,  α_0 = -0.05, σ_0 = 0.1)

        new_ε = draw_ε_crown(tree, edge_trees, lefts)
        ε = new_ε
    end


    return new_ε, σ, α
end

function update_edges_quad_hp!(tree::Tree, edge_trees::Array{EdgeTreeRates,1}, σ::Float64, α::Float64, ε::Float64, rates::Array{Float64,1}, parents::Array{Int64,1}; it_rates=1, scales = 1_000, with_ε = true, ε_version = 2, update_hp = true)
    lefts = n_left(tree)
    update_edges_quad_hp!(tree, edge_trees, σ, α, ε, rates, lefts, parents , it_rates = it_rates)
end
