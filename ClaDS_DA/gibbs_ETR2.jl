function draw_ε_crown(S::Float64, n_extinct::Int64, n_cond::Int64 ; n_it = 10) where {T<:Number}

    n_extinct_eff = n_extinct + 1  # because we want a flat prior on the nonloged value
    δ = sqrt((n_extinct_eff + n_cond - S)^2 + 4*n_extinct_eff*S)
    best_ε = log(max((n_extinct_eff + n_cond - S + δ)/(2*S),(n_extinct_eff + n_cond - S - δ)/(2*S)))
    max_f = best_ε * n_extinct_eff +
        log(exp(best_ε)+1) * n_cond -
        S * exp(best_ε)

    function fx(x ; ne = n_extinct_eff, nc = n_cond, s = S)
         exp(x * ne + log(exp(x)+1) * nc - s * exp(x) - max_f)
     end


    λ = best_ε
    fλ = fx(λ)

    for j in 1:n_it
        u = fλ * rand()
        reject = true

        min_eff = best_ε - 10
        max_eff = best_ε + 10

        while reject
            λ = rand() * (max_eff - min_eff) + min_eff
            fλ = fx(λ)
            reject = fλ < u

            if reject
                if λ < best_ε
                    min_eff = λ
                else
                    max_eff = λ
                end
            end
        end
    end

    if !isnan(λ)
        return exp(λ)
    else
        return 0.001
    end
end

function draw_ε_crown(tree::Tree, edge_trees::Array{EdgeTreeRates2,1}, lefts::Array{Int64,1})
    S = extract_S(tree, edge_trees)
    n_extinct = extract_nextinct(tree, edge_trees)
    crown_edges = [1, 1 + lefts[1]]

    n_cond = 0
    if tree.offsprings[1].n_nodes > 1
        n_cond += 1
    end
    if tree.offsprings[2].n_nodes > 1
        n_cond += 1
    end
    bl = [tree.offsprings[1].branch_length, tree.offsprings[2].branch_length]
    for edge in 1:2
        if edge_trees[crown_edges[edge]].tree.n_nodes > 1
            ltt = LTT(edge_trees[crown_edges[edge]].tree)
            s = edge_trees[crown_edges[edge]].scale
            b = bl[edge]*s
            for e in 1:length(ltt[2])
                if ltt[1][e] >= b
                    break
                end
                if ltt[2][e] == 1
                    n_cond += 1
                end
            end
        end
    end

    return draw_ε_crown(S, n_extinct, n_cond)
end

function draw_λi_quad!(rates::Array{Float64,1}, edge_trees::Array{EdgeTreeRates2,1},
    σ::Float64, α::Float64, ε::Float64, tree::Tree, lefts::Array{Int64,1}; n_it = 5, bounds = 10)

    log_α = log(α)
    function aux(subtree, id, parent_rate, sampled_rates)

        if subtree.n_nodes < 2
            former_λ = subtree.attributes[1]
            x0_eff = parent_rate + log_α

            n = 0
            bl_eff = 0.
            if id>0
                bl_eff = edge_trees[id].effective_bl * (1 + ε)
                #s = sum(extract_branch_lengths(edge_trees[id].tree) .*extract_rates(edge_trees[id].tree) / edge_trees[id].scale)
                n = edge_trees[id].tree.n_nodes - edge_trees[id].tip_number
            end

            λ = slicing(x0_eff, σ, bl_eff, n, former_λ)
            pushfirst!(sampled_rates, λ)


            return  n, bl_eff*λ

        else

            tip_λ = subtree.attributes[1]
            former_λ = subtree.attributes[1]
            if id>0
                tip_λ *= edge_trees[id].tip_rate
            end
            tip_λ = log(tip_λ)

            n_right, bl_right = aux(subtree.offsprings[2], id + 1 + lefts[id+1], tip_λ, sampled_rates)
            n_left, bl_left = aux(subtree.offsprings[1], id + 1, tip_λ, sampled_rates)

            x0_eff = parent_rate + log_α
            n = n_left + n_right
            bl_eff = (bl_right + bl_left) / former_λ

            λ = 0.

            if id>0
                bl_eff += edge_trees[id].effective_bl * (1 + ε)
                n += edge_trees[id].tree.n_nodes - edge_trees[id].tip_number

                λ = slicing(x0_eff, σ, bl_eff, n, former_λ)
            else
                λ = slicing0(bl_eff, n, former_λ)
            end

            ratio = λ/former_λ

            for k in 1:(subtree.n_nodes - 1)
                sampled_rates[k] *= ratio
            end

            pushfirst!(sampled_rates, λ)
            return  n, bl_eff*λ
        end
    end

    r = Array{Float64,1}(undef,0)
    aux(tree, 0, 0., r);
    rates[1] = r[1]

    for i in 1:length(edge_trees)
        rates[i+1] = r[i+1]
        edge_trees[i].stem_rate[1] = rates[i+1]
    end

    update_rates!(tree,rates)
end

function draw_λi_quad!(rates::Array{Float64,1}, edge_trees::Array{EdgeTreeRates2,1}, σ::Float64, α::Float64, ε::Float64, tree::Tree; n_it = 5, bounds = 10)
    lefts = n_left(tree)
    draw_λi_quad!(rates, edge_trees, σ, α, ε, tree, lefts, n_it = n_it, bounds = bounds)
end

function draw_λ0_slicing!(rates::Array{Float64,1}, edge_trees::Array{EdgeTreeRates2,1}, σ::Float64, α::Float64, ε::Float64, lefts; n_it = 5, bounds = 20)

    parent_λ = NaN

    daughter_λ = map(log,[rates[2],rates[2+lefts[1]]])
    x0_eff = (-2*log(α)+daughter_λ[1]+daughter_λ[2])/2
    σ_eff = σ / sqrt(2)

    law = Normal(x0_eff, σ_eff)
    new_λ = exp(rand(law,1)[1])
    rates[1] = new_λ

    return new_λ
end

function draw_λ0!(tree::Tree, ε::Float64, rates::Array{Float64,1},
    edge_trees::Array{EdgeTreeRates2,1})

    S = extract_S(tree, edge_trees)
    S *= 1+ε
    n = n_tip(tree, edge_trees)
    n_events = 2 * n - 2 - n_extant_tips(tree,edge_trees)
    former_λ0 = rates[1]
    S /= former_λ0
    λ0 = slicing0(S, n_events, former_λ0)

    #rates[1] = λ0
    ratio = λ0 / former_λ0
    rates .*= ratio

    update_rates!(tree, ratio)
    for i in 1:length(edge_trees)
        edge_trees[i].stem_rate[1] *= ratio
    end
end

function draw_σ2(relative_rates, α ; α0 = 1., β0 = 0.01)
    n = length(relative_rates)
    m = mean(relative_rates)
    α_n = α0 + n/2
    log_α = m #log(m)
    β_n = β0

    for rate in relative_rates
        β_n += (rate - log_α)^2 /2
    end

    if isnan(β_n)
        println("inv gam $α0 $β0; $α_n , $β_n ;    ")
        print(relative_rates)
    end
    sqrt(rand(InverseGamma(α_n,β_n)))
end
