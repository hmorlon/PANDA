function newton_LambertW(x ; rel_tol=0.0001)
    w0 = 0
    if x<10
        w1 = 1
    else
        w1 = max(1,log(x)-log(log(x))+ log(log(x))/log(x))
    end

    while abs((w0-w1)/min(abs(w0),abs(w1))) > rel_tol
        w0 = w1
        w1 = w1 - (w1 - x * exp(-w1))/(1+w1)
    end
    return w1
end

function draw_σ(relative_rates, α ; α0 = 1., β0 = 0.01)
    n = length(relative_rates)
    m = mean(relative_rates)
    α_n = α0 + n/2
    log_α = log(α)#m #log(m)
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

function draw_α(relative_rates, σ ; α_0 = 0., σ_0 = 1)
    # know posterior for the mean of a normal with normal prior, see for example
    # www.ams.sunysb.edu/~zhu/ams570/Bayesian_Normal.pdf

    # so here we have a normal prior on log(α) (with default parameters (α_0,σ_0)=(0,1))
    n = length(relative_rates)
    #α_n = (α_0 * σ^2 + σ_0 * sum(map(x -> x^2, relative_rates)))/(σ^2 + n * σ_0^2)
    σ_n = sqrt((σ^2 * σ_0^2)/(σ^2 + n * σ_0^2))
    α_n = (α_0 * σ^2 + σ_0^2 * sum(relative_rates))/(σ^2 + n * σ_0^2)


    a = exp(rand(Normal(α_n, σ_n)))
    if a==0
        println("$α_n , $σ_n, $α_0 ; $σ ; $relative_rates")
    end
    return a
end

function draw_ε_crown(branch_lengths::Array{T,1}, rates::Array{T,1}, n_extinct::Int64, n_cond::Int64 ; n_it = 10) where {T<:Number}
    n = length(branch_lengths)
    S = 0.
    for i in 1:n
        S += branch_lengths[i]*rates[i]
    end
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

    return exp(λ)

end

function draw_ε_crown(tree::Tree, edge_trees::Union{Array{EdgeTree,1},Array{EdgeTreeRates,1}}, lefts::Array{Int64,1}; scale_ε = 1000)
    enhanced_branch_lengths = extract_branch_lengths(tree, edge_trees)
    enhanced_rates = extract_rates(tree, edge_trees)
    n_extinct = n_tip(tree, edge_trees) - n_extant_tips(tree, edge_trees)
    crown_edges = [1, 1 + lefts[1]]
    bl = [tree.offsprings[1].branch_length,tree.offsprings[2].branch_length]
    n_cond = 0
    if tree.offsprings[1].n_nodes > 1
        n_cond += 1
    end
    if tree.offsprings[2].n_nodes > 1
        n_cond += 1
    end
    for edge in 1:2
        if edge_trees[crown_edges[edge]].tree.n_nodes > 1
            ltt = LTT(edge_trees[crown_edges[edge]].tree)
            for e in 1:length(ltt[2])
                if ltt[1][e] >= bl[edge]
                    break
                end
                if ltt[2][e] == 1
                    n_cond += 1
                end
            end
        end
    end
    return draw_ε_crown(enhanced_branch_lengths[2:end], enhanced_rates[2:end], n_extinct, n_cond)
end

function draw_λ0_slicing!(tree::Tree, ε::Float64, edge_trees::Union{Array{EdgeTree,1},Array{EdgeTreeRates,1}}, enhanced_branch_lengths, enhanced_rates ; n_it = 5)

    branch_lengths = deepcopy(enhanced_branch_lengths)[2:end]
    rates = deepcopy(enhanced_rates)[2:end]
    former_λ0 = tree.attributes[1]
    n_events = 2 * n_tip(tree, edge_trees) - 1 - n_extant_tips(tree, edge_trees) - 1
    n = length(branch_lengths)
    S = 0.

    for i in 1:n
        S += branch_lengths[i] * rates[i]
    end
    S *= (1 + ε)/former_λ0

    best_λ0 = log(n_events/S)
    max_λ0 = best_λ0+20
    min_λ0 = best_λ0-20
    extant_λ0 = max_λ0 - min_λ0

    max_f = n_events * (best_λ0 - 1)

    function fx(x ; ne = n_events, s = S)
         exp(x * ne - s * exp(x) - max_f)
     end

    λ = best_λ0
    fλ = fx(λ)

    if -1e10 < max_f < 1e10
        for j in 1:n_it
            min_eff = min_λ0
            max_eff = max_λ0

            u = fλ * rand()
            reject = true
            if u < Inf
                while reject
                    λ = rand() * (max_eff - min_eff) + min_eff
                    if max_eff - min_eff <= 0
                        println(" oh oh $j, $(max_eff - min_eff), $fλ, $u, $max_f, $min_eff, $max_eff, $best_λ0, $min_λ0, $max_λ0")
                        break
                    end
                    fλ = fx(λ)
                    reject = fλ < u
                    if reject
                        if λ < best_λ0
                            min_eff = λ
                        else
                            max_eff = λ
                        end
                    end

                end
            else
                λ = log(former_λ0)
            end
        end
    else
        new_λ = former_λ0
    end

    new_λ0 = exp(λ)

    λ0_ratio = new_λ0 / former_λ0

    #rates .*= λ0_ratio#new_λ0

    for i in 1:length(edge_trees)
        #tip_tree, graft_tip, nalive = edge_trees[i]
        update_rates!(edge_trees[i].tree, λ0_ratio)
        #edge_trees[i] = (tip_tree, graft_tip, nalive)
    end

    update_rates!(tree, λ0_ratio)
    return tree, extract_rates(tree)
end

function draw_λ0_slicing_new!(tree::Tree, ε::Float64, edge_trees::Union{Array{EdgeTree,1},Array{EdgeTreeRates,1}}, enhanced_branch_lengths, enhanced_rates ; n_it = 5)

    branch_lengths = deepcopy(enhanced_branch_lengths)[2:end]
    rates = deepcopy(enhanced_rates)[2:end]
    former_λ0 = tree.attributes[1]
    n_events = 2 * n_tip(tree, edge_trees) - 1 - n_extant_tips(tree, edge_trees) - 1
    n = length(branch_lengths)
    S = 0.

    for i in 1:n
        S += branch_lengths[i] * rates[i]
    end
    S *= (1 + ε)/former_λ0

    best_λ0 = log(n_events/S)
    max_λ0 = best_λ0+20
    min_λ0 = best_λ0-20
    extant_λ0 = max_λ0 - min_λ0

    max_f = n_events * (best_λ0 - 1)

    function fx(x ; ne = n_events, s = S)
         exp(x * ne - s * exp(x) - max_f)
     end

    λ = best_λ0
    fλ = fx(λ)

    if -1e10 < max_f < 1e10
        for j in 1:n_it
            min_eff = min_λ0
            max_eff = max_λ0

            u = fλ * rand()
            reject = true
            if u < Inf
                while reject
                    λ = rand() * (max_eff - min_eff) + min_eff
                    if max_eff - min_eff <= 0
                        println(" oh oh $j, $(max_eff - min_eff), $fλ, $u, $max_f, $min_eff, $max_eff, $best_λ0, $min_λ0, $max_λ0")
                        break
                    end
                    fλ = fx(λ)
                    reject = fλ < u
                    if reject
                        if λ < best_λ0
                            min_eff = λ
                        else
                            max_eff = λ
                        end
                    end

                end
            else
                λ = log(former_λ0)
            end
        end
    else
        new_λ = former_λ0
    end

    new_λ0 = exp(λ)

    λ0_ratio = new_λ0 / former_λ0

    #rates .*= λ0_ratio#new_λ0

    for i in 1:length(edge_trees)
        #tip_tree, graft_tip, nalive = edge_trees[i]
        update_rates!(edge_trees[i].tree, λ0_ratio)
        #edge_trees[i] = (tip_tree, graft_tip, nalive)
    end

    update_rates!(tree, λ0_ratio)
    return tree, λ0_ratio
end

function draw_λ0!(tree::Tree, ε, edge_trees)
    enhanced_tree = graft_edge_trees(tree, edge_trees)
    enhanced_branch_lengths = extract_branch_lengths(enhanced_tree)
    enhanced_rates = extract_rates(enhanced_tree)

    draw_λ0_slicing!(tree, ε, edge_trees, enhanced_branch_lengths, enhanced_rates)
end


function draw_λi_slicing!(rates, edge_trees::Union{Array{EdgeTree,1},Array{EdgeTreeRates,1}}, i, σ::Float64, α::Float64, ε::Float64, tree, lefts; n_it = 5, bounds = 20)
    if i==0
        parent_λ = NaN
        daughter_λ = map(log,get_daughter_rates(tree, i, edge_trees, rates, lefts))

        bl_eff = 0.

        x0_eff = (-2*log(α)+daughter_λ[1]+daughter_λ[2])/2
        σ_eff = σ / sqrt(2)
        n = 0.
    else
        parent_λ = log(get_parent_rate(tree, i, edge_trees, rates))
        daughter_λ = map(log,get_daughter_rates(tree, i, edge_trees, rates, lefts))
        branch_length = edge_trees[i].tree.branch_length
        bl_eff = branch_length * (1 + ε)

        if isnan(daughter_λ[1])
            x0_eff = parent_λ+log(α)
            σ_eff = σ
            n = 0.
        else
            x0_eff = (-log(α)+parent_λ+daughter_λ[1]+daughter_λ[2])/3
            σ_eff = σ / sqrt(3)
            n = 1.
        end
    end

    #n = 0.
    a = x0_eff + n * σ_eff^2
    t_eff = bl_eff * σ_eff^2

    if (exp(a) * t_eff) < 10000
        best_λi = a - newton_LambertW(t_eff * exp(a), rel_tol=0.000001)
    else
        lt = log(t_eff) + a
        best_λi = a - lt + log(lt) - (log(lt)/lt)
    end
    max_λi = best_λi + bounds * σ_eff
    min_λi = best_λi - bounds * σ_eff
    extant_λi = max_λi - min_λi

    max_f = n * best_λi - bl_eff * exp(best_λi) - (best_λi - x0_eff)^2/(2 * σ_eff^2)
    j = 1.

    function fx(x ; bl = bl_eff, ne = n, x0 = x0_eff, s = 2*σ_eff^2)
         exp(x * ne - bl * exp(x) - (x-x0)^2/s - max_f)
     end

    λ = best_λi
    fλ = fx(λ)

    for j in 1:n_it
        min_eff = min_λi
        max_eff = max_λi

        u = fλ * rand()

        reject = true

        if u < 1
            n_iter = 0
            while reject
                n_iter += 1
                λ = rand() * (max_eff - min_eff) + min_eff
                if max_eff - min_eff <= 0
                    println("oh oh $(max_eff - min_eff), $fλ, $max_f")
                    break
                end
                fλ = fx(λ)
                reject = (fλ < u)
                if reject
                    if λ < best_λi
                        min_eff = λ
                    else
                        max_eff = λ
                    end
                end

            end
        else
            λ = best_λi
        end
    end

    new_λ = exp(λ)
    if isnan(new_λ) || new_λ <= 0
        println("it's here $i : $(rates[i + 1]) , $new_λ , $best_λi, $a, $t_eff, $x0_eff, $bl_eff ; $(log(α)), $α")
        new_λ = rates[i + 1]#exp(best_λi)
    end
    if i>0
        edge_trees[i].tree.attributes[1] = new_λ
    end
    rates[i + 1] = new_λ



    return new_λ
end

function draw_λi_rel!(rates::Array{Float64,1}, edge_trees::Union{Array{EdgeTree,1},Array{EdgeTreeRates,1}}, i::Int64,
    σ::Float64, α::Float64, ε::Float64, tree::Tree, lefts::Array{Int64,1}, parents::Array{Int64,1}; n_it = 5, bounds = 10)
    former_λ = rates[i+1]

    if i==0
        parent_λ = NaN
        daughter_λ = map(log,get_daughter_rates(tree, i, edge_trees, rates, lefts))

        bl_eff = 0.

        x0_eff = (-2*log(α)+daughter_λ[1]+daughter_λ[2])/2
        σ_eff = σ / sqrt(2)
        n = 0.
    else
        parent_λ = log(get_parent_rate(tree, edge_trees, rates, parents[i]))
        tip_rate = extract_tip_rates(edge_trees[i].tree)[edge_trees[i].tip_id]
        daughter_λ = map(log,get_daughter_rates(tree, i, edge_trees, rates, lefts, with_edge_tree = false)) .- log(tip_rate) .+ log(former_λ)
        branch_length = extract_branch_lengths(edge_trees[i].tree)

        bl_eff = sum(branch_length .* extract_rates(edge_trees[i].tree)) * (1 + ε) / former_λ

        if isnan(daughter_λ[1])
            x0_eff = parent_λ+log(α)
            σ_eff = σ
            n = edge_trees[i].tree.n_nodes - n_extant_tips(edge_trees[i].tree)
        else
            x0_eff = (-log(α)+parent_λ+daughter_λ[1]+daughter_λ[2])/3
            σ_eff = σ / sqrt(3)
            n = edge_trees[i].tree.n_nodes - n_extant_tips(edge_trees[i].tree) + 1. #- 1.
        end
    end

    a = x0_eff + n * σ_eff^2
    t_eff = bl_eff * σ_eff^2

    if (exp(a) * t_eff) < 10000
        best_λi = a - newton_LambertW(t_eff * exp(a), rel_tol=0.000001)
    else
        lt = log(t_eff) + a
        best_λi = a - lt + log(lt) - (log(lt)/lt)
    end
    max_λi = best_λi + bounds * σ_eff
    min_λi = best_λi - bounds * σ_eff
    extant_λi = max_λi - min_λi

    max_f = n * best_λi - bl_eff * exp(best_λi) - (best_λi - x0_eff)^2/(2 * σ_eff^2)
    j = 1.
    function fx(x ; bl = bl_eff, ne = n, x0 = x0_eff, s = 2*σ_eff^2)
         exp(x * ne - bl * exp(x) - (x-x0)^2/s - max_f)
     end

    λ = best_λi
    fλ = fx(λ)

    if -1e10 < max_f < 1e10
        for j in 1:n_it
            min_eff = min_λi
            max_eff = max_λi

            u = fλ * rand()
            reject = true
            if u < 1
                n_iter = 0
                while reject
                    n_iter += 1
                    λ = rand() * (max_eff - min_eff) + min_eff
                    fλ = fx(λ)
                    reject = (fλ < u)
                    if reject
                        if λ < best_λi
                            min_eff = λ
                        else
                            max_eff = λ
                        end
                    end

                end
            else
                λ = best_λi
            end
        end
    else
        λ = log(former_λ)
    end

    new_λ = exp(λ)
    if isnan(new_λ) || new_λ <= 0.
        new_λ = rates[i + 1]
    end
    λ_ratio = new_λ / former_λ
    rates[i + 1] = new_λ

    if i>0
        update_rates!(edge_trees[i].tree,λ_ratio)
    end

    return new_λ
end


function slicing(x0_eff, σ_eff, bl_eff, n, former_λ ; n_it = 5, bounds = 5)
    a = x0_eff + n * σ_eff^2
    t_eff = bl_eff * σ_eff^2

    if (exp(a) * t_eff) < 10000
        best_λi = a - newton_LambertW(t_eff * exp(a), rel_tol=0.000001)
    else
        lt = log(t_eff) + a
        best_λi = a - lt + log(lt) - (log(lt)/lt)
    end
    max_λi = best_λi + bounds * σ_eff
    min_λi = best_λi - bounds * σ_eff
    extant_λi = max_λi - min_λi

    max_f = n * best_λi - bl_eff * exp(best_λi) - (best_λi - x0_eff)^2/(2 * σ_eff^2)
    j = 1.
    function fx(x ; bl = bl_eff, ne = n, x0 = x0_eff, s = 2*σ_eff^2)
         exp(x * ne - bl * exp(x) - (x-x0)^2/s - max_f)
     end

    λ = best_λi
    fλ = fx(λ)

    if -1e10 < max_f < 1e10
        for j in 1:n_it
            min_eff = min_λi
            max_eff = max_λi

            u = fλ * rand()
            reject = true
            if u < 1
                n_iter = 0
                while reject
                    n_iter += 1
                    λ = rand() * (max_eff - min_eff) + min_eff
                    fλ = fx(λ)
                    reject = (fλ < u)
                    if reject
                        if λ < best_λi
                            min_eff = λ
                        else
                            max_eff = λ
                        end
                    end

                end
            else
                λ = best_λi
            end
        end
    else
        λ = log(former_λ)
        #new_λ = former_λ
    end

    new_λ = exp(λ)
    if isnan(new_λ) || new_λ <= 0.
        new_λ = former_λ#rates[i + 1]
    end

    return new_λ
end

function slicing0(S, n_events, former_λ0 ; n_it = 5, bounds = 5)

    best_λ0 = log(n_events/S)
    max_λ0 = best_λ0+20
    min_λ0 = best_λ0-20
    extant_λ0 = max_λ0 - min_λ0

    max_f = n_events * (best_λ0 - 1)

    function fx(x ; ne = n_events, s = S)
         exp(x * ne - s * exp(x) - max_f)
     end

    λ = best_λ0
    fλ = fx(λ)

    if -1e10 < max_f < 1e10
        for j in 1:n_it
            min_eff = min_λ0
            max_eff = max_λ0

            u = fλ * rand()
            reject = true
            if u < Inf
                while reject
                    λ = rand() * (max_eff - min_eff) + min_eff
                    if max_eff - min_eff <= 0
                        println(" oh oh $j, $(max_eff - min_eff), $fλ, $u, $max_f, $min_eff, $max_eff, $best_λ0, $min_λ0, $max_λ0")
                        break
                    end
                    fλ = fx(λ)
                    reject = fλ < u
                    if reject
                        if λ < best_λ0
                            min_eff = λ
                        else
                            max_eff = λ
                        end
                    end

                end
            else
                λ = log(former_λ0)
            end
        end
    else
        new_λ = former_λ0
    end

    new_λ0 = exp(λ)

    if isnan(new_λ0) || new_λ0 <= 0.
        new_λ0 = former_λ0#rates[i + 1]
    end

    return new_λ0
end

function draw_λi_quad!(rates::Array{Float64,1}, edge_trees::Array{EdgeTreeRates,1},
    σ::Float64, α::Float64, ε::Float64, tree::Tree, lefts::Array{Int64,1}, parents::Array{Int64,1}; n_it = 5, bounds = 10)

    function aux(subtree, id, parent_rate, sampled_rates)
        if subtree.n_nodes < 2
            former_λ = subtree.attributes[1]
            x0_eff = parent_rate+log(α)
            σ_eff = σ
            n = 0
            bl_eff = 0.
            if id>0
                if edge_trees[id].tree.n_nodes > 1
                    n += edge_trees[id].tree.n_nodes - n_extant_tips(edge_trees[id].tree)
                    branch_length = extract_branch_lengths(edge_trees[id].tree)
                    bl_eff = sum(branch_length .* extract_rates(edge_trees[id].tree)) * (1 + ε)
                else
                    bl_eff = subtree.branch_length *former_λ * (1 + ε)
                end
            end
            λ = slicing(x0_eff, σ_eff, bl_eff/former_λ, n, former_λ)
            pushfirst!(sampled_rates, λ)
            return  n, bl_eff*λ/former_λ
        else
            tip_λ = subtree.attributes[1]
            if id>0
                edge_tree_tip_rates = extract_tip_rates(edge_trees[id].tree, return_extinct = false)
                tip_λ = edge_tree_tip_rates[edge_trees[id].tip_id]
                #tip_λ = edge_trees[id].rate #??
            end
            tip_λ = log(tip_λ)
            n_right, bl_right = aux(subtree.offsprings[2], id + 1 + lefts[id+1], tip_λ, sampled_rates)
            n_left, bl_left = aux(subtree.offsprings[1], id + 1, tip_λ, sampled_rates)
            former_λ = subtree.attributes[1]

            x0_eff = parent_rate + log(α)
            σ_eff = σ
            n = n_left + n_right
            bl_eff = bl_right + bl_left
            λ = 0.
            if id>0
                if edge_trees[id].tree.n_nodes > 1
                    n += edge_trees[id].tree.n_nodes - n_extant_tips(edge_trees[id].tree) + 1.
                    branch_length = extract_branch_lengths(edge_trees[id].tree)
                    bl_eff += sum(branch_length .* extract_rates(edge_trees[id].tree)) * (1 + ε)
                else
                    bl_eff += subtree.branch_length * (1 + ε) *former_λ
                    n += 1
                end
                λ = slicing(x0_eff, σ_eff, bl_eff/former_λ, n, former_λ)

            else
                λ = slicing0(bl_eff/former_λ, n, former_λ)

            end

            for k in 1:(subtree.n_nodes - 1)
                sampled_rates[k] *= λ/former_λ
            end
            pushfirst!(sampled_rates, λ)
            return  n, bl_eff*λ/former_λ
        end
    end

    r = Array{Float64,1}(undef,0)
    aux(tree, 0, 0., r);
    rates[1] = r[1]

    for i in 1:length(edge_trees)
        λ_ratio = r[i+1] / edge_trees[i].tree.attributes[1]
        update_rates!(edge_trees[i].tree,λ_ratio)
        rates[i+1] = r[i+1]
    end

    update_rates!(tree,rates)
end

function draw_λi_slicing!(rates, edge_trees, i, σ::Float64, α::Float64, ε::Float64, tree; n_it = 5, bounds = 10)
    lefts = n_left(tree)
    draw_λi_slicing!(rates, edge_trees, i, σ, α, ε, tree, lefts, n_it = n_it, bounds = bounds)
end

function draw_λi_quad!(rates, edge_trees::Union{Array{EdgeTree,1},Array{EdgeTreeRates,1}}, σ::Float64, α::Float64, ε::Float64, tree; n_it = 5, bounds = 10)
    lefts = n_left(tree)
    draw_λi_quad!(rates, edge_trees, σ, α, ε, tree, lefts, n_it = n_it, bounds = bounds)
end

function draw_λi_rel!(rates::Array{Float64,1}, edge_trees::Union{Array{EdgeTree,1},Array{EdgeTreeRates,1}}, i::Int64, σ::Float64,
    α::Float64, ε::Float64, tree::Tree; n_it = 5, bounds = 10)
    lefts = n_left(tree)
    parents = get_parent_edges(tree)
    draw_λi_rel!(rates, edge_trees, i, σ, α, ε, tree, lefts, parents, n_it = n_it, bounds = bounds)
end
