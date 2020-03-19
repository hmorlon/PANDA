function gelman_est(mcmc::MCMClist, Npar::Int64; thin = 1, burn = 0)
    Nit = length(mcmc.chain1[1])
    ini = Int64(floor(1 + Nit * burn))

    it = [ini]
    i = ini + thin
    while i<= Nit
        push!(it,i)
        i += thin
    end
    Nit = length(it)

    # compute means
    m1 = zeros(Npar)
    m2 = zeros(Npar)
    m3 = zeros(Npar)

    for i in 1:Npar
        if ((i == 1 )|| (i == 3))
            for j in it
                m1[i] += sum(mcmc.chain1[i][j])
                m2[i] += sum(mcmc.chain2[i][j])
                m3[i] += sum(mcmc.chain3[i][j])
            end
        else
            for j in it
                m1[i] += sum(log(mcmc.chain1[i][j]))
                m2[i] += sum(log(mcmc.chain2[i][j]))
                m3[i] += sum(log(mcmc.chain3[i][j]))
            end
        end
    end

    m1 ./= Nit
    m2 ./= Nit
    m3 ./= Nit

    m = m1 .+ m2 .+ m3
    m ./= 3

    # compute vars
    σ1 = zeros(Npar)
    σ2 = zeros(Npar)
    σ3 = zeros(Npar)

    for i in 1:Npar
        if ((i == 1 )|| (i == 3))
            for j in it
                σ1[i] += (mcmc.chain1[i][j] - m1[i])^2
                σ2[i] += (mcmc.chain2[i][j] - m2[i])^2
                σ3[i] += (mcmc.chain3[i][j] - m3[i])^2
            end
        else
            for j in it
                σ1[i] += (log(mcmc.chain1[i][j]) - m1[i])^2
                σ2[i] += (log(mcmc.chain2[i][j]) - m2[i])^2
                σ3[i] += (log(mcmc.chain3[i][j]) - m3[i])^2
            end
        end

    end

    σ1 ./= (Nit-1)
    σ2 ./= (Nit-1)
    σ3 ./= (Nit-1)


    # between chains variance
    B = zeros(Npar)

    for i in 1:Npar
        B[i] = (m1[i] - m[i])^2 + (m2[i] - m[i])^2 + (m3[i] - m[i])^2
    end

    B .* 4/6


    # within chains variance
    W = σ1 .+ σ2 .+ σ3
    W ./= 3
    Wd = W .* (Nit-1)/Nit
    #println((Nit-1)/Nit)
    #println(minimum(W))
    #println(B ./ W)

    # pooled variance
    V = Wd .+ B

    # PSRF
    PSRF = V ./ W

    #return PSRF
    id = 0
    maxi = 0.

    for i in 1:Npar
        if PSRF[i] > maxi
            id = i
            maxi = PSRF[i]
        end
    end

    return id, sqrt(maxi)
end

function gelman_est_old(mcmc::Array{Array{Array{Float64,1},1},1}, Npar::Int64 ; thin = 1, burn = 0)
    mcmcList = MCMClist(mcmc)
    gelman_est(mcmcList, Npar ; thin = thin, burn = burn)
end

function gelman_est(mcmc::Array{Array{Array{Float64,1},1},1}, Npar::Int64 ; thin = 1, burn = 0)
    nChains = length(mcmc)
    Nit = length(mcmc[1][1])
    ini = Int64(floor(1 + Nit * burn))

    it = [ini]
    i = ini + thin
    while i<= Nit
        push!(it,i)
        i += thin
    end
    Nit = length(it)

    # compute means
    M = [zeros(Npar) for i in 1:nChains]

    for n in 1:nChains
        for i in 1:Npar
            if ((i == 1 )|| (i == 3))
                for j in it
                    M[n][i] += sum(mcmc[n][i][j])
                end
            else
                for j in it
                    M[n][i] += sum(log(mcmc[n][i][j]))
                end
            end
        end
    end

    M ./= Nit

    m = zeros(Npar)
    for n in 1:nChains
        m .+= M[n]
    end

    m ./= nChains

    # compute vars
    Σ = [zeros(Npar) for i in 1:nChains]

    for n in 1:nChains
        for i in 1:Npar
            if ((i == 1 )|| (i == 3))
                for j in it
                    Σ[n][i] += (mcmc[n][i][j] - M[n][i])^2
                end
            else
                for j in it
                    Σ[n][i] += (log(mcmc[n][i][j]) - M[n][i])^2
                end
            end

        end
    end

    Σ ./= (Nit-1)

    # between chains variance
    B = zeros(Npar)

    for n in 1:nChains
        for i in 1:Npar
            B[i] += (M[n][i] - m[i])^2
        end
    end

    B .* (nChains + 1)/(nChains * (nChains - 1))


    # within chains variance

    W = zeros(Npar)
    for n in 1:nChains
        W .+= Σ[n]
    end
    W ./= nChains

    Wd = W .* (Nit-1)/Nit


    # pooled variance
    V = Wd .+ B

    # PSRF
    PSRF = V ./ W

    #return PSRF
    id = 0
    maxi = 0.

    for i in 1:Npar
        if PSRF[i] > maxi
            id = i
            maxi = PSRF[i]
        end
    end

    return id, sqrt(maxi)
end

function gelman_est2(mcmc::Array{Array{Array{Float64,1},1},1}, Npar::Int64 ; α = 0.05, thin = 1, burn = 0)
    Nit = length(mcmc[1][1])
    ini = Int64(floor(1 + Nit * burn))

    it = [ini]
    i = ini + thin
    while i<= Nit
        push!(it,i)
        i += thin
    end
    Nit = length(it)
end
