#=
load packages and other stuffs
=#

#=
first a function to construct the tree from a list of parents and offsprings
=#

function build_tree(branches::Array{Int64,2}, branch_lengths::Array{Float64,2}, attributes::Array{Array{T2,1},1} ;
    extinct = [], prune_extinct=true, return_void = false, stem_age=0., root_attributes = Array{T2,1}(undef,0)) where {T2<:Number}
    n_attributes = length(attributes)
    if length(extinct) == 0
        dead_branches = fill(false, size(branches)[1])
    else
        dead_branches = extinct
    end

    if size(branches)[1] == 0
        if return_void
            return Tree()
        else
            return Tree(Array{Tree,1}(undef,0),stem_age,root_attributes)
        end
    end
    if size(branches)[1] == 1
        if ! dead_branches[1]
            return Tree(Array{Tree,1}(undef,0),branch_lengths[1],root_attributes, true)
        elseif prune_extinct
            return Tree()
        else
            return Tree(Array{Tree,1}(undef,0),branch_lengths[1],root_attributes, false)
        end
    end

    root = -1
    for i in branches[:,1]
        has_no_parent = true
        for j in branches[:,2]
            if i == j
                has_no_parent = false
                break
            end
        end
        if has_no_parent
            root = i
            break
        end
    end

    max_node = maximum(branches[:,1:2])+1
    offsprings = fill(-1,max_node,2)
    parent_edges = fill(-1,max_node)
    daughter_edges = fill(-1,max_node)

    for i in 1:size(branches)[1]
        individual = branches[i,1]+1
        parent_edges[individual] = i
        daughter_edges[branches[i,2]+1] = i
        if offsprings[individual,1] < 0
            offsprings[individual,1] = branches[i,2]
        else
            offsprings[individual,2] = branches[i,2]
        end

        if dead_branches[i]
            individual = branches[i,2]+1
            offsprings[individual,:]= [-2 -2]
        end
    end

    if prune_extinct
        function aux(node)
            if offsprings[node+1,1]==-1
                return Tree(Array{Tree,1}(undef,0),branch_lengths[daughter_edges[node + 1]], [attributes[p][daughter_edges[node + 1]] for p in 1:n_attributes])
            elseif offsprings[node+1,1]==-2
                return Tree()
            else
                if node != root
                    add_time = branch_lengths[daughter_edges[node + 1]]
                    add_attribute = [attributes[p][daughter_edges[node + 1]] for p in 1:n_attributes]
                else
                    add_time = stem_age
                    add_attribute = root_attributes #repeat([Array{Float64}(undef,0)],n_attributes)
                end

                tree_left = aux(offsprings[node+1,1])
                tree_right = aux(offsprings[node+1,2])

                if tree_left.n_nodes + tree_right.n_nodes == 0
                    return Tree()
                elseif tree_left.n_nodes == 0
                    tree_right = Tree(tree_right.offsprings, tree_right.branch_length + add_time,
                        add_attribute, tree_right.n_nodes)

                    return tree_right
                elseif tree_right.n_nodes == 0
                    tree_left = Tree(tree_left.offsprings, tree_left.branch_length + add_time,
                        add_attribute, tree_left.n_nodes)
                    return tree_left
                else
                    return Tree([tree_left, tree_right], add_time, add_attribute)
                end
            end
        end

        return aux(root)

    else
        function aux_all(node)
            if offsprings[node+1,1]==-1
                return Tree(Array{Tree,1}(undef,0),branch_lengths[daughter_edges[node + 1]], [attributes[p][daughter_edges[node + 1]] for p in 1:n_attributes], true)
            elseif offsprings[node+1,1]==-2
                return Tree(Array{Tree,1}(undef,0),branch_lengths[daughter_edges[node + 1]], [attributes[p][daughter_edges[node + 1]] for p in 1:n_attributes], false)
            else
                if node != root
                    add_time = branch_lengths[daughter_edges[node + 1]]
                    add_attribute = [attributes[p][daughter_edges[node + 1]] for p in 1:n_attributes]
                else
                    add_time = stem_age
                    add_attribute = root_attributes#repeat([Array{Float64}(undef,0)],n_attributes)
                end
                tree_left = aux_all(offsprings[node+1,1])
                tree_right = aux_all(offsprings[node+1,2])
                return Tree([tree_left, tree_right], add_time, add_attribute)
            end
        end

        return aux_all(root)

    end
end

function build_tree(branches::Array{Array{T1,1}}, branch_lengths, attributes::Array{Array{T2,1},1} ;
    extinct = [], prune_extinct=true, return_void = false, stem_age=0., root_attributes = Array{T2,1}(undef,0)) where {T1<:Number, T2<:Number}

    n_attributes = length(attributes)
    if length(extinct) == 0
        dead_branches = fill(false, size(branches)[1])
    else
        dead_branches = extinct
    end

    if size(branches)[1] == 0
        if return_void
            return Tree()
        else
            return Tree(Array{Tree,1}(undef,0),stem_age,root_attributes)
        end
    end
    if size(branches)[1] == 1
        if ! dead_branches[1]
            return Tree(Array{Tree,1}(undef,0),branch_lengths[1],root_attributes, true)
        elseif prune_extinct
            return Tree()
        else
            return Tree(Array{Tree,1}(undef,0),branch_lengths[1],root_attributes, false)
        end
    end

    root = -1
    for b1 in branches
        has_no_parent = true
        for b2 in branches
            if b1[1] == b2[2]
                has_no_parent = false
                break
            end
        end
        if has_no_parent
            root = b1[1]
            break
        end
    end

    max_node = -1
    for b in branches
        max_node = max(b[1], b[2], max_node)
    end
    max_node += 1
    offsprings = fill(-1,max_node,2)
    parent_edges = fill(-1,max_node)
    daughter_edges = fill(-1,max_node)

    for i in 1:length(branches)
        individual = branches[i][1]+1
        parent_edges[individual] = i
        daughter_edges[branches[i][2]+1] = i
        if offsprings[individual,1] < 0
            offsprings[individual,1] = branches[i][2]
        else
            offsprings[individual,2] = branches[i][2]
        end

        if dead_branches[i]
            individual = branches[i][2]+1
            offsprings[individual,:]= [-2 -2]
        end
    end

    if prune_extinct
        function aux(node)
            if offsprings[node+1,1]==-1
                return Tree(Array{Tree,1}(undef,0),branch_lengths[daughter_edges[node + 1]], [attributes[p][daughter_edges[node + 1]] for p in 1:n_attributes])
            elseif offsprings[node+1,1]==-2
                return Tree()
            else
                if node != root
                    add_time = branch_lengths[daughter_edges[node + 1]]
                    add_attribute = [attributes[p][daughter_edges[node + 1]] for p in 1:n_attributes]
                else
                    add_time = stem_age
                    add_attribute = root_attributes #repeat([Array{Float64}(undef,0)],n_attributes)
                end

                tree_left = aux(offsprings[node+1,1])
                tree_right = aux(offsprings[node+1,2])

                if tree_left.n_nodes + tree_right.n_nodes == 0
                    return Tree()
                elseif tree_left.n_nodes == 0
                    tree_right = Tree(tree_right.offsprings, tree_right.branch_length + add_time,
                        add_attribute, tree_right.n_nodes)

                    return tree_right
                elseif tree_right.n_nodes == 0
                    tree_left = Tree(tree_left.offsprings, tree_left.branch_length + add_time,
                        add_attribute, tree_left.n_nodes)
                    return tree_left
                else
                    return Tree([tree_left, tree_right], add_time, add_attribute)
                end
            end
        end

        return aux(root)

    else
        function aux_all(node)
            if offsprings[node+1,1]==-1
                return Tree(Array{Tree,1}(undef,0),branch_lengths[daughter_edges[node + 1]], [attributes[p][daughter_edges[node + 1]] for p in 1:n_attributes], true)
            elseif offsprings[node+1,1]==-2
                return Tree(Array{Tree,1}(undef,0),branch_lengths[daughter_edges[node + 1]], [attributes[p][daughter_edges[node + 1]] for p in 1:n_attributes], false)
            else
                if node != root
                    add_time = branch_lengths[daughter_edges[node + 1]]
                    add_attribute = [attributes[p][daughter_edges[node + 1]] for p in 1:n_attributes]
                else
                    add_time = stem_age
                    add_attribute = root_attributes#repeat([Array{Float64}(undef,0)],n_attributes)
                end
                tree_left = aux_all(offsprings[node+1,1])
                tree_right = aux_all(offsprings[node+1,2])
                return Tree([tree_left, tree_right], add_time, add_attribute)
            end
        end

        return aux_all(root)

    end

end

function build_tree(branches::Array{Int64,2}, branch_lengths, attributes::Array{Array{T2,1},1}, tip_labels::Array{String,1} ;
    extinct = [], prune_extinct=true, return_void = false, stem_age=0., root_attributes = Array{T2,1}(undef,0)) where {T2<:Number}

    n_attributes = length(attributes)
    if length(extinct) == 0
        dead_branches = fill(false, size(branches)[1])
    else
        dead_branches = extinct
    end

    if size(branches)[1] == 0
        if return_void
            return Tree()
        else
            return Tree(Array{Tree,1}(undef,0),stem_age,root_attributes)
        end
    end
    if size(branches)[1] == 1
        if ! dead_branches[1]
            return Tree(Array{Tree,1}(undef,0),branch_lengths[1],root_attributes, true)
        elseif prune_extinct
            return Tree()
        else
            return Tree(Array{Tree,1}(undef,0),branch_lengths[1],root_attributes, false)
        end
    end

    root = -1
    for i in branches[:,1]
        has_no_parent = true
        for j in branches[:,2]
            if i == j
                has_no_parent = false
                break
            end
        end
        if has_no_parent
            root = i
            break
        end
    end

    max_node = maximum(branches[:,1:2])+1
    offsprings = fill(-1,max_node,2)
    parent_edges = fill(-1,max_node)
    daughter_edges = fill(-1,max_node)

    for i in 1:size(branches)[1]
        individual = branches[i,1]+1
        parent_edges[individual] = i
        daughter_edges[branches[i,2]+1] = i
        if offsprings[individual,1] < 0
            offsprings[individual,1] = branches[i,2]
        else
            offsprings[individual,2] = branches[i,2]
        end

        if dead_branches[i]
            individual = branches[i,2]+1
            offsprings[individual,:]= [-2 -2]
        end
    end

    if prune_extinct
        function aux(node)
            if offsprings[node+1,1]==-1
                return Tree(Array{Tree,1}(undef,0),branch_lengths[daughter_edges[node + 1]], [attributes[p][daughter_edges[node + 1]] for p in 1:n_attributes], tip_labels[node])
            elseif offsprings[node+1,1]==-2
                return Tree()
            else
                if node != root
                    add_time = branch_lengths[daughter_edges[node + 1]]
                    add_attribute = [attributes[p][daughter_edges[node + 1]] for p in 1:n_attributes]
                else
                    add_time = stem_age
                    add_attribute = root_attributes #repeat([Array{Float64}(undef,0)],n_attributes)
                end

                tree_left = aux(offsprings[node+1,1])
                tree_right = aux(offsprings[node+1,2])

                if tree_left.n_nodes + tree_right.n_nodes == 0
                    return Tree()
                elseif tree_left.n_nodes == 0
                    tree_right = Tree(tree_right.offsprings, tree_right.branch_length + add_time,
                        add_attribute, tree_right.n_nodes, tree_right.label)

                    return tree_right
                elseif tree_right.n_nodes == 0
                    tree_left = Tree(tree_left.offsprings, tree_left.branch_length + add_time,
                        add_attribute, tree_left.n_nodes, tree_left.label)
                    return tree_left
                else
                    return Tree([tree_left, tree_right], add_time, add_attribute)
                end
            end
        end

        return aux(root)

    else
        function aux_all(node)
            if offsprings[node+1,1]==-1
                return Tree(Array{Tree,1}(undef,0),branch_lengths[daughter_edges[node + 1]], [attributes[p][daughter_edges[node + 1]] for p in 1:n_attributes], true, tip_labels[node])
            elseif offsprings[node+1,1]==-2
                return Tree(Array{Tree,1}(undef,0),branch_lengths[daughter_edges[node + 1]], [attributes[p][daughter_edges[node + 1]] for p in 1:n_attributes], false)
            else
                if node != root
                    add_time = branch_lengths[daughter_edges[node + 1]]
                    add_attribute = [attributes[p][daughter_edges[node + 1]] for p in 1:n_attributes]
                else
                    add_time = stem_age
                    add_attribute = root_attributes#repeat([Array{Float64}(undef,0)],n_attributes)
                end
                tree_left = aux_all(offsprings[node+1,1])
                tree_right = aux_all(offsprings[node+1,2])
                return Tree([tree_left, tree_right], add_time, add_attribute)
            end
        end

        return aux_all(root)

    end
end

function ape2Tree(phylo)
    branches = phylo[:edge]
    branch_lengths = phylo[:edge_length]
    attribute = Array{Float64,1}([1:length(branch_lengths)...])
    tip_labels = phylo[:tip_label]
    return build_tree(branches, branch_lengths, [attribute], tip_labels, root_attributes=[0.])
end

function build_tree(branches, branch_lengths, attribute::Array{T,1} ;
    extinct = [], prune_extinct=true, return_void = false, stem_age=0., root_attributes = Array{T,1}(undef,0)) where {T<:Number}
    build_tree(branches, branch_lengths, [attribute], extinct=extinct, prune_extinct=prune_extinct,
        return_void = return_void, stem_age = stem_age, root_attributes= root_attributes)
end

function build_tree(branches, branch_lengths; extinct = [], prune_extinct=true, return_void = false, stem_age=0., root_attributes = Array{T,1}(undef,0))
    attributes = [1:length(branch_lengths)...]
    build_tree(branches, branch_lengths, [attributes], extinct=extinct,
        prune_extinct=prune_extinct, return_void = return_void, stem_age = stem_age, root_attributes= root_attributes)
end

function plot_ClaDS(tree::Tree ; id = 1, ln=true, lwd=3, show_labels = false, options="")
    plot_tree = Tree(tree.offsprings, 0., tree.attributes, tree.n_nodes)

    reval("""
        library(fields)
        library(ape)

        plot_ClaDS=function(phylo,rate1,rate2=NULL,same.scale=T,main=NULL,lwd=1,log=F, show_labels=F,...){
          Colors = colorRampPalette(c("steelblue2","paleturquoise3","palegreen2","yellow2","salmon1","darkorange", "red","red4"))( 100 )
          if(nrow(phylo[[1]]) <=1){
            plot(1000,xlim = c(0,1), ylim = c(0,2), axes = F, xlab = "", ylab = "")
            lines(c(0,1), c(1,1), lwd=lwd, col="steelblue2")
          }else{
            if(is.null(rate2)){
              if(log) rate1=log(rate1)
              if(isTRUE(all.equal(rep(as.numeric(rate1[1]),length(rate1)),as.numeric(rate1)))){
                col=rep(1,length(rate1))
                plot(phylo, edge.color = Colors[col], show.tip.label = show_labels,main=main,edge.width =lwd,...)
                if(log){
                  image.plot(z = c(exp(rate1[1]),2*exp(rate1[1])),col = Colors, horizontal=T,legend.only = T)
                }else{
                  image.plot(z = c(rate1[1],2*rate1[1]),col = Colors, horizontal=T,legend.only = T)
                }
              }else{
                col = round( (rate1 - min(rate1)) / diff(range(rate1))*99   )+1
                plot(phylo, edge.color = Colors[col], show.tip.label = show_labels,main=main,edge.width =lwd,...)
                if(log){
                  min=min(rate1)
                  max=max(rate1)
                  m10=floor(min/log(10))
                  M10=ceiling(max/log(10))
                  if((M10-m10)<4){
                    ticks=c(1,2,5)
                  }else{
                    ticks=1
                  }
                  ticks=as.vector(sapply(m10:M10,function(k){return(ticks*10^k)}))
                  lt=length(ticks[ticks>exp(min) & ticks<exp(max)])
                  if(lt<4){ticks=c(round(exp(min),max(0,-1*m10+(lt<2))),ticks,round(exp(max),max(0,-1*M10+1+(lt<2))))}
                  image.plot(z = c(min,max),col = Colors, horizontal=T,legend.only = T,axis.args=list( at=log(ticks), labels=ticks))
                }else{
                  image.plot(z = as.matrix(rate1),col = Colors, horizontal=T,legend.only = T)
                }
              }
            }else{
              if(log){
                rate1=log(rate1)
                rate2=log(rate2)
              }
              if(same.scale){
                min=min(min(rate1),min(rate2))
                max=max(max(rate1),max(rate2))
                par(mfrow=c(1,2))
                col = round(( (rate1 - min) / (max-min))*99   )+1
                plot(phylo, edge.color = Colors[col], show.tip.label = show_labels,edge.width =lwd,...)
                col = round(( (rate2 - min) / (max-min))*99   )+1
                plot(phylo, edge.color = Colors[col], show.tip.label = show_labels,edge.width =lwd,...)
                par(mfrow=c(1,1))
                if(log){
                  m10=floor(min/log(10))
                  M10=ceiling(max/log(10))
                  if((M10-m10)<4){
                    ticks=c(1,2,5)
                  }else{
                    ticks=1
                  }
                  ticks=as.vector(sapply(m10:M10,function(k){return(ticks*10^k)}))
                  lt=length(ticks[ticks>exp(min) & ticks<exp(max)])
                  if(lt<4){ticks=c(round(exp(min),max(0,-1*m10+(lt<2))),ticks,round(exp(max),max(0,-1*M10+1+(lt<2))))}
                  # ticks=seq(min,max,length.out = 5)
                  image.plot(z = c(min,max),col = Colors, horizontal=T,legend.only = T,axis.args=list( at=log(ticks), labels=ticks))
                }else{
                  image.plot(z = c(min,max),col = Colors, horizontal=T,legend.only = T)
                }
              }else{
                par(mfrow=c(1,2))
                if(isTRUE(all.equal(rep(rate1[1],length(rate1)),rate1))){
                  col=rep(1,length(rate1))
                  plot(phylo, edge.color = Colors[col], show.tip.label = show_labels,edge.width =lwd,...)
                  if(log){

                    image.plot(z = c(exp(rate1[1]),2*exp(rate1[1])),col = Colors, horizontal=T,legend.only = T)
                  }else{
                    image.plot(z = c(rate1[1],2*rate1[1]),col = Colors, horizontal=T,legend.only = T)
                  }
                }else{
                  col = round(( (rate1 - min(rate1)) / (max(rate1)-min(rate1)))*99   )+1
                  plot(phylo, edge.color = Colors[col], show.tip.label = show_labels,edge.width =lwd,...)
                  if(log){
                    min=min(rate1)
                    max=max(rate1)
                    m10=floor(min/log(10))
                    M10=ceiling(max/log(10))
                    if((M10-m10)<4){
                      ticks=c(1,2,5)
                    }else{
                      ticks=1
                    }
                    ticks=as.vector(sapply(m10:M10,function(k){return(ticks*10^k)}))
                    lt=length(ticks[ticks>exp(min) & ticks<exp(max)])
                    if(lt<4){ticks=c(round(exp(min),max(0,-1*m10+(lt<2))),ticks,round(exp(max),max(0,-1*M10+1+(lt<2))))}
                    image.plot(z = c(min,max),col = Colors, horizontal=T,legend.only = T,axis.args=list( at=log(ticks), labels=ticks))
                  }else{
                    image.plot(z = as.matrix(rate1),col = Colors, horizontal=T,legend.only = T)
                  }
                }
                if(isTRUE(all.equal(rep(rate2[1],length(rate2)),rate2))){
                  col=rep(1,length(rate2))
                  plot(phylo, edge.color = Colors[col], show.tip.label = show_labels,edge.width =lwd,...)
                  if(log){
                    image.plot(z = c(exp(rate2[1]),2*exp(rate2[1])),col = Colors, horizontal=T,legend.only = T)
                  }else{
                    image.plot(z = c(rate2[1],2*rate2[1]),col = Colors, horizontal=T,legend.only = T)
                  }
                }else{
                  col = round(( (rate2 - min(rate2)) / (max(rate2)-min(rate2)))*99   )+1
                  plot(phylo, edge.color = Colors[col], show.tip.label = show_labels,edge.width =lwd,...)
                  if(log){
                    min=min(rate2)
                    max=max(rate2)
                    m10=floor(min/log(10))
                    M10=ceiling(max/log(10))
                    if((M10-m10)<4){
                      ticks=c(1,2,5)
                    }else{
                      ticks=1
                    }
                    ticks=as.vector(sapply(m10:M10,function(k){return(ticks*10^k)}))
                    lt=length(ticks[ticks>exp(min) & ticks<exp(max)])
                    if(lt<4){ticks=c(round(exp(min),max(0,-1*m10+(lt<2))),ticks,round(exp(max),max(0,-1*M10+1+(lt<2))))}
                    image.plot(z = c(min,max),col = Colors, horizontal=T,legend.only = T,axis.args=list( at=log(ticks), labels=ticks))
                  }else{
                    image.plot(z = as.matrix(rate2),col = Colors, horizontal=T,legend.only = T)
                  }
                }
              }
              par(mfrow=c(1,1))
            }
          }

        }

    #source("/Users/maliet/code/plot_ClaDS.R")
    """)
    edges, branch_lengths, rates, tip_labels = make_ape(plot_tree, id = 1)
    ntip = (tree.n_nodes + 1)/2
    @rput edges
    @rput branch_lengths
    @rput rates
    @rput ntip
    @rput ln
    @rput lwd
    @rput show_labels
    @rput tip_labels

    reval("""
        tree = list(edge = edges, Nnode = ntip - 1, edge.lengths = branch_lengths, tip.labels = tip_labels)
        class(tree) = "phylo"
        plot_ClaDS(tree, rates, log=ln, lwd=lwd, show_labels = show_labels $options)
    """)
end

function plot_ClaDS(tree::Tree, rates ; id = 1, ln=true, lwd=3)
    plot_tree = Tree(tree.offsprings, 0., tree.attributes, tree.n_nodes)

    reval("""
    source("/Users/maliet/ownCloud/My_folder/ClaDS_Julia/plot_ClaDS.R")
    """)
    edges, branch_lengths, new_rates = make_ape(plot_tree, id = 1)
    ntip = (tree.n_nodes + 1)/2

    @rput edges
    @rput branch_lengths
    @rput rates
    @rput ntip
    @rput ln
    @rput lwd

    reval("""
        tree = list(edge = edges, Nnode = ntip - 1, edge.lengths = branch_lengths, tip.labels = 1:ntip)
        class(tree) = "phylo"
        plot_ClaDS(tree, rates, log=ln, lwd=lwd)
    """)
end

function plot_ClaDS(tree::Tree, rates1, rates2 ; id = 1, ln=true, lwd=3)
    plot_tree = Tree(tree.offsprings, 0., tree.attributes, tree.n_nodes)

    reval("""
    source("/Users/maliet/ownCloud/My_folder/ClaDS_Julia_list_of_lists/plot_ClaDS_2trees.R")
    """)
    edges, branch_lengths, new_rates = make_ape(plot_tree, id = 1)
    ntip = (tree.n_nodes + 1)/2

    @rput edges
    @rput branch_lengths
    @rput rates1
    @rput rates2
    @rput ntip
    @rput ln
    @rput lwd

    reval("""
        tree = list(edge = edges, Nnode = ntip - 1, edge.lengths = branch_lengths, tip.labels = 1:ntip)
        class(tree) = "phylo"
        leg = plot_ClaDS_noLeg(tree, rates1, rates2, log=ln, lwd=lwd)
        leg = plot_ClaDS_noLeg(tree, rates2, rates1, log=ln, lwd=lwd)
        image.plot(z = leg[[1]],col = leg[[2]], horizontal=F,legend.only = T,axis.args=leg[[5]], legend.mar=4.5)
    """)
end

function sim_ClaDS2_ntips_aux(n,σ,α,ε,λ0 ; return_if_extinct = false, make_tree = true, prune_extinct = true, tree_only = true)

    # accesory functions that will be called latter
    #println(ε)

    bernouilli=Bernoulli(ε/(1+ε))
    function dead()
        return (rand(bernouilli)[1] == 1)
    end

    lambda_law = LogNormal(log(α),σ)
    function new_rates(λ)
        λ * rand(lambda_law, 2)
    end

    while true
        current_node = 3
        alive = [1,2]
        rates = new_rates(λ0)
        branches = [[[0;1]]; [[0;2]]]
        branch_lengths = [0., 0.]
        dead_branches = [false,false]
        node_time = 0
        times = [0.]
        time_diffs = [0.]
        n_lineages = [1]

        time_int = randexp()/(sum(rates[alive])*(1+ε))
        node_time += time_int

        while 0 < length(alive) < n
            individual= splice!(alive, sample(eachindex(alive),Weights(rates[alive])))

            is_dead=dead()
            push!(times, node_time)
            push!(time_diffs, time_int)

            for node in alive
                branch_lengths[node] += time_int
            end
            if is_dead
                branch_lengths[individual] += time_int#node_time - branch_lengths[individual]
                #println(branch_lengths[individual] * rates[individual])
                dead_branches[individual]=true
                push!(n_lineages,length(alive))
            else
                sample_rates = new_rates(rates[individual])
                push!(dead_branches,false,false)
                push!(rates, sample_rates[1], sample_rates[2])
                #push!(branch_lengths, node_time, node_time);
                push!(branch_lengths, 0., 0.);
                push!(alive,current_node,current_node+1)
                push!(n_lineages,length(alive))
                push!(branches,[individual;current_node], [individual;current_node+1])

                current_node = current_node+2
                branch_lengths[individual] += time_int#node_time - branch_lengths[individual]
                #println(branch_lengths[individual] * rates[individual])

            end

            time_int = randexp()/(sum(rates[alive])*(1+ε))
            node_time += time_int
            if time_int==0
                println("zeror ! $(sum(rates[alive]))")
            end
        end

        #println(ε)

        push!(times, node_time)
        push!(time_diffs, time_int)
        push!(n_lineages,length(alive))

        for i in alive
            branch_lengths[i] += time_int#node_time - branch_lengths[i]
        end
        if length(alive) > 0 || return_if_extinct
            if make_tree
                if tree_only
                    return build_tree(branches, branch_lengths, rates, extinct = dead_branches, prune_extinct = prune_extinct, root_attributes=[λ0])
                else
                    t = build_tree(branches, branch_lengths, rates, extinct = dead_branches, prune_extinct = prune_extinct, root_attributes=[λ0])
                    return t, times, n_lineages, time_diffs
                end
            else
                return branches, branch_lengths, rates, dead_branches, times, n_lineages
            end
        end
    end
end


function sim_ClaDS2_ntips_bu(n,σ,α,ε,λ0 ; prune_extinct = true, n_times = 5., tree_only = true, verbose=false)
    r = 0.001 * λ0
    while true
        keep_going = true
        aux = []
        while keep_going
            aux = sim_ClaDS2_ntips_aux(n*n_times,σ,α,ε,λ0 ; return_if_extinct = true, make_tree = true, prune_extinct = false, tree_only = false)
            keep_going = maximum(aux[3])<n
        end

        n_lineages = aux[3]
        times = aux[2]
        tree = aux[1]
        time_diffs = aux[4]
        intervals = [times[i] for i in 1:(length(times)-1) if n_lineages[i] == n]
        interval_widths = [time_diffs[i+1] for i in 1:(length(times)-1) if n_lineages[i] == n]
        tot_time = [0; cumsum(interval_widths)]

        u = rand()
        if verbose
            print(" $r ; $(r*tot_time[end]) , $u     ")
            println(maximum(extract_rates(tree)))
        end
        if u < (r*tot_time[end])
            u = rand() * tot_time[end]
            end_time = 0.
            for i in 1:(length(tot_time)-1)
                if tot_time[i+1] > u > tot_time[i]
                    end_time = u - tot_time[i] + intervals[i]
                    break
                end
            end
            ct = cut_tree(tree, end_time, prune_extinct = prune_extinct)
            if n_extant_tips(ct) <= 2*n
                if tree_only
                    return ct
                else
                    return tree, cut_tree(tree, end_time, prune_extinct = prune_extinct), times, n_lineages, end_time
                end
            end
        else
            r *= 2
        end
    end

end

function sim_ClaDS2_ntips(n,σ,α,ε,λ0 ; return_if_extinct = false,
    make_tree = true, prune_extinct = true, tree_only = true, sed = 0.001,
    max_time = 5, sc=true)

    # accesory functions that will be called latter
    #println(ε)

    bernouilli=Bernoulli(ε/(1+ε))
    function dead()
        return (rand(bernouilli)[1] == 1)
    end

    lambda_law = LogNormal(log(α),σ)
    function new_rates(λ)
        λ * rand(lambda_law, 2)
    end
    N = n*max_time

    s = λ0*sed

    while true
        #println()
        current_node = 3
        alive = [1,2]
        rates = new_rates(λ0)
        branches = [[[0;1]]; [[0;2]]]
        branch_lengths = [0., 0.]
        dead_branches = [false,false]
        node_time = 0
        times = [0.]
        time_diffs = [0.]
        n_lineages = [1]
        scale = 1.
        #println(s)
        while 0 < length(alive) < N
            time_int = randexp()/(sum(rates[alive])*(1+ε))
            if time_int==0
                println("zeror ! $(sum(rates[alive]))")
                break
            end
            if length(alive) == n
                #print(" $scale")
                samp = s * time_int
                node_time += time_int
                u = rand()
                if u < samp*scale
                    time_int *= rand()
                    push!(times, node_time)
                    push!(time_diffs, time_int)
                    push!(n_lineages,length(alive))

                    for i in alive
                        branch_lengths[i] += time_int#node_time - branch_lengths[i]
                    end
                    if make_tree
                        t = build_tree(branches, branch_lengths, rates, extinct = dead_branches, prune_extinct = prune_extinct, root_attributes=[λ0])
                        if tree_only
                            return t
                        else
                            return t, times, n_lineages, time_diffs
                        end
                    else
                        return branches, branch_lengths, rates, dead_branches, times, n_lineages
                    end
                elseif sc

                    scale /= (1. - samp)
                end
            else
                node_time += time_int
            end
            individual= splice!(alive, sample(eachindex(alive),Weights(rates[alive])))

            is_dead=dead()
            push!(times, node_time)
            push!(time_diffs, time_int)

            for node in alive
                branch_lengths[node] += time_int
            end
            if is_dead
                branch_lengths[individual] += time_int#node_time - branch_lengths[individual]
                #println(branch_lengths[individual] * rates[individual])
                dead_branches[individual]=true
                push!(n_lineages,length(alive))
            else
                sample_rates = new_rates(rates[individual])
                push!(dead_branches,false,false)
                push!(rates, sample_rates[1], sample_rates[2])
                #push!(branch_lengths, node_time, node_time);
                push!(branch_lengths, 0., 0.);
                push!(alive,current_node,current_node+1)
                push!(n_lineages,length(alive))
                push!(branches,[individual;current_node], [individual;current_node+1])

                current_node = current_node+2
                branch_lengths[individual] += time_int#node_time - branch_lengths[individual]
                #println(branch_lengths[individual] * rates[individual])

            end


        end

        #println(ε)
        if length(alive) > 0
            s *=2
        end
    end
end

function sim_ClaDS2_time(root_age,σ,α,ε,λ0 ; return_if_extinct = true, max_node_number = Inf, make_tree = true, prune_extinct = false, return_if_max = true)

    # accesory functions that will be called latter
    bernouilli=Bernoulli(ε/(1+ε))
    function dead()
        return (rand(bernouilli)[1] == 1)
    end

    lambda_law = LogNormal(log(α),σ)
    function new_rates(λ::Float64)
        λ * rand(lambda_law, 2)
    end

    function aux_ext(time, rate, n_max)
        #println(n_max)
        if n_max == 0
            return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false), 0
        end
        node_time = randexp()/(rate*(1+ε))
        if node_time > time
                #offsprings, branch_length, attributes, n_nodes
                return Tree(Array{Tree,1}(undef,0), time, [rate], 1, true), n_max
        else
                is_dead = dead()
                if is_dead
                    return Tree(Array{Tree,1}(undef,0), node_time, [rate], 1, false), n_max + 1
                else
                    offspring_rates = new_rates(rate)
                    left_time = time - node_time
                    #print("$depth ;")
                    left_tree = aux_ext(left_time, offspring_rates[1], n_max - 1)
                    #print(left_tree[2])
                    if left_tree[1].n_nodes == -1
                        return left_tree
                    end
                    right_tree = aux_ext(left_time, offspring_rates[2], left_tree[2])
                    if right_tree[1].n_nodes == -1
                        return right_tree
                    end
                    return Tree([left_tree[1], right_tree[1]], node_time, [rate]), right_tree[2]
                end
        end
    end
    tree = aux_ext(root_age, λ0, max_node_number)[1]


    if prune_extinct
        if (tree.n_nodes == -1) & return_if_max
            return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false)
        else
            return prune_extinct_lineages(tree)
        end
    else
        if (tree.n_nodes == -1) & return_if_max
            return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false)
        else
            return tree
        end
    end
end

function sim_ClaDS2_time(root_age,σ,α,ε,λ0,u,sf,lf ; return_if_extinct = true, max_node_number = Inf, make_tree = true, prune_extinct = false, return_if_max = true)
    # accesory functions that will be called latter
    bernouilli=Bernoulli(ε/(1+ε))
    function dead()
        return (rand(bernouilli)[1] == 1)
    end

    lambda_law = LogNormal(log(α),σ)
    function new_rates(λ::Float64)
        X = rand(lambda_law, 2)
        λ * X
    end
    #println()
    function aux_ext(time, rate, n_max, s)
        #println(n_max)
        if n_max == 0
            #println("max")
            return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false), 0, -In
        end
        node_time = randexp()/(rate*(1+ε))
        if node_time <= 0
            return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false), 0
        end
        if node_time > time
            s += lf
            #print(" $s $u ;")
            if u<s
                return Tree(Array{Tree,1}(undef,0), time, [rate], 1, true), n_max + 1, s
            else
                #print(" !")
                return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false), 0, s
            end
        else
                is_dead = dead()
                if is_dead
                    #print("dead")
                    return Tree(Array{Tree,1}(undef,0), node_time, [rate], 1, false), n_max + 1, s
                else
                    offspring_rates = new_rates(rate)
                    left_time = time - node_time
                    if time == left_time
                        return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false), 0
                    end
                    #print("$depth ;")
                    left_tree = aux_ext(left_time, offspring_rates[1], n_max - 1, s)
                    #print(left_tree[2])
                    if left_tree[1].n_nodes == -1
                        return left_tree
                    end
                    right_tree = aux_ext(left_time, offspring_rates[2], left_tree[2], left_tree[3])
                    if right_tree[1].n_nodes == -1
                        return right_tree
                    end
                    return Tree([left_tree[1], right_tree[1]], node_time, [rate]), right_tree[2], right_tree[3]
                end
        end
    end
    tree = aux_ext(root_age, λ0, max_node_number, sf)


    if prune_extinct
        if (tree[1].n_nodes == -1) & return_if_max
            return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false), -Inf
        else
            return prune_extinct_lineages(tree[1]), tree[3]
        end
    else
        if (tree[1].n_nodes == -1) & return_if_max
            return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false), -Inf
        else
            return tree[1], tree[3]
        end
    end
end

function sim_ClaDS2_time_rates(root_age,σ,α,ε,λ0 ; return_if_extinct = true, max_node_number = 1_000,
    make_tree = true, prune_extinct = false, return_if_max = true, return_na = false)

    # accesory functions that will be called latter
    bernouilli=Bernoulli(ε/(1+ε))
    function dead()
        return (rand(bernouilli)[1] == 1)
    end

    lambda_law = LogNormal(log(α),σ)
    function new_rates(λ::Float64)
        λ * rand(lambda_law, 2)
    end
    if return_na
        function aux_na(time, rate, n_max, rates, live)
            #println(n_max)
            if n_max == 0
                return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false), 0
            end
            node_time = randexp()/(rate*(1+ε))
            if node_time <= 0
                return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false), 0
            end
            if node_time > time
                #offsprings, branch_length, attributes, n_nodes
                push!(rates, rate)
                push!(live, true)
                return Tree(Array{Tree,1}(undef,0), time, [rate], 1, true), n_max
            else
                is_dead = dead()
                if is_dead
                    #print("dead")
                    push!(live, false)
                    return Tree(Array{Tree,1}(undef,0), node_time, [rate], 1, false), n_max + 1
                else
                    offspring_rates = new_rates(rate)
                    left_time = time - node_time
                    if time == left_time
                        return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false), 0
                    end
                    #print("$depth ;")
                    left_tree = aux_na(left_time, offspring_rates[1], n_max - 1, rates, live)
                    #print(left_tree[2])
                    if left_tree[1].n_nodes == -1
                        return left_tree
                    end
                    right_tree = aux_na(left_time, offspring_rates[2], left_tree[2], rates, live)
                    if right_tree[1].n_nodes == -1
                        return right_tree
                    end
                    return Tree([left_tree[1], right_tree[1]], node_time, [rate]), right_tree[2]
                end
            end
        end

        rates = Array{Float64,1}(undef,0)
        alive = Array{Bool,1}(undef,0)

        tree = aux_na(root_age, λ0, max_node_number, rates, alive)[1]
        if (tree.n_nodes == -1) & return_if_max
            return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false), Array{Float64,1}(undef,0), Array{Bool,1}(undef,0)
        else
            return tree, rates, alive
        end
    else
        function aux_ext(time, rate, n_max, rates)
            #println(n_max)
            if n_max == 0
                return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false), 0
            end
            node_time = randexp()/(rate*(1+ε))
            if node_time > time
                #offsprings, branch_length, attributes, n_nodes
                push!(rates, rate)
                return Tree(Array{Tree,1}(undef,0), time, [rate], 1, true), n_max
            else
                is_dead = dead()
                if is_dead
                    #print("dead")
                    return Tree(Array{Tree,1}(undef,0), node_time, [rate], 1, false), n_max + 1
                else
                    offspring_rates = new_rates(rate)
                    left_time = time - node_time
                    #print("$depth ;")
                    left_tree = aux_ext(left_time, offspring_rates[1], n_max - 1, rates)
                    #print(left_tree[2])
                    if left_tree[1].n_nodes == -1
                        return left_tree
                    end
                    right_tree = aux_ext(left_time, offspring_rates[2], left_tree[2], rates)
                    if right_tree[1].n_nodes == -1
                        return right_tree
                    end
                    return Tree([left_tree[1], right_tree[1]], node_time, [rate]), right_tree[2]
                end
            end
        end

        rates = Array{Float64,1}(undef,0)

        tree = aux_ext(root_age, λ0, max_node_number, rates)[1]
        if (tree.n_nodes == -1) & return_if_max
            return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false), Array{Float64,1}(undef,0)
        else
            return tree, rates
        end
    end

end

function sim_ClaDS2_time_rates_tip(root_age,σ,α,ε,λ0,u,sf,fl ; return_if_extinct = true, max_node_number = 1_000,
    make_tree = true, prune_extinct = false, return_if_max = true, return_na = false)

    # accesory functions that will be called latter
    #print(" $(ε/(1+ε))")
    bernouilli=Bernoulli(ε/(1+ε))
    function dead()
        return (rand(bernouilli)[1] == 1)
    end

    lambda_law = LogNormal(log(α),σ)
    function new_rates(λ::Float64)
        X = (rand(lambda_law, 2))
        λ * X
    end
    if return_na
    else
        function aux_ext(time, rate, n_max, rates, s ,n)
            if n_max == 0
                return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false), 0,0,0
            end
            node_time = randexp()/(rate*(1+ε))
            #println(node_time)
            if node_time <= 0
                return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false), 0,0,0
            end#print(" $s $n $u;")# $time $node_time $rate;")
            if node_time > time
                #print("zzzz")
                new_s = s
                #offsprings, branch_length, attributes, n_nodes
                push!(rates, rate)
                if n>0
                    new_s = s + log(1-sf)+log(n+1)-log(n)
                else
                    new_s = s + log(sf)
                end
                if new_s >= s
                    #print("a")
                    return Tree(Array{Tree,1}(undef,0), time, [rate], 1, true), n_max, new_s , n+1
                elseif u < new_s
                    #print("b")
                    return Tree(Array{Tree,1}(undef,0), time, [rate], 1, true), n_max, new_s , n+1
                else
                    #print("c")
                    return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false), 0,new_s ,n+1
                end
            else
                is_dead = dead()
                if is_dead
                    return Tree(Array{Tree,1}(undef,0), node_time, [rate], 1, false), n_max + 1, s, n
                else
                    offspring_rates = new_rates(rate)
                    left_time = time - node_time
                    if time == left_time
                        return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false), 0,0,0
                    end
                    #print("$depth ;")
                    left_tree = aux_ext(left_time, offspring_rates[1], n_max - 1, rates, s, n)
                    #print(left_tree[2])
                    if left_tree[1].n_nodes == -1
                        return left_tree
                    end
                    right_tree = aux_ext(left_time, offspring_rates[2], left_tree[2], rates, left_tree[3], left_tree[4])
                    if right_tree[1].n_nodes == -1
                        return right_tree
                    end
                    return Tree([left_tree[1], right_tree[1]], node_time, [rate]), right_tree[2], right_tree[3], right_tree[4]
                end
            end
        end

        rates = Array{Float64,1}(undef,0)

        tree = aux_ext(root_age, λ0, max_node_number, rates,fl,0)
        if ((tree[1].n_nodes == -1) & return_if_max) | (tree[4] == 0) | (u > tree[3])
            return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false), Array{Float64,1}(undef,0)
        else
            return tree[1], rates
        end
    end

end

function sim_ClaDS2_time_unsampled(root_age,σ,α,ε,λ0 ;
    sampling_proba = 1., return_if_sampled = false, return_if_extinct = true, max_node_number = 1_000, not_sampled = true,
    make_tree = false, prune_extinct = false)

    while true
        tree = sim_ClaDS2_time(root_age,σ,α,ε,λ0 ,
            return_if_extinct = true, max_node_number = max_node_number,
            prune_extinct = false)
        if tree.n_nodes == -1
            if return_if_sampled
                return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false), false, Inf
            end
        else
            n_alive = n_extant_tips(tree)
            if not_sampled
                if n_alive == 0
                    log_sampling_proba = 0
                else
                    log_sampling_proba = n_alive * log(1-sampling_proba)
                end
            else
                if n_alive == 0
                    log_sampling_proba = -Inf
                elseif n_alive == 1
                    log_sampling_proba = log(sampling_proba)
                else
                    log_sampling_proba = log(n_alive) + log(sampling_proba) + (n_alive - 1) * log(1-sampling_proba)
                end
            end

            u = log(rand())
            unsampled = u < log_sampling_proba
            if unsampled | return_if_sampled
                return tree, unsampled, n_alive
            end
        end
    end
end

function sim_ClaDS2_time_unsampled_rates(root_age,σ,α,ε,λ0 ;
    sampling_proba = 1., return_if_sampled = false, return_if_extinct = true, max_node_number = 1_000, not_sampled = true,
    make_tree = false, prune_extinct = false)

    while true
        tree, rates = sim_ClaDS2_time_rates(root_age,σ,α,ε,λ0 ,
            return_if_extinct = true, max_node_number = max_node_number,
            prune_extinct = false)
        if tree.n_nodes == -1
            if return_if_sampled
                return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false), false, Inf, rates
            end
        else
            n_alive = length(rates)
            if not_sampled
                if n_alive == 0
                    log_sampling_proba = 0
                else
                    log_sampling_proba = n_alive * log(1-sampling_proba)
                end
            else
                if n_alive == 0
                    log_sampling_proba = -Inf
                elseif n_alive == 1
                    log_sampling_proba = log(sampling_proba)
                else
                    log_sampling_proba = log(n_alive) + log(sampling_proba) + (n_alive - 1) * log(1-sampling_proba)
                end
            end

            u = log(rand())
            unsampled = u < log_sampling_proba
            if unsampled | return_if_sampled
                return tree, unsampled, n_alive, rates
            end
        end
    end
end

function sim_ClaDS2_time_unsampled_rates(root_age,σ,α,ε,λ0, u,fl ;
    sampling_proba = 1., return_if_sampled = false, return_if_extinct = true, max_node_number = 1_000, not_sampled = true,
    make_tree = false, prune_extinct = false)

    while true
        tree, rates = sim_ClaDS2_time_rates_tip(root_age,σ,α,ε,λ0,u,sampling_proba,fl,
            return_if_extinct = true, max_node_number = max_node_number,
            prune_extinct = false)
        if tree.n_nodes == -1
            if return_if_sampled
                return Tree(Array{Tree,1}(undef,0), -1., [-1], -1, false), false, Inf, rates
            end
        else
            n_alive = length(rates)
            return tree, true, n_alive, rates
        end
    end
end
