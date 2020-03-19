function prune_extinct_lineages(tree::Tree)
    function combine(left, right, bl, att)
        if ! left.extant
            if ! right.extant
                return Tree()
            else
                return Tree(right.offsprings, right.branch_length +bl,
                    att, true, right.label)
            end
        elseif ! right.extant
            return Tree(left.offsprings, left.branch_length + bl,
                att, true, left.label)
        else
            return Tree([left, right], bl, att, true, "")
        end
    end

    function aux(subtree)
        if length(subtree.offsprings) == 0
            return subtree
        else
            combine(aux(subtree.offsprings[1]), aux(subtree.offsprings[2]), subtree.branch_length, subtree.attributes)
        end
    end
    aux(tree)
end

function sample_tips_tipRates(tree::Tree, f::Float64 ; root_age = true)
    function combine(left, right, bl, att)
        extant = (left.extant || right.extant)
        return Tree([left,right], bl, att, extant, "")
    end
    function aux(subtree)
        u = 0.5
        if subtree.extant
            if length(subtree.offsprings) == 0
                if f < 1.
                    u = rand()
                end
                if u < f
                    return subtree
                else
                    return Tree()
                end
            else
                combine(aux(subtree.offsprings[1]), aux(subtree.offsprings[2]), subtree.branch_length, subtree.attributes)
            end
        else
            return subtree
        end
    end

    if root_age
        left_tree = Tree()
        while (!left_tree.extant)
            left_tree = aux(tree.offsprings[1])
        end
        right_tree = Tree()
        while (!right_tree.extant)
            right_tree = aux(tree.offsprings[2])
        end
        new_tree = Tree([left_tree,right_tree], tree.branch_length, tree.attributes, tree.label)
        tip_rates = extract_tip_rates_light(new_tree, return_extinct = false)
        return (prune_extinct_lineages(new_tree), tip_rates)
    else
        new_tree = aux(tree)
        tip_rates = extract_tip_rates_light(new_tree, return_extinct = false)
        return (prune_extinct_lineages(new_tree), tip_rates)
    end
end
