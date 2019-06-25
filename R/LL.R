LL <-
function(mu,symbiont_tree,sequences, n, nodes, edges, el, eig_val, eig_vect, ivp, propinv, N){-sum(vapply(1:N, function(k) L_computation(k,mu,symbiont_tree,sequences,n, nodes, edges, el, eig_val, eig_vect, ivp, propinv),numeric(1)))}
