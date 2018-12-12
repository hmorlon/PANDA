LL <-
function(mu,symbiont_tree,nodes,Nd,n,eigQ,ivp,sequences,PI){-sum(sapply(1:Nd, function(k) L_computation(k,mu,symbiont_tree,nodes,n,eigQ,ivp,sequences,PI)))}
