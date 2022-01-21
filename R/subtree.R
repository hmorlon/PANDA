#take a list of species and return the monophyletic group containing these species (i.e the species and all the descendants of the MRCA), along with the root length, contained in subtree$root.edge

#source("./tools/extract.clade.ln.R")

subtree<-function(tree,species_list)

{
# find the MRCA of the species in the list
node<-unique(tree$edge[which(tree$edge[,2] %in% which(tree$tip.label %in% species_list)),1])
subtree<-extract.clade.ln(tree,min(node))

while(sum(!(species_list %in% subtree$tip.label))>0)
{
	node<-unique(tree$edge[which(tree$edge[,2] %in% node),1])
	subtree<-extract.clade.ln(tree,min(node))
	#print(node)
	}

root_length<-tree$edge.length[which(tree$edge[,2] == min(node))]
subtree$root.edge<-root_length

return(subtree)}
