##This function creates a geography object to pass to trait model scripts from the output of a BioGeoBEARS stochastic mapping exercise
##
##INPUTS
#full.phylo is the phylogeny that was used to build the stochastic map in biogeobears
#trimmed.phylo is the phylogeny of the subset of species for which the geography object will be built (if NULL, this is the full phylo)
#ana.events is one of the "ana.events.tables" that results from a stochastic mapping exercise in BioGeoBEARS
#clado.events is one of the "clado.events.tables" that results from a stochastic mapping exercise in BioGeoBEARS
#stratified is a logical indicating whether a stratified analysis was used (default is FALSE)


CreateGeoObject_BioGeoBEARS<-function(full.phylo,trimmed.phylo=NULL,ana.events,clado.events,stratified=FALSE){

	if(stratified){
		if(is.null(trimmed.phylo)){
			clado_events_tables<-list()
			clado_events_tables[[1]]<-clado.events
			smap<-stratified_BGB_to_tables(full.phylo,clado_events_tables,1)
			x<-CreateBioGeoB_Object_subclade(anc.phylo=full.phylo,subclade.phylo=full.phylo,ana.events=smap$ana.int,clado.events=smap$clado.int,nat.only=FALSE)
			return(x)
			}else{
			clado_events_tables<-list()
			clado_events_tables[[1]]<-clado.events
			smap<-stratified_BGB_to_tables(full.phylo,clado_events_tables,1)
			x<-CreateBioGeoB_Object_subclade(anc.phylo=full.phylo,subclade.phylo=trimmed.phylo,ana.events=smap$ana.int,clado.events=smap$clado.int,nat.only=FALSE)
			return(x)		
			}		
	} else{
		if(is.null(trimmed.phylo)){
			x<-CreateBioGeoB_Object_subclade(anc.phylo=full.phylo,subclade.phylo=full.phylo,ana.events=ana.events,clado.events=clado.events,nat.only=FALSE)
			return(x)
		} else{
			x<-CreateBioGeoB_Object_subclade(anc.phylo=full.phylo,subclade.phylo=trimmed.phylo,ana.events=ana.events,clado.events=clado.events,nat.only=FALSE)
			return(x)
		}
	}
}