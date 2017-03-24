#code from Nick Matzke's code on BioGeoBEARS wiki
.events_txt_into_events_table<-function(branch_events_txt)
	{
	words = strsplit(branch_events_txt, split=";")[[1]]
	
	events_table_for_branch = t(sapply(X=words, FUN=.event_txt_into_events_row))
	row.names(events_table_for_branch) = NULL
	events_table_for_branch
	
	events_table_for_branch = BioGeoBEARS::adf2(events_table_for_branch)
	events_table_for_branch
	names(events_table_for_branch) = c("nodenum_at_top_of_branch", "trynum", "brlen", "current_rangenum_1based", "new_rangenum_1based", "current_rangetxt", "new_rangetxt", "abs_event_time", "event_time", "event_type", "event_txt", "new_area_num_1based", "lost_area_num_1based", "dispersal_to", "extirpation_from")
	
	return(events_table_for_branch)
	}
