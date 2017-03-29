#code from Nick Matzke's code on BioGeoBEARS wiki
.events_txt_list_into_events_table<-function(events_txt_list, trtable=NULL, recalc_abs_ages=TRUE)
	{
	
	if (is.null(events_txt_list))
		{
		errortxt = paste("\nWARNING in events_txt_list_into_events_table(): your events_txt_list has NO events!\n\nThis means your tree has NO d/e/a events across the whole tree.\nThis is *expected* e.g. if you inferred d=e=0 under DEC+J. Input a list of '' or NA to avoid this error.\n\n", sep="")
		cat(errortxt)
		errortxt2 = paste("events_txt_list_into_events_table() is returning NULL which will might cause issues later.\n\n", sep="")
		cat(errortxt2)
		return(NULL)
		}
	
	# Convert NAs to "none"
	events_txt_list[is.na(events_txt_list)] = "none"
	
	# Remove lines with no events or NA:
	noneTF = events_txt_list == "none"
	keepTF = (noneTF == FALSE)
	events_txt_list = events_txt_list[keepTF]
	
	
	# If no anagenetic events, return NULL
	if (length(events_txt_list) == 0)
		{
		events_table = NULL
		return(events_table)
		}


	# Include the trtable, if that is input
	if (length(trtable) > 0)
		{
		trtable_subset = NULL
		}

	
	# Convert the events text back into a table:
	tmptable = NULL
	for (i in 1:length(events_txt_list))
		{
		#print(events_txt_list)
		tmptable_rows = .events_txt_into_events_table(events_txt_list[i])
		rownums_in_trtable = as.numeric(tmptable_rows$nodenum_at_top_of_branch)
		#print(tmptable_rows)
		num_newrows = nrow(tmptable_rows)
		tmptable = rbind(tmptable, tmptable_rows)
		} # END for (i in 1:length(events_txt_list))
	events_table = BioGeoBEARS::dfnums_to_numeric(BioGeoBEARS::adf2(tmptable))
	names(events_table) = c("nodenum_at_top_of_branch", "trynum", "brlen", "current_rangenum_1based", "new_rangenum_1based", "current_rangetxt", "new_rangetxt", "abs_event_time", "event_time", "event_type", "event_txt", "new_area_num_1based", "lost_area_num_1based", "dispersal_to", "extirpation_from")	
	
	return(events_table)
	}