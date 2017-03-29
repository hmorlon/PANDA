#code from Nick Matzke's code on BioGeoBEARS wiki
.event_txt_into_events_row<-function(word)
	{

	split_key_item<-function(word2)
		{
		output_pair = c("", "")
		words3 = strsplit(word2, split=":")[[1]]
		numwords = length(words3)
		
		output_pair[1:numwords] = words3[1:numwords]
		
		return(output_pair)
		}


	words2 = strsplit(word, split=",")[[1]]
	output = sapply(X=words2, FUN=split_key_item)
	tmprow = matrix(data=output[2,], nrow=1)
	return(tmprow)
	}
