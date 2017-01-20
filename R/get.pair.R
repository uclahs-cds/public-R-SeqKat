get.pair <- function(x) {
	dict <- c('A' = 'T', 'C' = 'G', 'G' = 'C', 'T' = 'A', 'N' = 'N');
	neu <- strsplit(x, '');
	sapply(neu, function(z) paste(rev(dict[z]), collapse=''))
	};