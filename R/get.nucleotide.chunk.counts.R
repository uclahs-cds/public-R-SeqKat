library(Rcpp);

get.nucleotide.chunk.counts  <-  function (key, chr, upstream = 1, downstream = 1, start = 1, end = -1){
	cget_nucleotide_chunk_counts(key, chr, upstream, downstream, start, end);
	};
