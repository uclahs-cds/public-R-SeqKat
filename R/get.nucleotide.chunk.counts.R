library(Rcpp);

#' Get Nucleotide Chunk Counts
#'
#' Short description
#'
#' Detailed description
#' 
#' @param key List of specify trinucleotides to count
#' @param chr Chromosome
#' @param upstream Length upstream to read
#' @param downstream Length downstream to read
#' @param start Starting position
#' @param end Ending position
#'
#' @examples
#' \dontrun{
#' 
#' }
#'
#' @author \email{Fouad.Yousif@oicr.on.ca}

get.nucleotide.chunk.counts  <-  function (key, chr, upstream = 1, downstream = 1, start = 1, end = -1){
	cget_nucleotide_chunk_counts(key, chr, upstream, downstream, start, end);
	};
