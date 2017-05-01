library(Rcpp);

#' Get Context
#'
#' Get the 3 context positions around the mutation
#'
#' Detailed description
#' 
#' @param file Reference files directory
#' @param start The position of the mutation gene

#' @return The trinucleotide context.
#'
#' @examples
#' \dontrun{
#' get.context(file.path(referencie.genome.dir, 'chr1.fa'), c(158297133, 161176181))
#' }
#'
#' @author \email{Fouad.Yousif@oicr.on.ca}


get.context  <-  function (file,start){
	cpp_get_context(file, start, length(start));
	};
