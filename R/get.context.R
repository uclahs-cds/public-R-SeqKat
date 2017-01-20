### GET.CONTEXT ####################################################################################
# Description: get the 3 context positions around the mutation one
# Input variables:
# - file: reference files directory
# - start: the position of the mutation gene
# Output variables: the 3 base ATCG

### PREAMBLE #######################################################################################
# loading libraries
library(Rcpp);

# tranform the 24 reference fastas in to the form that could be used, and change the last two names
# ------------------------prepare-------------------------------------------------------------------#
#chrs <- paste0('chr', c(as.character(1:22), 'X', 'Y'));
#setwd('/oicr/data/genomes/homo_sapiens_mc/UCSC/hg19_random/Genomic/references/fasta');
#sapply(
#	X   = chrs, 
#	FUN = function(chr) {
#		line.nums <- as.numeric(system(sprintf("wc -l %s.fa | cut -d' ' -f1", chr), intern = TRUE));
#		system(
#			sprintf(
#				"tail -n %d %s.fa | tr -d '\\n' > ~/isilon/scratch/ffan/hg19_ref/%s.fa", 
#				line.nums - 1, chr, chr
#				), 
#			intern = TRUE
#			);
#		}
#	);
#setwd("~/isilon/scratch/ffan/hg19_ref");
#file.copy("chrX.fa", "chr23.fa");
#file.copy("chrY.fa", "chr24.fa");
#
# --------------------------------------------------------------------------------------------------#


### FUNCTION: get.context(file, start) ##############################################################

get.context  <-  function (file,start){
	cpp_get_context(file, start, length(start));
	};

# test if the base is given correctly 
# get.context(file.path(referencie.genome.dir, 'chr1.fa'), c(158297133, 161176181))

