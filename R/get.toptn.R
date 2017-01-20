### GET.TOPTN  #######################################################################
# somatic.subset is the .bed file for a single chr
get.toptn <- function(somatic.subset, chr, start.bp, end.bp, ref.dir){
	dict <- c('A' = 'T', 'C' = 'G', 'G' = 'C', 'T' = 'A', 'N' = 'N');
	
	index <- which(somatic.subset$POS %in% c(start.bp, end.bp));
	index <- index[1]:index[2]
	mut.index <- somatic.subset$POS[index];
	mut.tn <- toupper(get.context(
			file.path(ref.dir,paste(paste('chr',chr,sep=''),'fa',sep='.')),
			mut.index
			));
	mut.snv <- toupper(somatic.subset$ALT[index]);
	mut.tn.adj <- ifelse(substr(mut.tn,2,2) %in% c('A','G'), get.pair(mut.tn),mut.tn);
	mut.snv.adj <- ifelse(substr(mut.tn,2,2) %in% c('A','G'),dict[mut.snv],mut.snv);
	changes <- paste(mut.tn.adj, mut.snv.adj, sep ='>');
	tmp <- table(changes);
	return(paste(paste(dimnames(tmp)[[1]],as.vector(tmp),sep = ':'), collapse = ', '));
	}