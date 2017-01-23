library(Rcpp);
library(foreach);

seqkat <- function(sigcutoff, mutdistance, segnum, ref.dir, bed.dir = "./") {
	setwd(bed.dir);

	# save all samples in a list 
	somatic.list <- list.files(
		path = bed.dir,
		pattern = '_snvs.bed'
		);

	output.list <- sapply(
		somatic.list,
		function(x) {
			somatic.directory <- paste(strsplit(x,'_')[[1]][strsplit(x,'_')[[1]]!='snvs.bed'],collapse='_');
			dir.create(somatic.directory);
			somatic.directory;
			}
		);

	for ( i in 1: length(somatic.list)){
		somatic.file <- somatic.list[i];
		print(somatic.file);
		#print(getwd());
		somatic<-read.delim(
			somatic.file,
			#"BTCA-SG_1e2dcbcc-771c-43c5-8c8d-e0eb77cb3494_snvs.bed",
			header = TRUE,  ####### 
			stringsAsFactors = FALSE
			);
		#print("\ntest");
		names(somatic) <- c("CHR","POS","REF","ALT");
		#Rename chrX and chrY
		somatic$CHR <- gsub("chrX","chr23",somatic$CHR);
		somatic$CHR <- gsub("chrY","chr24",somatic$CHR);
		# output.name
		output.name <- output.list[i];
		exprobntcx <- get.exprobntcx(somatic, ref.dir);
		chrs <- sort(na.omit(as.numeric(unlist(strsplit(unique(somatic$CHR),'[^0-9]+')))));
		sapply(
			chrs,
			function(x){
				combine.table(
					 final.score(
						 test.table =  test.kataegis(
							 chromosome.num = x, 
							 somatic, 
							 unit = 2,
							 exprobntcx = exprobntcx,
							 output.name = output.name,
							 ref.dir = ref.dir
							 ), 
						 cutoff = sigcutoff, 
						 somatic,
						 output.name = output.name
						 ),
					somatic, 
					mutdistance, 
					segnum,
					output.name = output.name
					);
				}
			);
		}
	}
