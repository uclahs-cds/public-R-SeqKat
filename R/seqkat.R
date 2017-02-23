library(Rcpp);
library(foreach);

seqkat <- function(
	sigcutoff = 5,
	mutdistance = 3.2,
	segnum = 4,
	ref.dir = NULL,
	bed.dir = "./",
	chromosome.length.file = NULL
	) {
	if (is.null(ref.dir)) {
		stop("Please supply a path to the reference genome with the ref.dir argument.")
		}

	if (is.null(chromosome.length.file)) {
		warning("No chromosome.length.file provided, using hg19 lengths by default.");
		chromosome.length.file <- paste0(path.package("SeqKat"),"/inst/extdata/length_hg19_chr.txt")
		}
	else {
		chr.length <- read.table(file=chromosome.length.file, header = TRUE, stringsAsFactors = FALSE);
		if (is.null(chr.length$num)) {
			stop("The supplied chromosome.length.file is missing the 'num' column.");
			}
		if (is.null(chr.length$length)) {
			stop("The supplied chromosome.length.file is missing the 'length' column.");
			}
		if (!all.equal(c(as.character(1:24), "sum.f", "sum.m"), chr.length$num)) {
			stop("The supplied chromosome.length.file is missing required rows.");
			}
		}

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
		somatic<-read.delim(
			somatic.file,
			header = TRUE,
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
							 ref.dir = ref.dir,
							 chromosome.length.file = chromosome.length.file
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
