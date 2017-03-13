library(Rcpp);
library(foreach);

seqkat <- function(
	sigcutoff = 5,
	mutdistance = 3.2,
	segnum = 4,
	ref.dir  = NULL,
	bed.file = "./",
	output.dir = "./",
	chromosome = "all",
	chromosome.length.file = NULL,
	trinucleotide.count.file = NULL
	) {

	# Validate sigcutoff, mutdistance, segnum
	stopifnot(is.numeric(sigcutoff), is.numeric(mutdistance), is.numeric(segnum));

	# Validate the bed.file
	if (!file.exists(bed.file)) {
		stop("The bed.file provided does not exist.");
		}

	# Validate the chromosome
	if ((chromosome != "all") & (!grepl('^(\\d*)$', chromosome))) {
		stop("The provided value for chromosome should be of the form (1, 2, ..., 23, 24), without a prepended 'chr'.");
		}
	else if (grepl('^(\\d*)$', chromosome)) {
		stopifnot(as.numeric(chromosome) >= 1, as.numeric(chromosome) <= 24, as.numeric(chromosome)%%1 == 0);
		}

	# Validate the ref.dir
	if (is.null(ref.dir)) {
		stop("Please supply a path to the reference genome with the ref.dir argument.")
		}
	if (!all(file.exists(sapply(1:24, function(x) { paste0(ref.dir, "/chr", x, ".fa") })))) {
		stop("The ref.dir directory must contain separate fasta files for each chromosome of the form (chr1.fa, chr2.fa, ..., chr23.fa, chr24.fa).")
		}

	# Validate the chromosome.length.file
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

	# Validate the trinucleotide.count.file
	if (is.null(trinucleotide.count.file)) {
		warning("No trinucleotide.count.file provided, using hg19 counts by default. This file can be generated using the get.trinucleotide.counts() function in this package.");
		trinucleotide.count.file <- paste0(path.package("SeqKat"),"/inst/extdata/tn_count.txt")
		}

	if (!file.exists(output.dir)) {
		dir.create(output.dir);
		}
	setwd(output.dir);

	somatic.directory <- paste(strsplit(basename(bed.file),'_')[[1]][strsplit(basename(bed.file),'_')[[1]]!='snvs.bed'],collapse='_');
	dir.create(somatic.directory);
	somatic.file <- bed.file;
	
	print(somatic.file);
	somatic<-read.delim(
		somatic.file,
		header = TRUE,
		stringsAsFactors = FALSE
		);

	names(somatic) <- c("CHR","POS","REF","ALT");

	#Rename chrX and chrY
	somatic$CHR <- gsub("chrX","chr23",somatic$CHR);
	somatic$CHR <- gsub("chrY","chr24",somatic$CHR);

	output.name <- somatic.directory;
	exprobntcx <- get.exprobntcx(somatic, ref.dir, trinucleotide.count.file);
	if (chromosome == "all") {
		chrs <- sort(na.omit(as.numeric(unlist(strsplit(unique(somatic$CHR),'[^0-9]+')))));
		}
	else {
		chrs <- c(as.numeric(chromosome));
		}
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
