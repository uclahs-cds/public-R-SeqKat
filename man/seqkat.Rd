\name{seqkat}
\alias{seqkat}
\title{SeqKat: Kataegis detection from SNV BED files}
\usage{
seqkat(sigcutoff, mutdistance, segnum, ref.dir, bed.dir = "./")
}
\arguments{
\item{sigcutoff}{the minimum hypermutation score used to classify the windows in the sliding binomial test as significant windows. The score is calculated per window as follows: -log10(binomial test p-value). Recommended value: 5)}

\item{mutdistance}{the maximum intermutational distance allowed for SNVs to be grouped in the same kataegic event. Recommended value: 3.2))}

\item{segnum}{minimum mutation count. The minimum number of mutations required within a cluster to be identified as kataegic. Recommended value: 4}

\item{ref.dir}{Path to a directory containing the reference genome.}

\item{bed.file}{Path to a directory containing the SNV data in BED format. The BED file should contain the following information: Chromosome, Position, Reference allele, Alternate allele}

\item{output.dir}{Path to a directory where output will be created.}

\item{chromosome}{The chromosome to be analysed. This can be (1, 2, ..., 23, 24) or "all" to run sequentially on all chromosomes.}

\item{chromomsome.length.file}{A tab separated file containing the lengths of all chromosomes in the reference genome.}

\item{trinucleotide.count.file}{A tab seprarated file containing a count of all trinucleotides present in the reference genome. This can be generated with the get.trinucleotide.counts() function in this package.}

}
\description{
Analyzes SNV data in BED format to detect Kataegic events
}
\details{
	The default paramters in SeqKat have been optimized using Alexanrov's "Signatures of mutational processes in human cancer" dataset. SeqKat accepts a BED file and outputs the results in TXT format. A file per chromosome is generated if a kataegic event is detected, otherwise no file is generated. SeqKat reports two scores per kataegic event, a hypermutation score and an APOBEC mediated kataegic score.
}

\examples{
\dontrun{
seqkat(5, 3.2, 4, ref.dir="~/path/to/hg19", bed.file="~/path/to/snv/data/file.bed", output.dir = ".", chromosome = "all", chromosome.length.file = "lengths.tsv", trinucleotide.count.file = "counts.tsv")
}
}
\author{
\email{Fouad.Yousif@oicr.on.ca}
}

