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

\item{bed.dir}{Path to a directory containing the SNV data in BED format. The BED file should contain the following information: Chromosome, Position, Reference allele, Alternate allele}
}
\description{
Analyzes SNV data in BED format to detect Kataegic events
}
\details{
	The default paramters in SeqKat have been optimized using Alexanrov's "Signatures of mutational processes in human cancer" dataset. SeqKat accepts a BED file and outputs the results in TXT format. A file per chromosome is generated if a kataegic event is detected, otherwise no file is generated. SeqKat reports two scores per kataegic event, a hypermutation score and an APOBEC mediated kataegic score.
}

\examples{
seqkat(5, 3.2, 4, ref.dir="~/path/to/hg19", bed.dir="~/path/to/snv/data/")

}
\author{
\email{Fouad.Yousif@oicr.on.ca}
}

