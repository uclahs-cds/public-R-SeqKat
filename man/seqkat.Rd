\name{seqkat}
\alias{seqkat}
\title{seqkat}
\usage{
seqkat(sigcutoff, mutdistance, segnum, ref.dir, bed.dir = "./")
}
\arguments{
\item{sigcutoff}{The sigcutoff.}

\item{mutdistance}{The mutdistance.}

\item{segnum}{The segnum.}

\item{ref.dir}{Path to a directory containing the reference genome.}

\item{bed.dir}{Path to a directory containing the SNV data in BED format.}
}
\description{
Analyzes SNV data in BED format to detect Kataegic events
}
\examples{
seqkat(5, 3.2, 4, ref.dir="~/path/to/hg19", bed.dir="~/path/to/snv/data/")

}
\author{
\email{Fouad.Yousif@oicr.on.ca}
}

