#' miRspeciesCodes
#'
#' A data frame containing the code of the species as reported in miRBase
#' 22 release. See also \url{http://www.mirbase.org}
#' 
#' @docType data
#' 
#' @usage data(miRspeciesCodes)
#' 
#' @format A data frame with 271 rows and 2 columns
#' \describe{
#'     \item{code}{is the unique identifier of the species as reported in
#'     miRBase db}
#'     \item{species}{is the species}
#'
#' }
#'
#' @examples
#' data(miRspeciesCodes)
#'
"miRspeciesCodes"

#' gtf
#'
#' A data frame containing a short version of the already formatted gencode v19
#' based-genome annotations. This data frame can be generated with
#' \code{\link{formatGTF}}.
#'
#' @docType data
#' 
#' @usage data(gtf)
#'
#' @examples
#' data(gtf)
#'
"gtf"


#' gwasTraits
#'
#' A data frame containing the traits/disease extracted on the 31th October 2018
#' from the GWAS catalog. See also \url{ https://www.ebi.ac.uk/gwas/} for more
#' detail.
#'
#' @docType data
#' 
#' @usage data(gwasTraits)
#'
#' @format A data frame with 653 rows and 1 columns
#' \describe{
#'     \item{id}{is the trait/disease as reported in the GWAS catalog}
#' }
#'
#' @examples
#' data(gwasTraits)
#'
"gwasTraits"



#' iupac
#'
#' A data frame containing IUPAC codes and the correspoding regular expression.
#' 
#' @docType data
#' 
#' @usage data(iupac)
#' 
#' @format A data frame with 18 rows and 4 columns
#' \describe{
#'     \item{code}{is the IUPAC nucleotide code}
#'     \item{base}{is the base}
#'     \item{regexDNA}{is the regular expression of the corresponding
#'     IUPAC code}
#'     \item{regexRNA}{is the regular expression of the corresponding
#'     IUPAC code}
#' }
#'
#' @examples
#' data(iupac)
#'
"iupac"


#' memeDB
#'
#' A dataframe containing the file paths of the files containing the RBP motifs
#' in meme format available in MEME db. See also \url{http://meme-suite.org/doc/download.html}.
#' File paths extracted from RNA folder of motif_databases.12.19.tgz
#' (Motif Databases updated 28 Oct 2019)
#'
#' @docType data
#' 
#' @usage data(memeDB)
#'
#' @format A character vector with 25 rows ans 2 columns
#' \describe{
#'     \item{path}{is the file path }
#'     \item{index}{index of the file path }
#' }
#'
#' @examples
#' data(memeDB)
#'
"memeDB"

#' attractSpecies
#'
#' A data frame containing the species for which RBP motifs are available
#' in ATtRACT db. See also \url{http://attract.cnic.es}.
#' 
#' @docType data
#' 
#' @usage data(attractSpecies)
#' 
#' @format A data frame with 37 rows and 1 columns
#' \describe{
#'     \item{species}{is the species}
#' }
#'
#' @examples
#' data(attractSpecies)
#'
"attractSpecies"

#' ahRepeatMasker
#'
#' A data frame containing the AnnotationHub id related to the repeatMasker db,
#' specific for each species and genome.
#' 
#' @docType data
#' 
#' @usage data(ahRepeatMasker)
#' 
#' @format A data frame with 84 rows and 3 columns
#' \describe{
#'     \item{id}{is the AnnotationHub id}
#'     \item{species}{is the species}
#'     \item{genome}{is the genome assembly}
#' }
#'
#' @examples
#' data(ahRepeatMasker)
#'

"ahRepeatMasker"

#' ahChainFile
#'
#' A data frame containing the AnnotationHub id related to the *.chain.gz file,
#' specific for each species and genome. The file name reflects the assembly
#' conversion data contained within in the format <db1>To<Db2>.over.chain.gz.
#' For example, a file named mm10ToHg19.over.chain.gz file contains the
#' liftOver data needed to convert mm10 (Mouse GRCm38) coordinates to hg19
#' (Human GRCh37).
#' 
#' @docType data
#' 
#' @usage data(ahChainFiles)
#' 
#' @format A data frame with 1113 rows and 2 columns
#' \describe{
#'     \item{id}{is the AnnotationHub id}
#'     \item{chainFile}{is the *.chain.gz file}
#' }
#'
#' @examples
#' data(ahChainFiles)
#'
"ahChainFiles"


#' backSplicedJunctions
#'
#' A data frame containing genomic coordinates (genome assembly hg19) of
#' circRNAs detected by three detection tools (CircMarker(cm), MapSplice2 (m)
#' and NCLscan (n)) in the human left ventricle tissues of 3 controls, 3
#' patients with dilated cardiomyopathies (DCM) and 3 patients with hypertrophic
#' cardiomyopathies (HCM). This data frame was generated with
#' \code{\link{getBackSplicedJunctions}}.
#' 
#' @docType data
#' 
#' @usage data(backSplicedJunctions)
#' 
#' @format A data frame with 63521 rows and 16 columns
#' \describe{
#'     \item{id}{Unique identifier}
#'     \item{gene}{is the gene name whose exon coordinates overlap that of
#'     the given back-spliced junctions}
#'     \item{strand}{is the strand where the gene is transcribed}
#'     \item{chrom}{the chromosome from which the circRNA is derived}
#'     \item{startUpBSE}{ is the 5' coordinate of the upstream back-spliced
#'     exon in the transcript. This corresponds to the back-spliced junction /
#'     acceptor site}
#'     \item{endDownBSE}{is the 3' coordinate of the downstream back-spliced
#'     exon in the transcript. This corresponds to the back-spliced junction /
#'     donor site}
#'     \item{tool}{ are the tools that identified the back-spliced
#'      junctions}
#'     \item{C1}{Number of occurences of each circRNA in control 1}
#'     \item{C2}{Number of occurences of each circRNA in control 2}
#'     \item{C3}{Number of occurences of each circRNA in control 3}
#'     \item{D1}{Number of occurences of each circRNA in DCM 1}
#'     \item{D2}{Number of occurences of each circRNA in DCM 2}
#'     \item{D3}{Number of occurences of each circRNA in DCM 3}
#'     \item{H1}{Number of occurences of each circRNA in HCM 1}
#'     \item{H2}{Number of occurences of each circRNA in HCM 2}
#'     \item{H3}{Number of occurences of each circRNA in HCM 3}
#' }
#'
#' @examples
#' data(backSplicedJunctions)
#'
"backSplicedJunctions"


#' mergedBSJunctions
#'
#' A data frame containing genomic coordinates (genome assembly hg19) detected
#' in human LV tissues of controls and diseased hearts. This data frame was
#' generated with \code{\link{mergedBSJunctions}} which grouped circRNAs
#' commonly identified by the three tools (CircMarker(cm), MapSplice2 (m) and
#' NCLscan (n)) used for circRNAs detection.
#'
#' @docType data
#' 
#' @usage data(mergedBSJunctions)
#'
#' @format A data frame with 41558 rows and 16 columns
#' \describe{
#'     \item{id}{Unique identifier}
#'     \item{gene}{is the gene name whose exon coordinates overlap that of
#'     the given back-spliced junctions}
#'     \item{strand}{is the strand where the gene is transcribed}
#'     \item{chrom}{the chromosome from which the circRNA is derived}
#'     \item{startUpBSE}{ is the 5' coordinate of the upstream back-spliced
#'     exon in the transcript. This corresponds to the back-spliced junction /
#'     acceptor site}
#'     \item{endDownBSE}{is the 3' coordinate of the downstream back-spliced
#'     exon in the transcript. This corresponds to the back-spliced junction /
#'     donor site}
#'     \item{tool}{ are the tools that identified the back-spliced
#'      junctions}
#'     \item{C1}{Number of occurences of each circRNA in control 1}
#'     \item{C2}{Number of occurences of each circRNA in control 2}
#'     \item{C3}{Number of occurences of each circRNA in control 3}
#'     \item{D1}{Number of occurences of each circRNA in DCM 1}
#'     \item{D2}{Number of occurences of each circRNA in DCM 2}
#'     \item{D3}{Number of occurences of each circRNA in DCM 3}
#'     \item{H1}{Number of occurences of each circRNA in HCM 1}
#'     \item{H2}{Number of occurences of each circRNA in HCM 2}
#'     \item{H3}{Number of occurences of each circRNA in HCM 3}
#' }
#'
#' @examples
#' data(mergedBSJunctions)
#'
"mergedBSJunctions"

