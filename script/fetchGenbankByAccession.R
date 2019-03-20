
## Extract Genbank data using sccession number.
## To use this script, please install the RCurl package.
## Usage: fetchGenbankByAccession(accession=c("AY986491", "AY242058"), format="genbank")

library(RCurl)

fetchGenbankByAccession <- function(accession=c("AY986491", "AY242058"),
                                    format="genbank") 
{
  if (format=="genbank") {
    genbank <- sapply(accession, function(x){
      url <- paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=",
                   x, "&rettype=gb&retmode=text", sep = "")
      urlRes <- getURL(url)
      return(urlRes)
    })
  } else if (format=="fasta") {
    genbank <- sapply(accession, function(x){
      url <- paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=",
                   x, "&rettype=fasta&retmode=text", sep = "")
      urlRes <- getURL(url)
      return(urlRes)
    })
  }
  return(genbank)
}


