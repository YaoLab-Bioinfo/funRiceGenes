

fetchGenbankByAccession <- function(accession=c("AF045770", "AY242058"),
                                    format="genbank") 
{
  if (format=="genbank") {
    genbank <- sapply(accession, function(x){
      url <- paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=",
                   x, "&rettype=gb&retmode=text", sep = "")
      urlRes <- getURL(url)
      return(urlRes)
    })
  } else if (format=="fasta") {
    genbank <- sapply(accession, function(x){
      url <- paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=",
                   x, "&rettype=fasta&retmode=text", sep = "")
      urlRes <- getURL(url)
      return(urlRes)
    })
  }
  return(genbank)
}


