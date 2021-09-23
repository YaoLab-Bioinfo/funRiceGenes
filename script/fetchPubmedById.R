
## Extract information of a publication using the PubMed ID.
## To use this script, please install the RCurl and XML packages.
## Usage: fetchPubmedByID(id="33089895")

fetchPubmedByID <- function(id="33089895") {
  url <- "https://pubmed.ncbi.nlm.nih.gov/"
  finalUrl <- paste(url, id, '/?format=pubmed', sep='')
  urlRes <- readLines(finalUrl)
  
  journal <- urlRes[grepl("^TA\\s", urlRes, ignore.case = FALSE)]
  journal <- gsub(".+-\\s", "", journal)
  title.start <- which(grepl("^TI\\s", urlRes, ignore.case = FALSE))
  title.end <- which(grepl("^PG\\s", urlRes, ignore.case = FALSE))
  title <- urlRes[title.start:(title.end - 1)]
  title <- gsub(".+-\\s", "", title)
  title <- gsub("^\\s+", "", title)
  title <- paste(title, collapse = "", sep="")
  year <- urlRes[grepl("^DEP\\s", urlRes, ignore.case = FALSE)]
  year <- gsub(".+-\\s", "", year)
  year <- substr(year, 1, 4)
  abstract.start <- which(grepl("^AB\\s", urlRes, ignore.case = FALSE))
  abstract.end <- which(grepl("^CI\\s", urlRes, ignore.case = FALSE))
  abstract <- urlRes[abstract.start:(abstract.end - 1)]
  abstract <- gsub(".+-\\s", "", abstract)
  abstract <- gsub("^\\s+", "", abstract)
  abstract <- paste(abstract, collapse = "", sep="")
  affiliation.start <- which(grepl("^AD\\s", urlRes, ignore.case = FALSE))[1]
  affiliation.end <- affiliation.start
  for (i in (affiliation.start+1):(affiliation.start+10)) {
    if (grepl("^\\s+", urlRes[i])) {
      affiliation.end <- i
    } else {
      break
    }
  }
  affiliation <- urlRes[affiliation.start:affiliation.end]
  affiliation <- gsub(".+-\\s", "", affiliation)
  affiliation <- gsub("^\\s+", "", affiliation)
  affiliation <- paste(affiliation, collapse = "", sep="")
  
  return(c(journal, title, year, affiliation, abstract))
}

