
## Extract information of a publication using the PubMed ID.
## To use this script, please install the RCurl and XML packages.
## Usage: fetchPubmedByID(id="29051769")

library(RCurl); library(XML)

fetchPubmedByID <- function(id="29051769") {
  url <- "https://www.ncbi.nlm.nih.gov/pubmed"
  finalUrl <- paste(url,'?term=', id, '&report=xml&format=text', sep='')
  urlRes <- getURL(finalUrl)
  if (grepl("title", urlRes, ignore.case=TRUE)) {
    urlRes <- gsub("&lt;","<",urlRes)
    urlRes <- gsub("&gt;",">",urlRes)
    xmlData <- xmlTreeParse(urlRes, useInternalNodes = TRUE)
    journal <- xpathSApply(xmlData, "//PubmedArticle/MedlineCitation/MedlineJournalInfo/MedlineTA", xmlValue)
    title <- xpathSApply(xmlData, "//PubmedArticle/MedlineCitation/Article/ArticleTitle", xmlValue)
    year <- xpathSApply(xmlData, "//PubmedArticle/MedlineCitation/Article/ArticleDate/Year", xmlValue)
    #    month <- xpathSApply(xmlData, "//PubmedArticle/MedlineCitation/Article/ArticleDate/Month", xmlValue)
    abstract <- xpathSApply(xmlData, "//PubmedArticle/MedlineCitation/Article/Abstract/AbstractText", xmlValue)
    abstract <- paste(abstract,sep="",collapse="")
    
    affiliation <- xpathSApply(xmlData, "//PubmedArticle/MedlineCitation/Article/AuthorList/Author/AffiliationInfo/Affiliation", xmlValue)
    affiliation <- affiliation[1]
    if (class(affiliation)!="character") {
      affiliation <- "Fail"
    }
    if (is.list(year)) {
      year <- xpathSApply(xmlData, "//PubmedArticle/MedlineCitation/Article/Journal/JournalIssue/PubDate/Year", xmlValue)
      #      month <- xpathSApply(xmlData, "//PubmedArticle/MedlineCitation/Article/Journal/JournalIssue/PubDate/Month", xmlValue)
    }
    
    return(c(journal, title, year, affiliation, abstract))
  } else {
    return(c("", "", "", "", ""))
  }
}

