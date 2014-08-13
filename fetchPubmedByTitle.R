
fetchPubmedByTitle <- function(title="")
{
  x <- URLencode(title, reserved=TRUE)
  x <- gsub(",", "%2C", x)
  url <- "http://www.ncbi.nlm.nih.gov/pubmed"
  finalUrl <- paste(url,'?term=(', x,'%5BTitle%5D', ')&report=xml&format=text', sep='')
  urlRes <- getURL(finalUrl)
  if (grepl("title", urlRes, ignore.case=TRUE)) {
    urlRes <- gsub("&lt;","<",urlRes)
    urlRes <- gsub("&gt;",">",urlRes)
    xmlData <- xmlTreeParse(urlRes, useInternalNodes = TRUE)
    journal <- xpathSApply(xmlData, "//PubmedArticle/MedlineCitation/MedlineJournalInfo/MedlineTA", xmlValue)
    title <- xpathSApply(xmlData, "//PubmedArticle/MedlineCitation/Article/ArticleTitle", xmlValue)
    year <- xpathSApply(xmlData, "//PubmedArticle/MedlineCitation/Article/ArticleDate/Year", xmlValue)
    month <- xpathSApply(xmlData, "//PubmedArticle/MedlineCitation/Article/ArticleDate/Month", xmlValue)
    abstract <- xpathSApply(xmlData, "//PubmedArticle/MedlineCitation/Article/Abstract/AbstractText", xmlValue)
    abstract <- paste(abstract,sep="",collapse="")
    
    affiliation <- xpathSApply(xmlData, "//PubmedArticle/MedlineCitation/Article/AuthorList/Author/Affiliation", xmlValue)
    affiliation <- affiliation[1]
    if (is.list(year)) {
      year <- xpathSApply(xmlData, "//PubmedArticle/MedlineCitation/Article/Journal/JournalIssue/PubDate/Year", xmlValue)
      month <- xpathSApply(xmlData, "//PubmedArticle/MedlineCitation/Article/Journal/JournalIssue/PubDate/Month", xmlValue)
    }
    
    return(c(journal, title, year, month, affiliation, abstract))
  } else {
    return(c("", title, "", "", "", ""))
  }
}

