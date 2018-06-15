
library(RCurl); library(XML); library(stringr); library(plyr); 
mypasswd = "asdfghjkl;'"
system('git config --global user.name "venyao"')
system('git config --global user.email "ywhzau@gmail.com"')

fetchPubmedById <- function(id="") {
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

load("RiceNet_V2_GS.RData")
names(ricenet) <- c("Gene1", "Gene2")
load("rapmsu.rda")
ref.info <- read.table("reference.table", head=T, as.is=T, sep="\t", 
                       quote="", comment="")

fam.gene.info <- read.table("famInfo.table", head=T, sep="\t", as.is=T)

fam.gene.msu <- 1:nrow(fam.gene.info)
names(fam.gene.msu) <- fam.gene.info$MSU
fam.gene.msu <- fam.gene.msu[names(fam.gene.msu)!="None"]
fam.gene.msu.new <- sapply(names(fam.gene.msu), function(x) {
  x.name <- unlist(strsplit(x, split='|', fixed=TRUE))
  if (length(x.name)==1) {
    return(fam.gene.msu[x])
  } else {
    y <- rep(fam.gene.msu[x], length(x.name))
    names(y) <- x.name
    return(y)
  }
})
fam.gene.msu.final <- unlist(unname(fam.gene.msu.new))

fam.gene.rap <- 1:nrow(fam.gene.info)
names(fam.gene.rap) <- fam.gene.info$RAPdb
fam.gene.rap <- fam.gene.rap[names(fam.gene.rap)!="None"]
fam.gene.rap.new <- sapply(names(fam.gene.rap), function(x) {
  x.name <- unlist(strsplit(x, split='|', fixed=TRUE))
  if (length(x.name)==1) {
    return(fam.gene.rap[x])
  } else {
    y <- rep(fam.gene.rap[x], length(x.name))
    names(y) <- x.name
    return(y)
  }
})
fam.gene.rap.final <- unlist(unname(fam.gene.rap.new))


gene.info <- read.table("geneInfo.table", head=T, sep="\t", as.is=T, quote="")
gene.keyword <- read.table("geneKeyword.table", head=T, 
                           sep="\t", as.is=T, quote="", comment="")
all.key <- unique(gene.keyword$Keyword)
all.sym <- sapply(gene.info$Symbol, function(x) {
  symb <- unlist(strsplit(x, split="|", fixed=TRUE))
  return(symb)
})
all.sym <- unlist(all.sym)
names(all.sym) <- NULL

gene.msu <- 1:nrow(gene.info)
names(gene.msu) <- gene.info$MSU
gene.msu <- gene.msu[names(gene.msu)!="None"]
gene.msu.new <- sapply(names(gene.msu), function(x) {
  x.name <- unlist(strsplit(x, split='|', fixed=TRUE))
  if (length(x.name)==1) {
    return(gene.msu[x])
  } else {
    y <- rep(gene.msu[x], length(x.name))
    names(y) <- x.name
    return(y)
  }
})
gene.msu.final <- unlist(unname(gene.msu.new))

gene.rap <- 1:nrow(gene.info)
names(gene.rap) <- gene.info$RAPdb
gene.rap <- gene.rap[names(gene.rap)!="None"]
gene.rap.new <- sapply(names(gene.rap), function(x) {
  x.name <- unlist(strsplit(x, split='|', fixed=TRUE))
  if (length(x.name)==1) {
    return(gene.rap[x])
  } else {
    y <- rep(gene.rap[x], length(x.name))
    names(y) <- x.name
    return(y)
  }
})
gene.rap.final <- unlist(unname(gene.rap.new))

#### publication
fetchRefByKey <- function(keyword="") {
  keyword <- gsub("^\\s+", "", keyword)
  keyword <- gsub(" +$", "", keyword)
  keyword <- tolower(keyword)
  if (nchar(keyword)<1) {
    return(NULL)
  } else {
    datRes <- ref.info[grepl(keyword, ref.info$Title)|grepl(keyword, ref.info$Abstract, ignore.case=TRUE), ]
    if (nrow(datRes)>0) {
      return(datRes)
    } else {
      return(NULL)
    }
  }
}

####  MSU
fetchInfoByMsu <- function(locus="") {
  locus <- gsub("^\\s+", "", locus)
  locus <- gsub(" +$", "", locus)
  locus.line <- gene.msu.final[locus]
  if ( (length(locus.line)==1) && !is.na(locus.line) ) {
    path <- gene.info$path[locus.line]
    gene.fl <- paste(path, "gene.info", sep="/")
    if (file.exists(gene.fl)) {
      dat <- read.table(gene.fl, head=T, sep="\t", as.is=T)
      
      msu <- unlist(strsplit(dat$MSU, split='|', fixed=TRUE))
      msu.new <- sapply(msu, function(x){
        if (x!="None") {
          y <- paste("http://rice.plantbiology.msu.edu/cgi-bin/ORF_infopage.cgi?orf=", 
                     x, sep="")
          y <- paste('<a href="', y, '" target="_blank">', x, '</a>', sep="")
          y <- HTML(y)
          return(y)
        } else {
          return("None")
        }
      })
      msu.new <- paste(unname(msu.new), sep="", collapse="|")
      dat$MSU <- msu.new
      
      rap <- unlist(strsplit(dat$RAPdb, split='|', fixed=TRUE))
      rap.new <- sapply(rap, function(x){
        if (x!="None") {
          y <- paste("http://rapdb.dna.affrc.go.jp/viewer/gbrowse_details/irgsp1?name=", 
                     x, sep="")
          y <- paste('<a href="', y, '" target="_blank">', x, '</a>', sep="")
          y <- HTML(y)
          return(y)
        } else {
          return("None")
        }
      })
      rap.new <- paste(unname(rap.new), sep="", collapse="|")
      dat$RAPdb <- rap.new
      
      return(dat)
    }
  } else {
    return(NULL)
  }
}

fetchFamInfoByMsu <- function(locus="") {
  locus <- gsub("^\\s+", "", locus)
  locus <- gsub(" +$", "", locus)
  locus.line <- fam.gene.msu.final[locus]
  if ( (length(locus.line)==1) && !is.na(locus.line) ) {
      dat <- fam.gene.info[locus.line, ]
      
      msu <- unlist(strsplit(dat$MSU, split='|', fixed=TRUE))
      msu.new <- sapply(msu, function(x){
        if (x!="None") {
          y <- paste("http://rice.plantbiology.msu.edu/cgi-bin/ORF_infopage.cgi?orf=", 
                     x, sep="")
          y <- paste('<a href="', y, '" target="_blank">', x, '</a>', sep="")
          y <- HTML(y)
          return(y)
        } else {
          return("None")
        }
      })
      msu.new <- paste(unname(msu.new), sep="", collapse="|")
      dat$MSU <- msu.new
      
      rap <- unlist(strsplit(dat$RAPdb, split='|', fixed=TRUE))
      rap.new <- sapply(rap, function(x){
        if (x!="None") {
          y <- paste("http://rapdb.dna.affrc.go.jp/viewer/gbrowse_details/irgsp1?name=", 
                     x, sep="")
          y <- paste('<a href="', y, '" target="_blank">', x, '</a>', sep="")
          y <- HTML(y)
          return(y)
        } else {
          return("None")
        }
      })
      rap.new <- paste(unname(rap.new), sep="", collapse="|")
      dat$RAPdb <- rap.new
      
      dat$path <- NULL
      return(dat)
  } else {
    return(NULL)
  }
}

fetchRefByMsu <- function(locus="") {
  locus <- gsub("^\\s+", "", locus)
  locus <- gsub(" +$", "", locus)
  locus.line <- gene.msu.final[locus]
  if ( (length(locus.line)==1) && !is.na(locus.line) ) {
    path <- gene.info$path[locus.line]
    ref.fl <- paste(path, "reference.info", sep="/")
    if (file.exists(ref.fl)) {
      dat <- read.table(ref.fl, head=T, sep="\t", as.is=T, 
                        quote="", comment="")
      dat$Publication <- NULL
      dat$Affiliation <- NULL
      for (i in 1:nrow(dat)) {
        title <- dat$Title[i]
        dat$Title[i] <- paste("https://www.ncbi.nlm.nih.gov/pubmed", 
                              '?term=(', dat$Title[i],'%5BTitle%5D', 
                              ')', sep='')
        dat$Title[i] <- paste('<a href="', dat$Title[i], '" target="_blank">', title, 
                              '</a>', sep="")
        dat$Title[i] <- HTML(dat$Title[i])
      }
      return(dat)
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
}

fetchFamRefByMsu <- function(locus="") {
  locus <- gsub("^\\s+", "", locus)
  locus <- gsub(" +$", "", locus)
  locus.line <- fam.gene.msu.final[locus]
  if ( (length(locus.line)==1) && !is.na(locus.line) ) {
    path <- fam.gene.info$path[locus.line]
    ref.fl <- paste(path, "family.ref", sep="/")
    if (file.exists(ref.fl)) {
      dat <- read.table(ref.fl, head=T, sep="\t", as.is=T, 
                        quote="", comment="")
      dat$Affiliation <- NULL
      for (i in 1:nrow(dat)) {
        title <- dat$Title[i]
        dat$Title[i] <- paste("https://www.ncbi.nlm.nih.gov/pubmed", 
                              '?term=(', dat$Title[i],'%5BTitle%5D', 
                              ')', sep='')
        dat$Title[i] <- paste('<a href="', dat$Title[i], '" target="_blank">', title, 
                              '</a>', sep="")
        dat$Title[i] <- HTML(dat$Title[i])
      }
      return(dat)
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
}

fetchAccByMsu <- function(locus="") {
  locus <- gsub("^\\s+", "", locus)
  locus <- gsub(" +$", "", locus)
  locus.line <- gene.msu.final[locus]
  if ( (length(locus.line)==1) && !is.na(locus.line) ) {
    path <- gene.info$path[locus.line]
    acc.fls <- list.files(path, patter="*acc-*", full=T)
    if (length(acc.fls)>0) {
      acc.fls <- gsub(".+-", "", acc.fls)
      acc.fls <- sapply(acc.fls, function(x){
        y <- paste("https://www.ncbi.nlm.nih.gov/nuccore/", x, sep="")
        y <- paste('<a href="', y, '" target="_blank">', x, '</a>', sep="")
        y <- HTML(y)
        return(y)
      })
      dat <- data.frame(Accession=acc.fls, stringsAsFactors=FALSE)
      return(dat)
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
}

fetchExpByMsu <- function(locus="") {
  locus <- gsub("^\\s+", "", locus)
  locus <- gsub(" +$", "", locus)
  locus.line <- gene.msu.final[locus]
  if ( (length(locus.line)==1) && !is.na(locus.line) ) {
    path <- gene.info$path[locus.line]
    exp.fl <- paste(path, "expression.info", sep="/")
    if (file.exists(exp.fl)) {
      dat <- read.table(exp.fl, head=T, sep="\t", as.is=T)
      
      return(dat)
    }
  } else {
    return(NULL)
  }
}

fetchTextByMsu <- function(locus="") {
  locus <- gsub("^\\s+", "", locus)
  locus <- gsub(" +$", "", locus)
  locus.line <- gene.msu.final[locus]
  if ( (length(locus.line)==1) && !is.na(locus.line) ) {
    path <- gene.info$path[locus.line]
    text.fls <- list.files(path, patter="^pub.text.mining$", full=T)
    if (length(text.fls)>0) {
      text.con <- readLines(text.fls)
      text.con <- text.con[!grepl("^symbol:", text.con)]
      clone <- sapply(seq(2, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("^cloning:\\s+", "", y)
        return(y)
      })
      tdna <- sapply(seq(4, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("^tdna:\\s+", "", y)
        return(y)
      })
      tos <- sapply(seq(5, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("^tos:\\s+", "", y)
        return(y)
      })
      homo <- sapply(seq(6, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("^homol:\\s+", "", y)
        return(y)
      })
      rnai <- sapply(seq(7, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("^rnai:\\s+", "", y)
        return(y)
      })
      ove <- sapply(seq(8, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("^ove:\\s+", "", y)
        return(y)
      })
      rt <- sapply(seq(9, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("^rt:\\s+", "", y)
        return(y)
      })
      north <- sapply(seq(10, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("^N:\\s+", "", y)
        return(y)
      })
      south <- sapply(seq(11, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("^S:\\s+", "", y)
        return(y)
      })
      west <- sapply(seq(12, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("^W:\\s+", "", y)
        return(y)
      })
      access <- sapply(seq(3, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("^access:\\s+", "", y)
        return(y)
      })
      locus <- sapply(seq(1, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("^locus:\\s+", "", y)
        return(y)
      })
      
      dat <- data.frame(Cloning=clone, TDNA=tdna, Tos17=tos, Homolog=homo, 
                        RNAi=rnai, Overexp=ove, RTPCR=rt, Northern=north,
                        Southern=south, Western=west, Locus=locus, 
                        Accession=access,
                        stringsAsFactors=FALSE)
      return(dat)
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
}

fetchKeyByMsu <- function(locus="") {
  locus <- gsub("^\\s+", "", locus)
  locus <- gsub(" +$", "", locus)
  locus.line <- gene.msu.final[locus]
  if ( (length(locus.line)==1) && !is.na(locus.line) ) {
    path <- gene.info$path[locus.line]
    key.fl <- paste(path, "Keyword.trait", sep="/")
    if (file.exists(key.fl)) {
      dat <- read.table(key.fl, head=T, sep="\t", as.is=T, 
                        quote="", comment="")
      for (i in 1:nrow(dat)) {
        title <- dat$Title[i]
        dat$Title[i] <- paste("https://www.ncbi.nlm.nih.gov/pubmed", 
                              '?term=(', dat$Title[i],'%5BTitle%5D', 
                              ')', sep='')
        dat$Title[i] <- paste('<a href="', dat$Title[i], '" target="_blank">', title, 
                              '</a>', sep="")
        dat$Title[i] <- HTML(dat$Title[i])
      }
      return(dat)
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
}

fetchConneByMsu <- function(locus="") {
  locus <- gsub("^\\s+", "", locus)
  locus <- gsub(" +$", "", locus)
  locus.line <- gene.msu.final[locus]
  if ( (length(locus.line)==1) && !is.na(locus.line) ) {
    path <- gene.info$path[locus.line]
    conne.fl <- paste(path, "Connection", sep="/")
    if (file.exists(conne.fl)) {
      dat <- read.table(conne.fl, head=F, sep="\t", as.is=T, 
                        quote="", comment="")
      names(dat) <- c("Symbol1", "Symbol2", "Title", "Evidence")
      for (i in 1:nrow(dat)) {
        title <- dat$Title[i]
        dat$Title[i] <- paste("https://www.ncbi.nlm.nih.gov/pubmed", 
                              '?term=(', dat$Title[i],'%5BTitle%5D', 
                              ')', sep='')
        dat$Title[i] <- paste('<a href="', dat$Title[i], '" target="_blank">', title, 
                              '</a>', sep="")
        dat$Title[i] <- HTML(dat$Title[i])
      }
      return(dat)
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
}

#### RAPdb
fetchInfoByRap <- function(locus="") {
  locus <- gsub("^\\s+", "", locus)
  locus <- gsub(" +$", "", locus)
  locus.line <- gene.rap.final[locus]
  if ( (length(locus.line)==1) && !is.na(locus.line) ) {
    path <- gene.info$path[locus.line]
    gene.fl <- paste(path, "gene.info", sep="/")
    if (file.exists(gene.fl)) {
      dat <- read.table(gene.fl, head=T, sep="\t", as.is=T)
      msu <- unlist(strsplit(dat$MSU, split='|', fixed=TRUE))
      msu.new <- sapply(msu, function(x){
        if (x!="None") {
          y <- paste("http://rice.plantbiology.msu.edu/cgi-bin/ORF_infopage.cgi?orf=", 
                     x, sep="")
          y <- paste('<a href="', y, '" target="_blank">', x, '</a>', sep="")
          y <- HTML(y)
          return(y)
        } else {
          return("None")
        }
      })
      msu.new <- paste(unname(msu.new), sep="", collapse="|")
      dat$MSU <- msu.new
      
      rap <- unlist(strsplit(dat$RAPdb, split='|', fixed=TRUE))
      rap.new <- sapply(rap, function(x){
        if (x!="None") {
          y <- paste("http://rapdb.dna.affrc.go.jp/viewer/gbrowse_details/irgsp1?name=", 
                     x, sep="")
          y <- paste('<a href="', y, '" target="_blank">', x, '</a>', sep="")
          y <- HTML(y)
          return(y)
        } else {
          return("None")
        }
      })
      rap.new <- paste(unname(rap.new), sep="", collapse="|")
      dat$RAPdb <- rap.new
      
      return(dat)
    }
  } else {
    return(NULL)
  }
}

fetchFamInfoByRap <- function(locus="") {
  locus <- gsub("^\\s+", "", locus)
  locus <- gsub(" +$", "", locus)
  locus.line <- fam.gene.rap.final[locus]
  if ( (length(locus.line)==1) && !is.na(locus.line) ) {
      dat <- fam.gene.info[locus.line, ]
      
      msu <- unlist(strsplit(dat$MSU, split='|', fixed=TRUE))
      msu.new <- sapply(msu, function(x){
        if (x!="None") {
          y <- paste("http://rice.plantbiology.msu.edu/cgi-bin/ORF_infopage.cgi?orf=", 
                     x, sep="")
          y <- paste('<a href="', y, '" target="_blank">', x, '</a>', sep="")
          y <- HTML(y)
          return(y)
        } else {
          return("None")
        }
      })
      msu.new <- paste(unname(msu.new), sep="", collapse="|")
      dat$MSU <- msu.new
      
      rap <- unlist(strsplit(dat$RAPdb, split='|', fixed=TRUE))
      rap.new <- sapply(rap, function(x){
        if (x!="None") {
          y <- paste("http://rapdb.dna.affrc.go.jp/viewer/gbrowse_details/irgsp1?name=", 
                     x, sep="")
          y <- paste('<a href="', y, '" target="_blank">', x, '</a>', sep="")
          y <- HTML(y)
          return(y)
        } else {
          return("None")
        }
      })
      rap.new <- paste(unname(rap.new), sep="", collapse="|")
      dat$RAPdb <- rap.new
      dat$path <- NULL
      
      return(dat)
  } else {
    return(NULL)
  }
}

fetchRefByRap <- function(locus="") {
  locus <- gsub("^\\s+", "", locus)
  locus <- gsub(" +$", "", locus)
  locus.line <- gene.rap.final[locus]
  if ( (length(locus.line)==1) && !is.na(locus.line) ) {
    path <- gene.info$path[locus.line]
    ref.fl <- paste(path, "reference.info", sep="/")
    if (file.exists(ref.fl)) {
      dat <- read.table(ref.fl, head=T, sep="\t", as.is=T, 
                        quote="", comment="")
      dat$Publication <- NULL
      dat$Affiliation <- NULL
      for (i in 1:nrow(dat)) {
        title <- dat$Title[i]
        dat$Title[i] <- paste("https://www.ncbi.nlm.nih.gov/pubmed", 
                              '?term=(', dat$Title[i],'%5BTitle%5D', 
                              ')', sep='')
        dat$Title[i] <- paste('<a href="', dat$Title[i], '" target="_blank">', title, 
                              '</a>', sep="")
        dat$Title[i] <- HTML(dat$Title[i])
      }
      return(dat)
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
}

fetchFamRefByRap <- function(locus="") {
  locus <- gsub("^\\s+", "", locus)
  locus <- gsub(" +$", "", locus)
  locus.line <- fam.gene.rap.final[locus]
  if ( (length(locus.line)==1) && !is.na(locus.line) ) {
    path <- fam.gene.info$path[locus.line]
    ref.fl <- paste(path, "family.ref", sep="/")
    if (file.exists(ref.fl)) {
      dat <- read.table(ref.fl, head=T, sep="\t", as.is=T, 
                        quote="", comment="")
      dat$Affiliation <- NULL
      for (i in 1:nrow(dat)) {
        title <- dat$Title[i]
        dat$Title[i] <- paste("https://www.ncbi.nlm.nih.gov/pubmed", 
                              '?term=(', dat$Title[i],'%5BTitle%5D', 
                              ')', sep='')
        dat$Title[i] <- paste('<a href="', dat$Title[i], '" target="_blank">', title, 
                              '</a>', sep="")
        dat$Title[i] <- HTML(dat$Title[i])
      }
      return(dat)
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
}

fetchAccByRap <- function(locus="") {
  locus <- gsub("^\\s+", "", locus)
  locus <- gsub(" +$", "", locus)
  locus.line <- gene.rap.final[locus]
  if ( (length(locus.line)==1) && !is.na(locus.line) ) {
    path <- gene.info$path[locus.line]
    acc.fls <- list.files(path, patter="*acc-*", full=T)
    if (length(acc.fls)>0) {
      acc.fls <- gsub(".+-", "", acc.fls)
      acc.fls <- sapply(acc.fls, function(x){
        y <- paste("https://www.ncbi.nlm.nih.gov/nuccore/", x, sep="")
        y <- paste('<a href="', y, '" target="_blank">', x, '</a>', sep="")
        y <- HTML(y)
        return(y)
      })
      dat <- data.frame(Accession=acc.fls, stringsAsFactors=FALSE)
      return(dat)
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
}

fetchExpByRap <- function(locus="") {
  locus <- gsub("^\\s+", "", locus)
  locus <- gsub(" +$", "", locus)
  locus.line <- gene.rap.final[locus]
  if ( (length(locus.line)==1) && !is.na(locus.line) ) {
    path <- gene.info$path[locus.line]
    exp.fl <- paste(path, "expression.info", sep="/")
    if (file.exists(exp.fl)) {
      dat <- read.table(exp.fl, head=T, sep="\t", as.is=T)
      
      return(dat)
    }
  } else {
    return(NULL)
  }
}

fetchTextByRap <- function(locus="") {
  locus <- gsub("^\\s+", "", locus)
  locus <- gsub(" +$", "", locus)
  locus.line <- gene.rap.final[locus]
  if ( (length(locus.line)==1) && !is.na(locus.line) ) {
    path <- gene.info$path[locus.line]
    text.fls <- list.files(path, patter="^pub.text.mining$", full=T)
    if (length(text.fls)>0) {
      text.con <- readLines(text.fls)
      text.con <- text.con[!grepl("^symbol:", text.con)]
      clone <- sapply(seq(2, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("^cloning:\\s+", "", y)
        return(y)
      })
      tdna <- sapply(seq(4, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("^tdna:\\s+", "", y)
        return(y)
      })
      tos <- sapply(seq(5, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("^tos:\\s+", "", y)
        return(y)
      })
      homo <- sapply(seq(6, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("^homol:\\s+", "", y)
        return(y)
      })
      rnai <- sapply(seq(7, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("^rnai:\\s+", "", y)
        return(y)
      })
      ove <- sapply(seq(8, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("^ove:\\s+", "", y)
        return(y)
      })
      rt <- sapply(seq(9, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("^rt:\\s+", "", y)
        return(y)
      })
      north <- sapply(seq(10, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("^N:\\s+", "", y)
        return(y)
      })
      south <- sapply(seq(11, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("^S:\\s+", "", y)
        return(y)
      })
      west <- sapply(seq(12, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("^W:\\s+", "", y)
        return(y)
      })
      access <- sapply(seq(3, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("^access:\\s+", "", y)
        return(y)
      })
      locus <- sapply(seq(1, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("^locus:\\s+", "", y)
        return(y)
      })
      
      dat <- data.frame(Cloning=clone, TDNA=tdna, Tos17=tos, Homolog=homo, 
                        RNAi=rnai, Overexp=ove, RTPCR=rt, Northern=north,
                        Southern=south, Western=west, Locus=locus, 
                        Accession=access,
                        stringsAsFactors=FALSE)
      return(dat)
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
}

fetchKeyByRap <- function(locus="") {
  locus <- gsub("^\\s+", "", locus)
  locus <- gsub(" +$", "", locus)
  locus.line <- gene.rap.final[locus]
  if ( (length(locus.line)==1) && !is.na(locus.line) ) {
    path <- gene.info$path[locus.line]
    key.fl <- paste(path, "Keyword.trait", sep="/")
    if (file.exists(key.fl)) {
      dat <- read.table(key.fl, head=T, sep="\t", as.is=T, 
                        quote="", comment="")
      for (i in 1:nrow(dat)) {
        title <- dat$Title[i]
        dat$Title[i] <- paste("https://www.ncbi.nlm.nih.gov/pubmed", 
                              '?term=(', dat$Title[i],'%5BTitle%5D', 
                              ')', sep='')
        dat$Title[i] <- paste('<a href="', dat$Title[i], '" target="_blank">', title, 
                              '</a>', sep="")
        dat$Title[i] <- HTML(dat$Title[i])
      }
      return(dat)
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
}

fetchConneByRap <- function(locus="") {
  locus <- gsub("^\\s+", "", locus)
  locus <- gsub(" +$", "", locus)
  locus.line <- gene.rap.final[locus]
  if ( (length(locus.line)==1) && !is.na(locus.line) ) {
    path <- gene.info$path[locus.line]
    conne.fl <- paste(path, "Connection", sep="/")
    if (file.exists(conne.fl)) {
      dat <- read.table(conne.fl, head=F, sep="\t", as.is=T, 
                        quote="", comment="")
      names(dat) <- c("Symbol1", "Symbol2", "Title", "Evidence")
      for (i in 1:nrow(dat)) {
        title <- dat$Title[i]
        dat$Title[i] <- paste("https://www.ncbi.nlm.nih.gov/pubmed", 
                              '?term=(', dat$Title[i],'%5BTitle%5D', 
                              ')', sep='')
        dat$Title[i] <- paste('<a href="', dat$Title[i], '" target="_blank">', title, 
                              '</a>', sep="")
        dat$Title[i] <- HTML(dat$Title[i])
      }
      return(dat)
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
}

findDirBySym <- function(symbol="") {
  line.tar <- sapply(1:nrow(gene.info), function(x){
    sym.line <- gene.info$Symbol[x]
    sym.line <- unlist(strsplit(sym.line, split='|', fixed=TRUE))
    sym.line <- tolower(sym.line)
    if (any(sym.line==symbol)) {
      return(x)
    } else {
      return(NULL)
    }
  })
  return(unlist(unname(line.tar)))
}

findDirByFamSym <- function(symbol="") {
  line.tar <- sapply(1:nrow(fam.gene.info), function(x){
    sym.line <- fam.gene.info$Symbol[x]
    sym.line <- unlist(strsplit(sym.line, split='|', fixed=TRUE))
    sym.line <- tolower(sym.line)
    if (any(sym.line==symbol)) {
      return(x)
    } else {
      return(NULL)
    }
  })
  return(unlist(unname(line.tar)))
}

#### Symbol
fetchInfoBySym <- function(symbol="") {
  symbol <- gsub("^\\s+", "", symbol)
  symbol <- gsub(" +$", "", symbol)
  locus.line <- findDirBySym(tolower(symbol))
  if ( (length(locus.line)==1) && !is.na(locus.line) ) {
    path <- gene.info$path[locus.line]
    gene.fl <- paste(path, "gene.info", sep="/")
    if (file.exists(gene.fl)) {
      dat <- read.table(gene.fl, head=T, sep="\t", as.is=T)
      
      msu <- unlist(strsplit(dat$MSU, split='|', fixed=TRUE))
      msu.new <- sapply(msu, function(x){
        if (x!="None") {
          y <- paste("http://rice.plantbiology.msu.edu/cgi-bin/ORF_infopage.cgi?orf=", 
                     x, sep="")
          y <- paste('<a href="', y, '" target="_blank">', x, '</a>', sep="")
          y <- HTML(y)
          return(y)
        } else {
          return("None")
        }
      })
      msu.new <- paste(unname(msu.new), sep="", collapse="|")
      dat$MSU <- msu.new
      
      rap <- unlist(strsplit(dat$RAPdb, split='|', fixed=TRUE))
      rap.new <- sapply(rap, function(x){
        if (x!="None") {
          y <- paste("http://rapdb.dna.affrc.go.jp/viewer/gbrowse_details/irgsp1?name=", 
                     x, sep="")
          y <- paste('<a href="', y, '" target="_blank">', x, '</a>', sep="")
          y <- HTML(y)
          return(y)
        } else {
          return("None")
        }
      })
      rap.new <- paste(unname(rap.new), sep="", collapse="|")
      dat$RAPdb <- rap.new
      
      return(dat)
    }
  } else if (length(locus.line)>1) {
    dat <- gene.info[locus.line, ]
    dat$path <- NULL
    for (i in 1:nrow(dat)) {    
      msu <- unlist(strsplit(dat$MSU[i], split='|', fixed=TRUE))
      msu.new <- sapply(msu, function(x){
        if (x!="None") {
          y <- paste("http://rice.plantbiology.msu.edu/cgi-bin/ORF_infopage.cgi?orf=", 
                     x, sep="")
          y <- paste('<a href="', y, '" target="_blank">', x, '</a>', sep="")
          y <- HTML(y)
          return(y)
        } else {
          return("None")
        }
      })
      msu.new <- paste(unname(msu.new), sep="", collapse="|")
      dat$MSU[i] <- msu.new
      
      rap <- unlist(strsplit(dat$RAPdb[i], split='|', fixed=TRUE))
      rap.new <- sapply(rap, function(x){
        if (x!="None") {
          y <- paste("http://rapdb.dna.affrc.go.jp/viewer/gbrowse_details/irgsp1?name=", 
                     x, sep="")
          y <- paste('<a href="', y, '" target="_blank">', x, '</a>', sep="")
          y <- HTML(y)
          return(y)
        } else {
          return("None")
        }
      })
      rap.new <- paste(unname(rap.new), sep="", collapse="|")
      dat$RAPdb[i] <- rap.new
    }
    return(dat)
  } else {
    return(NULL)
  }
}

fetchFamInfoBySym <- function(symbol="") {
  symbol <- gsub("^\\s+", "", symbol)
  symbol <- gsub(" +$", "", symbol)
  locus.line <- findDirByFamSym(tolower(symbol))
  if ( (length(locus.line)==1) && !is.na(locus.line) ) {
    dat <- fam.gene.info[locus.line, ]
    dat$path <- NULL
    for (i in 1:nrow(dat)) {    
      msu <- unlist(strsplit(dat$MSU[i], split='|', fixed=TRUE))
      msu.new <- sapply(msu, function(x){
        if (x!="None") {
          y <- paste("http://rice.plantbiology.msu.edu/cgi-bin/ORF_infopage.cgi?orf=", 
                     x, sep="")
          y <- paste('<a href="', y, '" target="_blank">', x, '</a>', sep="")
          y <- HTML(y)
          return(y)
        } else {
          return("None")
        }
      })
      msu.new <- paste(unname(msu.new), sep="", collapse="|")
      dat$MSU[i] <- msu.new
      
      rap <- unlist(strsplit(dat$RAPdb[i], split='|', fixed=TRUE))
      rap.new <- sapply(rap, function(x){
        if (x!="None") {
          y <- paste("http://rapdb.dna.affrc.go.jp/viewer/gbrowse_details/irgsp1?name=", 
                     x, sep="")
          y <- paste('<a href="', y, '" target="_blank">', x, '</a>', sep="")
          y <- HTML(y)
          return(y)
        } else {
          return("None")
        }
      })
      rap.new <- paste(unname(rap.new), sep="", collapse="|")
      dat$RAPdb[i] <- rap.new
    }
    return(dat)
  } else {
    return(NULL)
  }
}

fetchRefBySym <- function(symbol="") {
  symbol <- gsub("^\\s+", "", symbol)
  symbol <- gsub(" +$", "", symbol)
  locus.line <- findDirBySym(tolower(symbol))
  if ( (length(locus.line)==1) && !is.na(locus.line) ) {
    path <- gene.info$path[locus.line]
    ref.fl <- paste(path, "reference.info", sep="/")
    if (file.exists(ref.fl)) {
      dat <- read.table(ref.fl, head=T, sep="\t", as.is=T, 
                        quote="", comment="")
      dat$Publication <- NULL
      dat$Affiliation <- NULL
      for (i in 1:nrow(dat)) {
        title <- dat$Title[i]
        dat$Title[i] <- paste("https://www.ncbi.nlm.nih.gov/pubmed", 
                              '?term=(', dat$Title[i],'%5BTitle%5D', 
                              ')', sep='')
        dat$Title[i] <- paste('<a href="', dat$Title[i], '" target="_blank">', title, 
                              '</a>', sep="")
        dat$Title[i] <- HTML(dat$Title[i])
      }
      return(dat)
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
}

fetchFamRefBySym <- function(symbol="") {
  symbol <- gsub("^\\s+", "", symbol)
  symbol <- gsub(" +$", "", symbol)
  locus.line <- findDirByFamSym(tolower(symbol))
  if ( (length(locus.line)==1) && !is.na(locus.line) ) {
    path <- fam.gene.info$path[locus.line]
    ref.fl <- paste(path, "family.ref", sep="/")
    if (file.exists(ref.fl)) {
      dat <- read.table(ref.fl, head=T, sep="\t", as.is=T, 
                        quote="", comment="")
      dat$Affiliation <- NULL
      for (i in 1:nrow(dat)) {
        title <- dat$Title[i]
        dat$Title[i] <- paste("https://www.ncbi.nlm.nih.gov/pubmed", 
                              '?term=(', dat$Title[i],'%5BTitle%5D', 
                              ')', sep='')
        dat$Title[i] <- paste('<a href="', dat$Title[i], '" target="_blank">', title, 
                              '</a>', sep="")
        dat$Title[i] <- HTML(dat$Title[i])
      }
      return(dat)
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
}

fetchAccBySym <- function(symbol="") {
  symbol <- gsub("^\\s+", "", symbol)
  symbol <- gsub(" +$", "", symbol)
  locus.line <- findDirBySym(tolower(symbol))
  if ( (length(locus.line)==1) && !is.na(locus.line) ) {
    path <- gene.info$path[locus.line]
    acc.fls <- list.files(path, patter="*acc-*", full=T)
    if (length(acc.fls)>0) {
      acc.fls <- gsub(".+-", "", acc.fls)
      acc.fls <- sapply(acc.fls, function(x){
        y <- paste("https://www.ncbi.nlm.nih.gov/nuccore/", x, sep="")
        y <- paste('<a href="', y, '" target="_blank">', x, '</a>', sep="")
        y <- HTML(y)
        return(y)
      })
      dat <- data.frame(Accession=acc.fls, stringsAsFactors=FALSE)
      return(dat)
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
}

fetchExpBySym <- function(symbol="") {
  symbol <- gsub("^\\s+", "", symbol)
  symbol <- gsub(" +$", "", symbol)
  locus.line <- findDirBySym(tolower(symbol))
  if ( (length(locus.line)==1) && !is.na(locus.line) ) {
    path <- gene.info$path[locus.line]
    exp.fl <- paste(path, "expression.info", sep="/")
    if (file.exists(exp.fl)) {
      dat <- read.table(exp.fl, head=T, sep="\t", as.is=T)
      
      return(dat)
    }
  } else {
    return(NULL)
  }
}

fetchTextBySym <- function(symbol="") {
  symbol <- gsub("^\\s+", "", symbol)
  symbol <- gsub(" +$", "", symbol)
  locus.line <- findDirBySym(tolower(symbol))
  if ( (length(locus.line)==1) && !is.na(locus.line) ) {
    path <- gene.info$path[locus.line]
    text.fls <- list.files(path, patter="^pub.text.mining$", full=T)
    if (length(text.fls)>0) {
      text.con <- readLines(text.fls)
      text.con <- text.con[!grepl("^symbol:", text.con)]
      clone <- sapply(seq(2, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("^cloning:\\s+", "", y)
        return(y)
      })
      tdna <- sapply(seq(4, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("^tdna:\\s+", "", y)
        return(y)
      })
      tos <- sapply(seq(5, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("^tos:\\s+", "", y)
        return(y)
      })
      homo <- sapply(seq(6, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("^homol:\\s+", "", y)
        return(y)
      })
      rnai <- sapply(seq(7, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("^rnai:\\s+", "", y)
        return(y)
      })
      ove <- sapply(seq(8, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("^ove:\\s+", "", y)
        return(y)
      })
      rt <- sapply(seq(9, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("^rt:\\s+", "", y)
        return(y)
      })
      north <- sapply(seq(10, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("^N:\\s+", "", y)
        return(y)
      })
      south <- sapply(seq(11, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("^S:\\s+", "", y)
        return(y)
      })
      west <- sapply(seq(12, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("^W:\\s+", "", y)
        return(y)
      })
      access <- sapply(seq(3, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("^access:\\s+", "", y)
        return(y)
      })
      locus <- sapply(seq(1, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("^locus:\\s+", "", y)
        return(y)
      })
      
      dat <- data.frame(Cloning=clone, TDNA=tdna, Tos17=tos, Homolog=homo, 
                        RNAi=rnai, Overexp=ove, RTPCR=rt, Northern=north,
                        Southern=south, Western=west, Locus=locus, 
                        Accession=access,
                        stringsAsFactors=FALSE)
      return(dat)
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
}

fetchKeyBySym <- function(symbol="") {
  symbol <- gsub("^\\s+", "", symbol)
  symbol <- gsub(" +$", "", symbol)
  locus.line <- findDirBySym(tolower(symbol))
  if ( (length(locus.line)==1) && !is.na(locus.line) ) {
    path <- gene.info$path[locus.line]
    key.fl <- paste(path, "Keyword.trait", sep="/")
    if (file.exists(key.fl)) {
      dat <- read.table(key.fl, head=T, sep="\t", as.is=T, 
                        quote="", comment="")
      for (i in 1:nrow(dat)) {
        title <- dat$Title[i]
        dat$Title[i] <- paste("https://www.ncbi.nlm.nih.gov/pubmed", 
                              '?term=(', dat$Title[i],'%5BTitle%5D', 
                              ')', sep='')
        dat$Title[i] <- paste('<a href="', dat$Title[i], '" target="_blank">', title, 
                              '</a>', sep="")
        dat$Title[i] <- HTML(dat$Title[i])
      }
      return(dat)
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
}

fetchConneBySym <- function(symbol="") {
  symbol <- gsub("^\\s+", "", symbol)
  symbol <- gsub(" +$", "", symbol)
  locus.line <- findDirBySym(tolower(symbol))
  if ( (length(locus.line)==1) && !is.na(locus.line) ) {
    path <- gene.info$path[locus.line]
    conne.fl <- paste(path, "Connection", sep="/")
    if (file.exists(conne.fl)) {
      dat <- read.table(conne.fl, head=F, sep="\t", as.is=T, 
                        quote="", comment="")
      names(dat) <- c("Symbol1", "Symbol2", "Title", "Evidence")
      for (i in 1:nrow(dat)) {
        title <- dat$Title[i]
        dat$Title[i] <- paste("https://www.ncbi.nlm.nih.gov/pubmed", 
                              '?term=(', dat$Title[i],'%5BTitle%5D', 
                              ')', sep='')
        dat$Title[i] <- paste('<a href="', dat$Title[i], '" target="_blank">', title, 
                              '</a>', sep="")
        dat$Title[i] <- HTML(dat$Title[i])
      }
      return(dat)
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
}

#### Keyword
fetchInfoByKey <- function(keyword="") {
  keyword <- gsub("^\\s+", "", keyword)
  keyword <- gsub(" +$", "", keyword)
  all.key <- unique(gene.keyword$Keyword)
  if (tolower(keyword) %in% all.key) {
    dat <- gene.keyword[gene.keyword$Keyword==tolower(keyword), ]
    dat$Keyword <- NULL
    dat$path <- NULL
    for (i in 1:nrow(dat)) {    
      msu <- unlist(strsplit(dat$MSU[i], split='|', fixed=TRUE))
      msu.new <- sapply(msu, function(x){
        y <- paste("http://rice.plantbiology.msu.edu/cgi-bin/ORF_infopage.cgi?orf=", 
                   x, sep="")
        y <- paste('<a href="', y, '" target="_blank">', x, '</a>', sep="")
        y <- HTML(y)
      })
      msu.new <- paste(unname(msu.new), sep="", collapse="|")
      dat$MSU[i] <- msu.new
      
      rap <- unlist(strsplit(dat$RAPdb[i], split='|', fixed=TRUE))
      rap.new <- sapply(rap, function(x){
        y <- paste("http://rapdb.dna.affrc.go.jp/viewer/gbrowse_details/irgsp1?name=", 
                   x, sep="")
        y <- paste('<a href="', y, '" target="_blank">', x, '</a>', sep="")
        y <- HTML(y)
      })
      rap.new <- paste(unname(rap.new), sep="", collapse="|")
      dat$RAPdb[i] <- rap.new
      
      title <- dat$Title[i]
      dat$Title[i] <- paste("https://www.ncbi.nlm.nih.gov/pubmed", 
                            '?term=(', dat$Title[i],'%5BTitle%5D', 
                            ')', sep='')
      dat$Title[i] <- paste('<a href="', dat$Title[i], '" target="_blank">', title, 
                            '</a>', sep="")
      dat$Title[i] <- HTML(dat$Title[i])
    }
    
    return(dat)
  } else {
    return(NULL)
  }
}

#### submit new information to this database
write.gene <- function(df) {
  symbol <- df$Symbol
  symbol.low <- tolower(symbol)
  if (grepl("^os", symbol.low)) {
    symbol.low <- gsub("^os", "", symbol.low)
    symbol.head <- substr(symbol.low, 1, 1)
    if (symbol.head %in% letters[1:26]) {
      dir.to <- paste("data/Gene/Abstract/OS", 
                      toupper(symbol.head), sep="/")
    } else {
      dir.to <- "data/Gene/Abstract/OS/0-9"
    }
  } else {
    symbol.head <- substr(symbol.low, 1, 1)
    if (symbol.head %in% letters[1:26]) {
      dir.to <- paste("data/Gene/Abstract", 
                      toupper(symbol.head), sep="/")
    } else {
      dir.to <- "data/Gene/Abstract/0-9"
    }
  }
  
  symbol <- gsub('|', "~", symbol, fixed=TRUE)
  dir.to <- paste(dir.to, symbol, sep="/")
  if (!file.exists(dir.to)) {
    dir.create(dir.to)
    out.fl <- paste(dir.to, "gene.info", sep="/")
    write.table(df, file=out.fl, sep="\t", quote=F, row.names=F)
    
    df$path <- dir.to
    gene.info.new <- rbind(gene.info, df)
    gene.info.new <- gene.info.new[order(gene.info.new$Symbol), ]
    write.table(gene.info.new, file="geneInfo.table", 
                sep="\t", quote=F, row.names=F)
  }
}

write.pub <- function(df) {
  if (df$Symbol=="") {
    df.sub <- df[, c("Title", "Year", "Journal", "Affiliation", "Abstract", "Symbol")]
    df.sub$Symbol <- "None"
    names(df.sub)[6] <- "Gene"
    ref.info <- rbind(ref.info, df.sub)
    ref.info.new <- ddply(ref.info, .(Title, Year, Journal, Affiliation, Abstract), function(df){
      symbol <- paste(df$Gene, sep=",", collapse=",")
      return(symbol)
    })
    names(ref.info.new)[6] <- "Gene"
    write.table(ref.info.new, file="reference.table", sep="\t", quote=F, row.names=F)
  } else {
    symbol <- df$Symbol
    symbol <- gsub("^\\s+", "", symbol)
    symbol <- gsub(" +$", "", symbol)
    locus.line <- findDirBySym(tolower(symbol))
    df$Symbol <- gene.info$Symbol[locus.line]
    if (length(locus.line)==1) {
      path <- gene.info$path[locus.line]
      out.fl <- paste(path, "reference.info", sep="/")
      df.sub <- df[, c("Title", "Year", "Journal", "Affiliation", "Abstract", "Symbol")]
      if (file.exists(out.fl)) {
        df$Publication <- NA
        df <- df[, c("Publication", "Title", "Year", "Journal", "Affiliation", "Abstract")]
        df.tmp <- read.table(out.fl, sep="\t", quote="", head=T, as.is=T, comment="")
        df.new <- unique(rbind(df.tmp, df))
        write.table(df.new, file=out.fl, sep="\t", quote=F, row.names=F)
      } else {
        df$Publication <- NA
        df <- df[, c("Publication", "Title", "Year", "Journal", "Affiliation", "Abstract")]
        write.table(df, file=out.fl, sep="\t", quote=F, row.names=F)
      }
      
      names(df.sub)[6] <- "Gene"
      ref.info <- rbind(ref.info, df.sub)
      ref.info <- unique(ref.info)
      ref.info.new <- ddply(ref.info, .(Title, Year, Journal, Affiliation, Abstract), function(df){
        symbol <- paste(df$Gene, sep=",", collapse=",")
        return(symbol)
      })
      names(ref.info.new)[6] <- "Gene"
      write.table(ref.info.new, file="reference.table", sep="\t", quote=F, row.names=F)
    }
  }
}

write.acc <- function(df) {
  symbol <- df$Symbol
  symbol <- gsub("^\\s+", "", symbol)
  symbol <- gsub(" +$", "", symbol)
  locus.line <- findDirBySym(tolower(symbol))
  if (length(locus.line)==1) {
    path <- gene.info$path[locus.line]
    out.fl <- paste("acc-", df$Accession, sep="")
    out.fl <- paste(path, out.fl, sep="/")
    if (!file.exists(out.fl)) {
      file.create(out.fl)
    }
  }
}

write.key <- function(df) {
  symbol <- df$Symbol[1]
  symbol <- gsub("^\\s+", "", symbol)
  symbol <- gsub(" +$", "", symbol)
  locus.line <- findDirBySym(tolower(symbol))
  df$Symbol <- gene.info$Symbol[locus.line]
  if (length(locus.line)==1) {
    path <- gene.info$path[locus.line]
    out.fl <- paste(path, "Keyword.trait", sep="/")
    if (file.exists(out.fl)) {
      df.tmp <- read.table(out.fl, sep="\t", quote="", head=T, as.is=T, comment="")
      df.new <- unique(rbind(df.tmp, df))
      write.table(df.new, file=out.fl, sep="\t", quote=F, row.names=F)
      
      df.new$RAPdb <- gene.info$RAPdb[locus.line]
      df.new$MSU <- gene.info$MSU[locus.line]
      df.new$path <- gene.info$path[locus.line]
      df.new <- df.new[, c("Symbol","RAPdb","MSU","Keyword","Title", "path")]
      gene.keyword <- gene.keyword[gene.keyword$path!=path, ]
      gene.keyword.new <- rbind(gene.keyword, df.new)
      gene.keyword.new <- gene.keyword.new[order(gene.keyword.new$Symbol), ]
      write.table(gene.keyword.new, file="geneKeyword.table", 
                  sep="\t", quote=F, row.names=F)
    } else {
      write.table(df, file=out.fl, sep="\t", quote=F, row.names=F)
      
      df$RAPdb <- gene.info$RAPdb[locus.line]
      df$MSU <- gene.info$MSU[locus.line]
      df$path <- gene.info$path[locus.line]
      df <- df[, c("Symbol","RAPdb","MSU","Keyword","Title", "path")]
      gene.keyword.new <- rbind(gene.keyword, df)
      gene.keyword.new <- gene.keyword.new[order(gene.keyword.new$Symbol), ]
      write.table(gene.keyword.new, file="geneKeyword.table", 
                  sep="\t", quote=F, row.names=F)
    }
  }
}

scanAndWriteKey <- function(df) {
  title <- df$Title; abstract <- df$Abstract; symbol <- df$Symbol
  all.sent <- unlist(strsplit(abstract, split="\\."))
  all.sent <- c(title, all.sent)
  all.sent <- all.sent[grepl(symbol, all.sent, ignore.case=TRUE)]
  lstRes <- lapply(all.key, function(key) {
    if (any(grepl(key, all.sent, ignore.case=TRUE))) {
      return(cbind(key, all.sent[grepl(key, all.sent, ignore.case=TRUE)]))
    }
  })
  dfRes <- do.call(rbind, lstRes)
  dfRes <- data.frame(dfRes, stringsAsFactors=FALSE)
  if (nrow(dfRes)>0) {
    names(dfRes) <- c("Keyword", "Evidence")
    dfRes$Symbol <- symbol
    dfRes$Title <- title
    dfRes <- dfRes[, c("Symbol",  "Keyword",	"Title",	"Evidence")]
    line.right <- vector()
    for (i in 1:nrow(dfRes)) {
      evid.words <- tolower(unlist(strsplit(dfRes$Evidence[i], split="\\s+")))
      if (tolower(dfRes$Symbol[i])%in%evid.words && grepl("\\s+", dfRes$Keyword[i])) {
        line.right <- c(line.right, i)
      }
      if (tolower(dfRes$Symbol[i])%in%evid.words && tolower(dfRes$Keyword[i])%in%evid.words) {
        line.right <- c(line.right, i)
      }
    }
    if (length(line.right)>0) {
      dfRes <- dfRes[line.right, ]
      write.key(dfRes)
    }
  }
}

write.exp <- function(df) {
  symbol <- df$Symbol
  symbol <- gsub("^\\s+", "", symbol)
  symbol <- gsub(" +$", "", symbol)
  locus.line <- findDirBySym(tolower(symbol))
  df$Symbol <- gene.info$Symbol[locus.line]
  if (length(locus.line)==1) {
    path <- gene.info$path[locus.line]
    out.fl <- paste(path, "expression.info", sep="/")
    if (file.exists(out.fl)) {
      df.tmp <- read.table(out.fl, sep="\t", quote="", head=T, as.is=T, comment="")
      df.new <- unique(rbind(df.tmp, df))
      write.table(df.new, file=out.fl, sep="\t", quote=F, row.names=F)
    } else {
      write.table(df, file=out.fl, sep="\t", quote=F, row.names=F)
    }
  }
}

scanAndWriteExp <- function(df) {
  title <- df$Title; abstract <- df$Abstract; symbol <- df$Symbol
  all.sent <- unlist(strsplit(abstract, split="\\."))
  all.sent <- c(title, all.sent)
  all.sent <- all.sent[grepl(symbol, all.sent, ignore.case=TRUE)]
  all.exp <- c("expression", "overexpression", "rnai")
  lstRes <- lapply(all.exp, function(exp) {
    if (any(grepl(exp, all.sent, ignore.case=TRUE))) {
      return(cbind(exp, all.sent[grepl(exp, all.sent, ignore.case=TRUE)]))
    }
  })
  dfRes <- do.call(rbind, lstRes)
  dfRes <- data.frame(dfRes, stringsAsFactors=FALSE)
  if (nrow(dfRes)>0) {
    names(dfRes) <- c("Keyword", "Evidence")
    dfRes$Symbol <- symbol
    dfRes$Title <- title
    df.new <- NULL
    df.new$Symbol <- symbol
    if (any(dfRes$Keyword=="expression")) {
      df.new$Expression <- paste(dfRes$Evidence[dfRes$Keyword=="expression"], sep=" | ", collapse=" | ")
    }
    if (any(dfRes$Keyword=="overexpression")) {
      df.new$Overexpression <- paste(dfRes$Evidence[dfRes$Keyword=="overexpression"], sep=" | ", collapse=" | ")
    }
    if (any(dfRes$Keyword=="rnai")) {
      df.new$RNAi <- paste(dfRes$Evidence[dfRes$Keyword=="rnai"], sep=" | ", collapse=" | ")
    }
    write.exp(df.new)
  }
}

write.con <- function(df) {
  symbol.1 <- df$Symbol1
  symbol.2 <- df$Symbol2
  locus.line.1 <- findDirBySym(tolower(symbol.1))
  locus.line.2 <- findDirBySym(tolower(symbol.2))
  df$Symbol1 <- gene.info$Symbol[locus.line.1]
  df$Symbol2 <- gene.info$Symbol[locus.line.2]
  
  if (length(locus.line.1)==1 && length(locus.line.2)==1 && locus.line.1!=locus.line.2) {
    path.1 <- gene.info$path[locus.line.1]
    out.fl.1 <- paste(path.1, "Connection", sep="/")
    path.2 <- gene.info$path[locus.line.2]
    out.fl.2 <- paste(path.2, "Connection", sep="/")
    if (file.exists(out.fl.1)) {
      df.tmp <- read.table(out.fl.1, sep="\t", quote="", head=F, as.is=T, comment="")
      names(df) <- names(df.tmp)
      df.new <- unique(rbind(df.tmp, df))
      df.new$v1n <- pmin(df.new$V1, df.new$V2)
      df.new$v2n <- pmax(df.new$V1, df.new$V2)
      df.new <- df.new[, c("v1n", "v2n", "V3", "V4")]
      df.new <- unique(df.new)
      write.table(df.new, file=out.fl.1, sep="\t", quote=F, row.names=F, col.names=F)
    } else {
      names(df) <- paste("V", 1:4, sep="")
      df$v1n <- pmin(df$V1, df$V2)
      df$v2n <- pmax(df$V1, df$V2)
      df <- df[, c("v1n", "v2n", "V3", "V4")]
      df <- unique(df)
      write.table(df, file=out.fl.1, sep="\t", quote=F, row.names=F,
                  col.names=F)
    }
    if (file.exists(out.fl.2)) {
      df.tmp <- read.table(out.fl.2, sep="\t", quote="", head=F, as.is=T, comment="")
      names(df) <- names(df.tmp)
      df.new <- unique(rbind(df.tmp, df))
      df.new$v1n <- pmin(df.new$V1, df.new$V2)
      df.new$v2n <- pmax(df.new$V1, df.new$V2)
      df.new <- df.new[, c("v1n", "v2n", "V3", "V4")]
      df.new <- unique(df.new)
      write.table(df.new, file=out.fl.2, sep="\t", quote=F, row.names=F, col.names=F)
    } else {
      names(df) <- paste("V", 1:4, sep="")
      df$v1n <- pmin(df$V1, df$V2)
      df$v2n <- pmax(df$V1, df$V2)
      df <- df[, c("v1n", "v2n", "V3", "V4")]
      df <- unique(df)
      write.table(df, file=out.fl.2, sep="\t", quote=F, row.names=F,
                  col.names=F)
    }
  }
}

scanAndWriteCon <- function(df) {
  title <- df$Title; abstract <- df$Abstract; symbol <- df$Symbol
  all.sent <- unlist(strsplit(abstract, split="\\."))
  all.sent <- c(title, all.sent)
  all.sent <- all.sent[grepl(symbol, all.sent, ignore.case=TRUE)]
  all.sym <- setdiff(all.sym, symbol)
  if (length(all.sent)>0) {
    lstRes <- lapply(all.sym, function(sym) {
      if (any(grepl(sym, all.sent, ignore.case=TRUE))) {
        return(cbind(sym, all.sent[grepl(sym, all.sent, ignore.case=TRUE)]))
      }
    })
    dfRes <- do.call(rbind, lstRes)
    dfRes <- data.frame(dfRes, stringsAsFactors=FALSE)
    if (nrow(dfRes)>0) {
      names(dfRes) <- c("Symbol2", "Evidence")
      dfRes$Symbol1 <- symbol
      dfRes$Title <- title
      dfRes <- dfRes[, c("Symbol1", "Symbol2", "Title", "Evidence")]
      for (i in 1:nrow(dfRes)) {
        evid.words <- tolower(unlist(strsplit(dfRes$Evidence[i], split="\\s+")))
        evid.words <- gsub(",$", "", evid.words)
        evid.words <- gsub("\\.$", "", evid.words)
        evid.words <- gsub("\\)$", "", evid.words)
        evid.words <- gsub("^\\(", "", evid.words)
        if (tolower(dfRes$Symbol1[i])%in%evid.words && tolower(dfRes$Symbol2[i])%in%evid.words) {
          write.con(dfRes[i, ])
        }
      }
    }
  } 
}

gene.edit <- function(df){
  gene.tar <- gene.info[gene.info$Symbol==df$oldsym&gene.info$RAPdb==df$oldrap&
                          gene.info$MSU==df$oldmsu,]
  df.dat <- df[, c("newsym", "newrap", "newmsu")]
  names(df.dat) <- c("Symbol", "RAPdb", "MSU")
  tar.fl <- paste(gene.tar$path, "gene.info", sep="/")
  if (file.exists(tar.fl)) {
    write.table(df.dat, file=tar.fl, sep="\t", quote=F, row.names=F)
    
    gene.info <- gene.info[gene.info$path!=gene.tar$path, ]
    df.dat$path <- gene.tar$path
    gene.info.new <- rbind(gene.info, df.dat)
    gene.info.new <- gene.info.new[order(gene.info.new$Symbol), ]
    write.table(gene.info.new, file="geneInfo.table", 
                sep="\t", quote=F, row.names=F)
  }
  
  key.fl <- paste(gene.tar$path, "Keyword.trait", sep="/")
  if (file.exists(key.fl)) {
    key.con <- read.table(key.fl, sep="\t", head=T, as.is=T, quote="", comment="")
    key.con$Symbol <- df$newsym
    write.table(key.con, file=key.fl, sep="\t", quote=F, row.names=F)
    
    gene.keyword <- gene.keyword[gene.keyword$path!=gene.tar$path, ]
    key.con$path <- gene.tar$path
    key.con$RAPdb <- df.dat$RAPdb
    key.con$MSU <- df.dat$MSU
    key.con <- key.con[, c("Symbol","RAPdb","MSU","Keyword","Title", "path")]
    gene.keyword.new <- rbind(gene.keyword, key.con)
    gene.keyword.new <- gene.keyword.new[order(gene.keyword.new$Symbol), ]
    write.table(gene.keyword.new, file="geneKeyword.table", 
                sep="\t", quote=F, row.names=F)
  }
  
  cone.fl <- paste(gene.tar$path, "Connection", sep="/")
  if (file.exists(cone.fl)) {
    cone.con <- read.table(cone.fl, sep="\t", head=F, as.is=T, quote="", comment="")
    for (i in 1:nrow(cone.con)) {
      if (cone.con$V1[i]==df$oldsym) {
        cone.con$V1[i] <- df$newsym
        cone.2.sym <- cone.con$V2[i]
        cone.2.tar <- gene.info[gene.info$Symbol==cone.2.sym, ]
        if (nrow(cone.2.tar)==1) {
          cone.2.fl <- paste(cone.2.tar$path, "Connection", sep="/")
          cone.2.con <- read.table(cone.2.fl, sep="\t", head=F, as.is=T, quote="", comment="")
          for (j in 1:nrow(cone.2.con)) {
            if (cone.2.con$V1[j]==df$oldsym) {cone.2.con$V1[j] <- df$newsym}
            if (cone.2.con$V2[j]==df$oldsym) {cone.2.con$V2[j] <- df$newsym}
          }
          write.table(cone.2.con, file=cone.2.fl, sep="\t", quote=F, row.names=F, col.names=F)
        }
        
      } else if (cone.con$V2[i]==df$oldsym) {
        cone.con$V2[i] <- df$newsym
        cone.2.sym <- cone.con$V1[i]
        cone.2.tar <- gene.info[gene.info$Symbol==cone.2.sym, ]
        if (nrow(cone.2.tar)==1) {
          cone.2.fl <- paste(cone.2.tar$path, "Connection", sep="/")
          cone.2.con <- read.table(cone.2.fl, sep="\t", head=F, as.is=T, quote="", comment="")
          for (j in 1:nrow(cone.2.con)) {
            if (cone.2.con$V1[j]==df$oldsym) {cone.2.con$V1[j] <- df$newsym}
            if (cone.2.con$V2[j]==df$oldsym) {cone.2.con$V2[j] <- df$newsym}
          }
          write.table(cone.2.con, file=cone.2.fl, sep="\t", quote=F, row.names=F, col.names=F)
        }
      }
    }
    write.table(cone.con, file=cone.fl, sep="\t", quote=F, row.names=F, col.names=F)
  }
  
}

updateGeneInfo <- function() {
  
  gene.info <<- read.table("geneInfo.table", head=T, sep="\t", as.is=T, quote="")
  gene.keyword <<- read.table("geneKeyword.table", head=T, 
                             sep="\t", as.is=T, quote="", comment="")
  ref.info <<- read.table("reference.table", head=T, as.is=T, sep="\t", 
                         quote="", comment="")
  
  gene.msu <- 1:nrow(gene.info)
  names(gene.msu) <- gene.info$MSU
  gene.msu <- gene.msu[names(gene.msu)!="None"]
  gene.msu.new <- sapply(names(gene.msu), function(x) {
    x.name <- unlist(strsplit(x, split='|', fixed=TRUE))
    if (length(x.name)==1) {
      return(gene.msu[x])
    } else {
      y <- rep(gene.msu[x], length(x.name))
      names(y) <- x.name
      return(y)
    }
  })
  gene.msu.final <<- unlist(unname(gene.msu.new))
  
  gene.rap <- 1:nrow(gene.info)
  names(gene.rap) <- gene.info$RAPdb
  gene.rap <- gene.rap[names(gene.rap)!="None"]
  gene.rap.new <- sapply(names(gene.rap), function(x) {
    x.name <- unlist(strsplit(x, split='|', fixed=TRUE))
    if (length(x.name)==1) {
      return(gene.rap[x])
    } else {
      y <- rep(gene.rap[x], length(x.name))
      names(y) <- x.name
      return(y)
    }
  })
  gene.rap.final <<- unlist(unname(gene.rap.new))
}

fetchInfoByChoice <- function(query="", text="") {
  if (query=="MSU Locus") {
    return(fetchInfoByMsu(text))
  } else if (query=="RAPdb Locus") {
    return(fetchInfoByRap(text))
  } else if (query=="Gene Symbol") {
    return(fetchInfoBySym(text))
  }
}

fetchRefByChoice <- function(query="", text="") {
  if (query=="MSU Locus") {
    return(fetchRefByMsu(text))
  } else if (query=="RAPdb Locus") {
    return(fetchRefByRap(text))
  } else if (query=="Gene Symbol") {
    return(fetchRefBySym(text))
  }
}

fetchAccByChoice <- function(query="", text="") {
  if (query=="MSU Locus") {
    return(fetchAccByMsu(text))
  } else if (query=="RAPdb Locus") {
    return(fetchAccByRap(text))
  } else if (query=="Gene Symbol") {
    return(fetchAccBySym(text))
  }
}

fetchExpByChoice <- function(query="", text="") {
  if (query=="MSU Locus") {
    return(fetchExpByMsu(text))
  } else if (query=="RAPdb Locus") {
    return(fetchExpByRap(text))
  } else if (query=="Gene Symbol") {
    return(fetchExpBySym(text))
  }
}

fetchTextByChoice <- function(query="", text="") {
  if (query=="MSU Locus") {
    return(fetchTextByMsu(text))
  } else if (query=="RAPdb Locus") {
    return(fetchTextByRap(text))
  } else if (query=="Gene Symbol") {
    return(fetchTextBySym(text))
  }
}

fetchKeyByChoice <- function(query="", text="") {
  if (query=="MSU Locus") {
    return(fetchKeyByMsu(text))
  } else if (query=="RAPdb Locus") {
    return(fetchKeyByRap(text))
  } else if (query=="Gene Symbol") {
    return(fetchKeyBySym(text))
  }
}

fetchConneByChoice <- function(query="", text="") {
  if (query=="MSU Locus") {
    return(fetchConneByMsu(text))
  } else if (query=="RAPdb Locus") {
    return(fetchConneByRap(text))
  } else if (query=="Gene Symbol") {
    return(fetchConneBySym(text))
  }
}

query.intext <- c("LOC_Os06g40780", "Os05g0158500", "Ghd7")
names(query.intext) <- c("MSU Locus", "RAPdb Locus", "Gene Symbol")

fetchFamInfoByChoice <- function(query="", text="") {
  if (query=="MSU Locus") {
    return(fetchFamInfoByMsu(text))
  } else if (query=="RAPdb Locus") {
    return(fetchFamInfoByRap(text))
  } else if (query=="Gene Symbol") {
    return(fetchFamInfoBySym(text))
  }
}

fetchFamRefByChoice <- function(query="", text="") {
  if (query=="MSU Locus") {
    return(fetchFamRefByMsu(text))
  } else if (query=="RAPdb Locus") {
    return(fetchFamRefByRap(text))
  } else if (query=="Gene Symbol") {
    return(fetchFamRefBySym(text))
  }
}

query.intext.fam <- c("LOC_Os10g38470", "Os02g0677300", "RCN1")
names(query.intext.fam) <- c("MSU Locus", "RAPdb Locus", "Gene Symbol")

convMSU <- function(locus="Os02g0677300") {
  if (is.null(locus) || is.na(locus)) {
    return(NULL)
  } else {
    locus <- unlist(strsplit(locus, split="\\s+"))
    datRes <- rapmsu[rapmsu$rap %in% locus,]
    names(datRes) <- c("RAPdb", "MSU")
    return(datRes)
  }
}

convRap <- function(locus="LOC_Os03g57940") {
  if (is.null(locus) || is.na(locus)) {
    return(NULL)
  } else {
    locus <- unlist(strsplit(locus, split="\\s+"))
    datRes <- rapmsu[rapmsu$msu %in% locus,]
    names(datRes) <- c("RAPdb", "MSU")
    return(datRes)
  }
}

convID <- function(query="", text="") {
  if (query=="RAPdb to MSU") {
    return(convMSU(text))
  } else if (query=="MSU to RAPdb") {
    return(convRap(text))
  }
}

query.intext.conv <- c("Os02g0677300", "LOC_Os03g57940")
names(query.intext.conv) <- c("RAPdb to MSU", "MSU to RAPdb")

query.IJ.conv <- c("LOC_Os03g57940", "Os02g0677300", "MH01g0011000", "ZS01g0017300")
names(query.IJ.conv) <- c("MSU Nipponbare", "RAPdb Nipponbare", "Minghui 63", "Zhenshan 97")

genes.NMZ <- read.table("NIP-MH63-ZS97.txt", head=T, as.is=T, sep="\t")
names(genes.NMZ) <- c("Nipponbare.MSU", "Nipponbare.RAPdb", "Minghui 63", "Zhenshan 97")

geneID <- function(query="MSU Nipponbare", text="LOC_Os03g57940") {
  if (is.null(text) || is.na(text)) {
    return(NULL)
  } else {
    text <- unlist(strsplit(text, split="\\s+"))
    text <- setdiff(text, "None")
    if (query=="MSU Nipponbare") {
      dat.res <- genes.NMZ[genes.NMZ[,1] %in% text, ]
    } else if (query=="RAPdb Nipponbare") {
      dat.res <- genes.NMZ[genes.NMZ[,2] %in% text, ]
    } else if (query=="Minghui 63") {
      dat.res <- genes.NMZ[genes.NMZ[,3] %in% text, ]
    } else if (query=="Zhenshan 97") {
      dat.res <- genes.NMZ[genes.NMZ[,4] %in% text, ]
    }
    
    dat.res$Nipponbare.MSU <- sapply(dat.res$Nipponbare.MSU, function(x){
      if (x!="None") {
        y <- paste("http://rice.plantbiology.msu.edu/cgi-bin/ORF_infopage.cgi?orf=", x, sep="")
        y <- paste('<a href="', y, '" target="_blank">', x, '</a>', sep="")
        y <- HTML(y)
        return(y)
      } else {
        return("None")
      }
    })
    
    dat.res$Nipponbare.RAPdb <- sapply(dat.res$Nipponbare.RAPdb, function(x){
      if (x!="None") {
        y <- paste("http://rapdb.dna.affrc.go.jp/viewer/gbrowse_details/irgsp1?name=", x, sep="")
        y <- paste('<a href="', y, '" target="_blank">', x, '</a>', sep="")
        y <- HTML(y)
        return(y)
      } else {
        return("None")
      }
    })
    
    dat.res$`Minghui 63` <- sapply(dat.res$`Minghui 63`, function(x){
      if (x!="None") {
        y <- paste("http://rice.hzau.edu.cn/cgi-bin/rice/gene?org=MH63&locus=", x, sep="")
        y <- paste('<a href="', y, '" target="_blank">', x, '</a>', sep="")
        y <- HTML(y)
        return(y)
      } else {
        return("None")
      }
    })
    
    dat.res$`Zhenshan 97` <- sapply(dat.res$`Zhenshan 97`, function(x){
      if (x!="None") {
        y <- paste("http://rice.hzau.edu.cn/cgi-bin/rice/gene?org=ZS97&locus=", x, sep="")
        y <- paste('<a href="', y, '" target="_blank">', x, '</a>', sep="")
        y <- HTML(y)
        return(y)
      } else {
        return("None")
      }
    })
    
    return(dat.res)
  }
}

save.image <- function(symbol="", phenofig="", expfig="") {
  symbol <- gsub("^\\s+", "", symbol)
  symbol <- gsub(" +$", "", symbol)
  locus.line <- findDirBySym(tolower(symbol))
  if (length(locus.line)==1) {
    path <- gene.info$path[locus.line]
    phenofig.source <- phenofig$datapath
    expfig.source <- expfig$datapath
    phenofig.suffix <- gsub(".+\\.", "", basename(phenofig$type))
    expfig.suffix <- gsub(".+\\.", "", basename(expfig$type))
    phenofig.target <- paste(path, "/", symbol, ".pheno.", phenofig.suffix, sep="")
    expfig.target <- paste(path, "/", symbol, ".exp.", expfig.suffix, sep="")
    file.copy(phenofig.source, phenofig.target)
    file.copy(expfig.source, expfig.target)
  }
}

fetctRiceNetByMSU <- function(locus="") {
  ricenet[ricenet$Gene1 %in% locus | ricenet$Gene2 %in% locus, ]
}

fetctRiceNetByRAPdb <- function(locus="") {
  locus <- gsub("^\\s+", "", locus)
  locus <- gsub(" +$", "", locus)
  locus.line <- gene.rap.final[locus]
  if ( (length(locus.line)==1) && !is.na(locus.line) ) {
    path <- gene.info$path[locus.line]
    gene.fl <- paste(path, "gene.info", sep="/")
    if (file.exists(gene.fl)) { 
      dat <- read.table(gene.fl, head=T, sep="\t", as.is=T)
      msu <- unlist(strsplit(dat$MSU, split='|', fixed=TRUE))
      return(ricenet[ricenet$Gene1 %in% msu | ricenet$Gene2 %in% msu, ])
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
}

fetctRiceNetBySym <- function(symbol="") {
  symbol <- gsub("^\\s+", "", symbol)
  symbol <- gsub(" +$", "", symbol)
  locus.line <- findDirBySym(tolower(symbol))
  if ( (length(locus.line)==1) && !is.na(locus.line) ) {
    path <- gene.info$path[locus.line]
    gene.fl <- paste(path, "gene.info", sep="/")
    if (file.exists(gene.fl)) {
      dat <- read.table(gene.fl, head=T, sep="\t", as.is=T)
      msu <- unlist(strsplit(dat$MSU, split='|', fixed=TRUE))
      return(ricenet[ricenet$Gene1 %in% msu | ricenet$Gene2 %in% msu, ])
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
}

fetctRiceNetByChoice <- function(query="", text="") {
  if (query=="MSU Locus") {
    return(fetctRiceNetByMSU(text))
  } else if (query=="RAPdb Locus") {
    return(fetctRiceNetByRAPdb(text))
  } else if (query=="Gene Symbol") {
    return(fetctRiceNetBySym(text))
  }
}

#### Shiny
shinyServer(function(input, output, session) {
  
  output$inText <- renderUI({
    textInput("inText", label=NULL, 
              value=query.intext[input$query])
  })
  
  output$mytable1 = renderDataTable({
    fetchInfoByChoice(input$query, input$inText)
  }, options = list(lengthMenu = c(1, 2), pageLength = 1, searching = FALSE, autoWidth = FALSE), escape = FALSE
  )
  
  output$mytable2 = renderDataTable({
    fetchRefByChoice(input$query, input$inText)
  }, options = list(lengthMenu = c(1, 2, 4), pageLength = 1,
                    searching = FALSE, autoWidth = FALSE), escape = FALSE
  )
  
  output$mytable3 = renderDataTable({
    fetchAccByChoice(input$query, input$inText)
  }, options = list(lengthMenu = c(2, 4, 6), searching = FALSE,
                    pageLength = 2, autoWidth = FALSE), escape = FALSE)

  output$mytable4 = renderDataTable({
    fetchExpByChoice(input$query, input$inText)
  }, options = list(lengthMenu = c(1, 2, 4), searching = FALSE,
                    pageLength = 1, autoWidth = FALSE), escape = FALSE)
  
  output$mytable5 = renderDataTable({
    fetchKeyByChoice(input$query, input$inText)
  }, options = list(lengthMenu = c(1, 2, 4), searching = FALSE,
                    pageLength = 1, autoWidth = FALSE), escape = FALSE)
  
  output$mytable6 = renderDataTable({
    fetchConneByChoice(input$query, input$inText)
  }, options = list(lengthMenu = c(1, 2, 4), searching = FALSE,
                    pageLength = 1, autoWidth = FALSE), escape = FALSE)
  
  output$mytable7 = renderDataTable({
    fetchInfoByKey(input$keyword)
  }, options = list(lengthMenu = c(2, 4, 6), searching = FALSE,
                    pageLength = 2, autoWidth = FALSE), escape = FALSE)
  
  output$inTextfam <- renderUI({
    textInput("inTextfam", label=NULL, 
              value=query.intext.fam[input$queryfam])
  })
  
  output$mytable8 = renderDataTable({
    fetchFamInfoByChoice(input$queryfam, input$inTextfam)
  }, options = list(lengthMenu = c(1, 2), pageLength = 1, searching = FALSE, autoWidth = FALSE), escape = FALSE
  )
  
  output$mytable9 = renderDataTable({
    fetchFamRefByChoice(input$queryfam, input$inTextfam)
  }, options = list(lengthMenu = c(1, 2, 4), pageLength = 1,
                    searching = FALSE, autoWidth = FALSE), escape = FALSE
  )

  output$inTextconv <- renderUI({
    textInput("inTextconv", label=NULL, 
            value=query.intext.conv[input$queryconv])
  })

  output$mytable10 = renderDataTable({
    convID(input$queryconv, input$inTextconv)
  }, options = list(lengthMenu = c(1, 4), pageLength = 4,
                  searching = FALSE, autoWidth = FALSE), escape = FALSE
  )
  
  output$inIJconv <- renderUI({
    textInput("inIJconv", label=NULL, 
              value=query.IJ.conv[input$IJconv])
  })
  
  output$mytable12 = renderDataTable({
    geneID(input$IJconv, input$inIJconv)
  }, options = list(lengthMenu = c(1, 4), pageLength = 4,
                    searching = FALSE, autoWidth = FALSE), escape = FALSE
  )

  output$mytable11 = renderDataTable({
    fetchRefByKey(input$publication)
  }, options = list(lengthMenu = c(1, 2, 4, 6), searching = FALSE,
                  pageLength = 1, autoWidth = FALSE), escape = FALSE)
  
  output$mytable13 = renderDataTable({
    fetctRiceNetByChoice(input$query, input$inText)
  }, options = list(lengthMenu = c(2, 4, 8), pageLength = 4,
                    searching = FALSE, autoWidth = FALSE), escape = FALSE
  )
  
  gene.info.NM <- apply(gene.info, 1, function(x){
    x.msu <- unlist(strsplit(x[3], split="\\|"))
    res <- lapply(1:length(x.msu), function(y){
      return(c(x[1:2], x.msu[y]))
    })
    return(do.call(rbind, res))
  })
  gene.info.NM <- as.data.frame(do.call(rbind, gene.info.NM), stringsAsFactors=FALSE)
  
  gene.info.NP <- apply(gene.info, 1, function(x){
    x.rap <- unlist(strsplit(x[2], split="\\|"))
    res <- lapply(1:length(x.rap), function(y){
      return(c(x[1], x.rap[y], x[3]))
    })
    return(do.call(rbind, res))
  })
  gene.info.NP <- as.data.frame(do.call(rbind, gene.info.NP), stringsAsFactors=FALSE)
  
  fam.gene.info.NM <- apply(fam.gene.info, 1, function(x){
    x.msu <- unlist(strsplit(x[3], split="\\|"))
    res <- lapply(1:length(x.msu), function(y){
      return(c(x[1:2], x.msu[y]))
    })
    return(do.call(rbind, res))
  })
  fam.gene.info.NM <- as.data.frame(do.call(rbind, fam.gene.info.NM), stringsAsFactors=FALSE)
  
  fam.gene.info.NP <- apply(fam.gene.info, 1, function(x){
    x.rap <- unlist(strsplit(x[2], split="\\|"))
    res <- lapply(1:length(x.rap), function(y){
      return(c(x[1], x.rap[y], x[3]))
    })
    return(do.call(rbind, res))
  })
  fam.gene.info.NP <- as.data.frame(do.call(rbind, fam.gene.info.NP), stringsAsFactors=FALSE)
  
  output$dMsuInfo <- downloadHandler(
    filename <- function(){paste('locusInfo.csv')},
    content <- function(file) {
      in.locus <- unlist(strsplit(input$msuarea, split="\\n"))
      in.locus <- gsub("^\\s+", "", in.locus)
      in.locus <- gsub("\\s+$", "", in.locus)
      dat.tmp.1 <- gene.info[gene.info$Symbol %in% unique(gene.info.NM$Symbol[gene.info.NM$MSU %in% in.locus]), ]
      dat.tmp.2 <- fam.gene.info[fam.gene.info$Symbol %in% unique(fam.gene.info.NM$Symbol[fam.gene.info.NM$MSU %in% in.locus]), ]
      dat.tmp.2 <- dat.tmp.2[!dat.tmp.2$MSU %in% dat.tmp.1$MSU, ]
      dat.tmp <- rbind(dat.tmp.1, dat.tmp.2)
      dat.tmp$path <- NULL
      
      write.csv(dat.tmp, file, row.names=FALSE)
    }, contentType = "text/csv"
  )
  
  output$dRapInfo <- downloadHandler(
    filename <- function(){paste('locusInfo.csv')},
    content <- function(file) {
      in.locus <- unlist(strsplit(input$raparea, split="\\n"))
      in.locus <- gsub("^\\s+", "", in.locus)
      in.locus <- gsub("\\s+$", "", in.locus)
      dat.tmp.1 <- gene.info[gene.info$Symbol %in% unique(gene.info.NP$Symbol[gene.info.NP$RAP %in% in.locus]), ]
      dat.tmp.2 <- fam.gene.info[fam.gene.info$Symbol %in% unique(fam.gene.info.NP$Symbol[fam.gene.info.NP$RAP %in% in.locus]), ]
      dat.tmp.2 <- dat.tmp.2[!dat.tmp.2$RAP
                             %in% dat.tmp.1$RAP, ]
      dat.tmp <- rbind(dat.tmp.1, dat.tmp.2)
      dat.tmp$path <- NULL
      
      write.csv(dat.tmp, file, row.names=FALSE)
    }, contentType = "text/csv"
  )
  
  output$dMsuKey <- downloadHandler(
    filename <- function(){paste('locusKeyword.txt')},
    content <- function(file) {
      in.locus <- unlist(strsplit(input$msuarea, split="\\n"))
      in.locus <- gsub("^\\s+", "", in.locus)
      in.locus <- gsub("\\s+$", "", in.locus)
      dat.tmp <- gene.keyword[gene.keyword$Symbol %in% unique(gene.info.NM$Symbol[gene.info.NM$MSU %in% in.locus]), ]
      
      dat.tmp$path <- NULL
      write.table(dat.tmp, file, row.names=FALSE, sep="\t", quote=F)
    }, contentType = "text/txt"
  )
  
  output$dRapKey <- downloadHandler(
    filename <- function(){paste('locusKeyword.txt')},
    content <- function(file) {
      in.locus <- unlist(strsplit(input$raparea, split="\\n"))
      in.locus <- gsub("^\\s+", "", in.locus)
      in.locus <- gsub("\\s+$", "", in.locus)
      dat.tmp <- gene.keyword[gene.keyword$Symbol %in% unique(gene.info.NP$Symbol[gene.info.NP$RAP %in% in.locus]), ]
      
      dat.tmp$path <- NULL
      write.table(dat.tmp, file, row.names=FALSE, sep="\t", quote=F)
    }, contentType = "text/txt"
  )
  
  output$dMsuPub <- downloadHandler(
    filename <- function(){paste('locusPublication.txt')},
    content <- function(file) {
      in.locus <- unlist(strsplit(input$msuarea, split="\\n"))
      in.locus <- gsub("^\\s+", "", in.locus)
      in.locus <- gsub("\\s+$", "", in.locus)
      in.sym <- unique(gene.info.NM$Symbol[gene.info.NM$MSU %in% in.locus])
      ref.info.n <- apply(ref.info, 1, function(x){
        x.gene <- unlist(strsplit(x[6], split=","))
        res <- lapply(1:length(x.gene), function(y){
          return(c(x[1:5], x.gene[y]))
        })
        return(do.call(rbind, res))
      })
      ref.info.n <- do.call(rbind, ref.info.n)
      
      dat.tmp <- ref.info.n[ref.info.n[, "Gene"] %in% in.sym, ]
      
      write.table(dat.tmp, file, row.names=FALSE, sep="\t", quote=F)
    }, contentType = "text/txt"
  )
  
  output$dRapPub <- downloadHandler(
    filename <- function(){paste('locusPublication.txt')},
    content <- function(file) {
      in.locus <- unlist(strsplit(input$raparea, split="\\n"))
      in.locus <- gsub("^\\s+", "", in.locus)
      in.locus <- gsub("\\s+$", "", in.locus)
      in.sym <- unique(gene.info.NP$Symbol[gene.info.NP$RAP %in% in.locus])
      ref.info.n <- apply(ref.info, 1, function(x){
        x.gene <- unlist(strsplit(x[6], split=","))
        res <- lapply(1:length(x.gene), function(y){
          return(c(x[1:5], x.gene[y]))
        })
        return(do.call(rbind, res))
      })
      ref.info.n <- do.call(rbind, ref.info.n)
      
      dat.tmp <- ref.info.n[ref.info.n[, "Gene"] %in% in.sym, ]
      
      write.table(dat.tmp, file, row.names=FALSE, sep="\t", quote=F)
    }, contentType = "text/txt"
  )
  
  output$dMsuRiceNet <- downloadHandler(
    filename <- function(){paste('locusRiceNet.txt')},
    content <- function(file) {
      in.locus <- unlist(strsplit(input$msuarea, split="\\n"))
      in.locus <- gsub("^\\s+", "", in.locus)
      in.locus <- gsub("\\s+$", "", in.locus)
      dat.tmp <- ricenet[ricenet$Gene1 %in% in.locus | ricenet$Gene2 %in% in.locus, ]
      
      write.table(dat.tmp, file, row.names=FALSE, sep="\t", quote=F)
    }, contentType = "text/txt"
  )
  
  # submit new reference
  observe({
    if (input$submit2>0) {
      isolate({
        if (input$key2==mypasswd) {
		symsub2 <- gsub("^\\s+", "", input$symsub2)
		symsub2 <- gsub("\\s+$", "", symsub2)
		yearsub2 <- gsub("^\\s+", "", input$yearsub2)
		yearsub2 <- gsub("\\s+$", "", yearsub2)
		jousub2 <- gsub("^\\s+", "", input$jousub2)
		jousub2 <- gsub("\\s+$", "", jousub2)

          df.pub <- data.frame(Symbol=symsub2, Title=input$tilsub2, Year=yearsub2,
                               Journal=jousub2, Affiliation=input$afisub2, Abstract=input$abssub2,
                               stringsAsFactors=FALSE)
          write.pub(df.pub)
          scanAndWriteKey(df.pub)
          scanAndWriteCon(df.pub)
          scanAndWriteExp(df.pub)
          updateGeneInfo()
          
          git.info <- paste("add new reference for ", symsub2, sep="")
          system("git checkout master")
          system("git add -A")
          system(paste('git commit -m ', '"', git.info, '"', sep=""))
          
          js_string <- 'alert("New reference successfully added!");'
          session$sendCustomMessage(type='jsCode', list(value = js_string))
        } else {
          js_string <- 'alert("Authorization Required!");'
          session$sendCustomMessage(type='jsCode', list(value = js_string))
        }
      })
    } else {NULL}
  })
  
  # submit new accession
  observe({
    if (input$submit3>0) {
      isolate({
        if (input$key3==mypasswd) {
		symsub3 <- gsub("^\\s+", "", input$symsub3)
		symsub3 <- gsub("\\s+$", "", symsub3)
		accsub3 <- gsub("^\\s+", "", input$accsub3)
		accsub3 <- gsub("^\\s+", "", accsub3)

          df.acc <- data.frame(Symbol=symsub3, Accession=accsub3, 
                               stringsAsFactors=FALSE)
          write.acc(df.acc)
          updateGeneInfo()
          
          git.info <- paste("add new accessions for ", symsub3, sep="")
          system("git checkout master")
          system("git add -A")
          system(paste('git commit -m ', '"', git.info, '"', sep=""))
          
          js_string <- 'alert("Accessions successfully added!");'
          session$sendCustomMessage(type='jsCode', list(value = js_string))
        } else {
          js_string <- 'alert("Authorization Required!");'
          session$sendCustomMessage(type='jsCode', list(value = js_string))
        }
      })
    } else {NULL}
  })
  
  # submit new keywords
  observe({
    if (input$submit4>0) {
      isolate({
        if (input$key4==mypasswd) {
		symsub4 <- gsub("^\\s+", "", input$symsub4)
		symsub4 <- gsub("\\s+$", "", symsub4)
		keysub4 <- gsub("^\\s+", "", input$keysub4)
		keysub4 <- gsub("\\s+$", "", keysub4)

          df.key <- data.frame(Symbol=symsub4, Keyword=keysub4, 
                               Title=tilsub4, Evidence=evisub4,
                               stringsAsFactors=FALSE)
          write.key(df.key)
          updateGeneInfo()
          
          git.info <- paste("add new keywords for ", symsub4, sep="")
          system("git checkout master")
          system("git add -A")
          system(paste('git commit -m ', '"', git.info, '"', sep=""))
          
          js_string <- 'alert("Keywords successfully added!");'
          session$sendCustomMessage(type='jsCode', list(value = js_string))
        } else {
          js_string <- 'alert("Authorization Required!");'
          session$sendCustomMessage(type='jsCode', list(value = js_string))
        }
      })
    } else {NULL}
  })
  
  # submit new connections
  observe({
    if (input$submit5>0) {
      isolate({
        if (input$key5==mypasswd) {
		symsub5 <- gsub("^\\s+", "", input$symsub5)
		symsub5 <- gsub("\\s+$", "", symsub5)
		sym2sub5 <- gsub("^\\s+", "", input$sym2sub5)
		sym2sub5 <- gsub("\\s+$", "", sym2sub5)

          df.con <- data.frame(Symbol1=symsub5, Symbol2=sym2sub5, 
                               Title=input$tilsub5, Evidence=input$evisub5,
                               stringsAsFactors=FALSE)
          write.con(df.con)
          updateGeneInfo()
          
          git.info <- paste("add new connections between ", 
                            symsub5, " and ", sym2sub5, sep="")
          system("git checkout master")
          system("git add -A")
          system(paste('git commit -m ', '"', git.info, '"', sep=""))
          
          js_string <- 'alert("Info successfully added!");'
          session$sendCustomMessage(type='jsCode', list(value = js_string))
          
        } else {
          js_string <- 'alert("Authorization Required!");'
          session$sendCustomMessage(type='jsCode', list(value = js_string))
        }
      })
    } else {NULL}
  })
  
  #submit new figures
  observe({
    if (input$submit8>0) {
      isolate({
        if (input$key8==mypasswd) {
		symsub8 <- gsub("^\\s+", "", input$symsub8)
		symsub8 <- gsub("\\s+$", "", symsub8)

           save.image(symsub8, input$phenofig, input$expfig)
          
           git.info <- paste("add new figures for ", symsub8, sep="")
           system("git checkout master")
           system("git add -A")
           system(paste('git commit -m ', '"', git.info, '"', sep=""))
          
           js_string <- 'alert("Figures successfully added!");'
           session$sendCustomMessage(type='jsCode', list(value = js_string))
        } else {
           js_string <- 'alert("Authorization Required!");'
           session$sendCustomMessage(type='jsCode', list(value = js_string))
        }
      })
    } else {NULL}
  })

  # edit gene info
  observe({
    if (input$submit6>0) {
      isolate({
        if (input$key6==mypasswd) {
          oldsym <- gsub("^\\s+", "", input$oldsym)
          oldsym <- gsub("\\s+$", "", oldsym)
          newsym <- gsub("^\\s+", "", input$newsym)
          newsym <- gsub("\\s+$", "", newsym)
          oldmsu <- gsub("^\\s+", "", input$oldmsu)
          oldmsu <- gsub("\\s+$", "", oldmsu)
        
          newmsu <- gsub("^\\s+", "", input$newmsu)
          newmsu <- gsub("\\s+$", "", newmsu)
          oldrap <- gsub("^\\s+", "", input$oldrap)
          oldrap <- gsub("\\s+$", "", oldrap)
          newrap <- gsub("^\\s+", "", input$newrap)
          newrap <- gsub("\\s+$", "", newrap)
          
          df.con <- data.frame(oldsym=oldsym, newsym=newsym, 
                               oldmsu=oldmsu, newmsu=newmsu,
                               oldrap=oldrap, newrap=newrap,
                               stringsAsFactors=FALSE)
          gene.edit(df.con)
          updateGeneInfo()
          
          git.info <- gsub("\\|", " == ", newsym)
          system("git checkout master")
          system("git add -A")
          system(paste('git commit -m ', '"', git.info, '"', sep=""))
          
          js_string <- 'alert("Edit gene info successfully!");'
          session$sendCustomMessage(type='jsCode', list(value = js_string))
        } else {
          js_string <- 'alert("Authorization Required!");'
          session$sendCustomMessage(type='jsCode', list(value = js_string))
        }
      })
    } else {NULL}
  })
  
  # submit gene info
  observe({
     if (input$submit7>0) {
       isolate({
         if (input$key7==mypasswd) {
           symsub7 <- gsub("^\\s+", "", input$symsub7)
           symsub7 <- gsub("\\s+$", "", symsub7)
           msusub7 <- gsub("^\\s+", "", input$msusub7)
           msusub7 <- gsub("\\s+$", "", msusub7)
           rapsub7 <- gsub("^\\s+", "", input$rapsub7)
           rapsub7 <- gsub("\\s+$", "", rapsub7)
           
           if (symsub7!="") {
             df.gene <- data.frame(Symbol=symsub7, RAPdb=rapsub7, MSU=msusub7,
                                   stringsAsFactors=FALSE)
             locus.line <- findDirBySym(tolower(symsub7))
             if (length(locus.line)==0) {
               write.gene(df.gene)
               updateGeneInfo()
               
               git.info <- paste("add new gene: ", symsub7, sep="") 
               system("git checkout master")
               system("git add -A")
               system(paste('git commit -m ', '"', git.info, '"', sep=""))
          
               js_string <- 'alert("Add new gene successfully!");'
			   session$sendCustomMessage(type='jsCode', list(value = js_string))
             }
           }

	   pubmed7 <- gsub("^\\s+", "", input$pubmed7)
	   pubmed7 <- gsub("\\s+$", "", pubmed7)
           pubmedRes <- fetchPubmedById(pubmed7)
           if (all(pubmedRes!="")) {  
             df.pub <- data.frame(Symbol=symsub7, Title=pubmedRes[2], Year=pubmedRes[3],
                                  Journal=pubmedRes[1], Affiliation=pubmedRes[4], Abstract=pubmedRes[5],
                                  stringsAsFactors=FALSE)
             write.pub(df.pub)
             git.info <- "add new pub."
             system("git checkout master")
             system("git add -A")
             system(paste('git commit -m ', '"', git.info, '"', sep=""))
          
             if (symsub7=="") {
				js_string <- 'alert("Add new pub successfully!");'
				session$sendCustomMessage(type='jsCode', list(value = js_string))
			 }

             if (symsub7!="") {
               scanAndWriteKey(df.pub)
               scanAndWriteCon(df.pub)
               scanAndWriteExp(df.pub)
               
			   git.info <- paste("add new info for gene: ", symsub7, sep="") 
               system("git checkout master")
               system("git add -A")
               system(paste('git commit -m ', '"', git.info, '"', sep=""))

			   js_string <- 'alert("Add new info successfully!");'
	           session$sendCustomMessage(type='jsCode', list(value = js_string))
             }
             updateGeneInfo()
           }

         } else {
           js_string <- 'alert("Authorization Required!");'
           session$sendCustomMessage(type='jsCode', list(value = js_string))
         }
       })
     } else {NULL}
   })

   observe({
	  if (input$clear1>0) {
		isolate({
			updateTextInput(session, "symsub7", value="")
			updateTextInput(session, "msusub7", value="")
			updateTextInput(session, "rapsub7", value="")
			updateTextInput(session, "pubmed7", value="")
		})
      } else {NULL}
   })

   observe({
	  if (input$clear2>0) {
		isolate({
			updateTextInput(session, "oldsym", value="")
			updateTextInput(session, "newsym", value="")
			updateTextInput(session, "oldmsu", value="")
			updateTextInput(session, "newmsu", value="")
			updateTextInput(session, "oldrap", value="")
			updateTextInput(session, "newrap", value="")
		})
      } else {NULL}
   })
 
   observe({
	  if (input$clear3>0) {
		isolate({
			updateTextInput(session, "symsub4", value="")
			updateTextInput(session, "keysub4", value="")
			updateTextInput(session, "tilsub4", value="")
			updateTextInput(session, "evisub4", value="")
		})
      } else {NULL}
   })

   observe({
	  if (input$clear4>0) {
		isolate({
			updateTextInput(session, "symsub2", value="")
			updateTextInput(session, "tilsub2", value="")
			updateTextInput(session, "yearsub2", value="")
			updateTextInput(session, "jousub2", value="")
			updateTextInput(session, "afisub2", value="")
			updateTextInput(session, "abssub2", value="")
		})
      } else {NULL}
   })

   observe({
	  if (input$clear5>0) {
		isolate({
			updateTextInput(session, "symsub3", value="")
			updateTextInput(session, "accsub3", value="")
		})
      } else {NULL}
   })

   observe({
	  if (input$clear7>0) {
		isolate({
			updateTextInput(session, "symsub5", value="")
			updateTextInput(session, "sym2sub5", value="")
			updateTextInput(session, "tilsub5", value="")
			updateTextInput(session, "evisub5", value="")
		})
      } else {NULL}
   })

   observe({
	  if (input$clear8>0) {
		isolate({
			updateTextInput(session, "symsub8", value="")
		})
      } else {NULL}
   })

   observe({
	  if (input$clear6>0) {
		isolate({
			updateTextInput(session, "pubmed10", value="")
		})
      } else {NULL}
   })

   # submit gene family
   observe({
     if (input$submit10>0) {
       isolate({
         if (input$key10==mypasswd) {
           df.genefam <- read.table(input$genfamin$datapath, head=T, sep="\t", as.is=T)
           
           symbol <- df.genefam$Name[1]
           symbol.low <- tolower(symbol)
           if (grepl("^os", symbol.low)) {
             symbol.low <- gsub("^os", "", symbol.low)
             symbol.head <- substr(symbol.low, 1, 1)
             if (symbol.head %in% letters[1:26]) {
               dir.to <- paste("data/Family/OS", 
                               toupper(symbol.head), sep="/")
             } else {
               dir.to <- "data/Family/OS/0-9"
             }
           } else {
             symbol.head <- substr(symbol.low, 1, 1)
             if (symbol.head %in% letters[1:26]) {
               dir.to <- paste("data/Family/", 
                               toupper(symbol.head), sep="/")
             } else {
               dir.to <- "data/Family/0-9"
             }
           }
           
           symbol <- gsub("\\s+", "_", symbol)
           dir.to <- paste(dir.to, symbol, sep="/")
           dir.to <- gsub("/+", "/", dir.to)
           if (!file.exists(dir.to)) {
             dir.create(dir.to)
           }
           file.to <- paste(dir.to, "family.info", sep="/")
           write.table(df.genefam, file=file.to, sep="\t", quote=F, row.names=F)
           df.genefam <- df.genefam[, c("Symbol", "RAPdb", "MSU")]
           
           df.genefam$path <- dir.to
           fam.info.new <- rbind(fam.gene.info, df.genefam)
           fam.info.new <- fam.info.new[order(fam.info.new$Symbol), ]
           write.table(fam.info.new, file="famInfo.table", 
                       sep="\t", quote=F, row.names=F)
           fam.gene.info <<- read.table("famInfo.table", head=T, sep="\t", as.is=T)
           fam.gene.msu <- 1:nrow(fam.gene.info)
           names(fam.gene.msu) <- fam.gene.info$MSU
           fam.gene.msu <- fam.gene.msu[names(fam.gene.msu)!="None"]
           fam.gene.msu.new <- sapply(names(fam.gene.msu), function(x) {
             x.name <- unlist(strsplit(x, split='|', fixed=TRUE))
             if (length(x.name)==1) {
               return(fam.gene.msu[x])
             } else {
               y <- rep(fam.gene.msu[x], length(x.name))
               names(y) <- x.name
               return(y)
             }
           })
           fam.gene.msu.final <<- unlist(unname(fam.gene.msu.new))
           
           fam.gene.rap <- 1:nrow(fam.gene.info)
           names(fam.gene.rap) <- fam.gene.info$RAPdb
           fam.gene.rap <- fam.gene.rap[names(fam.gene.rap)!="None"]
           fam.gene.rap.new <- sapply(names(fam.gene.rap), function(x) {
             x.name <- unlist(strsplit(x, split='|', fixed=TRUE))
             if (length(x.name)==1) {
               return(fam.gene.rap[x])
             } else {
               y <- rep(fam.gene.rap[x], length(x.name))
               names(y) <- x.name
               return(y)
             }
           }) 
           fam.gene.rap.final <<- unlist(unname(fam.gene.rap.new))
          
	   pubmed10 <- gsub("^\\s+", "", input$pubmed10)
	   pubmed10 <- gsub("\\s+$", "", pubmed10)
           pubmedRes <- fetchPubmedById(input$pubmed10)
           if (all(pubmedRes!="")) {  
             df.pub <- data.frame(Journal=pubmedRes[1], Title=pubmedRes[2], Year=pubmedRes[3],
                                  Affiliation=pubmedRes[4], Abstract=pubmedRes[5],
                                  stringsAsFactors=FALSE)
             file.to <- paste(dir.to, "family.ref", sep="/")
             write.table(df.pub, file=file.to, sep="\t", quote=F, row.names=F)
           }
           
           git.info <- paste("add new family: ", symbol.low, sep="")
           system("git checkout master")
           system("git add -A")
           system(paste('git commit -m ', '"', git.info, '"', sep=""))
           
           js_string <- 'alert("New gene family successfully added!");'
           session$sendCustomMessage(type='jsCode', list(value = js_string))
           
         } else {
           js_string <- 'alert("Authorization Required!");'
           session$sendCustomMessage(type='jsCode', list(value = js_string))
         }
       })
     } else {NULL}
   })
   
})

