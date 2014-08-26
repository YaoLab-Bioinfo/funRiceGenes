
gene.info <- read.table("geneInfo.table", head=T, sep="\t", as.is=T)
gene.keyword <- read.table("geneKeyword.table", head=T, 
                           sep="\t", as.is=T, quote="", comment="")

gene.msu <- 1:nrow(gene.info)
names(gene.msu) <- gene.info$MSU
gene.msu <- gene.msu[names(gene.msu)!="None"]
gene.msu.new <- sapply(names(gene.msu), function(x) {
  x.name <- unlist(strsplit(x, split="\\|"))
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
  x.name <- unlist(strsplit(x, split="\\|"))
  if (length(x.name)==1) {
    return(gene.rap[x])
  } else {
    y <- rep(gene.rap[x], length(x.name))
    names(y) <- x.name
    return(y)
  }
})
gene.rap.final <- unlist(unname(gene.rap.new))


####  MSU
fetchInfoByMsu <- function(locus="") {
  locus <- gsub("^\\s+", "", locus)
  locus <- gsub("\\s+$", "", locus)
  locus.line <- gene.msu.final[locus]
  if (length(locus.line)==1) {
    path <- gene.info$path[locus.line]
    gene.fl <- paste(path, "gene.info", sep="/")
    if (file.exists(gene.fl)) {
      dat <- read.table(gene.fl, head=T, sep="\t", as.is=T)
      
      msu <- unlist(strsplit(dat$MSU, split="\\|"))
      msu.new <- sapply(msu, function(x){
        y <- paste("http://rice.plantbiology.msu.edu/cgi-bin/ORF_infopage.cgi?orf=", 
                   x, sep="")
        y <- paste('<a href="', y, '">', x, '</a>', sep="")
        y <- HTML(y)
      })
      msu.new <- paste(unname(msu.new), sep="", collapse="|")
      dat$MSU <- msu.new
      
      rap <- unlist(strsplit(dat$RAPdb, split="\\|"))
      rap.new <- sapply(rap, function(x){
        y <- paste("http://rapdb.dna.affrc.go.jp/viewer/gbrowse_details/irgsp1?name=", 
                   x, sep="")
        y <- paste('<a href="', y, '">', x, '</a>', sep="")
        y <- HTML(y)
      })
      rap.new <- paste(unname(rap.new), sep="", collapse="|")
      dat$RAPdb <- rap.new
      
      return(dat)
    }
  } else {
    return(NULL)
  }
}

fetchRefByMsu <- function(locus="") {
  locus <- gsub("^\\s+", "", locus)
  locus <- gsub("\\s+$", "", locus)
  locus.line <- gene.msu.final[locus]
  if (length(locus.line)==1) {
    path <- gene.info$path[locus.line]
    ref.fl <- paste(path, "reference.info", sep="/")
    if (file.exists(ref.fl)) {
      dat <- read.table(ref.fl, head=T, sep="\t", as.is=T, 
                        quote="", comment="")
      dat$Publication <- NULL
      dat$Abstract <- NULL
      for (i in 1:nrow(dat)) {
        title <- dat$Title[i]
        dat$Title[i] <- paste("http://www.ncbi.nlm.nih.gov/pubmed", 
                              '?term=(', dat$Title[i],'%5BTitle%5D', 
                              ')', sep='')
        dat$Title[i] <- paste('<a href="', dat$Title[i], '">', title, 
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
  locus <- gsub("\\s+$", "", locus)
  locus.line <- gene.msu.final[locus]
  if (length(locus.line)==1) {
    path <- gene.info$path[locus.line]
    acc.fls <- list.files(path, patter="*acc-*", full=T)
    if (length(acc.fls)>0) {
      acc.fls <- gsub(".+-", "", acc.fls)
      acc.fls <- sapply(acc.fls, function(x){
        y <- paste("http://www.ncbi.nlm.nih.gov/nuccore/", x, sep="")
        y <- paste('<a href="', y, '">', x, '</a>', sep="")
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

fetchTextByMsu <- function(locus="") {
  locus <- gsub("^\\s+", "", locus)
  locus <- gsub("\\s+$", "", locus)
  locus.line <- gene.msu.final[locus]
  if (length(locus.line)==1) {
    path <- gene.info$path[locus.line]
    text.fls <- list.files(path, patter="^pub.text.mining$", full=T)
    if (length(text.fls)>0) {
      text.con <- readLines(text.fls)
      clone <- sapply(seq(2, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("cloning:\\s+", "", y)
        return(y)
      })
      tdna <- sapply(seq(4, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("tdna:\\s+", "", y)
        return(y)
      })
      tos <- sapply(seq(5, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("tos:\\s+", "", y)
        return(y)
      })
      homo <- sapply(seq(6, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("homol:\\s+", "", y)
        return(y)
      })
      rnai <- sapply(seq(7, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("rnai:\\s+", "", y)
        return(y)
      })
      ove <- sapply(seq(8, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("ove:\\s+", "", y)
        return(y)
      })
      rt <- sapply(seq(9, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("rt:\\s+", "", y)
        return(y)
      })
      north <- sapply(seq(10, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("N:\\s+", "", y)
        return(y)
      })
      south <- sapply(seq(11, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("S:\\s+", "", y)
        return(y)
      })
      west <- sapply(seq(12, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("W:\\s+", "", y)
        return(y)
      })
      access <- sapply(seq(3, length(text.con), by=12), function(x){
        y <- text.con[x]
        return(y)
      })
      locus <- sapply(seq(1, length(text.con), by=12), function(x){
        y <- text.con[x]
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
  locus <- gsub("\\s+$", "", locus)
  locus.line <- gene.msu.final[locus]
  if (length(locus.line)==1) {
    path <- gene.info$path[locus.line]
    key.fl <- paste(path, "Keyword.trait", sep="/")
    if (file.exists(key.fl)) {
      dat <- read.table(key.fl, head=T, sep="\t", as.is=T, 
                        quote="", comment="")
      for (i in 1:nrow(dat)) {
        title <- dat$Title[i]
        dat$Title[i] <- paste("http://www.ncbi.nlm.nih.gov/pubmed", 
                              '?term=(', dat$Title[i],'%5BTitle%5D', 
                              ')', sep='')
        dat$Title[i] <- paste('<a href="', dat$Title[i], '">', title, 
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
  locus <- gsub("\\s+$", "", locus)
  locus.line <- gene.msu.final[locus]
  if (length(locus.line)==1) {
    path <- gene.info$path[locus.line]
    conne.fl <- paste(path, "Connection", sep="/")
    if (file.exists(conne.fl)) {
      dat <- read.table(conne.fl, head=F, sep="\t", as.is=T, 
                        quote="", comment="")
      names(dat) <- c("Symbol1", "Symbol2", "Title", "Evidence")
      for (i in 1:nrow(dat)) {
        title <- dat$Title[i]
        dat$Title[i] <- paste("http://www.ncbi.nlm.nih.gov/pubmed", 
                              '?term=(', dat$Title[i],'%5BTitle%5D', 
                              ')', sep='')
        dat$Title[i] <- paste('<a href="', dat$Title[i], '">', title, 
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
  locus <- gsub("\\s+$", "", locus)
  locus.line <- gene.rap.final[locus]
  if (length(locus.line)==1) {
    path <- gene.info$path[locus.line]
    gene.fl <- paste(path, "gene.info", sep="/")
    if (file.exists(gene.fl)) {
      dat <- read.table(gene.fl, head=T, sep="\t", as.is=T)
      msu <- unlist(strsplit(dat$MSU, split="\\|"))
      msu.new <- sapply(msu, function(x){
        y <- paste("http://rice.plantbiology.msu.edu/cgi-bin/ORF_infopage.cgi?orf=", 
                   x, sep="")
        y <- paste('<a href="', y, '">', x, '</a>', sep="")
        y <- HTML(y)
      })
      msu.new <- paste(unname(msu.new), sep="", collapse="|")
      dat$MSU <- msu.new
      
      rap <- unlist(strsplit(dat$RAPdb, split="\\|"))
      rap.new <- sapply(rap, function(x){
        y <- paste("http://rapdb.dna.affrc.go.jp/viewer/gbrowse_details/irgsp1?name=", 
                   x, sep="")
        y <- paste('<a href="', y, '">', x, '</a>', sep="")
        y <- HTML(y)
      })
      rap.new <- paste(unname(rap.new), sep="", collapse="|")
      dat$RAPdb <- rap.new
      
      return(dat)
    }
  } else {
    return(NULL)
  }
}

fetchRefByRap <- function(locus="") {
  locus <- gsub("^\\s+", "", locus)
  locus <- gsub("\\s+$", "", locus)
  locus.line <- gene.rap.final[locus]
  if (length(locus.line)==1) {
    path <- gene.info$path[locus.line]
    ref.fl <- paste(path, "reference.info", sep="/")
    if (file.exists(ref.fl)) {
      dat <- read.table(ref.fl, head=T, sep="\t", as.is=T, 
                        quote="", comment="")
      dat$Publication <- NULL
      dat$Abstract <- NULL
      for (i in 1:nrow(dat)) {
        title <- dat$Title[i]
        dat$Title[i] <- paste("http://www.ncbi.nlm.nih.gov/pubmed", 
                              '?term=(', dat$Title[i],'%5BTitle%5D', 
                              ')', sep='')
        dat$Title[i] <- paste('<a href="', dat$Title[i], '">', title, 
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
  locus <- gsub("\\s+$", "", locus)
  locus.line <- gene.rap.final[locus]
  if (length(locus.line)==1) {
    path <- gene.info$path[locus.line]
    acc.fls <- list.files(path, patter="*acc-*", full=T)
    if (length(acc.fls)>0) {
      acc.fls <- gsub(".+-", "", acc.fls)
      acc.fls <- sapply(acc.fls, function(x){
        y <- paste("http://www.ncbi.nlm.nih.gov/nuccore/", x, sep="")
        y <- paste('<a href="', y, '">', x, '</a>', sep="")
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

fetchTextByRap <- function(locus="") {
  locus <- gsub("^\\s+", "", locus)
  locus <- gsub("\\s+$", "", locus)
  locus.line <- gene.rap.final[locus]
  if (length(locus.line)==1) {
    path <- gene.info$path[locus.line]
    text.fls <- list.files(path, patter="^pub.text.mining$", full=T)
    if (length(text.fls)>0) {
      text.con <- readLines(text.fls)
      clone <- sapply(seq(2, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("cloning:\\s+", "", y)
        return(y)
      })
      tdna <- sapply(seq(4, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("tdna:\\s+", "", y)
        return(y)
      })
      tos <- sapply(seq(5, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("tos:\\s+", "", y)
        return(y)
      })
      homo <- sapply(seq(6, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("homol:\\s+", "", y)
        return(y)
      })
      rnai <- sapply(seq(7, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("rnai:\\s+", "", y)
        return(y)
      })
      ove <- sapply(seq(8, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("ove:\\s+", "", y)
        return(y)
      })
      rt <- sapply(seq(9, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("rt:\\s+", "", y)
        return(y)
      })
      north <- sapply(seq(10, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("N:\\s+", "", y)
        return(y)
      })
      south <- sapply(seq(11, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("S:\\s+", "", y)
        return(y)
      })
      west <- sapply(seq(12, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("W:\\s+", "", y)
        return(y)
      })
      access <- sapply(seq(3, length(text.con), by=12), function(x){
        y <- text.con[x]
        return(y)
      })
      locus <- sapply(seq(1, length(text.con), by=12), function(x){
        y <- text.con[x]
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
  locus <- gsub("\\s+$", "", locus)
  locus.line <- gene.rap.final[locus]
  if (length(locus.line)==1) {
    path <- gene.info$path[locus.line]
    key.fl <- paste(path, "Keyword.trait", sep="/")
    if (file.exists(key.fl)) {
      dat <- read.table(key.fl, head=T, sep="\t", as.is=T, 
                        quote="", comment="")
      for (i in 1:nrow(dat)) {
        title <- dat$Title[i]
        dat$Title[i] <- paste("http://www.ncbi.nlm.nih.gov/pubmed", 
                              '?term=(', dat$Title[i],'%5BTitle%5D', 
                              ')', sep='')
        dat$Title[i] <- paste('<a href="', dat$Title[i], '">', title, 
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
  locus <- gsub("\\s+$", "", locus)
  locus.line <- gene.rap.final[locus]
  if (length(locus.line)==1) {
    path <- gene.info$path[locus.line]
    conne.fl <- paste(path, "Connection", sep="/")
    if (file.exists(conne.fl)) {
      dat <- read.table(conne.fl, head=F, sep="\t", as.is=T, 
                        quote="", comment="")
      names(dat) <- c("Symbol1", "Symbol2", "Title", "Evidence")
      for (i in 1:nrow(dat)) {
        title <- dat$Title[i]
        dat$Title[i] <- paste("http://www.ncbi.nlm.nih.gov/pubmed", 
                              '?term=(', dat$Title[i],'%5BTitle%5D', 
                              ')', sep='')
        dat$Title[i] <- paste('<a href="', dat$Title[i], '">', title, 
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
    sym.line <- unlist(strsplit(sym.line, split="\\|"))
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
  symbol <- gsub("\\s+$", "", symbol)
  locus.line <- findDirBySym(tolower(symbol))
  if (length(locus.line)==1) {
    path <- gene.info$path[locus.line]
    gene.fl <- paste(path, "gene.info", sep="/")
    if (file.exists(gene.fl)) {
      dat <- read.table(gene.fl, head=T, sep="\t", as.is=T)
      
      msu <- unlist(strsplit(dat$MSU, split="\\|"))
      msu.new <- sapply(msu, function(x){
        y <- paste("http://rice.plantbiology.msu.edu/cgi-bin/ORF_infopage.cgi?orf=", 
                   x, sep="")
        y <- paste('<a href="', y, '">', x, '</a>', sep="")
        y <- HTML(y)
      })
      msu.new <- paste(unname(msu.new), sep="", collapse="|")
      dat$MSU <- msu.new
      
      rap <- unlist(strsplit(dat$RAPdb, split="\\|"))
      rap.new <- sapply(rap, function(x){
        y <- paste("http://rapdb.dna.affrc.go.jp/viewer/gbrowse_details/irgsp1?name=", 
                   x, sep="")
        y <- paste('<a href="', y, '">', x, '</a>', sep="")
        y <- HTML(y)
      })
      rap.new <- paste(unname(rap.new), sep="", collapse="|")
      dat$RAPdb <- rap.new
      
      return(dat)
    }
  } else if (length(locus.line)>1) {
    dat <- gene.info[locus.line, ]
    dat$path <- NULL
    for (i in 1:nrow(dat)) {    
      msu <- unlist(strsplit(dat$MSU[i], split="\\|"))
      msu.new <- sapply(msu, function(x){
        y <- paste("http://rice.plantbiology.msu.edu/cgi-bin/ORF_infopage.cgi?orf=", 
                   x, sep="")
        y <- paste('<a href="', y, '">', x, '</a>', sep="")
        y <- HTML(y)
      })
      msu.new <- paste(unname(msu.new), sep="", collapse="|")
      dat$MSU[i] <- msu.new
      
      rap <- unlist(strsplit(dat$RAPdb[i], split="\\|"))
      rap.new <- sapply(rap, function(x){
        y <- paste("http://rapdb.dna.affrc.go.jp/viewer/gbrowse_details/irgsp1?name=", 
                   x, sep="")
        y <- paste('<a href="', y, '">', x, '</a>', sep="")
        y <- HTML(y)
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
  symbol <- gsub("\\s+$", "", symbol)
  locus.line <- findDirBySym(tolower(symbol))
  if (length(locus.line)==1) {
    path <- gene.info$path[locus.line]
    ref.fl <- paste(path, "reference.info", sep="/")
    if (file.exists(ref.fl)) {
      dat <- read.table(ref.fl, head=T, sep="\t", as.is=T, 
                        quote="", comment="")
      dat$Publication <- NULL
      dat$Abstract <- NULL
      for (i in 1:nrow(dat)) {
        title <- dat$Title[i]
        dat$Title[i] <- paste("http://www.ncbi.nlm.nih.gov/pubmed", 
                              '?term=(', dat$Title[i],'%5BTitle%5D', 
                              ')', sep='')
        dat$Title[i] <- paste('<a href="', dat$Title[i], '">', title, 
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
  symbol <- gsub("\\s+$", "", symbol)
  locus.line <- findDirBySym(tolower(symbol))
  if (length(locus.line)==1) {
    path <- gene.info$path[locus.line]
    acc.fls <- list.files(path, patter="*acc-*", full=T)
    if (length(acc.fls)>0) {
      acc.fls <- gsub(".+-", "", acc.fls)
      acc.fls <- sapply(acc.fls, function(x){
        y <- paste("http://www.ncbi.nlm.nih.gov/nuccore/", x, sep="")
        y <- paste('<a href="', y, '">', x, '</a>', sep="")
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

fetchTextBySym <- function(symbol="") {
  symbol <- gsub("^\\s+", "", symbol)
  symbol <- gsub("\\s+$", "", symbol)
  locus.line <- findDirBySym(tolower(symbol))
  if (length(locus.line)==1) {
    path <- gene.info$path[locus.line]
    text.fls <- list.files(path, patter="^pub.text.mining$", full=T)
    if (length(text.fls)>0) {
      text.con <- readLines(text.fls)
      clone <- sapply(seq(2, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("cloning:\\s+", "", y)
        return(y)
      })
      tdna <- sapply(seq(4, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("tdna:\\s+", "", y)
        return(y)
      })
      tos <- sapply(seq(5, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("tos:\\s+", "", y)
        return(y)
      })
      homo <- sapply(seq(6, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("homol:\\s+", "", y)
        return(y)
      })
      rnai <- sapply(seq(7, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("rnai:\\s+", "", y)
        return(y)
      })
      ove <- sapply(seq(8, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("ove:\\s+", "", y)
        return(y)
      })
      rt <- sapply(seq(9, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("rt:\\s+", "", y)
        return(y)
      })
      north <- sapply(seq(10, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("N:\\s+", "", y)
        return(y)
      })
      south <- sapply(seq(11, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("S:\\s+", "", y)
        return(y)
      })
      west <- sapply(seq(12, length(text.con), by=12), function(x){
        y <- text.con[x]
        y <- gsub("W:\\s+", "", y)
        return(y)
      })
      access <- sapply(seq(3, length(text.con), by=12), function(x){
        y <- text.con[x]
        return(y)
      })
      locus <- sapply(seq(1, length(text.con), by=12), function(x){
        y <- text.con[x]
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
  symbol <- gsub("\\s+$", "", symbol)
  locus.line <- findDirBySym(tolower(symbol))
  if (length(locus.line)==1) {
    path <- gene.info$path[locus.line]
    key.fl <- paste(path, "Keyword.trait", sep="/")
    if (file.exists(key.fl)) {
      dat <- read.table(key.fl, head=T, sep="\t", as.is=T, 
                        quote="", comment="")
      for (i in 1:nrow(dat)) {
        title <- dat$Title[i]
        dat$Title[i] <- paste("http://www.ncbi.nlm.nih.gov/pubmed", 
                              '?term=(', dat$Title[i],'%5BTitle%5D', 
                              ')', sep='')
        dat$Title[i] <- paste('<a href="', dat$Title[i], '">', title, 
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
  symbol <- gsub("\\s+$", "", symbol)
  locus.line <- findDirBySym(tolower(symbol))
  if (length(locus.line)==1) {
    path <- gene.info$path[locus.line]
    conne.fl <- paste(path, "Connection", sep="/")
    if (file.exists(conne.fl)) {
      dat <- read.table(conne.fl, head=F, sep="\t", as.is=T, 
                        quote="", comment="")
      names(dat) <- c("Symbol1", "Symbol2", "Title", "Evidence")
      for (i in 1:nrow(dat)) {
        title <- dat$Title[i]
        dat$Title[i] <- paste("http://www.ncbi.nlm.nih.gov/pubmed", 
                              '?term=(', dat$Title[i],'%5BTitle%5D', 
                              ')', sep='')
        dat$Title[i] <- paste('<a href="', dat$Title[i], '">', title, 
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
  keyword <- gsub("\\s+$", "", keyword)
  all.key <- unique(gene.keyword$Keyword)
  if (tolower(keyword) %in% all.key) {
    dat <- gene.keyword[gene.keyword$Keyword==tolower(keyword), ]
    dat$Keyword <- NULL
    for (i in 1:nrow(dat)) {    
      msu <- unlist(strsplit(dat$MSU[i], split="\\|"))
      msu.new <- sapply(msu, function(x){
        y <- paste("http://rice.plantbiology.msu.edu/cgi-bin/ORF_infopage.cgi?orf=", 
                   x, sep="")
        y <- paste('<a href="', y, '">', x, '</a>', sep="")
        y <- HTML(y)
      })
      msu.new <- paste(unname(msu.new), sep="", collapse="|")
      dat$MSU[i] <- msu.new
      
      rap <- unlist(strsplit(dat$RAPdb[i], split="\\|"))
      rap.new <- sapply(rap, function(x){
        y <- paste("http://rapdb.dna.affrc.go.jp/viewer/gbrowse_details/irgsp1?name=", 
                   x, sep="")
        y <- paste('<a href="', y, '">', x, '</a>', sep="")
        y <- HTML(y)
      })
      rap.new <- paste(unname(rap.new), sep="", collapse="|")
      dat$RAPdb[i] <- rap.new
      
      title <- dat$Title[i]
      dat$Title[i] <- paste("http://www.ncbi.nlm.nih.gov/pubmed", 
                            '?term=(', dat$Title[i],'%5BTitle%5D', 
                            ')', sep='')
      dat$Title[i] <- paste('<a href="', dat$Title[i], '">', title, 
                            '</a>', sep="")
      dat$Title[i] <- HTML(dat$Title[i])
    }
    
    return(dat)
  } else {
    return(NULL)
  }
}

write.gene <- function(df) {
  symbol <- df$Symbol
  symbol.low <- tolower(symbol)
  if (grepl("^os", symbol.low)) {
    symbol.low <- gsub("^os", "", symbol.low)
    symbol.head <- substr(symbol.low, 1, 1)
    if (symbol.head %in% letters[1:24]) {
      dir.to <- paste("data/Gene/Abstract/OS", 
                      toupper(symbol.head), sep="/")
    } else {
      dir.to <- "data/Gene/Abstract/OS/0-9"
    }
  } else {
    symbol.head <- substr(symbol.low, 1, 1)
    if (symbol.head %in% letters[1:24]) {
      dir.to <- paste("data/Gene/Abstract", 
                      toupper(symbol.head), sep="/")
    } else {
      dir.to <- "data/Gene/Abstract/0-9"
    }
  }
  
  symbol <- gsub("\\|", "~", symbol)
  dir.to <- paste(dir.to, symbol, sep="/")
  if (!file.exists(dir.to)) {
    dir.create(dir.to)
    out.fl <- paste(dir.to, "gene.info", sep="/")
    write.table(df, file=out.fl, sep="\t", quote=F, row.names=F)
  }
}

write.pub <- function(df) {
  symbol <- df$Symbol
  symbol <- gsub("^\\s+", "", symbol)
  symbol <- gsub("\\s+$", "", symbol)
  locus.line <- findDirBySym(tolower(symbol))
  df$Symbol <- gene.info$Symbol[locus.line]
  if (length(locus.line)==1) {
    path <- gene.info$path[locus.line]
    out.fl <- paste(path, "reference.info", sep="/")
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
  }
}

write.acc <- function(df) {
  symbol <- df$Symbol
  symbol <- gsub("^\\s+", "", symbol)
  symbol <- gsub("\\s+$", "", symbol)
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
  symbol <- df$Symbol
  symbol <- gsub("^\\s+", "", symbol)
  symbol <- gsub("\\s+$", "", symbol)
  locus.line <- findDirBySym(tolower(symbol))
  df$Symbol <- gene.info$Symbol[locus.line]
  if (length(locus.line)==1) {
    path <- gene.info$path[locus.line]
    out.fl <- paste(path, "Keyword.trait", sep="/")
    if (file.exists(out.fl)) {
      df.tmp <- read.table(out.fl, sep="\t", quote="", head=T, as.is=T, comment="")
      df.new <- unique(rbind(df.tmp, df))
      write.table(df.new, file=out.fl, sep="\t", quote=F, row.names=F)
    } else {
      write.table(df, file=out.fl, sep="\t", quote=F, row.names=F)
    }
  }
}

write.con <- function(df) {
  symbol.1 <- df$Symbol1
  symbol.2 <- df$Symbol2
  locus.line.1 <- findDirBySym(tolower(symbol.1))
  locus.line.2 <- findDirBySym(tolower(symbol.2))
  df$Symbol1 <- gene.info$Symbol[locus.line.1]
  df$Symbol2 <- gene.info$Symbol[locus.line.2]
  
  if (length(locus.line.1)==1 && length(locus.line.2)==1) {
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
      df$v1n <- pmin(df$V1, df$V2)
      df$v2n <- pmax(df$V1, df$V2)
      df <- df[, c("v1n", "v2n", "V3", "V4")]
      df <- unique(df)
      write.table(df, file=out.fl.2, sep="\t", quote=F, row.names=F,
                  col.names=F)
    }
  }
}

gene.edit <- function(df){
  gene.tar <- gene.info[gene.info$Symbol==df$oldsym&gene.info$RAPdb==df$oldrap&
                          gene.info$MSU==df$oldmsu,]
  df.dat <- df[, c("newsym", "newrap", "newmsu")]
  names(df.dat) <- c("Symbol", "RAPdb", "MSU")
  tar.fl <- paste(gene.tar$path, "gene.info", sep="/")
  write.table(df.dat, file=tar.fl, sep="\t", quote=F, row.names=F)
  
  key.fl <- paste(gene.tar$path, "Keyword.trait", sep="/")
  key.con <- read.table(key.fl, sep="\t", head=T, as.is=T, quote="", comment="")
  key.con$Symbol <- df$newsym
  write.table(key.con, file=key.fl, sep="\t", quote=F, row.names=F)
  
  cone.fl <- paste(gene.tar$path, "Connection", sep="/")
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

updateGeneInfo <- function() {
  all.gene.fls <- list.files("data/Gene",
                             patter="^gene.info$", full=T, recur=T)
  all.gene.lst <- lapply(all.gene.fls, function(x){
    dat <- read.table(x, head=T, sep="\t", as.is=T, quote="", comment="")
    dir.cwd <- dirname(x)
    dat$path <- dir.cwd
    return(dat)
  })
  all.gene.df <- do.call(rbind, all.gene.lst)
  write.table(all.gene.df, file="geneInfo.table",
              sep="\t", quote=F, row.names=F)
}

updateKeyword <- function() {
  all.key.fls <- list.files("data/Gene", patter="^Keyword.trait$",
                            full=T, recur=T)
  all.key.lst <- lapply(all.key.fls, function(x){
    cwd.dir <- dirname(x)
    cwd.gene <- paste(cwd.dir, "gene.info", sep="/")
    gene.dat <- read.table(cwd.gene, head=T, sep="\t", as.is=T, quote="", comment="")
    key.dat <- read.table(x, head=T, sep="\t", as.is=T, quote="", comment="")
    key.dat$RAPdb <- gene.dat$RAPdb
    key.dat$MSU <- gene.dat$MSU
    key.dat <- key.dat[, c("Symbol", "RAPdb",  "MSU",	"Keyword",	"Title")]
    return(key.dat)
  })
  all.key.df <- do.call(rbind, all.key.lst)
  write.table(all.key.df, file="geneKeyword.table", sep="\t", quote=F, row.names=F)
}


#### Shiny
shinyServer(function(input, output) {
  
  output$mytable1 = renderDataTable({
    fetchInfoByMsu(input$msu)
  }, options = list(aLengthMenu = c(1, 2), bFilter = FALSE,
                    bPaginate = FALSE)
  )
  
  output$mytable2 = renderDataTable({
    fetchRefByMsu(input$msu)
  }, options = list(aLengthMenu = c(2, 3, 4), iDisplayLength = 2,
                    bFilter = FALSE, bAutoWidth = FALSE)
  )
  
  output$mytable3 = renderDataTable({
    fetchAccByMsu(input$msu)
  }, options = list(aLengthMenu = c(1, 2), bFilter = FALSE,
                    iDisplayLength = 2))
  
  output$mytable4 = renderDataTable({
    fetchTextByMsu(input$msu)
  }, options = list(aLengthMenu = c(1, 2), bFilter = FALSE,
                    iDisplayLength = 1, bAutoWidth = FALSE))
  
  output$mytable5 = renderDataTable({
    fetchKeyByMsu(input$msu)
  }, options = list(aLengthMenu = c(1, 2), bFilter = FALSE,
                    iDisplayLength = 1, bAutoWidth = FALSE))
  
  output$mytable6 = renderDataTable({
    fetchConneByMsu(input$msu)
  }, options = list(aLengthMenu = c(1, 2), bFilter = FALSE,
                    iDisplayLength = 1, bAutoWidth = FALSE))
  
  output$mytable7 = renderDataTable({
    fetchInfoByRap(input$rap)
  }, options = list(aLengthMenu = c(1, 2), bFilter = FALSE,
                    bPaginate = FALSE))
  
  output$mytable8 = renderDataTable({
    fetchRefByRap(input$rap)
  }, options = list(aLengthMenu = c(2, 3, 4), iDisplayLength = 2,
                    bFilter = FALSE, bAutoWidth = FALSE)
  )
  
  output$mytable9 = renderDataTable({
    fetchAccByRap(input$rap)
  }, options = list(aLengthMenu = c(1, 2), bFilter = FALSE,
                    iDisplayLength = 2))
  
  output$mytable10 = renderDataTable({
    fetchTextByRap(input$rap)
  }, options = list(aLengthMenu = c(1, 2), bFilter = FALSE,
                    iDisplayLength = 1, bAutoWidth = FALSE))
  
  output$mytable11 = renderDataTable({
    fetchKeyByRap(input$rap)
  }, options = list(aLengthMenu = c(1, 2), bFilter = FALSE,
                    iDisplayLength = 1, bAutoWidth = FALSE))
  
  output$mytable12 = renderDataTable({
    fetchConneByRap(input$rap)
  }, options = list(aLengthMenu = c(1, 2), bFilter = FALSE,
                    iDisplayLength = 1, bAutoWidth = FALSE))
  
  output$mytable13 = renderDataTable({
    fetchInfoBySym(input$symbol)
  }, options = list(aLengthMenu = c(1, 2), bFilter = FALSE,
                    bPaginate = FALSE))
  
  output$mytable14 = renderDataTable({
    fetchRefBySym(input$symbol)
  }, options = list(aLengthMenu = c(2, 3, 4), iDisplayLength = 2,
                    bFilter = FALSE, bAutoWidth = FALSE)
  )
  
  output$mytable15 = renderDataTable({
    fetchAccBySym(input$symbol)
  }, options = list(aLengthMenu = c(1, 2), bFilter = FALSE,
                    iDisplayLength = 2))
  
  output$mytable16 = renderDataTable({
    fetchTextBySym(input$symbol)
  }, options = list(aLengthMenu = c(1, 2), bFilter = FALSE,
                    iDisplayLength = 1, bAutoWidth = FALSE))
  
  output$mytable17 = renderDataTable({
    fetchKeyBySym(input$symbol)
  }, options = list(aLengthMenu = c(1, 2), bFilter = FALSE,
                    iDisplayLength = 1, bAutoWidth = FALSE))
  
  output$mytable18 = renderDataTable({
    fetchConneBySym(input$symbol)
  }, options = list(aLengthMenu = c(1, 2), bFilter = FALSE,
                    iDisplayLength = 1, bAutoWidth = FALSE))
  
  output$mytable19 = renderDataTable({
    fetchInfoByKey(input$keyword)
  }, options = list(aLengthMenu = c(2, 4, 6), bFilter = FALSE,
                    iDisplayLength = 2, bAutoWidth = FALSE))
  
  observe({
    if (input$submit1>0) {
      isolate({
        df.gene <- data.frame(Symbol=input$symsub1, MSU=input$msusub1, RAPdb=input$rapsub1,
                              stringsAsFactors=FALSE)
        write.gene(df.gene)
      })
    } else {NULL}
  })
  
  observe({
    if (input$submit2>0) {
      isolate({
        df.pub <- data.frame(Symbol=input$symsub2, Title=input$tilsub2, Year=input$yearsub2,
                             Journal=input$jousub2, Affiliation=input$afisub2, Abstract=input$abssub2,
                             stringsAsFactors=FALSE)
        write.pub(df.pub)
      })
    } else {NULL}
  })
  
  observe({
    if (input$submit3>0) {
      isolate({
        df.acc <- data.frame(Symbol=input$symsub3, Accession=input$accsub3, 
                             stringsAsFactors=FALSE)
        write.acc(df.acc)
      })
    } else {NULL}
  })
  
  observe({
    if (input$submit4>0) {
      isolate({
        df.key <- data.frame(Symbol=input$symsub4, Keyword=input$keysub4, 
                             Title=input$tilsub4, Evidence=input$evisub4,
                             stringsAsFactors=FALSE)
        write.key(df.key)
      })
    } else {NULL}
  })
  
  observe({
    if (input$submit5>0) {
      isolate({
        df.con <- data.frame(Symbol1=input$symsub5, Symbol2=input$sym2sub5, 
                             Title=input$tilsub5, Evidence=input$evisub5,
                             stringsAsFactors=FALSE)
        write.con(df.con)
      })
    } else {NULL}
  })
  
  observe({
    if (input$submit6>0) {
      isolate({
        df.con <- data.frame(oldsym=input$oldsym, newsym=input$newsym, 
                             oldmsu=input$oldmsu, newmsu=input$newmsu,
                             oldrap=input$oldrap, newrap=input$newrap,
                             stringsAsFactors=FALSE)
        gene.edit(df.con)
      })
    } else {NULL}
  })
  
  observe({
    if (input$submit7>0) {
      updateGeneInfo()
      updateKeyword()
    } else {NULL}
  })
  
})









