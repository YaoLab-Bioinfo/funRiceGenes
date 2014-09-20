
gene.lst <- read.table("geneInfo.table", head=T, 
                       as.is=T, sep="\t", quote="", comment="")
# for (j in 1:nrow(gene.lst)) {
for (j in 6:nrow(gene.lst)) {
  sym <- gene.lst$Symbol[j]
  sym <- gsub("\\|", ",", sym)
  msu <- gene.lst$MSU[j]
  rap <- gene.lst$RAPdb[j]
  path <- gene.lst$path[j]
  path <- paste("E:/git/RicENcode/RICENCODE", path, sep="/")
  
  ### msu
  msu <- unlist(strsplit(msu, split='|', fixed=TRUE))
  msu.new <- sapply(msu, function(x){
    if (x!="None") {
      y <- paste("http://rice.plantbiology.msu.edu/cgi-bin/ORF_infopage.cgi?orf=", 
                 x, sep="")
      y <- paste('[',x,']', '(', y, ')', sep="")
      return(y)
    } else {
      return(x)
    }
  })
  msu.new <- paste(unname(msu.new), sep="", collapse=",")
  
  ### rapdb
  rap <- unlist(strsplit(rap, split='|', fixed=TRUE))
  rap.new <- sapply(rap, function(x){
    if (x!="None") {
      y <- paste("http://rapdb.dna.affrc.go.jp/viewer/gbrowse_details/irgsp1?name=", 
                 x, sep="")
      y <- paste('[', x, ']', '(', y, ')', sep="")
      return(y)
    } else {
      return(x)
    }
  })
  rap.new <- paste(unname(rap.new), sep="", collapse=",")
  
  ### reference 
  ref.fl <- paste(path, "reference.info", sep="/")
  ref <- NULL
  if (file.exists(ref.fl)) {
    ref <- read.table(ref.fl, head=T, sep="\t", as.is=T, 
                      quote="", comment="")
    ref$Publication <- NULL
    ref$Affiliation <- NULL
    for (i in 1:nrow(ref)) {
      title <- ref$Title[i]
      ref$Title[i] <- paste("http://www.ncbi.nlm.nih.gov/pubmed", 
                            '?term=(', ref$Title[i],'%5BTitle%5D', 
                            ')', sep='')
      ref$Title[i] <- paste('[', title, ']', '(', ref$Title[i], ')', 
                             sep="")
    }
  }
  
  ### accession
  acc.fls <- list.files(path, patter="*acc-*", full=T)
  acc <- NULL
  if (length(acc.fls)>0) {
    acc.fls <- gsub(".+-", "", acc.fls)
    acc.fls <- sapply(acc.fls, function(x){
      y <- paste("http://www.ncbi.nlm.nih.gov/nuccore/", x, sep="")
      y <- paste('[', x, ']', '(', y, ')', sep="")
    })
    acc <- data.frame(Accession=acc.fls, stringsAsFactors=FALSE)
  }
  
  ### key
  key.fl <- paste(path, "Keyword.trait", sep="/")
  key <- NULL
  if (file.exists(key.fl)) {
    key <- read.table(key.fl, head=T, sep="\t", as.is=T, 
                      quote="", comment="")
    for (i in 1:nrow(key)) {
      title <- key$Title[i]
      key$Title[i] <- paste("http://www.ncbi.nlm.nih.gov/pubmed", 
                            '?term=(', key$Title[i],'%5BTitle%5D', 
                            ')', sep='')
      key$Title[i] <- paste('[', title, ']', '(', key$Title[i], ')', 
                            sep="")
    }
    key$Evidence <- gsub("\\|", ",", key$Evidence)
  }
  
  ### connection
  conne.fl <- paste(path, "Connection", sep="/")
  cone <- NULL
  if (file.exists(conne.fl)) {
    cone <- read.table(conne.fl, head=F, sep="\t", as.is=T, 
                      quote="", comment="")
    names(cone) <- c("Symbol1", "Symbol2", "Title", "Evidence")
    for (i in 1:nrow(cone)) {
      title <- cone$Title[i]
      cone$Title[i] <- paste("http://www.ncbi.nlm.nih.gov/pubmed", 
                            '?term=(', cone$Title[i],'%5BTitle%5D', 
                            ')', sep='')
      cone$Title[i] <- paste('[', title, ']', '(', cone$Title[i], ')', 
                            sep="")
    }
    cone$Evidence <- gsub("\\|", ",", cone$Evidence)
    cone$Symbol1 <- gsub("\\|", ",", cone$Symbol1)
    cone$Symbol2 <- gsub("\\|", ",", cone$Symbol2)
  }
  
  ###  output file
  md.cont <- ""
  md.cont[1] <- "---"
  md.cont[2] <- "layout: post"
  md.cont[3] <- paste('title: "', sym, '"', sep="")
  md.cont[4] <- 'description: ""'
  md.cont[5] <- "category: genes"
  if (!is.null(key)) {
    key.tmp <- paste(key$Keyword, sep="", collapse=", ")
    md.cont[6] <- paste('tags: [', key.tmp, ']', sep="")
  } else {
    md.cont[6] <- 'tags: '
  }
  md.cont[7] <- "---"
  md.cont[8] <- "{% include JB/setup %}"
  md.cont[9] <- ""
  md.cont[10] <- "## Information"
  md.sym <- paste("__Symbol__: ", sym, "  ", sep="")
  md.cont <- c(md.cont, md.sym)
  md.msu <- paste("__MSU__: ", msu.new, "  ", sep="")
  md.cont <- c(md.cont, md.msu)
  md.rap <- paste("__RAPdb__: ", rap.new, "  ", sep="")
  md.cont <- c(md.cont, md.rap)
  md.cont <- c(md.cont, "", "## Publication")
  if (!is.null(ref)) {
    for (i in 1:nrow(ref)) {
      md.ref <- paste(i, ". ", ref$Title[i], ", ", ref$Year[i], ", ",
                      ref$Journal[i], ".", sep="")
      md.cont <- c(md.cont, md.ref)
    }
  }
  md.cont <- c(md.cont, "", "## Genbank accession number")
  if (!is.null(acc)) {
    md.acc <- paste(acc$Accession, sep="", collapse=", ")
    md.cont <- c(md.cont, md.acc)
  }
  md.cont <- c(md.cont, "", "## Key message")
  if (!is.null(key)) {
    for (i in 1:nrow(key)) {
      md.key <- paste(i, ". __", key$Keyword[i], "__, ", key$Title[i],
                      ", ", key$Evidence[i], sep="")
      md.cont <- c(md.cont, md.key)
    }
  }
  md.cont <- c(md.cont, "", "## Connection")
  if (!is.null(cone)) {
    for (i in 1:nrow(cone)) {
      md.cone <- paste(i, ". __", cone$Symbol1[i], "__, __",cone$Symbol2[i],
                       "__, ", cone$Title[i], ", ", cone$Evidence[i], sep="")
      md.cont <- c(md.cont, md.cone)
    }
  }
  md.cont <- c(md.cont, "", "")
  
  gene.lst$Symbol[j] <- gsub("\\|", "~", gene.lst$Symbol[j])
  out.fl.name <- paste("2014-09-20-", gene.lst$Symbol[j], ".md", sep="")
  out.fl.name <- paste(path, out.fl.name, sep="/")
  
  writeLines(md.cont, con=out.fl.name)
}









