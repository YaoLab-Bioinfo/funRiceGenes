
setwd("E:/GIT/RICENCODE")
unlink("E:/GIT/funRiceGenes.github.io/_posts", recur=T, force=T)
gene.lst <- read.table("geneInfo.table", head=T, 
                       as.is=T, sep="\t", quote="", comment="")

dat.gene <- gene.lst
dat.gene$path <- NULL
write.table(dat.gene, file="E:/GIT/funRiceGenes.github.io/geneInfo.table.txt", 
            quote=F, sep="\t", row.names=F)

dat.fam <- read.table("famInfo.table", head=T, 
                      as.is=T, sep="\t", quote="", comment="")
dat.fam$Name <- basename(dat.fam$path)
dat.fam$path <- NULL
write.table(dat.fam, file="E:/GIT/funRiceGenes.github.io/famInfo.table.txt", 
            quote=F, sep="\t", row.names=F)

file.copy(from="reference.table", to="E:/GIT/funRiceGenes.github.io/reference.table.txt", overwrite = TRUE)

dat.key <- read.table("geneKeyword.table", head=T, as.is = T, sep="\t", quote="", comment="")
dat.key$path <- NULL
write.table(dat.key, file="E:/GIT/funRiceGenes.github.io/geneKeyword.table.txt", 
            quote=F, sep="\t", row.names=F)


for (j in 1:nrow(gene.lst)) {
  sym <- gene.lst$Symbol[j]
  sym <- gsub("\\|", ",", sym)
  msu <- gene.lst$MSU[j]
  rap <- gene.lst$RAPdb[j]
  path <- gene.lst$path[j]
  path <- paste("E:/GIT/RICENCODE", path, sep="/")
  
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
      title <- gsub("\\)", "", title)
      title <- gsub("\\(", "", title)
      
      ref$Title[i] <- gsub("\\)", "", ref$Title[i])
      ref$Title[i] <- gsub("\\(", "", ref$Title[i])
      
      ref$Title[i] <- paste("http://www.ncbi.nlm.nih.gov/pubmed", 
                            '?term=', ref$Title[i],'%5BTitle%5D', 
                            sep='')
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
  
  acc <- unique(acc)
  
  ### key
  key.fl <- paste(path, "Keyword.trait", sep="/")
  key <- NULL
  if (file.exists(key.fl)) {
    key <- read.table(key.fl, head=T, sep="\t", as.is=T, 
                      quote="", comment="")
    key$Evidence <- gsub("\\|", ",", key$Evidence)
    key <- key[, c("Keyword", "Evidence"), drop=F]
    key$Evidence <- gsub("^\\s+", "", key$Evidence)
    key$Evidence <- gsub("\\s+$", "", key$Evidence)
    key <- unique(key)
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
      title <- gsub("\\)", "", title)
      title <- gsub("\\(", "", title)
      cone$Title[i] <- gsub("\\)", "", cone$Title[i])
      cone$Title[i] <- gsub("\\(", "", cone$Title[i])
      
      cone$Title[i] <- paste("http://www.ncbi.nlm.nih.gov/pubmed", 
                            '?term=', cone$Title[i],'%5BTitle%5D', 
                            sep='')
      cone$Title[i] <- paste('[', title, ']', '(', cone$Title[i], ')', 
                            sep="")
    }
    cone$Evidence <- gsub("\\|", ",", cone$Evidence)
    cone$Symbol1 <- gsub("\\|", "~", cone$Symbol1)
    cone$Symbol2 <- gsub("\\|", "~", cone$Symbol2)
  }
  
  ### figures
  exp.fig.fl <- list.files(path, patter=".*\\.exp\\..+", full=T)
  pheno.fig.fl <- list.files(path, patter=".*\\.pheno\\..+", full=T)
  
  ###  output file
  md.cont <- ""
  md.cont[1] <- "---"
  md.cont[2] <- "layout: post"
  md.cont[3] <- paste('title: "', sym, '"', sep="")
  md.cont[4] <- 'description: ""'
  md.cont[5] <- "category: genes"
  if (!is.null(key)) {
    key.tmp <- paste(unique(key$Keyword), sep="", collapse=", ")
    md.cont[6] <- paste('tags: [', key.tmp, ']', sep="")
  } else {
    md.cont[6] <- 'tags: '
  }
  md.cont[7] <- "---"
  md.cont[8] <- ""
  md.cont[9] <- "* **Information**  "
  md.sym <- paste("    + Symbol: ", sym, "  ", sep="")
  md.cont <- c(md.cont, md.sym)
  md.msu <- paste("    + MSU: ", msu.new, "  ", sep="")
  md.cont <- c(md.cont, md.msu)
  md.rap <- paste("    + RAPdb: ", rap.new, "  ", sep="")
  md.cont <- c(md.cont, md.rap)
  
  md.cont <- c(md.cont, "", "* **Publication**  ")
  if (!is.null(ref)) {
    for (i in 1:nrow(ref)) {
      md.ref <- paste("    + ", ref$Title[i], ", ", ref$Year[i], ", ",
                      ref$Journal[i], ".", sep="")
      md.cont <- c(md.cont, md.ref)
    }
  }

  md.cont <- c(md.cont, "", "* **Genbank accession number**  ")
  if (!is.null(acc)) {
    for (i in 1:nrow(acc)) {
      md.acc <- paste("    + ", acc$Accession, sep="")
      md.cont <- c(md.cont, md.acc)
    }
  }

  md.cont <- c(md.cont, "", "* **Key message**  ")
  if (!is.null(key)) {
    key$Keyword <- NULL
    key <- unique(key)
    for (i in 1:nrow(key)) {
      md.key <- paste("    + ", key$Evidence[i], sep="")
      md.cont <- c(md.cont, md.key)
    }
  }

  md.cont <- c(md.cont, "", "* **Connection**  ")
  if (!is.null(cone)) {
    for (i in 1:nrow(cone)) {
      md.cone <- paste("    + __", cone$Symbol1[i], "__, __",cone$Symbol2[i],
                       "__, ", cone$Title[i], ", ", cone$Evidence[i], sep="")
      md.cont <- c(md.cont, md.cone)
    }
  }
  
  md.cont <- c(md.cont, "", "[//]: # * **Key figures**  ")
  if (length(pheno.fig.fl)==1) {
    md.cont <- c(md.cont, paste('[//]: # <img src="http://funRiceGenes.github.io/images/', 
                                basename(pheno.fig.fl), '" alt="phenotype"  style="width: 600px;"/>', sep=""))
    md.cont <- c(md.cont, "")
  }
  
  if (length(exp.fig.fl)==1) {
    md.cont <- c(md.cont, paste('[//]: # <img src="http://funRiceGenes.github.io/images/', 
                                basename(exp.fig.fl), '" alt="expression"  style="width: 600px;"/>', sep=""))
  }
  
  md.cont <- c(md.cont, "", "")
  
  gene.lst$Symbol[j] <- gsub("\\|", "~", gene.lst$Symbol[j])
  out.fl.name <- paste("2015-01-20-", gene.lst$Symbol[j], ".md", sep="")
  if (grepl("^os", gene.lst$Symbol[j], ignore.case=T)) {
    tmp.tag <- substr(gsub("^os", "", gene.lst$Symbol[j], ignore.case=T), 1,1)
	tmp.tag <- toupper(tmp.tag)
	if (tmp.tag %in% LETTERS[1:26]) {
	  path <- paste("E:/GIT/funRiceGenes.github.io/_posts/OS/", tmp.tag, sep="")
	} else {
	  path <- "E:/GIT/funRiceGenes.github.io/_posts/OS/0-9"
	}
  } else {
    tmp.tag <- substr(gene.lst$Symbol[j], 1,1)
	tmp.tag <- toupper(tmp.tag)
	if (tmp.tag %in% LETTERS[1:26]) {
	  path <- paste("E:/GIT/funRiceGenes.github.io/_posts/", tmp.tag, sep="")
	} else {
	  path <- "E:/GIT/funRiceGenes.github.io/_posts/0-9"
	}
  }
	  
  out.fl.name <- paste(path, out.fl.name, sep="/")
  if (length(exp.fig.fl)==1) {
#    file.copy(from=exp.fig.fl, to="E:/GIT/funRiceGenes.github.io/images")
  }
  
  if (length(pheno.fig.fl)==1) {
#    file.copy(from=pheno.fig.fl, to="E:/GIT/funRiceGenes.github.io/images")
  }
  
  dir.name <- dirname(out.fl.name)
  if (!file.exists(dir.name)) {
    dir.create(dir.name, recursive=TRUE)
  }
  writeLines(md.cont, con=out.fl.name)
}


####### gene family
fam.lst <- read.table("famInfo.table", head=T, 
                       as.is=T, sep="\t", quote="", comment="")
fam.path <- unique(fam.lst$path)
for (j in 1:length(fam.path)) {
  path <- fam.path[j]
  path <- paste("E:/GIT/RICENCODE", path, sep="/")
  
  info.fl <- paste(path, "family.info", sep="/")
  fam.info <- NULL
  if (file.exists(info.fl)) {
    fam.info <- read.table(info.fl, head=T, sep="\t", as.is=T, 
                      quote="", comment="")
    
    for (i in 1:nrow(fam.info)) {
      ### msu
      msu <- fam.info$MSU[i]
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
      fam.info$MSU[i] <- msu.new
      
      ### rapdb
      rap <- fam.info$RAPdb[i]
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
      fam.info$RAPdb[i] <- rap.new
    }
  }
  
  ### reference 
  ref.fl <- paste(path, "family.ref", sep="/")
  ref <- NULL
  if (file.exists(ref.fl)) {
    ref <- read.table(ref.fl, head=T, sep="\t", as.is=T, 
                      quote="", comment="")
    ref$Publication <- NULL
    ref$Affiliation <- NULL
    for (i in 1:nrow(ref)) {
      title <- ref$Title[i]
      title <- gsub("\\)", "", title)
      title <- gsub("\\(", "", title)
      ref$Title[i] <- gsub("\\)", "", ref$Title[i])
      ref$Title[i] <- gsub("\\(", "", ref$Title[i])
      
      ref$Title[i] <- paste("http://www.ncbi.nlm.nih.gov/pubmed", 
                            '?term=', ref$Title[i],'%5BTitle%5D', 
                            sep='')
      ref$Title[i] <- paste('[', title, ']', '(', ref$Title[i], ')', 
                            sep="")
    }
  }
  
  ###  output file
  name <- basename(path)
  
  md.cont <- ""
  md.cont[1] <- "---"
  md.cont[2] <- "layout: post"
  md.cont[3] <- paste('title: "', name, '"', sep="")
  md.cont[4] <- 'description: ""'
  md.cont[5] <- "category: gene family"
  md.cont[6] <- "---"
  md.cont[7] <- ""
  md.cont[8] <- "* **Information**  "
  if (!is.null(fam.info)) {
    for (i in 1:nrow(fam.info)) {
      md.info <- paste("    + ", gsub("\\|", ",", fam.info$Symbol[i]), ", ", fam.info$MSU[i], ", ",
                      fam.info$RAPdb[i], ".", sep="")
      md.cont <- c(md.cont, md.info)
    }
  }
  
  md.cont <- c(md.cont, "", "* **Publication**  ")
  if (!is.null(ref)) {
    for (i in 1:nrow(ref)) {
      md.ref <- paste( "    + ", ref$Title[i], ", ", ref$Year[i], ", ",
                      ref$Journal[i], ".", sep="")
      md.cont <- c(md.cont, md.ref)
    }
  }
  
  md.cont <- c(md.cont, "", "")
  
  out.fl.name <- paste("2015-01-20-", name, ".md", sep="")
  if (grepl("^os", name, ignore.case=T)) {
    tmp.tag <- substr(gsub("^os", "", name, ignore.case=T), 1,1)
    tmp.tag <- toupper(tmp.tag)
    if (tmp.tag %in% LETTERS[1:26]) {
      path <- paste("E:/GIT/funRiceGenes.github.io/_posts/FAM/OS/", tmp.tag, sep="")
    } else {
      path <- "E:/GIT/funRiceGenes.github.io/_posts/FAM/OS/0-9"
    }
  } else {
    tmp.tag <- substr(name, 1,1)
    tmp.tag <- toupper(tmp.tag)
    if (tmp.tag %in% LETTERS[1:26]) {
      path <- paste("E:/GIT/funRiceGenes.github.io/_posts/FAM/", tmp.tag, sep="")
    } else {
      path <- "E:/GIT/funRiceGenes.github.io/_posts/FAM/0-9"
    }
  }
  
  out.fl.name <- paste(path, out.fl.name, sep="/")
  dir.name <- dirname(out.fl.name)
  if (!file.exists(dir.name)) {
    dir.create(dir.name, recursive=TRUE)
  }

  writeLines(md.cont, con=out.fl.name)
}


pub.df <- read.table("reference.table", head=T, as.is=T, sep="\t", quote="", comment="")
pub.df <- pub.df[order(-pub.df$Year), ]
md.cont <- ""
md.cont[1] <- "---"
md.cont[2] <- "layout: page"
md.cont[3] <- "title: Literatrues"
md.cont[4] <- "group: navigation"
md.cont[5] <- "---"
md.cont[6] <- ""
for (i in 1:nrow(pub.df)) {
  title <- pub.df$Title[i]
  title <- gsub("\\)", "", title)
  title <- gsub("\\(", "", title)
  pub.df$Title[i] <- gsub("\\)", "", pub.df$Title[i])
  pub.df$Title[i] <- gsub("\\(", "", pub.df$Title[i])
  
  pub.df$Title[i] <- paste("http://www.ncbi.nlm.nih.gov/pubmed", 
                        '?term=', pub.df$Title[i],'%5BTitle%5D', 
                        sep='')
  pub.df$Title[i] <- paste('[', title, ']', '(', pub.df$Title[i], ')', 
                        sep="")
  md.ref <- paste(i, ". ", pub.df$Title[i], ", ", pub.df$Year[i], ", ",
                  pub.df$Journal[i], ".", sep="")
  md.cont <- c(md.cont, md.ref)
}
writeLines(md.cont, con="E:/GIT/funRiceGenes.github.io/docs/index.md")


meg <- readLines("git.log")
md.cont <- ""
md.cont[1] <- "---"
md.cont[2] <- "layout: page"
md.cont[3] <- "title: News"
md.cont[4] <- "comments: no"
md.cont[5] <- "thread: 616"
md.cont[6] <- "---"
md.cont[7] <- ""
for (i in seq(1, length(meg), by=6)) {
  meg.date <- meg[i+2]
  meg.con <- meg[i+4]
  meg.con <- gsub("^\\s+","", meg.con)
  meg.date.lst <- strsplit(meg.date, split="\\s+")
  meg.date.lst <- unlist(meg.date.lst)
  meg.date.vec <- meg.date.lst[c(6,3:4)]
  meg.date.str <- paste(meg.date.vec, sep="", collapse="/")
  meg.final <- paste("*", meg.date.str, meg.con, sep=" ")
  md.cont <- c(md.cont, meg.final)
}
md.cont <- gsub("\\|", "/", md.cont)
writeLines(md.cont, con="E:/GIT/funRiceGenes.github.io/news/index.md")


