

library(data.table)

## keyword file
fls <- list.files("E:/GIT/funRiceGenes/data/Gene", pattern = "^Keyword.trait$", full=T, recur=T)

fls.c <- lapply(fls, function(i) {
  i.dat <- fread(i, data.table=F)
  i.path <- dirname(i)
  i.gene <- basename(i.path)
  
  i.gene.fl <- gsub("Keyword.trait", "gene.info", i)
  i.gene.d <- fread(i.gene.fl, data.table=F)
  i.msu <- i.gene.d$MSU
  i.rap <- i.gene.d$RAPdb
  
  i.dat$MSU <- i.msu
  i.dat$RAPdb <- i.rap
  i.dat$path <- i.path
  
  i.dat <- i.dat[, c("Symbol",	"RAPdb",	"MSU",	"Keyword",	"Title",	"Evidence",	"path")]
  
  return(i.dat)
})

key.dat <- do.call(rbind, fls.c)
key.dat$path <- gsub("E:/GIT/funRiceGenes/", "", key.dat$path)

fwrite(key.dat, file="E:/GIT/funRiceGenes/geneKeyword.table", sep="\t")
fwrite(unique(key.dat[, 1:5]), file="E:/GIT/funRiceGenes/key.txt", sep="\t")


## reference file
refs <- list.files("E:/GIT/funRiceGenes/data/Gene", pattern = "^reference.info$", full=T, recur=T)
gene.info <- fread("E:/GIT/funRiceGenes/geneInfo.table", data.table=F)

refs.path <- gsub("E:/GIT/funRiceGenes/", "", dirname(refs))
setdiff(gene.info$path, refs.path)

writeLines(setdiff(gene.info$path, refs.path), "E:/genes_without_ref.txt")

