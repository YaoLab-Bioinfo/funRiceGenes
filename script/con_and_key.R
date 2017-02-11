

con.fls <- list.files("E:/GIT/RICENCODE/data/Gene", patter="Connection", full=T, recur=T)
key.fls <- list.files("E:/GIT/RICENCODE/data/Gene", patter="Keyword.trait", full=T, recur=T)

con.fls[1:5]
key.fls[1:5]

for (i in con.fls) {
  dat <- read.table(i, head=F, as.is=T, sep="\t", quote="", comment="")
  dat$V3 <- gsub("^\\s+", "", dat$V3)
  dat$V3 <- gsub('^"+', "", dat$V3)
  dat$V3 <- gsub("\\s+$", "", dat$V3)
  dat$V3 <- gsub('"+$', "", dat$V3)
  
  dat$V4 <- gsub("^\\s+", "", dat$V4)
  dat$V4 <- gsub('^"+', "", dat$V4)
  dat$V4 <- gsub("\\s+$", "", dat$V4)
  dat$V4 <- gsub('"+$', "", dat$V4)
  
  write.table(dat, file=i, sep="\t", quote=F, row.names=F, col.names=F)
}

for (i in key.fls) {
  dat <- read.table(i, head=T, as.is=T, sep="\t", quote="", comment="")
  dat$Title <- gsub("^\\s+", "", dat$Title)
  dat$Title <- gsub('^"+', "", dat$Title)
  dat$Title <- gsub("\\s+$", "", dat$Title)
  dat$Title <- gsub('"+$', "", dat$Title)
  
  dat$Evidence <- gsub("^\\s+", "", dat$Evidence)
  dat$Evidence <- gsub('^"+', "", dat$Evidence)
  dat$Evidence <- gsub("\\s+$", "", dat$Evidence)
  dat$Evidence <- gsub('"+$', "", dat$Evidence)
  
  write.table(dat, file=i, sep="\t", quote=F, row.names=F)
}

