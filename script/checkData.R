
dat <- read.table("E:/GIT/RICENCODE/geneInfo.table", 
                  head=T, as.is=T, sep="\t", quote="", comment="")
dim(dat)
dat[1, ]
which(is.na(dat$Symbol))
which(is.na(dat$MSU))
which(is.na(dat$RAPdb))
which(is.na(dat$path))

which(dat$Symbol=="")
which(dat$MSU=="")
which(dat$RAPdb=="")
which(dat$path=="")

which(dat$Symbol=="NA")
which(dat$MSU=="NA")
which(dat$RAPdb=="NA")
which(dat$path=="NA")

which(dat$Symbol=="none")
which(dat$MSU=="none")
which(dat$RAPdb=="none")
which(dat$path=="none")

dat[1, ]
sym.tab <- table(dat$Symbol)
rap.tab <- table(dat$RAPdb)
msu.tab <- table(dat$MSU)
path.tab <- table(dat$path)
sort(sym.tab, decreasing=T)[1:4]
sort(rap.tab, decreasing=T)[1:2]
sort(msu.tab, decreasing=T)[1:2]
sort(path.tab, decreasing=T)[1:2]

dat$MSU[dat$MSU!="None"&nchar(dat$MSU)!=14]
dat$RAPdb[dat$RAPdb!="None"&nchar(dat$RAPdb)!=12]





