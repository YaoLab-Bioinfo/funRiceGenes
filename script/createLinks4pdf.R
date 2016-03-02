
library(BBmisc)

setwd("~")
all.pdfs <- list.files("E:/坚果云/GIT/RICENCODE/data/Gene", 
                       patter="*.pdf$", full=T, recur=T)
all.pdfs <- gsub("-冲突-andrew", "", all.pdfs)
length(all.pdfs)
all.pdfs[1:5]
length(unique(basename(all.pdfs)))

ct.fls <- all.pdfs[grepl("冲突", all.pdfs)]
str(ct.fls)
ct.fls.rt <- gsub("-冲突-andrew", "", ct.fls)
str(ct.fls.rt)
file.remove(ct.fls.rt)
file.rename(ct.fls, ct.fls.rt)

getRelativePath(to="E:/坚果云/GIT/RICENCODE/data/Gene/Abstract/0-9/alpha1a/JBC-2001-Jiang-9322-9.pdf", 
               from="E:/坚果云/GIT/RICENCODE/data/Gene/Abstract/0-9/alpha1b/JBC-2001-Jiang-9322-9.pdf")

all.pdfs.df <- data.frame(all.pdfs, basename(all.pdfs), file.size(all.pdfs),
                          stringsAsFactors = FALSE)
library(plyr)
names(all.pdfs.df) <- c("file", "pdf", "size")
d_ply(all.pdfs.df, .(pdf), function(df){
  if (nrow(df)>=2) {
    froms <- df$file[df$size==0]
    tos <- df$file[df$size>0]
    if (length(froms)>0 && length(tos)>0) {
      for (i in froms) {
        file.remove(i)
        from.pt <- getRelativePath(to=tos[1], from=i)
        from.pt <- sub("\\.\\.\\/", "", from.pt)
        to.path <- dirname(i)
        setwd(to.path)
        file.symlink(from=from.pt, to=".")
        setwd("~")
      }
    }
  }
})


