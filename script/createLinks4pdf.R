
# library(BBmisc)

setwd("~")
all.pdfs <- list.files("E:/坚果云/GIT/RICENCODE/data/Gene", 
                       patter="*.pdf$", full=T, recur=T)
length(all.pdfs)
all.pdfs[1:5]
length(unique(basename(all.pdfs)))

# getRelativePath(to="E:/坚果云/GIT/RICENCODE/data/Gene/Abstract/0-9/alpha1a/J. Biol. Chem.-2001-Jiang-9322-9.pdf", 
#                from="E:/坚果云/GIT/RICENCODE/data/Gene/Abstract/0-9/alpha1b/J. Biol. Chem.-2001-Jiang-9322-9.pdf")

all.pdfs.df <- data.frame(all.pdfs, basename(all.pdfs), stringsAsFactors = FALSE)
library(plyr)
names(all.pdfs.df) <- c("file", "pdf")
d_ply(all.pdfs.df[11:nrow(all.pdfs.df), ], .(pdf), function(df){
  if (nrow(df)>=2) {
    for (i in 2:nrow(df)) {
      unlink(df$file[i])
      file.symlink(from=df$file[1], to=df$file[i])
    }
  }
})

