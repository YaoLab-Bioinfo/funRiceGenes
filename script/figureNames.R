
fl.lst <- list.files("E:/坚果云/GIT/RICENCODE/data/Gene", 
                     patter="*.png", full=T, recur=T)
length(fl.lst)
fl.lst[1:3]

fl.lst <- fl.lst[!grepl("exp.png", fl.lst)]
fl.lst <- fl.lst[!grepl("pheno.png", fl.lst)]

length(fl.lst)

fl.lst[1:3]

fl.lst.nm <- basename(dirname(fl.lst))
fl.lst.nm[1:5]

fl.lst.df <- data.frame(fl.lst.nm, fl.lst, stringsAsFactors = FALSE)

library(plyr)
head(fl.lst.df, 2)
names(fl.lst.df) <- c("sym", "old")
fl.lst.df.new <- ddply(fl.lst.df, .(sym), function(df){
  if (nrow(df)==1) {
    df$new <- paste(dirname(df$old), "/", df$sym, ".pheno.png", sep="")
  } else if (nrow(df)==2) {
    df$new[1] <- paste(dirname(df$old[1]), "/", df$sym[1], ".pheno.png", sep="")
    df$new[2] <- paste(dirname(df$old[2]), "/", df$sym[2], ".exp.png", sep="")
  }
  return(df)
})

for (i in 1:nrow(fl.lst.df.new)) {
  file.rename(from=fl.lst.df.new$old[i], to=fl.lst.df.new$new[i])
}


