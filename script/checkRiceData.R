
library(RCurl)

i.start <- 2872
i.end <- 2900

for (i in i.start:i.end) {
  url <- paste("http://www.ricedata.cn/gene/list/", i, ".htm", sep="")
  if (url.exists(url)) {
    des.fl <- paste("E:/", basename(url), sep="")
    download.file(url, destfile=des.fl)
  }
}

