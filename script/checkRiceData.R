
library(RCurl)

i.start <- 2900
i.end <- 70000

for (i in i.start:i.end) {
  url <- paste("http://www.ricedata.cn/gene/list/", i, ".htm", sep="")
  if (url.exists(url)) {
    des.fl <- paste("E:/", basename(url), sep="")
    download.file(url, destfile=des.fl)
  }
}

