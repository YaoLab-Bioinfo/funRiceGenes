
extractInfoFromPdf <- function(file="")
{
  pdf <- try(readPDF(control = list(text = "-layout"))(elem = list(uri = file),
                                                       language = "en",
                                                       id = "id1"))
  if (class(pdf)=="try-error") {
    return(NULL)
  } else {
    locusStat <- pdf[grepl("os[0-1][0-9]g[0-9]+.*",
                           pdf,ignore.case=T)]
    cloneStat <- pdf[grepl("map-based cloning|positional cloning",
                           pdf,ignore.case=T)]
    accessStat <- pdf[grepl("accession",
                            pdf,ignore.case=T)]
    tdnaStat <- pdf[grepl("t-DNA",
                          pdf,ignore.case=T)]
    tosStat <- pdf[grepl("tos17",
                         pdf,ignore.case=T)]
    homologStat <- pdf[grepl("homolog|ortholog",
                             pdf,ignore.case=T)]
    rnaiStat <- pdf[grepl("rnai",
                          pdf,ignore.case=T)]
    oveStat <- pdf[grepl("overexpres",
                         pdf,ignore.case=T)]
    rtStat <- pdf[grepl("rt-pcr",
                        pdf,ignore.case=T)]
    NthernStat <- pdf[grepl("northern",
                            pdf,ignore.case=T)]
    WthernStat <- pdf[grepl("western",
                            pdf,ignore.case=T)]
    SthernStat <- pdf[grepl("southern",
                            pdf,ignore.case=T)]
    resDf <- as.data.frame(t(locusStat, cloneStat, accessStat, tdnaStat, 
                             tosStat, homologStat, rnaiStat, oveStat,
                             rtStat, NthernStat, WthernStat, SthernStat))
    names(resDf) <- c("locus", "cloning", "accession", "TDNA", 
                      "Tos17", "homolog", "RNAi", "overexpression",
                      "RT-PCR", "Northern", "Western", "Southern")
    return(resDf)
  }
}

