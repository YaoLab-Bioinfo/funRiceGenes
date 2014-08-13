
geocoding <- function(address="") 
{
  address <- URLencode(address, reserved=TRUE)
  url <- "http://www.datasciencetoolkit.org/maps/api/geocode/json?sensor=false&address="
  url <- paste(url, address, sep="")
  urlRes <- getURL(url)
  urlRes <- gsub("\"","",urlRes)
  if (grepl("location:", urlRes)) {
    urlRes <- gsub(".+location:", "location:", urlRes)
    urlRes <- gsub("location_type.*", "", urlRes)
    lng <- gsub(".+lng:\\s+(-*\\d+\\.\\d+).+","\\1",urlRes)
    lat <- gsub(".+lat:\\s+(-*\\d+\\.\\d+).+","\\1",urlRes)
    return(c(lng, lat))
  } else {
    return(c(NA, NA))
  }
}

