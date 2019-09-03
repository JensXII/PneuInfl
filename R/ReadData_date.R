################################
### Create data for analyses ###
################################


#' Temperature data from WeatherData (http://euromomo.eu/methods/weather)
#'
#' @param CountryCode NUTS0 country code
#' @param start_year first year to be included. Must be >= 2000
#' @param end_year last week to be included
#' @return Return date, temp = mean messaured temperature, ptemp = predicted temperature from year \code{start_year} to \code{end_year}.
#' @examples
#' TemperatureData("DK", 2010, 2017) Danish temperature data from 2010 to and including 2017
TemperatureData <- function(CountryCode, start_year, end_year) {
  ET <- read.table(paste0("http://euromomo.eu/methods/weather/wdata_", CountryCode,".txt"), header = TRUE, sep = ";", dec = ".", stringsAsFactors = FALSE)[, c("date","pop3","NUTS3","temp")]
  ET$date <- as.Date(as.character(ET$date,"%Y%m%d"))
  ET[is.na(ET$pop3), "pop3"] <- 1
  ET <- aggregate(cbind(temp, pop3) ~ NUTS3 + date, data=ET, function(x) mean(x, na.rm = TRUE))
  ET$pop3.sum <- with(ET, tapply(pop3, date, sum))[as.character(ET$date)]
  ET$temp <- ET$temp * ET$pop3 / ET$pop3.sum
  ET <- aggregate(temp ~ date, data = ET, function(x) sum(x, na.rm = TRUE))
  ET <- merge(ET, data.frame(date = seq.Date(as.Date(paste0(start_year,"-01-01")), as.Date(paste0(end_year,"-12-31")), by = "day")), by = "date", all.y = TRUE)
  ET$sin52 <- sin((2*pi/(365.25)) * as.numeric(ET$date))
  ET$cos52 <- cos((2*pi/(365.25)) * as.numeric(ET$date))
  ET$ptemp <- predict(glm(temp ~ sin52 + cos52, data=ET[!(is.na(ET$temp) | is.infinite(ET$temp)),]), ET)
  ET[,c("cos52","sin52")] <- NULL
  return(ET <- ET[order(ET$date),])
}

#' Population by date.
#'
#' @param start start year.
#' @param end end year.
#' @return Return date, N = population from year \code{start} to \code{end}.
#' @examples
#' pop(2010, 2017)
pop <- function(start, end) {
  con <- RODBC::odbcConnect("DKMOMO")
  res <- RODBC::sqlQuery(con, paste0("
                          select date, sum(N) as N from
                          EuroMOMO.DKMOMO.DKpopDateRegionSexAgegrp with(nolock)
                          where (", start, " <= datepart(year, date)) and (datepart(year, date) <= ",end, ")
                          group by date
                          order by date
                         "))
  RODBC::odbcClose(con)
  rm(con)
  return(res)
}

#' Data for conditional probability by date.
#'
#' @param X name of dataframe with dependent cases.
#' @param Y name of dataframe with conditioned on cases.
#' @param s start days in conditional period.
#' @param e end days in conditional period.
#' @return Return date, N\code{Y, d, e} = number of conditioned on cases, N\code{X, Y, d, e} = number in intersection.
#' @examples
#' InflPneu("Pneu", "Infl", 0, 15) -> P(Pneu | Infl within 0 to 15 days before)
#' InflPneu("Infl", "Pneu", 1, 30) -> P(Infl | Pneu within 1 to 30 days before)
InflPneu <- function(X, Y, s, e) {
  res <- data.frame(date = seq.Date(as.Date("2010-09-15"), as.Date("2017-12-31"), by = "day"))
  res <- sqldf::sqldf(paste0("
                        select a.date, case when b.cprnr is NULL then 'NA' else b.cprnr end as cprnr, (a.date - b.prdate) as N from
                        res as a
                        left join ",
                        Y, "Data as b
                        on (", s, " <= (a.date - b.prdate)) & ((a.date - b.prdate) <= ", e, ")
                        order by a.date, cprnr, N
                      ")
          )
  res <- aggregate(N ~ date + cprnr, data = res, min, na.rm = FALSE, na.action = NULL)
  res$N <- ifelse(is.na(res$N), 0, 1)

  both <- sqldf::sqldf(paste0("
                        select a.prdate as date, a.cprnr, (a.prdate - b.prdate) as Nboth from ",
                        Y, "Data as a
                        join ",
                        X ,"Data as b
                        on (a.cprnr = b.cprnr) and (", s," <= (a.prdate - b.prdate)) & ((a.prdate - b.prdate) <= ", e, ")
                        order by a.prdate, a.cprnr, Nboth
                      ")
          )
  both <- aggregate(Nboth ~ date + cprnr, data = both, min)
  both$Nboth <- 1

  res <- merge(aggregate(N ~ date, data = res, sum), aggregate(Nboth ~ date, data = both, sum), by = "date", all.x = TRUE)
  res[is.na(res)] <- 0
  colnames(res) <- c("date", paste0("N", Y, ".prev"), paste0("N", X, Y, ".prev"))
  res <- res[with(res, order(date)),]
  return(res)
}

#' Number of cases by date.
#'
#' @param X Type ("Pneu" or "Infl") dataframe with cases.
#' @return Return date, N\code{X} = number of cases.
#' @examples
#' CaseNumber("Pneu")
#' CaseNumber("Infl")
CaseNumber <- function(X) {
  Number <- paste0("N", substr(deparse(substitute(X)),1,4))
  X$N <- 1
  res <- data.frame(date = seq.Date(as.Date("2010-09-15"), as.Date("2017-12-31"), by = "day"))
  res <- merge(res, X, by.x = "date", by.y = "prdate", all.x = TRUE)
  res <- aggregate(N ~ date, data = res, sum, na.rm = FALSE, na.action = NULL)
  res[is.na(res)] <- 0
  colnames(res) <- c("date", Number)
  res <- res[with(res, order(date)),]
  return(res)
}

#' Data for analyses.
#'
#' @param s start days in conditional period.
#' @param e end days in conditional period.
#' @return Return date, N\code{Y, d, e} = number of conditioned on cases, N\code{X, Y, d, e} = number in intersection.
#' @examples
#' CreateData(0, 15) -> data/PneuInfl015.data.RData
#' CreateData(1, 30) -> data/PneuInfl130.data.RData
CreateData <- function(s, e) {

  ### P(Pneu | Infl within d days before) ###
  Pneu <- merge(CaseNumber(PneuData), InflPneu("Pneu", "Infl", s, e), by = "date", all = TRUE)

  ### P(Infl | Pneu within d days before) ###
  Infl <- merge(CaseNumber(InflData), InflPneu("Infl", "Pneu", s, e), by = "date", all = TRUE)

  # Merge Pneu and Infl
  PneuInfl.data <- merge(Pneu, Infl, by = "date", all = TRUE)

  # Add population data
  PneuInfl.data <- merge(PneuInfl.data, pop(2010, 2017), by = "date", all.x = TRUE)

  # Add temperature data
  PneuInfl.data <- merge(PneuInfl.data, TemperatureData("DK", 2010, 2017), by = "date", all.x = TRUE)

  rm(Infl, Pneu)

  return(PneuInfl.data)
}
