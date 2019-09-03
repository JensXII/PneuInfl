################
### Raw data ###
################

### Influenza ###

# Samples
InflSamples <- read.table("S:/BigData/data/influenza/Influenza data 2010-2018.csv", header = TRUE, sep =";")
InflSamples$cprnr <- as.character(InflSamples$cprnr)
InflSamples$prdate <- as.Date(InflSamples$prdate)
InflSamples[is.na(InflSamples)] <- 0
InflSamples <- transform(InflSamples, a = as.numeric((a == 1) | (h1n1 == 1) | (h3n2 == 1)))
InflSamples <- transform(InflSamples,
                         h1n1 = as.numeric((h1n1 == 1) | ((a == 1) & (h3n2 == 0))),
                         h3n2 = as.numeric((h3n2 == 1) | ((a == 1) & (h1n1 == 0)))
)
InflSamples <- transform(InflSamples, pos = as.numeric((a == 1) | (b == 1)))
InflSamples <- InflSamples[with(InflSamples, order(cprnr, prdate)),]
InflSamples <- aggregate(cbind(pos, a, h3n2, h1n1, b) ~ cprnr + prdate, data = InflSamples, max)
save(InflSamples, file = "data/InflSamples.RData")

# Cases
InflSamples$season <- as.numeric(format(InflSamples$prdate, "%Y")) - (as.numeric(format(InflSamples$prdate, "%m")) <= 9)
InflData <- aggregate(prdate ~ season + cprnr, data = subset(InflSamples, pos == 1), min)[, c("cprnr", "prdate")]
save(InflData, file = "data/InflData.RData")
rm(InflSamples, InflData)

### IPD ###

# Cases
IPDData <- read.table("S:/BigData/data/Pneumococcus/pneumokokker_kun_data.csv", header = TRUE, sep =";")
IPDData$cprnr <- sprintf("%010.0f", IPDData$cprnr)
IPDData$prdate <- as.Date(IPDData$pr_vedato, format = "%d-%B-%y")
IPDData <- IPDData[!duplicated(IPDData[, c("cprnr", "prdate")]), c("cprnr", "prdate")]
save(IPDData, file = "data/IPDData.RData")
rm(IPDData)

### Streptococcus pneumoniae ###

# Samples
con <- RODBC::odbcConnect("DKMOMO")
PneuSamples <- RODBC::sqlQuery(con,
                      "select distinct replace(header_cprnr, '-', '') as cprnr, header_Prdate as prdate, concat(lar, '-', header_labnr) as LarLabnr,
                      TabSpecimen_Text,
                      TabMicroorganism_Banr, TabMicroorganism_Text,
                      Quantitative_Qtnr, TabAnalysis_Text from
                      IB_EpiMiba.mikro.vw_pneumokok_all with(nolock)
                      where TabSpecimen_Text != 'Kvalitetskontrolprøve'",
                      as.is = TRUE, stringsAsFactors = FALSE
)
RODBC::odbcClose(con)
rm(con)
save(PneuSamples, file = "data/PneuSamples.RData")

# Cases
PneuData <- merge(PneuSamples, read.table("S:/BigData/BigData/data/IPD_SampleType.csv", header = TRUE, sep =";"), by = "TabSpecimen_Text", all.x = TRUE)
PneuData <- aggregate(IPD ~ cprnr + prdate, data = PneuData, max)
PneuData$prdate <- as.Date(PneuData$prdate)
PneuData <- PneuData[with(PneuData, order(cprnr, prdate)),]

PneuData.Pneu <- PneuData
PneuData.Pneu$lag.prdate <- c(NA, PneuData.Pneu$prdate[-nrow(PneuData.Pneu)])
PneuData.Pneu$lag.prdate[which(!duplicated(PneuData.Pneu$cprnr))] <- NA
PneuData.Pneu <- subset(PneuData.Pneu, ((prdate - lag.prdate) > 30) | is.na(lag.prdate))[, c("cprnr", "prdate")]
PneuData.Pneu$Pneu <- 1

PneuData.IPD <- subset(PneuData, IPD == 1)
PneuData.IPD$lag.prdate <- c(NA, PneuData.IPD$prdate[-nrow(PneuData.IPD)])
PneuData.IPD$lag.prdate[which(!duplicated(PneuData.IPD$cprnr))] <- NA
PneuData.IPD <- subset(PneuData.IPD, ((prdate - lag.prdate) > 30) | is.na(lag.prdate))[, c("cprnr", "prdate", "IPD")]

PneuData <- merge(PneuData.Pneu, PneuData.IPD, by = c("cprnr", "prdate"), all = TRUE)
PneuData <- PneuData[with(PneuData, order(cprnr, prdate)),]
rm(PneuData.Pneu, PneuData.IPD)
PneuData[is.na(PneuData$IPD), "IPD"] <- 0
PneuData[is.na(PneuData$Pneu), "Pneu"] <- 0
save(PneuData, file = "data/PneuData.RData")
rm(PneuSamples, PneuData)
