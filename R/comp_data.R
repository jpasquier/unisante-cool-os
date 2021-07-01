library(readxl)
setwd("~/Projects/Consultations/Favre Lucie (COOL-OS)")
dta1 <- read_xls("data-raw/Fichier Matching pour J Pasquier.xls",
                 sheet = "COOL")
dta2 <- read_xlsx("data-raw/COOL DEXA_DO_20210607_ for stats.xlsx")
(ids <- dta2$`Subject ID`[!(dta2$`Subject ID` %in% dta1$`Subject ID`)])
dta1$`Subject ID`[!(dta1$`Subject ID` %in% dta2$`Subject ID`)]
paste(ids, collapse = ", ")
dta2[dta2$`Subject ID` %in% ids, c("DX_age", "DX_BMI")]
