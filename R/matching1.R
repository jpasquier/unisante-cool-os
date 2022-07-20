#                __  __    _  _____ ____ _   _ ___ _   _  ____
#               |  \/  |  / \|_   _/ ___| | | |_ _| \ | |/ ___|
#               | |\/| | / _ \ | || |   | |_| || ||  \| | |  _
#               | |  | |/ ___ \| || |___|  _  || || |\  | |_| |
#               |_|  |_/_/   \_\_| \____|_| |_|___|_| \_|\____|
#
#                   ____ ___   ___  _           ___  ____
#                  / ___/ _ \ / _ \| |         / _ \/ ___|
#                 | |  | | | | | | | |   _____| | | \___ \
#                 | |__| |_| | |_| | |__|_____| |_| |___) |
#                  \____\___/ \___/|_____|     \___/|____/
#

library(MatchIt)
library(ggplot2)
library(grid)
library(gridExtra)
library(readxl)
library(writexl)

set.seed(666)
options(width = 120)

# Working directory
setwd("~/Projects/Consultations/Favre Lucie (COOL-OS)")

# Data
osteolaus <- read_xls("data-raw/Controls (OsteoLaus) for Matching.xls",
                      sheet = "selection")
cool <- read_xlsx("data-raw/Cases (COOL) for Matching.xlsx")
cool$DX_age <- as.numeric(cool$DX_age)
cool$`DX_BMI***` <- as.numeric(cool$`DX_BMI***`)
coln <- c("id", "age", "bmi")
data <- rbind(cbind(grp = 0, setNames(osteolaus, coln)),
              cbind(grp = 1, setNames(cool, coln)))
if (any(is.na(data))) stop("missing values")
rm(osteolaus, cool, coln)

# Sélection des personnes ayant plus de 55 ans
data <- data[data$age >= 55, ]

# Matching
m <- matchit(grp ~ age + bmi, method = "optimal", data = data,
             distance = "mahalanobis")

# Export results - Summary, matching plots, match data
mdir <- paste0("results/matching_", format(Sys.time(), "%Y%m%d"))
if (!dir.exists(mdir)) dir.create(mdir, recursive = TRUE)
sink(file.path(mdir, "match_smy.txt"))
print(summary(m))
sink()
pdf(file.path(mdir, "match_figs.pdf"))
plot(m, interactive = FALSE)
plot(m, "ecdf", interactive = FALSE)
dev.off()
write_xlsx(match.data(m), file.path(mdir, "match_data.xlsx"))

# Export results - Boxplots
bps <- unlist(recursive = FALSE, lapply(1:2, function(k) {
  ttl <- c("Before matching", "After mahalanobis distance matching")[k]
  d <- data
  if (k == 2) {
    d <- d[d$id %in% match.data(m)$id, ]
  }
  d$grp = factor(d$grp, 0:1, c("OsteoLaus", "COOL"))
  lapply(c("Age", "BMI"), function(Y) {
    ggplot(d, aes_string(x = "grp", y = tolower(Y))) +
      geom_boxplot() +
      labs(subtitle = ttl, x = "", y = Y)
  })
}))
pdf(file.path(mdir, "boxplots.pdf"))
grid.arrange(bps[[1]], bps[[2]], bps[[3]], bps[[4]], nrow = 2)
dev.off()

# Export results - m object and session infos
saveRDS(m, file = file.path(mdir, "matching.rda"), compress = "xz")
sink(file.path(mdir, "sessionInfo.txt"))
print(sessionInfo(), locale = FALSE)
sink()
