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

library(readxl)
library(MatchIt)
library(ggplot2)
library(grid)
library(gridExtra)

set.seed(666)
options(width = 120)

# Working directory
setwd("~/Projects/Consultations/Favre Lucie (COOL-OS)")

# Data
fileName <- "data-raw/Fichier Matching pour J Pasquier.xls"
osteolaus <- read_xls(fileName, sheet = "osteoLaus")
cool <- read_xls(fileName, sheet = "COOL")
cool$DX_age <- as.numeric(cool$DX_age)
cool$`BMI (kg/m^2)` <- as.numeric(cool$`BMI (kg/m^2)`)
coln <- c("id", "age", "bmi")
data <- rbind(cbind(grp = 0, setNames(osteolaus, coln)),
              cbind(grp = 1, setNames(cool, coln)))
if (any(is.na(data))) stop("missing values")
rm(fileName, osteolaus, cool, coln)

# Apperçu des données
tmp <- data
tmp$grp <- factor(tmp$grp, 0:1, c("OsteoLaus", "COOL"))
bps <- list()
bps$age <- ggplot(tmp, aes(x = factor(grp), y = age)) +
  geom_boxplot() + labs(title = "Age", x = "")
bps$bmi <- ggplot(tmp, aes(x = factor(grp), y = bmi)) +
  geom_boxplot() + labs(title = "BMI", x = "")
rm(tmp)

# Sélection des personnes ayant plus de 55 ans
data <- data[data$age >= 55, ]

# Matching
mlist <- lapply(c(cs0 = 1, cs1 = 2, md = 3), function(l) {
  lapply(c(r11 = 1, r21 = 2, r31 = 3), function(r) {
    discard <- c("none", "both", "none")[l]
    distance <- c("glm", "glm", "mahalanobis")[l]
    m <- matchit(grp ~ age + bmi, method = "optimal", data = data,
                 ratio = r, distance = distance, discard = discard)
    v <- c("id", "age", "bmi")
    mm <- data[rownames(m$match.matrix), v]
    if (l %in% 1:2) {
      mm <- cbind(mm, ps = m$distance[rownames(m$match.matrix)])
    }
    names(mm) <- paste0(names(mm), "_exp")
    for (u in 1:r) {
      mm1 <- data[m$match.matrix[, u], v]
      if (l %in% 1:2) {
        mm1 <- cbind(mm1, ps = m$distance[m$match.matrix[, u]])
      }
      names(mm1) <- paste0(names(mm1), "_ctrl", u)
      mm <- cbind(mm, mm1)
    }
    list(m = m, mm = mm)
  })
})

# Included persons
for (s in names(mlist)) {
  for (r in names(mlist[[s]])) {
    inc_id <- na.omit(mlist[[s]][[r]]$mm[c("id_exp", "id_ctrl1")])
    inc_id <- c(inc_id[[1]], inc_id[[2]])
    data[[paste("inc", s, r, sep = "_")]] <- data$id %in% inc_id
  }
}
rm(s, r, inc_id)

# Nombre de contrôles à inclure si on veut considérer tous les appariements
ctrl_ids <- c()
for (s in names(mlist)) {
  for (r in names(mlist[[s]])) {
    mm <- mlist[[s]][[r]]$mm
    ctrl_ids <-
      c(ctrl_ids, na.omit(do.call(c, mm[grep("^id_ctrl", names(mm))])))
  }
}
ctrl_ids <- unique(ctrl_ids)
length(ctrl_ids)

# Export results - Summaries, matrices and matching plots
mdir <- "results/matching_20210607"
if (!dir.exists(mdir)) dir.create(mdir, recursive = TRUE)
for (s in names(mlist)) {
  for (r in names(mlist[[s]])) {
    m <- mlist[[s]][[r]]$m
    mm <- mlist[[s]][[r]]$mm
    u <- paste(s, r, sep = "_")
    sink(file.path(mdir, paste0("match_smy_", u, ".txt")))
    print(summary(m))
    sink()
    pdf(file.path(mdir, paste0("match_figs_", u, ".pdf")))
    plot(m, interactive = FALSE)
    plot(m, "ecdf", interactive = FALSE)
    if (s != "md") {
      plot(m, type = "jitter", interactive = FALSE)
      plot(m, type = "hist")
    }
    dev.off()
    write.table(mm, file.path(mdir, paste0("matching_matrix_", u, ".csv")),
                sep = ";", row.names = FALSE)
  }
}
rm(s, r, m, mm, u)

# Export results - Osteolaus ID
write.table(data.frame(`OsteoLaus ID` = ctrl_ids), row.names = FALSE,
            quote = FALSE, file = file.path(mdir, "matched_osteolaus_ids.csv"))

# Export results - Boxplots before matching
pdf(file.path(mdir, "boxplots_before_matching.pdf"))
print(bps$age)
print(bps$bmi)
dev.off()

# Export results - Boxplots after matching
pdf(file.path(mdir, "boxplots_after_matching.pdf"))
for (k in 1:3) {
  s <- c("cs0", "cs1", "md")[k]
  ttl <- c(
    "Propensity score matching with all COOL participants",
    "Propensity score matching with common support",
    "Mahalanobis distance matching"
  )[k]
  p <- unlist(recursive = FALSE, lapply(1:3, function(r) {
    tmp <- data[data[[paste0("inc_", s, "_r", r, "1")]], ]
    tmp$grp <- factor(tmp$grp, 0:1, c("OsteoLaus", "COOL"))
    p1 <- ggplot(tmp, aes(x = grp, y = age)) +
      geom_boxplot() +
      labs(subtitle = paste0("Ratio ", r, ":1"), x = "", y = "Age")
    p2 <- ggplot(tmp, aes(x = grp, y = bmi)) +
      geom_boxplot() +
      labs(subtitle = paste0("Ratio ", r, ":1"), x = "", y = "BMI")
    list(p1, p2)
  }))
  grid.arrange(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], p[[6]], nrow = 3,
               top = textGrob(ttl, gp = gpar(fontsize = 20)))
}
dev.off()
rm(s, ttl, p)

# Export results - Figs - Included participants
V <- grep("^inc_", names(data), value = TRUE)
tmp <- reshape(data, varying = V, v.names = "inc", timevar = "matching",
               times = sub("inc_", "", V), idvar = "id", direction = "long")
pdf(file.path(mdir, "included_participants.pdf"))
for (k in 1:3) {
  ttl <- paste(c(rep("Propensity score", 2), "Mahalanobis distance"),
               "matching")[k]
  sttl <- c("All COOL observations", "Common support")[c(1, 2, 1)[k]]
  fig <- ggplot(tmp[grepl(c("^cs0", "^cs1", "^md")[k], tmp$matching), ],
                aes(x = age, y = bmi, color = inc)) +
    geom_point() +
    facet_grid(matching ~ factor(grp, 0:1, c("OsteoLaus", "COOL"))) +
    labs(title = ttl, subtitle = sttl, x = "Age", y = "BMI",
         color = "Included")
  print(fig)
}
dev.off()
rm(V, tmp, k, ttl, sttl, fig)

# Export results - mlist object and session infos
save(mlist, file = file.path(mdir, "mlist.rda"), compress = "xz")
sink(file.path(mdir, "sessionInfo.txt"))
print(sessionInfo(), locale = FALSE)
sink()

# --------------------------------------------------------------------------- #

if (FALSE) {
fit <- glm(grp ~ scale(age) + scale(bmi), family = binomial, data = data)
summary(fit)
data$ps <- predict(fit, type = "response")
sapply(mlist, function(z) {
  sapply(z, function(w) {
           max(abs(data$ps - w$m$distance))
  })
})
boxplot(ps~factor(grp), data)
sort(data[data$grp == 0, "ps"], decreasing = TRUE)[1:10]
}
