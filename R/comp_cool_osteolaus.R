library(ggplot2)
library(pander)
library(parallel)
library(readxl)
library(writexl)

options(mc.cores = detectCores())

# Set working directory
setwd("~/Projects/Consultations/Favre Lucie (COOL-OS)")

# Output directory
outdir <- paste0("results/comp_cool_osteolaus_", format(Sys.Date(), "%Y%m%d"))
if (!dir.exists(outdir)) dir.create(outdir)

# Import data
dta <- c(
  cool = "data-raw/COOL DEXA_DO_2022_20221130_for_stat.xlsx",
  osteolaus = "data-raw/OsteoLaus (Controls)_2022_20221130_for_stat.xls"
)
dta <- lapply(dta, read_excel)

# Rename variables
dta <- lapply(dta, function(d) {
  names(d) <- gsub("-", "", gsub(" ", "_", names(d)))
  d
})

# Keep only matched observation
dta$cool <- dta$cool[dta$cool$DOF50Ymatch %in% "true", ]

# Keep only common variables
v <- lapply(dta, names)
v <- intersect(v[[1]], v[[2]])
dta <- lapply(dta, function(d) d[v])

# Recoding
dta <- lapply(dta, function(d) {
  for (j in which(sapply(d, class) == "character")) {
    if (all(grepl("^(-)?[0-9]+(\\.[0-9]+)?(E(-)?[0-9]+)?$", d[[j]]) |
              is.na(d[[j]]))) {
      d[[j]] <- as.numeric(d[[j]])
    }
  }
  for (j in which(sapply(d, class) == "numeric")) {
    if (all(d[[j]] %in% 0:1 | is.na(d[[j]]))) {
      d[[j]] <- as.logical(d[[j]])
    }
  }
  d
})

# Bind dataframes
dta <- rbind(cbind(dta$cool, Cohort = "COOL"),
             cbind(dta$osteolaus, Cohort = "OsteoLaus"))

# Import matching data
tarfile <- "results/matching_20220720.txz"
xlsxfile <- "matching_20220720/match_data.xlsx"
exdir <- paste0("tmp-", round(runif(1) * 10^6))
untar(tarfile, xlsxfile, exdir = exdir)
md <- read_xlsx(file.path(exdir, xlsxfile))
unlink(exdir, recursive = TRUE)

# Help funtion - identifiy paired data which are complete for a given variable
ids <- function(x) {
  d <- md[md$id %in% dta$Subject_ID[!is.na(dta[[x]])], ]
  d <- aggregate(id ~ subclass, d, function(z) length(unique(z)))
  s <- d$subclass[d$id == 2]
  md$id[md$subclass %in% s]
}

# Comparison of the numeric variables - Table
X <- names(dta)[sapply(dta, class) == "numeric"]
tab_num <- do.call(rbind, mclapply(X, function(x) {
  dta <- dta[dta$Subject_ID %in% ids(x), ]
  n <- nrow(dta) / 2
  m <- do.call(cbind, lapply(c("COOL", "OsteoLaus"), function(g) {
    u <- dta[dta$Cohort == g, x]
    r <- data.frame(mean = mean(u), SD = sd(u))
    names(r) <- paste0(names(r), " (", g, ")")
    return(r)
  }))
  pv <- t.test(as.formula(paste(x, "~ Cohort")), dta)$p.value
  m <- cbind(m, `t-test p-value` = pv)
  cbind(variable = x, n = n, m)
}))
rm(X)
write_xlsx(tab_num, file.path(outdir, "comparisons_numeric_variables.xlsx"))

# Comparison of the numeric variables - Figures
pdf(file.path(outdir, "comparisons_numeric_variables.pdf"))
for (x in names(dta)[sapply(dta, class) == "numeric"]) {
  subdta <- dta[dta$Subject_ID %in% ids(x), ]
  p <- ggplot(subdta[c(x, "Cohort")], aes(x = Cohort, y = !!sym(x))) +
    geom_boxplot() +
    labs(x = "", caption = paste("n (per group):", nrow(subdta) / 2))
  print(p)
}
dev.off()
rm(p, x)

# Comparison of the binary variables
X <- names(dta)[sapply(dta, class) == "logical"]
tab_bin <- do.call(rbind, mclapply(X, function(x) {
  dta <- dta[dta$Subject_ID %in% ids(x), ]
  n <- nrow(dta) / 2
  m <- do.call(cbind, lapply(c("COOL", "OsteoLaus"), function(g) {
    u <- dta[dta$Cohort == g, x]
    r <- data.frame(npos = sum(u), prop = mean(u))
    names(r) <- paste0(names(r), " (", g, ")")
    return(r)
  }))
  r <- evals("chisq.test(dta[[x]], dta$Cohort)", env = environment(),
             graph.dir = "/tmp/r_pander_graph_dir")[[1]]
  pv <- r$result$p.value
  w <- paste(r$msg$warnings, collapse = "; ")
  e <- paste(r$msg$errors, collapse = "; ")
  pv2 <- fisher.test(dta[[x]], dta$Cohort)$p.value
  cbind(variable = x, n = n, m, `χ² test p-value` = pv, `χ² test warning` = w,
        `χ² test error` = e, `Fisher test p-value` = pv2)
}))
rm(X)
tab_bin <- tab_bin[, apply(tab_bin != "", 2, any)]
write_xlsx(tab_bin, file.path(outdir, "comparisons_binary_variables.xlsx"))

# Session info
sink(file.path(outdir, "sessionInfo.txt"))
print(sessionInfo(), locale = FALSE)
sink()
