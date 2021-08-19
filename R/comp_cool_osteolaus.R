library(ggplot2)
library(pander)
library(parallel)
library(readxl)
library(writexl)

options(mc.cores = detectCores())

# Set working directory
setwd("~/Projects/Consultations/Favre Lucie (COOL-OS)")

# Output directory
outdir <- "results/comp_cool_osteolaus_20210819"
if (!dir.exists(outdir)) dir.create(outdir)

# Import data
dta <- read_xlsx("data-raw/COOL_OsteoLaus_47n2Matcheddata.xlsx")
dta <- as.data.frame(dta)
names(dta) <- gsub("-", "", gsub(" ", "_", names(dta)))

# Recoding
for (j in which(sapply(dta, class) == "character")) {
  if (all(grepl("^(-)?[0-9]+(\\.[0-9]+)?(E(-)?[0-9]+)?$", dta[[j]]) | 
            is.na(dta[[j]]))) {
    dta[[j]] <- as.numeric(dta[[j]])
  }
}
for (j in which(sapply(dta, class) == "numeric")) {
  if (all(dta[[j]] %in% 0:1 | is.na(dta[[j]]))) {
    dta[[j]] <- as.logical(dta[[j]])
  }
}

# Import matching matrix
tarfile <- "results/matching_20210607.txz"
csvfile <- "matching_20210607/matching_matrix_md_r11.csv"
exdir <- paste0("tmp-", round(runif(1) * 10^6))
untar(tarfile, csvfile, exdir = exdir)
mm <- read.csv(file.path(exdir, csvfile), sep = ";")
unlink(exdir, recursive = TRUE)

# Help funtion - identifiy paired data which are complete for a given variable
ids <- function(x) {
  id_exp <- dta[!is.na(dta[[x]]) & dta$Cohort == "COOL", "Subject_ID"]
  id_ctrl <- dta[!is.na(dta[[x]]) & dta$Cohort == "OsteoLaus", "Subject_ID"]
  i <- mm$id_exp %in% id_exp & mm$id_ctrl1 %in% id_ctrl
  c(mm$id_exp[i], mm$id_ctrl1[i])
}

# Subgroup DX_BMi <= 37 (COOL) for L1L4sTBS
dta$L1L4sTBS_BMI.le.37 <- ifelse(
  dta$Cohort == "OsteoLaus" | dta$DX_BMI <= 37, dta$L1L4sTBS, NA)

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
  p <- ggplot(subdta[c(x, "Cohort")], aes_string(x = "Cohort", y = x)) +
    geom_boxplot() +
    labs(x = "", caption = paste("n (per group):", nrow(subdta)))
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

