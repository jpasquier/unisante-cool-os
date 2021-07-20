library(ggplot2)
library(mice)
library(miceadds)
library(pander)
library(parallel)
library(readxl)
library(writexl)

options(mc.cores = detectCores())

source("https://raw.githubusercontent.com/stamats/MKmisc/master/R/mi.t.test.R")


# Set working directory
setwd("~/Projects/Consultations/Favre Lucie (COOL-OS)")

# Output directory
outdir <- "results/comp_cool_osteolaus_20210720"
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

# Imputed data
imp <- parlmice(dta, maxit = 30, cluster.seed = 666, n.core = 16,
                n.imp.core = 2)
pdf(file.path(outdir, "imputations_diagnostic_plots.pdf"))
print(plot(imp))
dev.off()
imp_list <- lapply(1:imp$m, function(i) complete(imp, i))
imp_list <- lapply(1:imp$m, function(i) {
  z <- complete(imp, i)
  z$Cohort <- factor(z$Cohort)
  z$femurTot_Tscore <- predict(lm(femurTot_Tscore ~ femurTot_DMO, dta), z)
  z$col_Tscore <- predict(lm(col_Tscore ~ col_DMO, dta), z)
  z
})

# Comparison of the numeric variables - Table
X <- names(dta)[sapply(dta, class) == "numeric"]
tab_num <- do.call(rbind, mclapply(X, function(x) {
  m <- do.call(cbind, lapply(c("COOL", "OsteoLaus"), function(g) {
    u <- dta[dta$Cohort == g, x]
    n <- sum(!is.na(u))
    u <- u[!is.na(u)]
    r <- data.frame(n = n, mean = mean(u), SD = sd(u))
    names(r) <- paste0(names(r), " (", g, ")")
    return(r)
  }))
  pv <- t.test(as.formula(paste(x, "~ Cohort")), dta)$p.value
  m <- cbind(m, `t-test p-value` = pv)
  #M_imp <- t(sapply(c("COOL", "OsteoLaus"), function(g) {
  #  n <- sum(dta$Cohort == g)
  #  se <- sapply(imp_list, function(d) sd(d[d$Cohort == g, x]) / sqrt(n))
  #  m <- sapply(imp_list, function(d) mean(d[d$Cohort == g, x]))
  #  se <- sqrt(mean(se^2) + (1 + 1 / length(m)) * var(m))
  #  m <- mean(m)
  #  c(n = n, m = m, s = se * sqrt(n), se = se)
  #}))
  m_imp <- mi.t.test(imp_list, x = x, y = "Cohort")
  m_imp <- as.data.frame(t(c(m_imp$estimate,
                             `t-test p-value` = m_imp$p.value)))
  names(m_imp) <- paste(names(m_imp), "(imp)")
  cbind(variable = x, m, m_imp)
}))
rm(X)
write_xlsx(tab_num, file.path(outdir, "comparisons_numeric_variables.xlsx"))

# Comparison of the numeric variables - Figures
pdf(file.path(outdir, "comparisons_numeric_variables.pdf"))
for (x in names(dta)[sapply(dta, class) == "numeric"]) {
  p <- ggplot(na.omit(dta[c(x, "Cohort")]), aes_string(x = "Cohort", y = x)) +
    geom_boxplot() +
    labs(x = "")
  print(p)
}
dev.off()
rm(p, x)

# Comparison of the binary variables
X <- names(dta)[sapply(dta, class) == "logical"]
tab_bin <- do.call(rbind, mclapply(X, function(x) {
  m <- do.call(cbind, lapply(c("COOL", "OsteoLaus"), function(g) {
    u <- dta[dta$Cohort == g, x]
    n <- sum(!is.na(u))
    u <- u[!is.na(u)]
    r <- data.frame(n = n, npos = sum(u), prop = mean(u))
    names(r) <- paste0(names(r), " (", g, ")")
    return(r)
  }))
  r <- evals("chisq.test(dta[[x]], dta$Cohort)", env = environment(),
             graph.dir = "/tmp/r_pander_graph_dir")[[1]]
  pv <- r$result$p.value
  w <- paste(r$msg$warnings, collapse = "; ")
  e <- paste(r$msg$errors, collapse = "; ")
  pv2 <- fisher.test(dta[[x]], dta$Cohort)$p.value
  m <- cbind(m, `χ² test p-value` = pv, `χ² test warning` = w,
             `χ² test error` = e, `Fisher test p-value` = pv2)
  m_imp <- do.call(cbind, lapply(c("COOL", "OsteoLaus"), function(g) {
    npos <- mean(sapply(imp_list, function(d) sum(d[d$Cohort == g, x])))
    prop <- npos / sum(dta$Cohort == g)
    r <- data.frame(npos = npos, prop = prop)
    names(r) <- paste0(names(r), " (", g, ")")
    return(r)
  }))
  chi2_imp <- do.call(rbind, lapply(1:imp$m, function(k) {
    d <- complete(imp, k)
    r <- evals("chisq.test(d[[x]], d$Cohort)", env = environment(),
               graph.dir = "/tmp/r_pander_graph_dir")[[1]]
    chi2 <- r$result
    w <- paste(r$msg$warnings, collapse = "; ")
    e <- paste(r$msg$errors, collapse = "; ")
    data.frame(X.squared = chi2$statistic, df = chi2$parameter,
               warning = w, error = e)
  }))
  chi2_imp_wrn <- paste0(unique(chi2_imp$warning), collapse = "; ")
  chi2_imp_err <- paste0(unique(chi2_imp$error), collapse = "; ")
  chi2_imp_pvl <- micombine.chisquare(
    dk = chi2_imp$X.squared, df = chi2_imp$df[1], display = F
  )[["p"]]
  m_imp <- cbind(m_imp, `χ² test p-value` = chi2_imp_pvl,
                 `χ² test warning` = chi2_imp_wrn,
                 `χ² test error` = chi2_imp_err)
  names(m_imp) <- paste(names(m_imp), "(imp)")
  cbind(variable = x, m, m_imp)
}))
rm(X)
tab_bin <- tab_bin[, apply(tab_bin != "", 2, any)]
write_xlsx(tab_bin, file.path(outdir, "comparisons_binary_variables.xlsx"))

# Session info
sink(file.path(outdir, "sessionInfo.txt"))
print(sessionInfo(), locale = FALSE)
sink()

