library(car)
library(ggplot2)
library(gridExtra)
library(parallel)
library(readxl)
library(writexl)

options(mc.cores = detectCores())

# Set working directory
setwd("~/Projects/Consultations/Favre Lucie (COOL-OS)")

# Output directory
outdir <- "results/analyses_cool_only_20210819"
if (!dir.exists(outdir)) dir.create(outdir)

# Import data
dta <- read_xlsx("data-raw/1. Final for stat COOL DEXA VAT 20210715.xlsx")
dta <- as.data.frame(dta)
names(dta) <- sub("^1", "One", gsub("-", "", gsub(" ", "_", names(dta))))

# Recoding
for (j in which(sapply(dta, class) == "character")) {
  if (all(grepl("^(-)?[0-9]+(\\.[0-9]+)?$", dta[[j]]) | is.na(dta[[j]]))) {
    dta[[j]] <- as.numeric(dta[[j]])
  }
}
dta$Time_postop <- as.numeric(sub("Y$", "", sub("_", ".", dta$Time_postop)))

# Outcomes and predictors
Y <- c("L1L4s_BMD", "femurTot_DMO", "col_DMO", "L1L4sTBS")
X <- c("BS_Age", "PO_BMI", "DX_age", "DX_BMI", "Time_postop", "deltaWeight",
       "OneY_TWL", "OneY_EBMIL", "ALMI", "FMI", "VitD1", "PTH", "P1NP",
       "Bcrosslaps", "glyc", "HOMA", "HbA1c", "VAT")

# Plots
pdf(file.path(outdir, "distributions.pdf"))
par(mfrow = c(2, 2))
for (v in c(Y, X)) {
  if (v == "L1L4sTBS") {
    dta1 <- subset(dta, DX_BMI <= 37)
    sttl <- "Only patients with DX_BMI <= 27"
  } else {
    dta1 <- dta
    sttl <- NULL
  }
  hist(dta1[[v]], main = paste("Histogram of", v), sub = sttl, xlab = v)
  boxplot(dta1[[v]], main = paste("Boxplot of", v), sub = sttl)
}
par(mfrow = c(1, 1))
dev.off()
rm(v)

# Univariate regressions - figures
figs_unireg <- mclapply(setNames(Y, Y), function(y) {
  if (y == "L1L4sTBS") {
    dta <- subset(dta, DX_BMI <= 37)
  }
  dta$HOMA_no_ex <- ifelse(dta$HOMA > 10, NA, dta$HOMA)
  J <- length(X)
  X2 <- c(X, "HOMA_no_ex")[c(1:(J - 1), J + 1, J)]
  figs <- mclapply(setNames(X2, X2), function(x) {
    fml <- as.formula(paste(y, "~", x))
    fit <- lm(fml, dta)
    n <- nrow(fit$model)
    b <- signif(coef(fit)[[2]], 3)
    ci <- signif(confint(fit)[2, ], 3)
    p <- coef(summary(fit))[2, "Pr(>|t|)" ]
    p <- if (p >= 0.001) paste0("p=", round(p, 3)) else "p<0.001"
    r2 <- paste0("R2=", round(summary(fit)$r.squared, 3))
    cap <- paste0("b = ", b, " (", ci[1], ",", ci[2], ")\n", "n = ", n, ", ",
                  p, ", ", r2)
    fig <- ggplot(na.omit(dta[c(x, y)]), aes_string(y = y, x = x)) +
      geom_point() +
      geom_smooth(method = lm, formula = y ~ x, se = FALSE) +
      labs(caption = cap) +
      theme(plot.caption = element_text(size = rel(.7), hjust = 0.5))
  })
})
pdf(file.path(outdir, "univariable_regressions.pdf"),
    width = 10.375, height = 14.625)
for (figs in figs_unireg) {
  print(do.call(grid.arrange, append(figs, list(ncol = 4))))
}
dev.off()
rm(figs)

# Univariable regressions - tables
uni_reg_tbl <- mclapply(setNames(Y,Y), function(y) {
  if (y == "L1L4sTBS") {
    dta <- subset(dta, DX_BMI <= 37)
  }
  tab <- do.call(rbind, mclapply(X, function(x) {
    fml <- paste(y, "~", x)
    fit <- lm(as.formula(fml), dta)
    tab <- cbind(coef(summary(fit)), confint(fit))[, c(1, 5:6, 4)]
    tab <- cbind(n = nrow(fit$model), tab)
    colnames(tab) <- paste(c("n", "beta", "lwr", "upr", "pval"), "(uni)")
    tab[2, , drop = FALSE]
  }))
  tab <- cbind(variable = rownames(tab), as.data.frame(tab))
  rownames(tab) <- NULL
  tab
})
write_xlsx(uni_reg_tbl, file.path(outdir, "univariable_regressions.xlsx"))

# Multivariable regressions
multi_reg_tbl <- mclapply(1:3, function(k) {
  m <- c("model20", "model21", "model22")[k]
  M <- c("Model 2.0", "Model 2.1", "Model 2.2")[k]
  X <- c("BS_Age", "DX_BMI", "Time_postop", "deltaWeight", "OneY_TWL", "PTH")
  if (k == 2) {
    X <- X[X != "deltaWeight"]
  } else if (k == 3) {
    X <- X[X != "OneY_TWL"]
  }
  tabs <- mclapply(setNames(Y,Y), function(y) {
    if (y == "L1L4sTBS") {
      dta <- subset(dta, DX_BMI <= 37)
    }
    fml <- paste(y, "~", paste(X, collapse = " + "))
    fit <- do.call("lm", list(formula = as.formula(fml), data = quote(dta)))
    tab <- cbind(coef(fit), confint(fit), coef(summary(fit))[, 4])
    colnames(tab) <- paste(c("beta", "lwr", "upr", "pval"), "(muti)")
    tab <- cbind(data.frame(variable = rownames(tab), n = nrow(fit$model)),
                 tab)
    tab$dummy_n <- 1:nrow(tab)
    vif <- vif(fit)
    vif <- data.frame(variable = names(vif), vif = vif)
    tab <- merge(tab, vif, by = "variable", all = TRUE)
    tab <- tab[order(tab$dummy_n), ]
    tab$dummy_n <- NULL
    # Diagnostic plots
    figdir <- file.path(outdir, "multivariable_regressions_diagnostic_plots")
    if (!dir.exists(figdir)) dir.create(figdir)
    pdf(file.path(figdir, paste0(y, "_", m, ".pdf")))
    par(mfrow = c(2, 2))
    for (j in 1:4) plot(fit, j)
    mtext(paste(y, "-", M), side = 3, line = -2, outer = TRUE)
    par(mfrow = c(1, 1))
    dev.off()
    tab
  })
  attr(tabs, "m") <- m
  tabs
})

# Multivariable regressions - Export tables
for (tabs in multi_reg_tbl) {
  f <- paste0("multivariable_regressions_", attr(tabs, "m"), ".xlsx")
  write_xlsx(tabs, file.path(outdir, f))
}
rm(f, tabs)

# Comparison of the tertiles
cmp_tert <- mclapply(setNames(Y, Y), function(y) {
  if (y == "L1L4sTBS") {
    dta <- subset(dta, DX_BMI <= 37)
  }
  dta <- dta[!is.na(dta[[y]]), ]
  quantile(dta[[y]], (0:3) / 3)
  dta$tert <- cut(dta[[y]], quantile(dta[[y]], (0:3) / 3),
                  include.lowest = TRUE, labels = paste0("T", 1:3))
  do.call(rbind, mclapply(X, function(x) {
    n <- function(z, na.rm) sum(!is.na(z))
    fml <- as.formula(paste(x, "~ tert"))
    r <- do.call(c, mclapply(c("n", "mean", "sd"), function(fct) {
      tab <- aggregate(fml, dta, get(fct))
      setNames(tab[[2]], paste0(fct, ".", tab[[1]]))
    }))
    r <- r[do.call(c, lapply(1:3, function(k) (0:2) * 3 + k))]
    pv <- anova(lm(fml, dta))$`Pr(>F)`[1]
    cbind(data.frame(variable = x), t(r), anova.pval = pv)
  }))
})
write_xlsx(cmp_tert, file.path(outdir, "comparison_tertiles.xlsx"))

# Comparison of the tertiles - Plots
cmp_tert_figs <- mclapply(setNames(Y, Y), function(y) {
  if (y == "L1L4sTBS") {
    dta <- subset(dta, DX_BMI <= 37)
  }
  dta <- dta[!is.na(dta[[y]]), ]
  quantile(dta[[y]], (0:3) / 3)
  dta$tert <- cut(dta[[y]], quantile(dta[[y]], (0:3) / 3),
                  include.lowest = TRUE, labels = paste0("T", 1:3))
  mclapply(setNames(X, X), function(x) {
    ggplot(na.omit(dta[c("tert", x)]), aes_string(x = "tert", y = x)) +
      geom_boxplot() +
      labs(title = paste("Tertiles of", y), x = "")
  })
})
pdf(file.path(outdir, "comparison_tertiles.pdf"))
for (x in X) {
  figs <- lapply(Y, function(y) cmp_tert_figs[[y]][[x]])
  print(do.call(grid.arrange, append(figs, list(ncol = 2))))
}
dev.off()
rm(x, figs)

# Session info
sink(file.path(outdir, "sessionInfo.txt"))
print(sessionInfo(), locale = FALSE)
sink()
