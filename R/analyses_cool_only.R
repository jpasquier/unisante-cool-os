library(car)
library(ggplot2)
library(gridExtra)
library(mice)
library(parallel)
library(readxl)
library(writexl)

options(mc.cores = detectCores())

# Set working directory
setwd("~/Projects/Consultations/Favre Lucie (COOL-OS)")

# Output directory
outdir <- "results/analyses_cool_only_20210701"
if (!dir.exists(outdir)) dir.create(outdir)

# Import data
dta <- read_xlsx("./data-raw/COOL DEXA_DO_20210607_ for stats.xlsx")
dta <- as.data.frame(dta)
names(dta) <- sub("^1", "One", gsub("-", "", gsub(" ", "_", names(dta))))

# Recoding
for (j in which(sapply(dta, class) == "character")) {
  if (all(grepl("^(-)?[0-9]+(\\.[0-9]+)?$", dta[[j]]) | is.na(dta[[j]]))) {
    dta[[j]] <- as.numeric(dta[[j]])
  }
}
dta$DO_is_included <- c(false = 0, true = 1)[dta$DO_is_included]
dta$Time_postop <- as.numeric(sub("Y$", "", sub("_", ".", dta$Time_postop)))

# Outcomes and predictors
Y <- c("L1L4s_BMD", "femurTot_DMO", "col_DMO")
X <- c("BS_Age", "PO_BMI", "DX_age", "DX_BMI", "Time_postop", "deltaWeight",
       "OneY_TWL", "OneY_EBMIL", "ALMI", "FMI", "VitD1", "PTH", "P1NP",
       "Bcrosslaps", "glyc", "HOMA", "HbA1c")

# Plots
pdf(file.path(outdir, "distributions.pdf"))
par(mfrow = c(2, 2))
for (v in c(Y, X)) {
  hist(dta[[v]], main = paste("Histogram of", v), xlab = v)
  boxplot(dta[[v]], main = paste("Boxplot of", v))
}
par(mfrow = c(1, 1))
dev.off()
rm(v)

# Imputed data
imp <- parlmice(dta[c(Y, X)], maxit = 30, cluster.seed = 666,
                n.core = 16, n.imp.core = 2)
pdf(file.path(outdir, "imputations_diagnostic_plots.pdf"))
print(plot(imp))
dev.off()

# Univariate regressions - figures
figs_unireg <- mclapply(setNames(Y, Y), function(y) {
  dta$HOMA_no_ex <- ifelse(dta$HOMA > 10, NA, dta$HOMA)
  J <- length(X)
  X2 <- c(X, "HOMA_no_ex")[c(1:(J - 1), J + 1, J)]
  figs <- mclapply(setNames(X2, X2), function(x) {
    fml <- as.formula(paste(y, "~", x))
    fit <- lm(fml, dta)
    #imp_fit <- with(imp, lm(as.formula(paste(y, "~", x))))
    n <- nrow(fit$model)
    b <- signif(coef(fit)[[2]], 3)
    #b_fit <- pool(imp_fit)$pooled$estimate
    ci <- signif(confint(fit)[2, ], 3)
    p <- coef(summary(fit))[2, "Pr(>|t|)" ]
    p <- if (p >= 0.001) paste0("p=", round(p, 3)) else "p<0.001"
    r2 <- paste0("R2=", round(summary(fit)$r.squared, 3))
    cap <- paste0("b = ", b, " (", ci[1], ",", ci[2], ")\n", "n = ", n, ", ",
                  p, ", ", r2)
    fig <- ggplot(na.omit(dta[c(x, y)]), aes_string(y = y, x = x)) +
      geom_point() +
      #geom_abline(intercept = b_fit[1], slope = b_fit[2], colour = "red") +
      geom_smooth(method = lm, formula = y ~ x, se = FALSE) +
      labs(caption = cap) +
      #theme(axis.title=element_text(size = rel(.75)))
      theme(plot.caption = element_text(size = rel(.7), hjust = 0.5))
  })
})
pdf(file.path(outdir, "univariable_regressions.pdf"),
    width = 10.375, height = 14.625)
for (figs in figs_unireg) {
  print(do.call(grid.arrange, append(figs, list(ncol = 4))))
}
dev.off()

# Univariate and multivariate regressions with imputed dataset
X <- c("BS_Age", "PO_BMI", "DX_age", "DX_BMI", "Time_postop", "deltaWeight",
       "OneY_TWL", "OneY_EBMIL", "ALMI", "FMI", "VitD1", "PTH", "P1NP",
       "Bcrosslaps", "glyc", "HOMA", "HbA1c")
reg_tbl <- mclapply(1:2, function(k) {
  mclapply(setNames(Y,Y), function(y) {
    if (k == 2) {
      fml0 <- paste(y, "~", paste(X, collapse = " + "))
      sel <- do.call(cbind, mclapply(1:imp$m, function(i) {
        m <- step(lm(as.formula(fml0), complete(imp, i)), trace = FALSE)
        X %in% names(coef(m))[coef(summary(m))[, 4] < .1]
      }))
      rownames(sel) <- X
      sel <- apply(sel, 1, mean)
      X <- names(sel)[sel >= 0.5]
    }
    tab_uni <- do.call(rbind, mclapply(X, function(x) {
      fml <- paste(y, "~", x)
      fit <- lm(as.formula(fml), dta)
      tab <- cbind(coef(summary(fit)), confint(fit))[, c(1, 5:6, 4)]
      tab <- cbind(n = nrow(fit$model), tab)
      colnames(tab) <- paste(c("n", "beta", "lwr", "upr", "pval"), "(uni)")
      imp_fit <- with(imp, lm(as.formula(fml)))
      imp_tab <- summary(pool(imp_fit), conf.int = TRUE)[, c(2, 7:8, 6)]
      colnames(imp_tab) <- paste(c("beta", "lwr", "upr", "pval"), "(uni imp)")
      cbind(tab, imp_tab)[2, ]
    }))
    tab_uni <- cbind(variable = rownames(tab_uni), tab_uni)
    rownames(tab_uni) <- NULL
    fml <- paste(y, "~", paste(X, collapse = " + "))
    imp_fit <- with(imp, lm(as.formula(fml)))
    tab_multi <- summary(pool(imp_fit), conf.int = TRUE)[, c(1:2, 7:8, 6)]
    colnames(tab_multi) <-
      c("variable", paste(c("beta", "lwr", "upr", "pval"), "(muti imp)"))
    tab_multi$dummy_n <- 1:nrow(tab_multi)
    vif <- sapply(1:imp$m, function(i) 
      vif(lm(as.formula(fml), complete(imp, i))))
    colnames(vif) <- paste0("vif_imp", 1:ncol(vif))
    vif <- cbind(variable = rownames(vif), as.data.frame(vif))
    tab <- merge(tab_uni, tab_multi, by = "variable", all = TRUE)
    tab <- merge(tab, vif, by = "variable", all = TRUE)
    tab <- tab[order(tab$dummy_n), ]
    tab$dummy_n <- NULL
    tab
  })
})
write_xlsx(reg_tbl[[1]], file.path(outdir, "regressions_all_predictors.xlsx"))
write_xlsx(reg_tbl[[2]],
           file.path(outdir, "regressions_selected_predictors.xlsx"))
