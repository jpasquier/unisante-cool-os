library(broom)
library(dplyr)
library(emmeans)
library(ggplot2)
library(ggpubr)
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

# Adjusted comparisons of the numeric variables
Y <- c("L1L4s_BMD", "L1L4sTBS", "femurTot_DMO", "col_DMO")
X <- list("DepTabac", "TTT_Hormmenop", c("DepTabac", "TTT_Hormmenop"))
adj_tab_num <- mclapply(setNames(Y, Y), function(y) {
  R <- mclapply(X, function(x) {
    subdta <- dta %>%
      filter(Subject_ID %in% ids(y)) %>%
      select(Cohort, any_of(c(x, y)))
    if (any(is.na(subdta))) stop("missing value(s)")
    fits <- lapply(1:2, function(k) {
      fml <- as.formula(paste0(
        y, "~", paste(c("Cohort", x), collapse = c("+", "*")[k])))
      lm(fml, subdta)
    })
    tab <- lapply(fits, function(fit) {
      tidy(fit, conf.int = TRUE) %>%
        select(term, estimate, conf.low, conf.high, p.value)
    }) %>%
      {full_join(.[[1]], .[[2]], by = "term", suffix = c("", ".I"))}
    na <- rep(NA, nrow(tab) - 1)
    tab <- cbind(outcome = c(y, na),
                 adjustment = c(paste(x, collapse = "/"), na), tab)
    sttl <- paste("Adjusted for", paste(x, collapse = " and "))
    d <- do.call(bind_rows, lapply(1:2, function(k) {
      fit <- fits[[k]]
      tidy(emmeans(fit, update(formula(fit), NULL ~ .))) %>%
        mutate(inter = k)
    })) %>%
      mutate(inter = factor(inter, 1:2,
        paste(c("without", "with"), "interaction")))
    if (length(x) == 2) {
      d$comp <- paste(d[[x[1]]], d[[x[2]]], sep = "/")
      x <- paste(x[1], x[2], sep = "/")
      names(d)[names(d) == "comp"] <- x
    }
    fig <- ggplot(d, aes(x = Cohort, y = estimate, color = !!sym(x),
                         group = !!sym(x))) +
      geom_point(position = position_dodge(width = .1)) +
      geom_line(position = position_dodge(width = .1), linetype = "dashed") +
      geom_errorbar(aes(ymin = estimate - std.error,
                        ymax = estimate + std.error),
                    position = position_dodge(width = .1), width = 0) +
      facet_grid(cols = vars(inter)) +
      labs(x = NULL, y = paste("Mean", y, "± SD"), subtitle = sttl) +
      theme_bw() +
      theme(legend.position = "bottom")
    list(tab = tab, fig = fig)
  })
  tab <- do.call(bind_rows, lapply(R, function(r) r$tab))
  fig <- ggarrange(plotlist = lapply(R, function(r) r$fig), ncol = 1) %>%
    annotate_figure(top = text_grob(y, face = "bold", size = 14))
  attr(tab, "fig") <- fig
  return(tab)
})
write_xlsx(adj_tab_num,
           file.path(outdir, "adjusted_comparisons_numeric_variables.xlsx"))
rm(X, Y)

# Adjusted comparisons of the numeric variables - Figures
pdf(file.path(outdir, "adjusted_comparisons_numeric_variables.pdf"))
for (z in adj_tab_num) print(attr(z, "fig"))
dev.off()
rm(z)

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
