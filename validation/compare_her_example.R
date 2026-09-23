# validation/compare_her_example.R
# Phase A step 1: compare the batch reproduction of Volponi's example_Release1 fit
# (validation/matlab_out/her_example/SMIMResults_repro.csv, produced by
# validation/matlab/repro_her_example.m) with the results shipped in her repository
# (SMIMResults.xls and example_Release1.mat). Pass criterion agreed at checkpoint 1:
# parameters within the reported standard errors and objective within 1 %; both reported.
#
#   Rscript validation/compare_her_example.R

suppressPackageStartupMessages(library(dplyr))
her <- "inputs/SabrinaVolponi-SMIMfit-058f87d"
repro <- read.csv("validation/matlab_out/her_example/SMIMResults_repro.csv")
orig <- readxl::read_excel(file.path(her, "SMIMResults.xls"))
info <- readLines("validation/matlab_out/her_example/run_info.txt")
pars <- c("U", "D", "Lambda", "Beta", "logT1", "logT2")
ses <- paste0("SE_", c("U", "D", "Lambda", "B", "logT1", "logT2"))
tab <- data.frame(parameter = pars,
                  original = as.numeric(orig[1, pars]), original_SE = as.numeric(orig[1, ses]),
                  reproduced = as.numeric(repro[1, pars]), reproduced_SE = as.numeric(repro[1, ses]))
tab$diff_in_SE <- abs(tab$reproduced - tab$original) / tab$original_SE
tab$rel_diff <- abs(tab$reproduced - tab$original) / abs(tab$original)
tab$pass <- tab$diff_in_SE <= 1 | (tab$original_SE == 0 & tab$rel_diff < 1e-6)
wm <- data.frame(quantity = c("WMAE (sum of weighted |resid|)", "mass recovery fraction", "Q"),
                 original = as.numeric(orig[1, c("WMAE", "MassRecoveryFractionQAv", "QEstimated")]),
                 reproduced = as.numeric(repro[1, c("WMAE", "MassRecoveryFractionQAv", "QEstimated")]))
wm$rel_diff <- abs(wm$reproduced - wm$original) / abs(wm$original)
cat("Run info:", paste(info, collapse = "; "), "\n\n")
cat("| parameter | original | SE (orig) | reproduced | SE (repro) | |diff| / SE | rel diff | pass |\n|---|---|---|---|---|---|---|---|\n")
for (i in seq_len(nrow(tab))) with(tab[i, ], cat(sprintf("| %s | %.6g | %.3g | %.6g | %.3g | %.3g | %.2e | %s |\n", parameter, original, original_SE, reproduced, reproduced_SE, diff_in_SE, rel_diff, if (pass) "yes" else "no")))
cat("\n| quantity | original | reproduced | rel diff |\n|---|---|---|---|\n")
for (i in seq_len(nrow(wm))) with(wm[i, ], cat(sprintf("| %s | %.6g | %.6g | %.2e |\n", quantity, original, reproduced, rel_diff)))
cat(sprintf("\nObjective (WMAE) within 1 %%: %s; all parameters within 1 SE: %s\n", wm$rel_diff[1] < 0.01, all(tab$pass)))
write.csv(tab, "validation/matlab_out/her_example/comparison.csv", row.names = FALSE)
