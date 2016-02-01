load('tests/lasso/roc.lasso.list.RData')
load('tests/plsr_imputation/roc.plsr.list.RData')
load('tests/tlasso_imputation/roc.tlasso.list.RData')
fetch_auc <- function(x, name){
  sapply(name, function(i) {
    x[[i]]$auc
  })
}

## ensure name pairing
auc_original <- fetch_auc(roc.lasso.list, names(roc.lasso.list))
auc_plsr     <- fetch_auc(roc.plsr.list, names(roc.plsr.list))
auc_tlasso   <- fetch_auc(roc.tlasso.list, names(roc.tlasso.list))

auc_medians <- format(as.numeric(sapply(list(auc_original, auc_plsr, auc_tlasso), median)), digits = 3)

plsr_lasso_wilcox <- wilcox.test(auc_original, auc_plsr)
tlasso_lasso_wilcox <- wilcox.test(auc_original, auc_tlasso)

pdf('tests/fig_log/imputation_schemes.pdf', 5, 7)
vioplot(auc_original, auc_plsr, auc_tlasso,
        names = c('Lasso', "PLSR", "Tilling Lasso"),
        col = 'lightblue', rectCol = 'grey', ylim = c(0, 1))
text(1:3, y = as.numeric(auc_medians), auc_medians)
# text(1.5, 0.5,
#      paste0('pval: ', format(plsr_lasso_wilcox$p.value, digits = 3)))
# text(2.5, 0.5,
#      paste0('pval: ', format(tlasso_lasso_wilcox$p.value, digits = 3)))
dev.off()
# pdf('tests/plsr_imputation/plsr_vs_lasso.pdf', 5, 7)
# vioplot(auc_original, auc_plsr,
#         names = c("Lasso", "PLSR"),
#         col = 'tomato', rectCol = 'grey', ylim = c(0, 1))
# text(x = 1:2, y = as.numeric(auc_medians[1:2]), auc_medians[1:2])
# text(x = 1.5, 0.5,
#      paste0("pval: ", format(plsr_lasso_wilcox$p.value, digits = 3))
# )
# dev.off()
#
#
# pdf('tests/tlasso_imputation/tlasso_vs_lasso.pdf', 5, 7)
# vioplot(auc_original, auc_tlasso,
#         names = c("Lasso", "tLasso"),
#         col = 'tomato', rectCol = 'grey', ylim = c(0, 1))
# text(x = 1:2, y = as.numeric(auc_medians[c(1, 3)]), auc_medians[c(1, 3)])
# text(x = 1.5, 0.5,
#      paste0("pval: ", format(tlasso_lasso_wilcox$p.value, digits = 3))
# )
dev.off()
