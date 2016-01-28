load('tests/fig_log/roc_list_original.Rdata')
load('tests/plsr_imputation/roc.list.plsr.Rdata')
load('tests/tlasso_imputation/roc.list.tlasso.Rdata')
fetch_auc <- function(x, name){
  sapply(name, function(i) {
    x[[i]]$auc
  })
}

## ensure name pairing
auc_original <- fetch_auc(roc.list.original, names(roc.list.original))
auc_plsr     <- fetch_auc(roc.list.plsr, names(roc.list.original))
auc_tlasso   <- fetch_auc(roc.list.tlasso, names(roc.list.original))

auc_medians <- format(as.numeric(sapply(list(auc_original, auc_plsr, auc_tlasso), median)), digits = 3)

plsr_lasso_wilcox <- wilcox.test(auc_original, auc_plsr)
tlasso_lasso_wilcox <- wilcox.test(auc_original, auc_tlasso)

pdf('tests/plsr_imputation/plsr_vs_lasso.pdf', 5, 7)
vioplot(auc_original, auc_plsr,
        names = c("Lasso", "PLSR"),
        col = 'tomato', rectCol = 'grey', ylim = c(0, 1))
text(x = 1:2, y = as.numeric(auc_medians[1:2]), auc_medians[1:2])
text(x = 1.5, 0.5,
     paste0("pval: ", format(plsr_lasso_wilcox$p.value, digits = 3))
)
dev.off()


pdf('tests/tlasso_imputation/tlasso_vs_lasso.pdf', 5, 7)
vioplot(auc_original, auc_tlasso,
        names = c("Lasso", "tLasso"),
        col = 'tomato', rectCol = 'grey', ylim = c(0, 1))
text(x = 1:2, y = as.numeric(auc_medians[c(1, 3)]), auc_medians[c(1, 3)])
text(x = 1.5, 0.5,
     paste0("pval: ", format(tlasso_lasso_wilcox$p.value, digits = 3))
)
dev.off()
