load('tests/fig_log/roc_list_original.Rdata')
load('tests/plsr_imputation/roc.list.plsr.Rdata')

fetch_auc <- function(x, name){
  sapply(name, function(i) {
    x[[i]]$auc
  })
}

## ensure name pairing
auc_original <- fetch_auc(roc.list.original, names(roc.list.original))
auc_plsr     <- fetch_auc(roc.list.plsr, names(roc.list.original))

auc_means <- format(as.numeric(sapply(list(auc_original, auc_plsr), mean)), digits = 3)

plsr_lasso_wilcox <- wilcox.test(auc_original, auc_plsr)

pdf('tests/plsr_imputation/plsr_vs_lasso.pdf', 5, 7)
vioplot(auc_original, auc_plsr,
        names = c("Lasso", "PLSR"),
        col = 'tomato', rectCol = 'grey', ylim = c(0, 1))
text(x = 1:2, y = as.numeric(auc_means), auc_means)
text(x = 1.5, 0.5,
     paste0("pval: ", format(plsr_lasso_wilcox$p.value, digits = 3))
)
dev.off()