# Studying differences in gene family sizes between Tfas and Tlei

library(ggplot2)
library(dplyr)

setwd("Documents/GitHub/Tillandsia-compgenomics/5. Gene-family-evo-tfas-tlei/")
setwd("/home/clara/Documents/GitHub/Tillandsia-compgenomics/Gene-family-evo-tfas-tlei/")
counts <- read.table("orthogroups_Tfas_Tlei_Acom.counts.with_functional_info.no_TEs.size_corrections.no_plastid-mito-ribo.blastandsearch.txt", sep = '\t')
colnames(counts) <- c("og_id", "Acom", "Tfas", "Tlei")
per_gene <- read.delim("orthogroups_Tfas_Tlei_Acom.per_gene.with_functional_info.no_TEs.size_corrections.no_plastid-mito-ribo.blastandsearch_noAcom.txt",
                       sep = "\t", header = F)

# Filter out unique orthogroups
counts_Tfas_Tlei <- counts[counts$Tfas != 0 & counts$Tlei != 0,]
ggplot(counts_Tfas_Tlei, aes(x=Tfas, y=Tlei)) + geom_point(size = .5) +
  geom_abline(intercept = 0, slope = 1) +
  labs(title = "Per-species gene counts in multi-copy orthogroups") +
  ylab(label = "T. leiboldiana") +
  xlab(label = "T. fasciculata") +
  theme_bw()

# Add colour gradient to show number of datapoints with the same count combination
counts_Tfas_Tlei_multi <- counts_Tfas_Tlei[!(counts_Tfas_Tlei$Tfas == 1 & counts_Tfas_Tlei$Tlei == 1),]
ggplot(counts_Tfas_Tlei_multi) + geom_hex(aes(Tfas, Tlei, fill = stat(log(count))), bins = 100) +
  labs(title = "Per-species gene counts in multi-copy orthogroups") +
  ylab(label = "T. leiboldiana") +
  xlab(label = "T. fasciculata") +
  scale_fill_continuous(type = "viridis") +
  theme_bw()

wilcox.test(counts_Tfas_Tlei_multi$Tfas, counts_Tfas_Tlei_multi$Tlei)
mean(counts_Tfas_Tlei_multi$Tfas)
mean(counts_Tfas_Tlei_multi$Tlei)
IQR(counts_Tfas_Tlei_multi$Tfas)
IQR(counts_Tfas_Tlei_multi$Tlei)

counts_more_Tlei <- counts_Tfas_Tlei_multi[(counts_Tfas_Tlei_multi$Tfas < counts_Tfas_Tlei_multi$Tlei),]
counts_more_Tfas <- counts_Tfas_Tlei_multi[(counts_Tfas_Tlei_multi$Tfas > counts_Tfas_Tlei_multi$Tlei),]

counts_multi_Tfas<- counts_Tfas_Tlei[counts_Tfas_Tlei$Tfas >1,]
counts_multi_Tlei<- counts_Tfas_Tlei[counts_Tfas_Tlei$Tlei >1,]
sum(counts_Tfas_Tlei_multi$Tfas)
sum(counts_Tfas_Tlei_multi$Tlei)

write.table(counts_Tfas_Tlei_multi, file = "orthogroup_selection_multicopy_for_GO_term_all.txt", sep = "\t", quote = F, row.names = F)
write.table(counts_more_Tfas, file = "orthogroup_selection_multicopy_larger_in_Tfas.txt", sep = "\t", quote = F, row.names = F)
write.table(counts_more_Tlei, file = "orthogroup_selection_multicopy_larger_in_Tlei.txt", sep = "\t", quote = F, row.names = F)

### Log-ratio test
counts_Tfas_Tlei$logratio <- log(counts_Tfas_Tlei$Tfas/counts_Tfas_Tlei$Tlei)
mean_logratio <- mean(counts_Tfas_Tlei$logratio) # 0.0203
counts_Tfas_Tlei$corr_logratio <- counts_Tfas_Tlei$logratio - mean_logratio

ggplot(counts_Tfas_Tlei, aes(x=corr_logratio)) + geom_density()
quantile(counts_Tfas_Tlei$corr_logratio)

top2percent_Tfas_larger <- counts_Tfas_Tlei %>%
   arrange(desc(corr_logratio)) %>%
   filter(corr_logratio >= quantile(corr_logratio, .98))
top2percent_Tlei_larger <- counts_Tfas_Tlei %>%
   arrange((corr_logratio)) %>%
   filter(corr_logratio <= quantile(corr_logratio, .02))

top2percent_Tfas_larger_multi <- counts_Tfas_Tlei_multi %>%
   arrange(desc(corr_logratio)) %>%
   filter(corr_logratio >= quantile(corr_logratio, .98))
top2percent_Tlei_larger_multi <- counts_Tfas_Tlei_multi %>%
   arrange((corr_logratio)) %>%
   filter(corr_logratio <= quantile(corr_logratio, .02))

write.table(top2percent_Tfas_larger_multi, file = "Top2percent_multicopy_orthogroups_LogRatio_larger_in_Tfas.txt", sep = "\t", quote = F, row.names = F)
write.table(top2percent_Tlei_larger_multi, file = "Top2percent_multicopy_orthogroups_LogRatio_larger_in_Tlei.txt", sep = "\t", quote = F, row.names = F)


# Because of the very high occurrence of 1:1 orthogroups, anything deviating from
# this will already be in the top 2 % (includes all duplications, also 2:1). This shows that any
# form of copy variance is in fact already "deviating" as far as we can tell.
