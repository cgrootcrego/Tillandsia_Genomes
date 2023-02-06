cast_check<-read.table("checklist_with_diff_completeness.txt", header=TRUE)
exon_numbers <-read.table("exon_numbers_both_sp.txt",header=FALSE)
colnames(exon_numbers)=c("og_id","Tfas_gene_id","Tfas_exon_nb","Tlei_gene_id","Tlei_exon_nb")
exon_numbers$exon_nb_diff <- exon_numbers$Tlei_exon_nb-exon_numbers$Tfas_exon_nb

cast_check$Tfas_exon_nb <- exon_numbers$Tfas_exon_nb
cast_check$Tlei_exon_nb <- exon_numbers$Tlei_exon_nb
cast_check$exon_nb_diff <- exon_numbers$exon_nb_diff
cast_check$diff_per_Tfas <- abs(cast_check$diff_per_Tfas)
cast_check$diff_per_Tlei <- abs(cast_check$diff_per_Tlei)

both_complete <- cast_check[which(cast_check$completeness=="both complete"),]
one_complete <- cast_check[which(cast_check$completeness=="one complete"),]

write.table(both_complete, "both_complete.txt", append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE)
write.table(one_complete, "one_complete.txt", append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE)
