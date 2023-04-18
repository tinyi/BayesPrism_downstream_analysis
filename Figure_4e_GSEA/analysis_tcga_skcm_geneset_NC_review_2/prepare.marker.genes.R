#work under /fscratch/tinyi/Ted/NC_review/analysis_embedding_phi/analysis_tcga_skcm_geneset_NC_review_2


#prepare marker genes
#marker from the melanoma scRNA-seq paper

mitf.axl.list <- as.list(read.table("mitf_axl.txt",header=T))

tcell.exclusion.list <- as.list(read.table("Tcell_exclusion.txt",header=T,sep="\t"))
tcell.exclusion.list  <- lapply(tcell.exclusion.list,function(i) i[i!=""])


merged.list <- c(mitf.axl.list, tcell.exclusion.list)

save(merged.list , file="marker.gene.rdata")





