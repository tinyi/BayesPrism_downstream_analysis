library(fgsea)
library(data.table)

library("parallel")
library("tidyverse")


do.fgsea <- function(stat.mat,
					 pathway.table,
					 pval.cut=0.001,
					 NES.cut=0, 
					 n.cores=10,
					 topN=NULL,
					 prefix=NULL){
					 	
	fgseaRes.list <- mclapply(1:ncol(stat.mat), function(idx) {
		fgseaRes <- fgsea(pathways= pathway.table, stats= stat.mat[, idx],eps=0)
		fgseaResTidy <- fgseaRes %>%
  						as_tibble() %>%
  						filter(padj < pval.cut) %>%
  						filter(abs(NES) > NES.cut) %>%
  						arrange(desc(NES))
 		
 		if(!is.null(topN)) {
 			NES.rank <- rank(fgseaResTidy$NES)
 			selected.idx.pos <- NES.rank>(length(NES.rank)-topN) & fgseaResTidy$NES>0
 			selected.idx.neg <- NES.rank<= topN  & fgseaResTidy$NES<0
 			selected.idx <- selected.idx.pos | selected.idx.neg
 			fgseaResTidy <- fgseaResTidy [selected.idx,]
 		}
 		fgseaResTidy
	},mc.cores= n.cores)

	names(fgseaRes.list) <- paste(prefix, colnames(stat.mat),sep="-")
	fgseaRes.list
}


plot.fgsea <- function(fgseaRes.list, pdf.name.merged, col.lim=10){
	for(i in 1:length(fgseaRes.list)){
		if(i<10) pdf.name <- paste("tmp.0",i,".pdf",sep="")
		else pdf.name <- paste("tmp.",i,".pdf",sep="")
				
		fgseaResTidy.i <- fgseaRes.list[[i]]
		
		minus.log10.padj <- -log10(fgseaResTidy.i$padj)
		minus.log10.padj[minus.log10.padj> col.lim] <- col.lim
		signed.minus.log10.padj <- minus.log10.padj
		signed.minus.log10.padj[fgseaResTidy.i$NES<0] <- - signed.minus.log10.padj[fgseaResTidy.i$NES<0]
		
		fgseaResTidy.i$signed.minus.log10.padj <- signed.minus.log10.padj
		
		p <- ggplot(fgseaResTidy.i, aes(reorder(pathway, NES), NES)) +
  		geom_col(aes(fill= signed.minus.log10.padj), colour="gray") +
  		#scale_color_gradientn(colours = c("blue", "white", "red"), 
  		#					  limits =c(-10, 10), breaks=c(-10,-5,0,5,10)) +\
  		scale_fill_gradientn("-log10(padj)",colours = c("blue", "white", "red"), 
  							 limits =c(-10, 10), breaks=c(-col.lim,-col.lim/2,0, col.lim/2, col.lim), 
  							 					 labels=c(col.lim,col.lim/2,"0",col.lim/2,col.lim))  + 
  		coord_flip() +
  		labs(x="Pathway", y="Normalized Enrichment Score",
       		title=paste( names(fgseaRes.list)[i], " Hallmark & GO pathways NES from GSEA",sep=": ")) + 
  		theme_minimal()
  		
		ggsave(filename = pdf.name, plot=p,width=16,height=9)

	}
	system(paste("pdfjam", " ./tmp*.pdf ",  "--nup 3x4 --landscape --outfile ", pdf.name.merged, sep="")) 
		#system(paste("pdfjam", " ./tmp*.pdf ",  "--nup ", ncol.plot, "x", nrow.plot, " --landscape  --outfile ", pdf.name.merged, sep="")) # needs to install pdfjam
		system("rm ./tmp*.pdf")	
	NULL
}


pathways.hallmark <- gmtPathways("/home/chut/data/MSigDB/v7.4/h.all.v7.4.symbols.gmt.txt")
pathways.GO <- gmtPathways("/home/chut/data/MSigDB/v7.4/c5.go.bp.v7.4.symbols.gmt.txt")
pathways.combined <- c(pathways.hallmark, pathways.GO)


#new filter
load("ted.cor.obj.noMerge.new.gbm8only.rdata")

#mask cor based on p

mask.out.cor.using.p <- function(cor.mat, p.mat, p.cut){
	stopifnot(prod(colnames(cor.mat)==colnames(p.mat))==1)
	stopifnot(prod(rownames(cor.mat)== rownames(p.mat))==1)


	masked.cor <- cor.mat
	
	for(i in 1:ncol(masked.cor)){
		masked.idx <- p.mat[,i]> p.cut & !is.na(p.mat[,i])
		masked.cor[masked.idx,i] <- NA
	}
	masked.cor
}


p.cut<-0.01

masked.cor.skcm <- mask.out.cor.using.p(three.tumor.list.noMerge[[1]]$cor, three.tumor.list.noMerge[[1]]$pval, p.cut)
masked.cor.hnscc <- mask.out.cor.using.p(three.tumor.list.noMerge[[2]]$cor, three.tumor.list.noMerge[[2]]$pval, p.cut)
masked.cor.gbm <- mask.out.cor.using.p(three.tumor.list.noMerge[[3]]$cor, three.tumor.list.noMerge[[3]]$pval, p.cut)



pval.cut<-0.01
topN<-10
fgsea.skcm<-do.fgsea(stat.mat = masked.cor.skcm, pathway.table= pathways.combined, prefix= "SKCM", pval.cut= pval.cut,topN= topN)
fgsea.hnscc<-do.fgsea(stat.mat = masked.cor.hnscc, pathway.table= pathways.combined, prefix= "HNSCC", pval.cut= pval.cut,topN= topN)
fgsea.gbm<-do.fgsea(stat.mat = masked.cor.gbm, pathway.table= pathways.combined, prefix= "GBM", pval.cut= pval.cut,topN= topN)


plot.fgsea(fgsea.skcm, pdf.name.merged="SKCM.fgsea.maskedCor.top10.p0.01.pdf")
plot.fgsea(fgsea.hnscc, pdf.name.merged="HSNCC.fgsea.maskedCor.top10.p0.01.pdf")
plot.fgsea(fgsea.gbm, pdf.name.merged="GBM.fgsea.maskedCor.top10.p0.01.pdf")



pval.cut<-0.01
topN<-20
fgsea.skcm<-do.fgsea(stat.mat = masked.cor.skcm, pathway.table= pathways.combined, prefix= "SKCM", pval.cut= pval.cut,topN= topN)
fgsea.hnscc<-do.fgsea(stat.mat = masked.cor.hnscc, pathway.table= pathways.combined, prefix= "HNSCC", pval.cut= pval.cut,topN= topN)
fgsea.gbm<-do.fgsea(stat.mat = masked.cor.gbm, pathway.table= pathways.combined, prefix= "GBM", pval.cut= pval.cut,topN= topN)


plot.fgsea(fgsea.skcm, pdf.name.merged="SKCM.fgsea.maskedCor.top20.p0.01.pdf")
plot.fgsea(fgsea.hnscc, pdf.name.merged="HSNCC.fgsea.maskedCor.top20.p0.01.pdf")
plot.fgsea(fgsea.gbm, pdf.name.merged="GBM.fgsea.maskedCor.top20.p0.01.pdf")























load("ted.cor.obj.noMerge.rdata")

# fgsea.ls <- list()
# for(i in 1:length(three.tumor.list.noMerge)){
	# dat.list.i <- three.tumor.list.noMerge[i]
	# tumor.name <- names(three.tumor.list.noMerge)[i]
	# fgsea.ls [[i]]<-do.fgsea(stat.mat =three.tumor.list.noMerge[[i]]$cor, pathway.table= pathways.combined, prefix= tumor.name)
# }	
# names(fgsea.ls) <- names(three.tumor.list.noMerge)

pval.cut<-1
NES.cut <- 0
fgsea.skcm<-do.fgsea(stat.mat =three.tumor.list.noMerge[[1]]$cor, pathway.table= pathways.combined, prefix= "SKCM", pval.cut= pval.cut,NES.cut=0)
fgsea.hnscc<-do.fgsea(stat.mat =three.tumor.list.noMerge[[2]]$cor, pathway.table= pathways.combined, prefix= "HNSCC", pval.cut= pval.cut,NES.cut=0)
fgsea.gbm<-do.fgsea(stat.mat =three.tumor.list.noMerge[[3]]$cor, pathway.table= pathways.combined, prefix= "GBM", pval.cut= pval.cut,NES.cut=0)

plot.fgsea(fgsea.skcm, pdf.name.merged="SKCM.fgsea.pdf")
plot.fgsea(fgsea.hnscc, pdf.name.merged="HSNCC.fgsea.pdf")
plot.fgsea(fgsea.gbm, pdf.name.merged="GBM.fgsea.pdf")



# # fgsea.skcm<-do.fgsea(stat.mat =three.tumor.list.noMerge[[1]]$cor, pathway.table= pathways.combined, prefix= "SKCM", pval.cut=0.001,NES.cut=1)
# fgsea.hnscc<-do.fgsea(stat.mat =three.tumor.list.noMerge[[2]]$cor, pathway.table= pathways.combined, prefix= "HNSCC", pval.cut=0.001,NES.cut=1)
# fgsea.gbm<-do.fgsea(stat.mat =three.tumor.list.noMerge[[3]]$cor, pathway.table= pathways.combined, prefix= "GBM", pval.cut=0.001,NES.cut=1)

# plot.fgsea(fgsea.skcm, pdf.name.merged="SKCM.fgsea.NESCut1.pdf")
# plot.fgsea(fgsea.hnscc, pdf.name.merged="HSNCC.fgsea.NESCut1.pdf")
# plot.fgsea(fgsea.gbm, pdf.name.merged="GBM.fgsea.NESCut1.pdf")



# fgsea.skcm<-do.fgsea(stat.mat =three.tumor.list.noMerge[[1]]$cor, pathway.table= pathways.combined, prefix= "SKCM", pval.cut=0.001,NES.cut=2)
# fgsea.hnscc<-do.fgsea(stat.mat =three.tumor.list.noMerge[[2]]$cor, pathway.table= pathways.combined, prefix= "HNSCC", pval.cut=0.001,NES.cut=2)
# fgsea.gbm<-do.fgsea(stat.mat =three.tumor.list.noMerge[[3]]$cor, pathway.table= pathways.combined, prefix= "GBM", pval.cut=0.001,NES.cut=2)

# plot.fgsea(fgsea.skcm, pdf.name.merged="SKCM.fgsea.NESCut2.pdf")
# plot.fgsea(fgsea.hnscc, pdf.name.merged="HSNCC.fgsea.NESCut2.pdf")
# plot.fgsea(fgsea.gbm, pdf.name.merged="GBM.fgsea.NESCut2.pdf")



pval.cut<-1
NES.cut <- 0
fgsea.skcm<-do.fgsea(stat.mat =three.tumor.list.noMerge[[1]]$cor, pathway.table= pathways.combined, prefix= "SKCM", pval.cut= pval.cut,topN=20)
fgsea.hnscc<-do.fgsea(stat.mat =three.tumor.list.noMerge[[2]]$cor, pathway.table= pathways.combined, prefix= "HNSCC", pval.cut= pval.cut,topN=20)
fgsea.gbm<-do.fgsea(stat.mat =three.tumor.list.noMerge[[3]]$cor, pathway.table= pathways.combined, prefix= "GBM", pval.cut= pval.cut,topN=20)

plot.fgsea(fgsea.skcm, pdf.name.merged="SKCM.fgsea.top20.pdf")
plot.fgsea(fgsea.hnscc, pdf.name.merged="HSNCC.fgsea.top20.pdf")
plot.fgsea(fgsea.gbm, pdf.name.merged="GBM.fgsea.top20.pdf")




pval.cut<-1
NES.cut <- 0
fgsea.skcm<-do.fgsea(stat.mat =three.tumor.list.noMerge[[1]]$cor, pathway.table= pathways.combined, prefix= "SKCM", pval.cut= pval.cut,topN=10)
fgsea.hnscc<-do.fgsea(stat.mat =three.tumor.list.noMerge[[2]]$cor, pathway.table= pathways.combined, prefix= "HNSCC", pval.cut= pval.cut,topN=10)
fgsea.gbm<-do.fgsea(stat.mat =three.tumor.list.noMerge[[3]]$cor, pathway.table= pathways.combined, prefix= "GBM", pval.cut= pval.cut,topN=10)

plot.fgsea(fgsea.skcm, pdf.name.merged="SKCM.fgsea.top10.pdf")
plot.fgsea(fgsea.hnscc, pdf.name.merged="HSNCC.fgsea.top10.pdf")
plot.fgsea(fgsea.gbm, pdf.name.merged="GBM.fgsea.top10.pdf")



pval.cut<-0.01
NES.cut <- 0
fgsea.skcm<-do.fgsea(stat.mat =three.tumor.list.noMerge[[1]]$cor, pathway.table= pathways.combined, prefix= "SKCM", pval.cut= pval.cut,topN=20)
fgsea.hnscc<-do.fgsea(stat.mat =three.tumor.list.noMerge[[2]]$cor, pathway.table= pathways.combined, prefix= "HNSCC", pval.cut= pval.cut,topN=20)
fgsea.gbm<-do.fgsea(stat.mat =three.tumor.list.noMerge[[3]]$cor, pathway.table= pathways.combined, prefix= "GBM", pval.cut= pval.cut,topN=20)

plot.fgsea(fgsea.skcm, pdf.name.merged="SKCM.fgsea.top20.p0.01.pdf")
plot.fgsea(fgsea.hnscc, pdf.name.merged="HSNCC.fgsea.top20.p0.01.pdf")
plot.fgsea(fgsea.gbm, pdf.name.merged="GBM.fgsea.top20.p0.01.pdf")


fgsea.skcm<-do.fgsea(stat.mat =three.tumor.list.noMerge[[1]]$cor, pathway.table= pathways.combined, prefix= "SKCM", pval.cut= pval.cut,topN=10)
fgsea.hnscc<-do.fgsea(stat.mat =three.tumor.list.noMerge[[2]]$cor, pathway.table= pathways.combined, prefix= "HNSCC", pval.cut= pval.cut,topN=10)
fgsea.gbm<-do.fgsea(stat.mat =three.tumor.list.noMerge[[3]]$cor, pathway.table= pathways.combined, prefix= "GBM", pval.cut= pval.cut,topN=10)

plot.fgsea(fgsea.skcm, pdf.name.merged="SKCM.fgsea.top10.p0.01.pdf")
plot.fgsea(fgsea.hnscc, pdf.name.merged="HSNCC.fgsea.top10.p0.01.pdf")
plot.fgsea(fgsea.gbm, pdf.name.merged="GBM.fgsea.top10.p0.01.pdf")





#try signed p 
# seems that too many terms are affected by outliers (because p is calculated from raw value, spearman is more robust)

convert.stat <- function(cor.mat, p.mat){
	stopifnot(prod(colnames(cor.mat)==colnames(p.mat))==1)
	stopifnot(prod(rownames(cor.mat)== rownames(p.mat))==1)

	minus.log10.p <- -log10(p.mat)

	signed.p <- matrix(NA,nrow=nrow(minus.log10.p),ncol=ncol(minus.log10.p), dimnames=dimnames(minus.log10.p))
	
	for(i in 1:ncol(signed.p)){
		pos.idx <- cor.mat[,i]>0 & !is.na(cor.mat[,i])
		neg.idx <- cor.mat[,i]<0 & !is.na(cor.mat[,i])
		signed.p[pos.idx,i] <-  minus.log10.p[pos.idx,i]
		signed.p[neg.idx,i] <-  - minus.log10.p[neg.idx,i] #flip the sign for negative cor
	}
	signed.p
}

signed.p.skcm <- convert.stat(three.tumor.list.noMerge[[1]]$cor, three.tumor.list.noMerge[[1]]$pval)
signed.p.hnscc <- convert.stat(three.tumor.list.noMerge[[2]]$cor, three.tumor.list.noMerge[[2]]$pval)
signed.p.gbm <- convert.stat(three.tumor.list.noMerge[[3]]$cor, three.tumor.list.noMerge[[3]]$pval)

pval.cut<-0.01
topN<-20
fgsea.skcm<-do.fgsea(stat.mat = signed.p.skcm, pathway.table= pathways.combined, prefix= "SKCM", pval.cut= pval.cut,topN= topN)
fgsea.hnscc<-do.fgsea(stat.mat = signed.p.hnscc, pathway.table= pathways.combined, prefix= "HNSCC", pval.cut= pval.cut,topN= topN)
fgsea.gbm<-do.fgsea(stat.mat = signed.p.gbm, pathway.table= pathways.combined, prefix= "GBM", pval.cut= pval.cut,topN= topN)


plot.fgsea(fgsea.skcm, pdf.name.merged="SKCM.fgsea.pstat.top20.p0.01.pdf")
plot.fgsea(fgsea.hnscc, pdf.name.merged="HSNCC.fgsea.pstat.top20.p0.01.pdf")
plot.fgsea(fgsea.gbm, pdf.name.merged="GBM.fgsea.pstat.top20.p0.01.pdf")




#mask cor based on p

mask.out.cor.using.p <- function(cor.mat, p.mat, p.cut){
	stopifnot(prod(colnames(cor.mat)==colnames(p.mat))==1)
	stopifnot(prod(rownames(cor.mat)== rownames(p.mat))==1)


	masked.cor <- cor.mat
	
	for(i in 1:ncol(masked.cor)){
		masked.idx <- p.mat[,i]> p.cut & !is.na(p.mat[,i])
		masked.cor[masked.idx,i] <- NA
	}
	masked.cor
}


p.cut<-0.01

masked.cor.skcm <- mask.out.cor.using.p(three.tumor.list.noMerge[[1]]$cor, three.tumor.list.noMerge[[1]]$pval, p.cut)
masked.cor.hnscc <- mask.out.cor.using.p(three.tumor.list.noMerge[[2]]$cor, three.tumor.list.noMerge[[2]]$pval, p.cut)
masked.cor.gbm <- mask.out.cor.using.p(three.tumor.list.noMerge[[3]]$cor, three.tumor.list.noMerge[[3]]$pval, p.cut)



pval.cut<-0.01
topN<-10
fgsea.skcm<-do.fgsea(stat.mat = masked.cor.skcm, pathway.table= pathways.combined, prefix= "SKCM", pval.cut= pval.cut,topN= topN)
fgsea.hnscc<-do.fgsea(stat.mat = masked.cor.hnscc, pathway.table= pathways.combined, prefix= "HNSCC", pval.cut= pval.cut,topN= topN)
fgsea.gbm<-do.fgsea(stat.mat = masked.cor.gbm, pathway.table= pathways.combined, prefix= "GBM", pval.cut= pval.cut,topN= topN)


plot.fgsea(fgsea.skcm, pdf.name.merged="SKCM.fgsea.maskedCor.top10.p0.01.pdf")
plot.fgsea(fgsea.hnscc, pdf.name.merged="HSNCC.fgsea.maskedCor.top10.p0.01.pdf")
plot.fgsea(fgsea.gbm, pdf.name.merged="GBM.fgsea.maskedCor.top10.p0.01.pdf")



pval.cut<-0.01
topN<-20
fgsea.skcm<-do.fgsea(stat.mat = masked.cor.skcm, pathway.table= pathways.combined, prefix= "SKCM", pval.cut= pval.cut,topN= topN)
fgsea.hnscc<-do.fgsea(stat.mat = masked.cor.hnscc, pathway.table= pathways.combined, prefix= "HNSCC", pval.cut= pval.cut,topN= topN)
fgsea.gbm<-do.fgsea(stat.mat = masked.cor.gbm, pathway.table= pathways.combined, prefix= "GBM", pval.cut= pval.cut,topN= topN)


plot.fgsea(fgsea.skcm, pdf.name.merged="SKCM.fgsea.maskedCor.top20.p0.01.pdf")
plot.fgsea(fgsea.hnscc, pdf.name.merged="HSNCC.fgsea.maskedCor.top20.p0.01.pdf")
plot.fgsea(fgsea.gbm, pdf.name.merged="GBM.fgsea.maskedCor.top20.p0.01.pdf")















#select (remove sementically redundant for visuliazation)

topN<-20
fgsea.skcm<-do.fgsea(stat.mat = masked.cor.skcm, pathway.table= pathways.combined, prefix= "SKCM", pval.cut= pval.cut,topN= topN)
fgsea.hnscc<-do.fgsea(stat.mat = masked.cor.hnscc, pathway.table= pathways.combined, prefix= "HNSCC", pval.cut= pval.cut,topN= topN)
fgsea.gbm<-do.fgsea(stat.mat = masked.cor.gbm, pathway.table= pathways.combined, prefix= "GBM", pval.cut= pval.cut,topN= topN)

save(fgsea.gbm, fgsea.hnscc, fgsea.skcm, file="fgsea.top20.p0.01.rdata")

gbm.macro <- fgsea.gbm[[1]]
> gbm.macro$pathway
 # [1] "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
 # [2] "HALLMARK_INTERFERON_GAMMA_RESPONSE"
 # [3] "GOBP_HUMORAL_IMMUNE_RESPONSE"
 # [4] "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
 # [5] "HALLMARK_INTERFERON_ALPHA_RESPONSE"
 # [6] "HALLMARK_INFLAMMATORY_RESPONSE"
 # [7] "HALLMARK_COAGULATION"
 # [8] "GOBP_REGULATION_OF_HUMORAL_IMMUNE_RESPONSE"
 # [9] "GOBP_INFLAMMATORY_RESPONSE"
# [10] "GOBP_ACUTE_PHASE_RESPONSE"
# [11] "GOBP_REGULATION_OF_COMPLEMENT_ACTIVATION"
# [12] "GOBP_POSITIVE_REGULATION_OF_INFLAMMATORY_RESPONSE"
# [13] "GOBP_COMPLEMENT_ACTIVATION"
# [14] "GOBP_ANTIMICROBIAL_HUMORAL_RESPONSE"
# [15] "GOBP_ACUTE_INFLAMMATORY_RESPONSE"
# [16] "GOBP_EXTERNAL_ENCAPSULATING_STRUCTURE_ORGANIZATION"
# [17] "HALLMARK_IL6_JAK_STAT3_SIGNALING"
# [18] "GOBP_MYELOID_LEUKOCYTE_MEDIATED_IMMUNITY"
# [19] "GOBP_DEFENSE_RESPONSE"
# [20] "GOBP_MYELOID_LEUKOCYTE_ACTIVATION"
# [21] "GOBP_HISTONE_METHYLATION"
# [22] "GOBP_MEIOTIC_CELL_CYCLE"
# [23] "GOBP_COVALENT_CHROMATIN_MODIFICATION"
# [24] "GOBP_SPINDLE_ORGANIZATION"
# [25] "GOBP_DNA_REPLICATION"
# [26] "GOBP_NUCLEOSOME_ORGANIZATION"
# [27] "GOBP_RECOMBINATIONAL_REPAIR"
# [28] "GOBP_CELL_CYCLE_CHECKPOINT"
# [29] "GOBP_CHROMATIN_ASSEMBLY_OR_DISASSEMBLY"
# [30] "GOBP_REGULATION_OF_GENE_EXPRESSION_EPIGENETIC"
# [31] "GOBP_NUCLEAR_CHROMOSOME_SEGREGATION"
# [32] "GOBP_MITOTIC_SISTER_CHROMATID_SEGREGATION"
# [33] "GOBP_CHROMOSOME_SEGREGATION"
# [34] "GOBP_CHROMATIN_ORGANIZATION"
# [35] "HALLMARK_E2F_TARGETS"
# [36] "GOBP_DNA_CONFORMATION_CHANGE"
# [37] "GOBP_SISTER_CHROMATID_SEGREGATION"
# [38] "HALLMARK_G2M_CHECKPOINT"
# [39] "GOBP_DNA_PACKAGING"
# [40] "GOBP_CHROMOSOME_ORGANIZATION"


gbm.macro$pathway[gbm.macro$NES<0]

selected.idx <- c(1,2,3,4,5,6,7,11,16,17,18,  40,39,38,35,30,28,27,26,25:21)
gbm.macro <- gbm.macro[selected.idx,]
fgsea.gbm[[1]] <- gbm.macro
plot.fgsea(fgsea.gbm, pdf.name.merged="GBM.fgsea.selectedMyeloid.pdf")


hnscc.macro <- fgsea.hnscc[[5]]
hnscc.macro$pathway
 # [1] "HALLMARK_INTERFERON_GAMMA_RESPONSE"
 # [2] "HALLMARK_INTERFERON_ALPHA_RESPONSE"
 # [3] "GOBP_INNATE_IMMUNE_RESPONSE"
 # [4] "GOBP_DEFENSE_RESPONSE_TO_OTHER_ORGANISM"
 # [5] "GOBP_DEFENSE_RESPONSE_TO_VIRUS"
 # [6] "GOBP_RESPONSE_TO_TYPE_I_INTERFERON"
 # [7] "GOBP_RESPONSE_TO_VIRUS"
 # [8] "GOBP_DEFENSE_RESPONSE"
 # [9] "GOBP_RESPONSE_TO_BIOTIC_STIMULUS"
# [10] "GOBP_RESPONSE_TO_INTERFERON_GAMMA"
# [11] "GOBP_NEGATIVE_REGULATION_OF_VIRAL_PROCESS"
# [12] "GOBP_POSITIVE_REGULATION_OF_IMMUNE_RESPONSE"
# [13] "GOBP_REGULATION_OF_RESPONSE_TO_BIOTIC_STIMULUS"
# [14] "GOBP_IMMUNE_EFFECTOR_PROCESS"
# [15] "GOBP_REGULATION_OF_INNATE_IMMUNE_RESPONSE"
# [16] "GOBP_CYTOKINE_MEDIATED_SIGNALING_PATHWAY"
# [17] "GOBP_REGULATION_OF_DEFENSE_RESPONSE"
# [18] "GOBP_REGULATION_OF_IMMUNE_RESPONSE"
# [19] "GOBP_POSITIVE_REGULATION_OF_RESPONSE_TO_BIOTIC_STIMULUS"
# [20] "GOBP_NEGATIVE_REGULATION_OF_VIRAL_GENOME_REPLICATION"
# [21] "GOBP_EPIDERMIS_DEVELOPMENT"
# [22] "GOBP_CELLULAR_RESPONSE_TO_INSULIN_STIMULUS"

selected.idx <- c(1,2,3,4,5,6,9,14,16,21,22)

hnscc.macro <- hnscc.macro[selected.idx,]
fgsea.hnscc.macro.edit <- fgsea.hnscc
fgsea.hnscc.macro.edit[[5]] <- hnscc.macro
plot.fgsea(fgsea.hnscc, pdf.name.merged="HNSCC.fgsea.selectedMyeloid.pdf")



skcm.macro <- fgsea.skcm[[1]]
skcm.macro $pathway
 # [1] "GOBP_RESPONSE_TO_TUMOR_NECROSIS_FACTOR"
 # [2] "GOBP_CYTOKINE_MEDIATED_SIGNALING_PATHWAY"
 # [3] "HALLMARK_INTERFERON_GAMMA_RESPONSE"
 # [4] "HALLMARK_HYPOXIA"
 # [5] "HALLMARK_COMPLEMENT"
 # [6] "GOBP_RESPONSE_TO_CYTOKINE"
 # [7] "HALLMARK_INTERFERON_ALPHA_RESPONSE"
 # [8] "GOBP_HUMORAL_IMMUNE_RESPONSE"
 # [9] "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
# [10] "HALLMARK_COAGULATION"
# [11] "GOBP_TUMOR_NECROSIS_FACTOR_MEDIATED_SIGNALING_PATHWAY"
# [12] "GOBP_REGULATION_OF_IMMUNE_EFFECTOR_PROCESS"
# [13] "GOBP_RESPONSE_TO_INTERLEUKIN_1"
# [14] "GOBP_LEUKOCYTE_MEDIATED_IMMUNITY"
# [15] "GOBP_MYELOID_DENDRITIC_CELL_ACTIVATION"
# [16] "GOBP_MYELOID_LEUKOCYTE_ACTIVATION"
# [17] "GOBP_EXOCYTOSIS"
# [18] "GOBP_REGULATION_OF_PEPTIDASE_ACTIVITY"
# [19] "GOBP_SECRETION"
# [20] "GOBP_DEFENSE_RESPONSE"
# [21] "GOBP_DNA_CONFORMATION_CHANGE"
# [22] "GOBP_ESTABLISHMENT_OF_RNA_LOCALIZATION"
# [23] "GOBP_METHYLATION"
# [24] "GOBP_REGULATION_OF_MRNA_PROCESSING"
# [25] "GOBP_MRNA_3_END_PROCESSING"
# [26] "GOBP_RIBOSOME_BIOGENESIS"
# [27] "GOBP_TRNA_METABOLIC_PROCESS"
# [28] "GOBP_RNA_3_END_PROCESSING"
# [29] "GOBP_CHROMOSOME_ORGANIZATION"
# [30] "GOBP_RNA_LOCALIZATION"
# [31] "GOBP_MACROMOLECULE_METHYLATION"
# [32] "GOBP_NCRNA_PROCESSING"
# [33] "GOBP_NCRNA_METABOLIC_PROCESS"
# [34] "GOBP_RIBONUCLEOPROTEIN_COMPLEX_SUBUNIT_ORGANIZATION"
# [35] "GOBP_RIBONUCLEOPROTEIN_COMPLEX_BIOGENESIS"
# [36] "GOBP_RNA_SPLICING"
# [37] "GOBP_RNA_SPLICING_VIA_TRANSESTERIFICATION_REACTIONS"
# [38] "GOBP_MRNA_METABOLIC_PROCESS"
# [39] "GOBP_RNA_PROCESSING"
# [40] "GOBP_MRNA_PROCESSING"
selected.idx <- - c(6,11,40,37,34,24,22,23)

skcm.macro <- skcm.macro[selected.idx,]

fgsea.skcm.macro.edit <- fgsea.skcm

fgsea.skcm.macro.edit[[1]] <- skcm.macro
plot.fgsea(fgsea.skcm.macro.edit, pdf.name.merged="SKCM.fgsea.selectedMyeloid.pdf")






















#write to table

fgsea.skcm.all<-do.fgsea(stat.mat = masked.cor.skcm, pathway.table= pathways.combined, prefix= "SKCM", pval.cut= pval.cut,topN= 9999)
fgsea.hnscc.all<-do.fgsea(stat.mat = masked.cor.hnscc, pathway.table= pathways.combined, prefix= "HNSCC", pval.cut= pval.cut,topN= 9999)
fgsea.gbm.all<-do.fgsea(stat.mat = masked.cor.gbm, pathway.table= pathways.combined, prefix= "GBM", pval.cut= pval.cut,topN= 9999)


write.to.table <- function(go.tab.list, prefix, topN){
	
	if(is.null(names(go.tab.list))) names(go.tab.list) <- paste(prefix, "-Tumor-",1:length(go.tab.list),sep="")
	
	lapply(1:length(go.tab.list), function(idx){
		write.table(go.tab.list[[idx]][1:min(nrow(go.tab.list[[idx]]), topN),1:7], file=paste(prefix, "_",names(go.tab.list)[idx],".txt",sep=""),
					sep="\t",quote=F,col.names=T,row.names=F)
	})
	NULL
}


write.to.table(fgsea.skcm.all, prefix="", topN=9999)
write.to.table(fgsea.hnscc.all, prefix="", topN=9999)
write.to.table(fgsea.gbm.all, prefix="", topN=9999)









