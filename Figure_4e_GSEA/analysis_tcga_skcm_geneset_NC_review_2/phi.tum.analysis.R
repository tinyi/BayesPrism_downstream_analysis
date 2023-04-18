library(GSVA)
library(gplots)

load("/fscratch/tinyi/Ted/NC_review/dat/skcm/tcga/tcga.skcm.rdata")
load("/fscratch/tinyi/Ted/NC_review/dat/bp.dat/tcga/Recls/pc/all/tcga.skcm.ted.pcGene.B_Macro_cd4t_cd8tRecls.rdata")
load("/fscratch/tinyi/Ted/NC_review/dat/bp.dat/tcga/NMF_EL/tcga.skcm.ebd.res.rdata")
load("marker.gene.rdata")


ted.res <- tcga.skcm.B_Macro_cd4t_cd8tRecls.pc.ted
ebd.res <- tcga.skcm.ebd.res


pdf.prefix <- "skcm"


phi.hat.tum <- t(ebd.res $ opt.phi.hat.tum)
colnames(phi.hat.tum) <- paste("Program-",1:ncol(phi.hat.tum),sep="")
rownames(phi.hat.tum) <- gene.annot$V5[match(rownames(phi.hat.tum), gene.annot$V4)]
all.score <- gsva(phi.hat.tum, gs= merged.list, verbose=FALSE, parallel.sz=6)


pdf("skcm.subtypesNMF.pdf",pointsize=8,useDingbats=FALSE )
my_palette <- colorRampPalette(c("green", "black", "red"))(n = 254)

heatmap.2(all.score, 
								col=my_palette, 
								density.info="none", 
								trace="none",
								dendrogram='none',
								symm=F, 
								symkey=F, 
								scale="none" , 
								Rowv= NULL, 
								Colv= NULL,
								margins=c(10,10),
								cellnote= round(all.score,2),
     							notecex=1,
     							notecol="white")
     		
dev.off()




#lets include GO and hallmark pathways 

library(fgsea)
library(data.table)


pathways.hallmark <- gmtPathways("/home/chut/data/MSigDB/v7.4/h.all.v7.4.symbols.gmt.txt")
pathways.GO <- gmtPathways("/home/chut/data/MSigDB/v7.4/c5.go.bp.v7.4.symbols.gmt.txt")


gene.files <- list.files(pattern=glob2rx("SKCM_*"),
						path="/fscratch/tinyi/Ted/NC_review/analysis_embedding_phi/GO_analysis/top10")

gene.files.full <- list.files(pattern=glob2rx("SKCM_*"),
						path="/fscratch/tinyi/Ted/NC_review/analysis_embedding_phi/GO_analysis/top10",
						full.name=T)

gene.set.names <- lapply(gene.files.full, function(gene.files.i) read.table(gene.files.i,header=T) [,"pathway"] )
names(gene.set.names) <- gsub(".txt","", gene.files)




selected.process <- lapply(1:length(gene.set.names),function(tum.idx){
	gene.set.tum.i <- gene.set.names[[tum.idx]]
	
	ret.list <- lapply(1:length(gene.set.tum.i),function(pathway.idx){
		gene.set.tum.i.pathway.j <- gene.set.tum.i[[pathway.idx]]
		if(grepl(glob2rx("*GO*"), gene.set.tum.i.pathway.j)) return(pathways.GO[[gene.set.tum.i.pathway.j]])
		if(grepl(glob2rx("*HALLMARK*"), gene.set.tum.i.pathway.j)) return(pathways.hallmark[[gene.set.tum.i.pathway.j]]) 
	})
	
	names(ret.list) <- gene.set.tum.i
	ret.list
}  )

names(selected.process) <- names(gene.set.names)


half.len <- length(selected.process)/2
selected.process.merged <- list()
for(i in 1: half.len) {
	selected.process.merged[[i]] <- c(selected.process[[i]],  selected.process[[i+ half.len]])
}

#names(selected.process.merged) <- paste(pdf.prefix,"_Tumor-",1:half.len,sep="")

 




process.plot  <- unlist(selected.process.merged,recursive=F)
process.plot <- process.plot[!is.na(names(process.plot))]

all.score <- gsva(phi.hat.tum, gs= process.plot, verbose=FALSE, parallel.sz=10)

order.vec <- -apply(all.score,1,which.max) + apply(all.score,1,max)

all.score.ordered <- all.score[order(order.vec,decreasing=TRUE),]
all.score.ordered <- unique(all.score.ordered)
pdf("skcm.subtypesNMF.genesets.pdf",pointsize=8,useDingbats=FALSE )
my_palette <- colorRampPalette(c("green", "black", "red"))(n = 254)


heatmap.2(all.score.ordered, 
								col=my_palette, 
								density.info="none", 
								trace="none",
								dendrogram='none',
								symm=F, 
								symkey=F, 
								scale="none" , 
								Rowv= NULL, 
								Colv= NULL,
								margins=c(10,10),
								cellnote= round(all.score.ordered,2),
     							notecex=0.5,
     							notecol="white",
     							cexRow=0.5)
     		
dev.off()


# # process.plot  <- c(merged.list, unlist(selected.process.merged,recursive=F))
# process.plot <- process.plot[!is.na(names(process.plot))]

# all.score <- gsva(phi.hat.tum, gs= process.plot, verbose=FALSE, parallel.sz=10)


# pdf("skcm.subtypesNMF.more.pdf",pointsize=8,useDingbats=FALSE )
# my_palette <- colorRampPalette(c("green", "black", "red"))(n = 254)

# heatmap.2(all.score, 
								# col=my_palette, 
								# density.info="none", 
								# trace="none",
								# dendrogram='none',
								# symm=F, 
								# symkey=F, 
								# scale="none" , 
								# Rowv= NULL, 
								# Colv= NULL,
								# margins=c(10,10),
								# cellnote= round(all.score,2),
     							# notecex=0.5,
     							# notecol="white",
     							# cexRow=0.3)
     		
# dev.off()




