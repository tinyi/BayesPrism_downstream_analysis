

#prepare marker genes
#table S1 from https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0056823

library(readxl)

read_excel_allsheets <- function(filename, tibble = FALSE) {
    # I prefer straight data.frames
    # but if you like tidyverse tibbles (the default with read_excel)
    # then just pass tibble = TRUE
    sheets <- readxl::excel_sheets(filename)
    x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
    if(!tibble) x <- lapply(x, as.data.frame)
    names(x) <- sheets
    x
}

xls.list <- read_excel_allsheets("pone.0056823.s008.xlsx")

names(xls.list)

selelcted.marker.sets <- c("Basal High vs. Others","Mesenchymal High vs. Others","Atypical High vs. Others","Classical High vs. Others")

tcga.list <- xls.list[selelcted.marker.sets]
names(tcga.list) <- c("Basal","Mesenchymal","Atypical","Classical")

tcga.list.filtered  <- lapply(tcga.list, function(tcga.df.i) tcga.df.i[tcga.df.i$"q-Value"<0.01,"Gene"] )
ct.tab <- table(unlist(tcga.list.filtered))
uniq.genes <- names(ct.tab[ct.tab==1])
tcga.list <- lapply(tcga.list.filtered,function(tcga.list.filtered.i) tcga.list.filtered.i[tcga.list.filtered.i %in% uniq.genes] )


sc.marker <- t(read.table("sc.subtype.txt",row.names=1,header=F,sep=" ",fill = TRUE))
sc.marker <- data.frame(sc.marker)

merged.list <- c(tcga.list, as.list(sc.marker))

save(merged.list , file="marker.gene.rdata")





