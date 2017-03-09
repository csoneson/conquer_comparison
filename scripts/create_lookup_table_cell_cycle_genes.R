## Create a lookup table for cell cycle genes
basedir <- "/home/Shared/data/seq/conquer/comparison"
cell_cycle_genes <- as.character(unlist(read.delim(paste0(basedir, "/data/mouse_cell_cycle_genes.txt"), header = FALSE)))

library(Biostrings)
cdna <- readDNAStringSet("/home/Shared/data/seq/conquer/database/reference-files/mus-musculus/Mus_musculus.GRCm38.cdna.all.fa")
id_to_symbol <- data.frame(t(sapply(names(cdna), function(nm) {
  s <- strsplit(nm, " ")[[1]]
  geneid <- gsub("^gene:", "", grep("^gene:", s, value = TRUE))
  symbol <- gsub("^gene_symbol:", "", grep("^gene_symbol:", s, value = TRUE))
  c(geneid = geneid, symbol = symbol)
})), stringsAsFactors = FALSE)
rownames(id_to_symbol) <- NULL

cell_cycle_ids <- setdiff(unique(id_to_symbol$geneid[match(cell_cycle_genes, id_to_symbol$symbol)]), NA)
saveRDS(list(Mus_musculus.GRCm38.84 = cell_cycle_ids,
             Mus_musculus = unique(cell_cycle_genes)), 
        file = paste0(basedir, "/data/cell_cycle_geneids.rds"))
