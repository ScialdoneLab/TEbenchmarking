library(splatter)

# old matrix
out_dir <- "splatter_old_mm"
num_genes <- 3173
num_cells <- 500

sim <- splatSimulate( 
	nGenes=num_genes, 
	batchCells=num_cells, 
	verbose = FALSE
)
dir.create(out_dir, showWarnings=FALSE)

selected_loci <- read.table(file="selected_loci_old.csv", header=F, row.names=NULL)$V1

write.table(colnames(sim), file= file.path(out_dir, "quants_mat_cols.txt"), quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(counts(sim), file= file.path(out_dir, "quants_mat.csv"), quote=FALSE, col.names=FALSE, row.names=FALSE, sep=",")
write.table(selected_loci, file= file.path(out_dir, "quants_mat_rows.txt"), quote=FALSE, col.names=FALSE, row.names=FALSE)


# young matrix

out_dir <- "splatter_young_mm"
num_genes <- 2100
num_cells <- 500

sim <- splatSimulate( 
	nGenes=num_genes, 
	batchCells=num_cells, 
	verbose = FALSE
)
dir.create(out_dir, showWarnings=FALSE)

selected_loci <- read.table(file="selected_loci_young.csv", header=F, row.names=NULL)$V1

write.table(colnames(sim), file= file.path(out_dir, "quants_mat_cols.txt"), quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(counts(sim), file= file.path(out_dir, "quants_mat.csv"), quote=FALSE, col.names=FALSE, row.names=FALSE, sep=",")
write.table(selected_loci, file= file.path(out_dir, "quants_mat_rows.txt"), quote=FALSE, col.names=FALSE, row.names=FALSE)


# matrix with genes

out_dir <- "splatter_all_wGenes_mm"
num_genes <- 7273
num_cells <- 500

sim <- splatSimulate( 
	nGenes=num_genes, 
	batchCells=num_cells, 
	verbose = FALSE
)
dir.create(out_dir, showWarnings=FALSE)

selected_loci <- c(read.table(file="selected_loci_old.csv", header=F, row.names=NULL)$V1,
					read.table(file="selected_loci_young.csv", header=F, row.names=NULL)$V1,
                    read.table("hvg_genes_ensembl_mm10.csv", header=F, row.names=NULL)$V1)

write.table(colnames(sim), file= file.path(out_dir, "quants_mat_cols.txt"), quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(counts(sim), file= file.path(out_dir, "quants_mat.csv"), quote=FALSE, col.names=FALSE, row.names=FALSE, sep=",")
write.table(selected_loci, file= file.path(out_dir, "quants_mat_rows.txt"), quote=FALSE, col.names=FALSE, row.names=FALSE)