library(argparse)
library(Seurat)
library(Matrix)
library(stringr)

args <- commandArgs(trailingOnly = FALSE)
scdir <- normalizePath(dirname(sub('--file=', '', args[grep('--file', args)])))

parser <- ArgumentParser()
parser$add_argument("--refmatrix", required = TRUE,
                    help = "ref_matrix.mtx")
parser$add_argument("--altmatrix", required = TRUE,
                    help = "alt_matrix.mtx")
parser$add_argument("--annofile", required = TRUE,
                    help = "final annotation file ;eg: final.vcf.txt")
parser$add_argument("--barcode", required = TRUE,
                    help = "barcodes file ")
parser$add_argument("--variants", required = TRUE,
                    help = "variants.txt")
parser$add_argument("--outdir", required = TRUE,
                    help = "oupput path")
Args <- parser$parse_args(args = commandArgs(trailingOnly = TRUE))

inmatrix <- readMM(Args$altmatrix)
inmatrix <- as.data.frame(as.matrix(inmatrix))
barcodes <- read.table(Args$barcode, header = F)
gene_file <- read.table(Args$annofile,header=T, sep="\t")
head(gene_file)
gene <- gene_file[,"SYMBOL"]
chr <- gene_file[,"chr"]
pos <- gene_file[,"pos"]
ref <- gene_file[,"ref"]
alt <- gene_file[,"Allele"]
genes <- paste0(gene, ":", chr, "-", pos, ":", ref, ">" ,alt)
genes<- data.frame(genes)
#vartions <- read.table("/PROJ2/FLOAT/zhaomeng/proj_zhm/YS2307047_ruijin_fusion/vatrix/pip_shell/demo_output/matrix//variants.txt", header = F)
colnames(inmatrix) <- barcodes$V1
row.names(inmatrix) <- genes$gene

ref_matrix <- readMM(Args$refmatrix)
refmatrix <- as.data.frame(as.matrix(ref_matrix))
colnames(refmatrix) <- barcodes$V1
row.names(refmatrix) <- genes$gene
coverage_matrix <- data.frame(Map(function(x, y) x + y, refmatrix, inmatrix))
row.names(coverage_matrix) <- genes$gene
inmatrix <- inmatrix[rowSums(inmatrix) > 2 ,]
final_alt_matrix <- inmatrix[apply(inmatrix > 0, 1, sum) > 2, ]
coverage_matrix <- coverage_matrix[rownames(final_alt_matrix), ]
length(rownames(coverage_matrix))
write.table(file=paste0(Args$outdir,"/final.alt.matrix"), data.frame("pos"= rownames(final_alt_matrix), final_alt_matrix),sep="\t",quote=F,col.names=T,row.names=F)
write.table(file=paste0(Args$outdir,"/final.all.matrix"), data.frame("pos"= rownames(coverage_matrix), coverage_matrix),sep="\t",quote=F,col.names=T,row.names=F)

