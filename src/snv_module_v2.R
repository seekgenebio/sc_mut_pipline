library(Seurat)
library(tidyverse)
library(argparse)
library(patchwork)
library(Matrix)
library(future.apply)

args <- commandArgs(trailingOnly = FALSE)
scdir <- normalizePath(dirname(sub('--file=', '', args[grep('--file', args)])))
parser <- ArgumentParser()
parser$add_argument("--snv_cover_mat", required = TRUE,
                    help = "Required. SNV coverage mat path, directory or file(in txt format).")
parser$add_argument("--snv_mut_mat", required = TRUE,
                    help = "Required. SNV mut mat path, directory or file(in txt format).")
parser$add_argument("--rna_mat", required = TRUE,
                    help = "Required. SINGLE sample rds path. ")
parser$add_argument("--name", required = TRUE,
                    help = "Required. Sample name.")
parser$add_argument("--outdir", default = "./",
                    help = "The output directory [default %(default)s]")
parser$add_argument("--hotspot",
                    help = "Hotspot list file. One snv position per line, 
                    the snv position must be the rowname of snv mat. 
                    Like: EGFR:7-55191822:T>G. List file WITHOUT header.")
parser$add_argument("--cellanno", 
                    default = "seurat_clusters",
                    help = "If do not have meta file, 
                    please specify one column names in @meta.data to group cells.")
parser$add_argument("--meta",help = "Meta data file. CSV format, and with HEADER(Barcode,CellType).")
parser$add_argument("--query",help = "Please specify tumor celltype label. 
                    Multiple celltypes use comma to speparate, 
                    eg: tumor_celltype1,tumor_celltype2, default is find SNV among all celltypes, 
                    this will take a long time.")
parser$add_argument("--memory",help = "memory size, default is 2G per cpu.", default = 2000)


Args <- parser$parse_args(args = commandArgs(trailingOnly = TRUE))
MYCOLOR <- c(
          "#6394ce", "#2a4c87", "#eed500", "#ed5858",
          "#f6cbc2", "#f5a2a2", "#3ca676", "#6cc9d8",
          "#ef4db0", "#992269", "#bcb34a", "#74acf3",
          "#3e275b", "#fbec7e", "#ec4d3d", "#ee807e",
          "#f7bdb5", "#dbdde6", "#f591e1", "#51678c",
          "#2fbcd3", "#80cfc3", "#fbefd1", "#edb8b5",
          "#5678a8", "#2fb290", "#a6b5cd", "#90d1c1",
          "#a4e0ea", "#837fd3", "#5dce8b", "#c5cdd9",
          "#f9e2d6", "#c64ea4", "#b2dfd6", "#dbdfe7",
          "#dff2ec", "#cce8f3", "#e74d51", "#f7c9c4",
          "#f29c81", "#c9e6e0", "#c1c5de", "#750000"
          )

snv_cover_mat <- Args$snv_cover_mat
snv_mut_mat <- Args$snv_mut_mat
rna_mat <- Args$rna_mat
samplename <- Args$name
outdir <- Args$outdir
hotspot <- Args$hotspot
cellanno <- Args$cellanno
meta <- Args$meta
query <- Args$query
memory <- Args$memory
cpu <- 8

print(Args)

if(! dir.exists(outdir)){
    dir.create(outdir, recursive = T)
}
options(future.globals.maxSize = memory * 1024^2)

# Function Field
generate_mat <- function(obj, mat, assay_name){
    # Function：生成突变/Cover矩阵的函数
    ## obj: Seurat 对象
    ## mat: 矩阵文件路径
    ## assay_name: 生成的Assay名称
    readed_mat <- as(as.matrix(data.table::fread(mat),rownames = 1), 'dgCMatrix')
    
    ## 向一个Seurat Object中添加assay要求barcode一致，这里是将barcode差集默认值都填充为0
    diff.barcodes <- setdiff(colnames(obj), colnames(readed_mat))
    # obj 的barcode有可能被添加了后缀
    if(length(diff.barcodes) == length(colnames(obj))){
        
        suffix <- gsub('.*_', '', colnames(obj)) %>% unique()
        if(length(suffix) == 1){
            warning(paste0('add ', suffix, ' suffix to SNV matrix...'))
            colnames(readed_mat) <- paste0(colnames(readed_mat),"_", suffix)
            diff.barcodes <- setdiff(colnames(obj), colnames(readed_mat))
            if(length(diff.barcodes) == length(colnames(obj))){
                stop('The intersection of cell barcode of rna data and snv data is None...')
            }
        }else{
            stop('The intersection of cell barcode of rna data and snv data is None... 
            And we could not infer suffix')
        }
    }
    # 缺少的barcode用0填充信息
    if(length(diff.barcodes) > 0) {
        tmp_mat <- Matrix::Matrix(
            data = 0L, 
            nrow = nrow(readed_mat), 
            ncol = length(diff.barcodes), 
            sparse = TRUE
        )
        rownames(tmp_mat) <- rownames(readed_mat)
        colnames(tmp_mat) <- diff.barcodes
        readed_mat <- cbind(readed_mat, tmp_mat)
    }
    ## readed_mat中按照将obj中的barcode顺序排序
    readed_mat <- readed_mat[, colnames(obj)]
    obj[[assay_name]] <- CreateAssayObject(readed_mat)
    return(obj)

}

get_mut_number_and_cover_number <- function(obj, pos, assay_name, min_cells = 3, column_name = 'mut_cells',...){
    # Function：获取突变的细胞和Cover的细胞数
    ## obj: Seurat 对象
    ## pos: 突变位点，是个字符型向量
    ## assay_name: 提取突变信息的Assay名称
    ## min_cells: 要求突变至少在min_cells个细胞中发生突变，或者这个位点有被覆盖到
    ## coulmn_name: 生成的突变位点和对应发生突变/有覆盖到的细胞数的 data.frame
    
    keep_pos <- intersect(pos, rownames(GetAssayData(obj, assay = assay_name)))
    if(length(keep_pos) != length(pos) & length(keep_pos) > 0){
        warning('Some pos are not in assays...')
    }else if(length(keep_pos) == 0){
        stop('No pos in assay...')
    }
    pos_umi_cells <- list()
    for(i in levels(Idents(obj))){
        # 获取每个突变的发生突变的细胞数，即UMI不等于0
        get_umi_cells <- rowSums(
            GetAssayData(obj, assay = assay_name, slot = 'counts')[keep_pos, WhichCells(obj, idents = i)] != 0
        )
        # 要求每个突变至少在min_cells个细胞中中发生了突变，默认：1
        pos_umi_cells[[i]] <- get_umi_cells[get_umi_cells >= min_cells]
        
    }
    lapply(names(pos_umi_cells), function(x) {
        tmp <- as.data.frame(pos_umi_cells[[x]]) %>% rownames_to_column('pos')
        colnames(tmp) <- c('pos', column_name)
        tmp <- tmp %>% mutate(Cluster = x)
        tmp
    }) %>%
    do.call('rbind',.) -> df
    return(df)
}
snv_fisher_test <- function(df, ident){
    # Function: Fisher test
    target_df <- df %>% dplyr::filter(Cluster == ident) %>% dplyr::filter(mut_cells != 0)
    other_df <- df %>% dplyr::filter(Cluster != ident)
    results <- c()
    for(i in 1:nrow(target_df)){
        pos <- target_df[i, 'pos']
        target_mut_number <- target_df[i, 'mut_cells']
        target_ref_number <- target_df[i, 'ref_cells']
        tmp_other <- other_df[other_df$pos == pos,] %>% distinct() %>%
                            summarise(total_mut = sum(mut_cells), total_ref = sum(ref_cells))
        if(is.null(tmp_other)){
            other_mut_number <- 0
            other_ref_number <- 0
            
        }else{
            other_mut_number <- tmp_other %>% pull(total_mut)
            other_ref_number <- tmp_other %>% pull(total_ref)
        }
        tmp_matrix <- matrix(
            c(target_mut_number, target_ref_number ,other_mut_number, other_ref_number),
            ncol=2
        )
        fisher_out = fisher.test(tmp_matrix, alternative = 'greater')
        results = rbind(results, data.frame(
            SNV = pos,
            p_val=fisher_out$p.value,
            ident1_cover =  target_ref_number + target_mut_number,
            ident1_ref = target_ref_number,
            ident1_mut = target_mut_number,
            ident2_cover = other_ref_number + other_mut_number,
            ident2_ref = other_ref_number,
            ident2_mut = other_mut_number,
            cluster = ident))
    }
    results$p_val_adj <- p.adjust(results$p_val)
    return(results)
}
fill_NA <- function(combind_df){
    # Function: 填充full_join 产生的 NA值
    combind_df[is.na(combind_df$pos), 'Cluster'] <- gsub('.*___','',combind_df[is.na(combind_df$pos), 'pos_cluster'])
    combind_df[is.na(combind_df$pos), 'mut_cells'] <- 0
    combind_df[is.na(combind_df$pos), 'ref_cells'] <- combind_df[is.na(combind_df$pos), 'cover_cells'] - combind_df[is.na(combind_df$pos), 'mut_cells']
    combind_df[is.na(combind_df$pos), 'pos'] <- gsub('___.*','',combind_df[is.na(combind_df$pos), 'pos_cluster'])

    combind_df[is.na(combind_df$pos_cluster), 'cover_cells'] <- 0
    combind_df[is.na(combind_df$pos_cluster), 'ref_cells'] <- 0
    combind_df[is.na(combind_df$pos_cluster), 'pos_cluster'] <- paste0(
        combind_df[is.na(combind_df$pos_cluster), 'pos'],
        '___',
        combind_df[is.na(combind_df$pos_cluster), 'Cluster']
    )
    return(combind_df)
}

# 流程开始运行
setwd(outdir)
# 读取单样本的rds文件。已经做好细胞过滤，细胞分群的rds文件。不能是整合的rds！因为细胞barcode会改变！！
obj <- readRDS(rna_mat)
obj$Sample <- samplename

if(!is.null(meta)){
    meta <- read.table(meta, sep = ',', header = T, check.names = F)
    meta <- meta %>% column_to_rownames('Barcode')
    cellanno <- 'CellType'
    if(length(meta$CellType) < ncol(obj)){
        obj <- subset(obj, cells = rownames(meta))
    }
    obj <- AddMetaData(obj, meta)
}

# 读取 SNV coverage 矩阵 
if(file.exists(snv_cover_mat)){
    obj <- generate_mat(obj, snv_cover_mat, assay_name = 'SNV_all')
}else{
    stop('SNV cover file not exists!')
}
# 读取 SNV mut 矩阵
if(file.exists(snv_mut_mat)){
    obj <- generate_mat(obj, snv_mut_mat, assay_name = 'SNV')
}else{
    stop('SNV mut file not exists!')
}

Idents(obj) <- cellanno
obj$CellType <- Idents(obj)

if (length(unique(obj@meta.data[, cellanno])) > length(MYCOLOR)){
    print(paste0('cluster number is ', length(unique(obj@meta.data[, cellanno]))))
    MYCOLOR <- scales::hue_pal()(length(unique(obj@meta.data[, cellanno])))
}

if(! is.null(hotspot)){
    hotspot <- read.table(hotspot, header = F,sep='\t', check.names = F, stringsAsFactors = F)
    dir.create('hotspots')
    mut_df <- get_mut_number_and_cover_number(
        obj, 
        hotspot$V1, 
        assay_name = 'SNV', 
        min_cells = 1, 
        column_name = 'mut_cells'
    )
    cover_df <- get_mut_number_and_cover_number(
        obj, 
        hotspot$V1, 
        assay_name = 'SNV_all', 
        min_cells = 0, 
        column_name = 'cover_cells'
    )

    combind_df <- mut_df %>%
    mutate(pos_cluster = paste0(pos,'___',Cluster)) %>%
    full_join(
        cover_df %>%
        mutate(pos_cluster = paste0(pos,'___',Cluster)) %>% 
        select(-Cluster) %>% 
        select(-pos),
        by = 'pos_cluster'
    ) %>%
    mutate(ref_cells = cover_cells - mut_cells)
    combind_df <- fill_NA(combind_df)
    future::plan('multisession', workers = cpu)
    cluster_number <- unique(combind_df$Cluster)
    system.time({
        y <- future_lapply(cluster_number, FUN = function(x) snv_fisher_test(combind_df, x))
    })
    all_results <- do.call('rbind',y)
    
    dir.create(paste0(outdir,'/hotspots'), recursive = T)
    write_delim(all_results, paste0('hotspots/',samplename, '_hotspot_fisher_test_results.xls'), delim = '\t')
    top20_sig_results <- all_results %>% arrange(p_val_adj) %>% dplyr::filter(p_val_adj < 0.05)
    for(i in top20_sig_results$pos){
        obj$tmp_mut <- 'Others'
        p1 <- DimPlot(obj, 
                    group.by = cellanno, 
                    label = F, 
                    repel = T, 
                    label.box = T, 
                    label.color = "white", 
                    cols = MYCOLOR) 
        
        gene <- unlist(strsplit(i,':'))[[1]]
        DefaultAssay(obj) <- 'RNA'
        p2 <- FeaturePlot(obj,features = gene)
        DefaultAssay(obj) <- 'SNV_all'
        dat_tmp <- FetchData(obj, i)
        cover.cell <- rownames(dat_tmp[dat_tmp[, i] != 0,drop = F, ])
        DefaultAssay(obj) <- 'SNV'
        dat_tmp <- FetchData(obj,i)
        mut.cell <- rownames(dat_tmp[dat_tmp[, i] != 0, drop = F, ])
        cover.cell = setdiff(cover.cell, mut.cell)
        obj@meta.data[cover.cell,'tmp_mut'] <- 'Cover'
        obj@meta.data[mut.cell,'tmp_mut'] <- 'Mut'
        obj$tmp_mut <- factor(obj$tmp_mut, levels = c('Cover','Mut','Others'))
        p3 <- DimPlot(
            obj, 
            group.by = 'tmp_mut',
            pt.size = 0.5,order = c('Mut', 'Cover', 'Others'),
            cols = c('grey90', MYCOLOR[c(3, 1)])) + 
        ggtitle(i) + 
        theme(plot.title = element_text(hjust = .5)) 
        
        p4 <- DimPlot(obj, 
                    cells.highlight = list('Mut' = mut.cell), 
                    cols.highlight = MYCOLOR[1]) + 
            ggtitle(i) + 
            theme(plot.title = element_text(hjust = .5))
        p <- (p1 | p2) / (p3 | p4)
        ggsave(
            paste0(
                'hotspots/',
                samplename, 
                '_',
                gsub(':|>', "_", i),
                '_featureplot.png'), 
            p,
            dpi = 300, 
            width = 12,
            height = 10)
    }
}
mutsites <- rownames(GetAssayData(obj, assay = 'SNV'))
print(paste0('Total SNV pos: ',length(mutsites), ' ...'))
mut_df <- get_mut_number_and_cover_number(
    obj, 
    mutsites, 
    assay_name = 'SNV', 
    min_cells = 1, 
    column_name = 'mut_cells'
)
cover_df <- get_mut_number_and_cover_number(
    obj, 
    mutsites, 
    assay_name = 'SNV_all', 
    min_cells = 0, 
    column_name = 'cover_cells'
)

combind_df <- mut_df %>%
mutate(pos_cluster = paste0(pos,'___',Cluster)) %>%
left_join(
    cover_df %>%
    mutate(pos_cluster = paste0(pos,'___',Cluster)) %>% 
    select(-Cluster) %>% 
    select(-pos),
    by = 'pos_cluster'
) %>%
mutate(ref_cells = cover_cells - mut_cells)
combind_df <- fill_NA(combind_df)

future::plan('multisession', workers = 8)
cluster_number <- unique(combind_df$Cluster)
print('fisher test')
head(combind_df)
system.time({
    y <- future_lapply(cluster_number, FUN = function(x) snv_fisher_test(combind_df, x))
})
all_results <- do.call('rbind', y)

dir.create(paste0(outdir,'/SNV_diff'), recursive = T)
all_results <- all_results %>% arrange(p_val_adj, p_val) 
write_delim(all_results, paste0('SNV_diff/',samplename, '_snv_markers.xls'), delim = '\t')
# 按照p_val_adj 排序，但是只要fisher 检验的p_val < 0.05即可，为了放宽阈值，能展示更多的突变
top10_sig_results <- all_results %>% 
group_by(cluster) %>%
dplyr::filter(p_val < 0.05) %>%
slice_min(n = 10, order_by = p_val_adj) 
sub_markers_list <- split(top10_sig_results, top10_sig_results$cluster)

for(ct in names(sub_markers_list)){
    ct <- gsub("[/\\s]","_",perl = T,ct)
    print(ct)
    dir.create(paste0('SNV_diff/',ct), recursive = T)
    # 根据p_val_adj筛选前10，可能会筛选出来多于10个突变，在这里再筛选一遍。 
    sub_markers_list[[ct]] <- head(sub_markers_list[[ct]], n = 10)
    for(i in sub_markers_list[[ct]]$SNV){
        obj$tmp_mut <- 'Others'
        DefaultAssay(obj) <- 'SNV_all'
        cells_all <- FetchData(obj, c(i))
        colnames(cells_all) <- 'count'
        cover.cell <- rownames(cells_all %>% dplyr::filter(count != 0 ))
        gene <- gsub(":.*","",i)
        
        DefaultAssay(obj) <- 'SNV'
        dat_tmp <- FetchData(obj,i)
        mut.cell <- rownames(dat_tmp[dat_tmp[,i] != 0, drop = F,])
        cover.cell <- setdiff(cover.cell,mut.cell)

        obj@meta.data[cover.cell,'tmp_mut'] <- 'Cover'
        obj@meta.data[mut.cell,'tmp_mut'] <- 'Mut'
        obj$tmp_mut <- factor(obj$tmp_mut ,levels = c('Cover','Mut','Others'))

	if(! gene %in% rownames(obj@assays$RNA@counts)){
	    print(paste0("skip plot ", gene, " expression.."))
        
        p1 <- DimPlot(
            obj,
            group.by = 'tmp_mut',
            order = c('Mut', 'Cover', 'Others'),
            cols = c('grey90', MYCOLOR[c(3,1)]),
            raster = F, 
            pt.size = .1
            ) + 
        ggtitle(i) + 
        theme(plot.title = element_text(hjust = .5))
        p3 <- DimPlot(obj, 
                    cells.highlight = list('Mut' = mut.cell), 
                    cols.highlight = MYCOLOR[1],
                    raster = F, 
                    pt.size = .1
            ) + 
            ggtitle(i) + 
            theme(plot.title = element_text(hjust = .5))

        p4 <-  DimPlot(
            obj,
            group.by = 'CellType',
            label = F,
            repel = T,
            label.box = T,
            label.color = "white",
            cols = MYCOLOR,
            raster = F, 
            pt.size = .1
            ) 
        p <- (p1 | plot_spacer()) / (p3 | p4)   

	}else{
        p1 <- DimPlot(
            obj, 
            group.by = 'tmp_mut',
            order = c('Mut', 'Cover', 'Others'),
            cols = c('grey90', MYCOLOR[c(3, 1)]),
            raster = F, 
            pt.size = .1
            ) + 
        ggtitle(i) + 
        theme(plot.title = element_text(hjust = .5))
        DefaultAssay(obj) <- 'RNA'
	    p2 <-  FeaturePlot(obj,features = gene, order = T)
        p3 <- DimPlot(
            obj, 
            cells.highlight = list('Mut' = mut.cell), 
            cols.highlight = MYCOLOR[1],
            raster = F, 
            pt.size = .1
            ) + 
            ggtitle(i) + 
            theme(plot.title = element_text(hjust = .5))

        p4 <- DimPlot(
            obj,
            group.by = 'CellType',
            label = F,
            repel = T,
            label.box = T,
            label.color = "white",
            cols = MYCOLOR,
            raster = F, 
            pt.size = .1) 

        p <- (p1 | p2) / (p3 | p4)   
	}
         
        ggsave(
            paste0('SNV_diff/',ct, '/',gsub(':|>',"_",i),".png"), 
            p, 
            dpi = 300, 
            width = 12, 
            height = 10
        )
    }
}


saveRDS(obj, paste0(samplename,"_add_snv.rds"))






