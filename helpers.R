
df_to_signature <- function(df) {
    
    colnames(df) <- c("pathway", "gene")
    
    pathways <- unique(df[[1]])
    
    sig <- list()
    
    for (i in seq_along(1:length(pathways))) {
        
        sig[[i]] <- df %>%
            filter(pathway == pathways[i]) %>%
            pull(gene)
    }
    
    names(sig) <- pathways
    
    return(sig)
}


deseq_format_test <- function(rds) {
    
    if ((class(rds[[1]]) == "DESeqDataSet") & (class(rds[[2]]) == "DESeqResults")) {
        
        return(TRUE)
    }
}

num_colors <- function(pal) {
    
    a <- rownames(brewer.pal.info) == pal
    return(brewer.pal.info$maxcolors[a])
}


deseq_transform <- function(res, p_co, lfc_co) {
    
    res %<>%
        as.data.frame() %>%
        rownames_to_column(var = "symbol") %>%
        as_tibble() %>%
        filter(!is.na(padj)) %>%
        mutate(significant = ifelse(padj <= p_co & log2FoldChange >= lfc_co, 
                                    "Up",
                                    ifelse(padj <= p_co & log2FoldChange <= -lfc_co, 
                                           "Down", 
                                           "Not Sig")))
    return(res)
}

get_rna_genes <- function(res, topn = 50) {
    
    genes <- res %>%
        as.data.frame() %>%
        rownames_to_column(var = "symbol") %>%
        as_tibble() %>%
        arrange(padj) %>%
        head(topn) %>%
        pull(symbol)
    
    return(genes)
}

parse_rna_genes <- function(gene_list) {
    
    gene_list <- gsub("[[:space:]]", "", gene_list)
    
    genes <- str_split(gene_list, "(,|;)")[[1]]
    
    return(genes[genes != ""])
}

get_mtx_dds <- function(dds, genes, raw = F) {
    
    if (raw) {
        vsd <- counts(dds)
    } else {
        vsd <- vst(dds, blind = FALSE)
        vsd <- as.matrix(vsd@assays@data[[1]])
    }
    
    mtx <- vsd[genes,,drop = FALSE]
    
    return(mtx)
}

get_nm_count_dds <- function(dds, genes, var) {
    
    df_list <- list()
    
    genes <- genes[genes %in% rownames(dds)]
    
    for (i in seq_along(1:length(genes))) {
        
        d <- plotCounts(dds, 
                        gene = genes[i], 
                        intgroup = var, 
                        returnData = TRUE)
        
        d %<>%
            rownames_to_column() %>%
            as_tibble() %>%
            mutate(symbol = genes[i])
        
        df_list[[i]] <- d
    }
    
    df <- do.call("rbind", df_list)
    
    return(df)
}


mtx_rescale <- function(mtx) {
    
    mtx2 <- mtx
    
    for (i in seq_along(1:nrow(mtx))) {
        mtx2[i,] <- (mtx[i,] - min(mtx[i,]))/(max(mtx[i,]) - min(mtx[i,])) * 2 - 1
    }
    return(mtx2)
}

mtx_logtrans <- function(mtx, b) {
    
    mtx2 <- mtx
    
    for (i in 1:nrow(mtx)) {
        mtx2[i,] <- logb(mtx[i,], b)
    }
    
    return(mtx2)
}

deseq_to_stat <- function(res) {
    
    stat <- res$stat
    names(stat) <- rownames(res)
    return(stat)
}
