
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


parse_rna_genes <- function(gene_list) {
    gene_list <- gsub("[[:space:]]", "", gene_list)
    genes <- str_split(gene_list, "(,|;)")[[1]]
    return(genes[genes != ""])
}


get_count_message <- function(mtx) {
    
    msg <- paste0(ncol(mtx), 
                  " HTSeq count files uploaded\n",
                  "Number of genes: ", nrow(mtx))
    return(msg)
}


filter_mt <- function(mtx, metadata) {
    metadata[match(colnames(mtx), metadata[["Sample"]]),]
}


cts_to_dds <- function(mtx, metadata, var = 1) {
    
    withProgress(message = "Loading Data..", value = 0.3, {
    
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = mtx,
                                  colData = metadata,
                                  design= as.formula(paste0("~", var)))
    
    incProgress(0.5, message = "Building DESeq2 Data Objects..")

    dds <- DESeq2::DESeq(dds)
    return(dds)
    })
}


assign_km_clu <- function(vsd, km_res) {
    vsd@colData$Kmeans <- LETTERS[km_res$cluster]
    return(vsd)
}


assign_km_clu_col <- function(coldata, km_res) {
    coldata$Kmeans <- LETTERS[km_res$cluster]
    return(coldata)
}


trubble <- function(cts) {
    tmp <- as.data.frame(
        rbind(
            cts[1:5, ],
            ... = rep("...", length(cts[1, ])),
            cts[(nrow(cts) - 4):(nrow(cts)), ]
        )
    )
    
    if (ncol(tmp) > 10) {
        tmp2 <- tmp[, 1:10]
    } else {
        tmp2 <- tmp
    }
    
    nr <- nrow(cts)
    nc <- ncol(cts)
    
    if (ncol(tmp) > 10) {
        output <- paste(
            "Your pre-processed data contains", nr, "genes and", nc, 
            "samples. Showing the first 10 samples:\n"
        )
    } else {
        output <- paste(
            "Your pre-processed data contains", nr, "genes and", 
            nc, "samples.\n"
        )
    }
    
    test <- as.matrix(tmp2)
    test <- rbind(colnames(tmp2), test)
    y <- sprintf(paste0("%",max(nchar(test)),"s"), test)
    y <- matrix(y, nrow = 12)
    
    gen <- c("", rownames(tmp2))
    gen <- gsub("\\s", " ", format(gen, width = max(nchar(gen))))
    
    if (ncol(tmp) > 10) {
        output2 <- paste("\n", ncol(tmp) - 10, "Samples not shown\n")
    } else {
        output2 <- NULL
    }
    
    cat(output, "\n")
    for(i in 1:nrow(y)) {
        cat(gen[i], y[i, ], "\n")
    }
    
    cat(output2)
}

