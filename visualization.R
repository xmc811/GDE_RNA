
vsd_km <- function(vsd, k) {
    pca <- DESeq2::plotPCA(vsd, 
                           intgroup = colnames(vsd@colData)[1],
                           returnData = TRUE)
    set.seed(42)
    km_res <- kmeans(pca[,1:2], k)
    return(km_res)
}

plot_pca_vsd_km <- function(vsd, km_res, pal) {
    vsd@colData$Kmeans <- LETTERS[km_res$cluster]
    plot_pca_vsd(vsd, "Kmeans", pal)
}


deseq_table <- function(res, p_co, lfc_co) {
    
    res %<>%
        deseq_transform(p_co, lfc_co)
    
    res %<>%
        filter(significant != "Not Sig") %>%
        select(symbol:padj) %>%
        arrange(padj, abs(log2FoldChange)) %>%
        datatable() %>%
        formatRound(columns = c(2:5), digits = 3) %>%
        formatSignif(columns = c(6:7), digits = 3)
    
    return(res)
}
