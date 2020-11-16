

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

deseq_cluster <- function(dds, genes, palette, dir) {
    
    withProgress(message = "Plotting...", value = 0, {
    
    mtx <- get_mtx_dds(dds, genes)
    
    mtx %<>%
        mtx_rescale()
    
    num <- num_colors(palette)
    colors <- brewer.pal(num, palette)
    
    if (dir) {
        colors <- rev(colors)
    }
    
    col_fun <- colorRamp2(seq(from = -1, 
                              to = 1, 
                              length.out = num), colors)
    
    Heatmap(mtx, 
            col = col_fun,
            rect_gp = gpar(col = "white", lwd = 2))
    })
}


deseq_box <- function(dds, genes, var, palette) {
    
    withProgress(message = "Plotting...", value = 0, {
    
    df <- get_nm_count_dds(dds, genes, var)
    
    ggplot(df) +
        geom_boxplot(aes(x = !!sym(var),
                         y = log10(count),
                         fill = !!sym(var))) +
        geom_point(aes(x = !!sym(var),
                       y = log10(count))) +
        facet_wrap(~symbol) +
        scale_fill_brewer(palette = palette)
    })
}


deseq_gsea <- function(res, pathways) {
    
    withProgress(message = "Plotting...", value = 0, {
    
    res <- deseq_to_stat(res)
    
    gsea_res <- fgsea(pathways = pathways, 
                      stats = res, 
                      nperm = 10000)
    
    gsea_res %>%
        mutate(pathway = str_remove(string = pathway, pattern = "HALLMARK_")) %>%
        mutate(color = -log10(padj) * ifelse(padj <= 1, 1, 0) * ifelse(NES > 0, 1, -1)) %>%
        ggplot() +
        geom_bar(aes(x = reorder(pathway, NES), y = NES, fill = color), stat = "identity") +
        scale_fill_gradient2(high = "#d7301f", 
                             mid = "#f0f0f0",
                             low = "#0570b0",
                             midpoint = 0) +
        coord_flip() +
        labs(x = "Pathway",
             fill = "-log10 adjusted p-value") +
        theme_bw()
    })
}
