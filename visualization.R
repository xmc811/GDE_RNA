
deseq_heatmap <- function(dds, var, palette, dir) {
    
    withProgress(message = "Plotting...", value = 0, {
        
        vsd <- vst(dds, blind = FALSE)
        
        incProgress(0.2)
        
        sampleDists <- dist(t(assay(vsd)))
        
        incProgress(0.1)
        
        sampleDistMatrix <- as.matrix(sampleDists)
        rownames(sampleDistMatrix) <- vsd[[var]]
        colnames(sampleDistMatrix) <- vsd[[var]]
        
        num <- num_colors(palette)
        
        colors <- brewer.pal(num,palette)
        
        if (dir) {
            colors <- rev(colors)
        }
        
        col_fun <- colorRamp2(seq(from = 0, 
                                  to = max(sampleDistMatrix), 
                                  length.out = num), colors)
        
        incProgress(0.4)
        
        Heatmap(sampleDistMatrix, 
                col = col_fun,
                rect_gp = gpar(col = "white", lwd = 2))
        
    })
}

deseq_pca <- function(dds, var, palette, dir) {
    
    vsd <- vst(dds, blind = FALSE)
    
    if (is.numeric(dds@colData[[var]])) {
        
        plotPCA(vst(dds, blind = FALSE), intgroup = var) +
            scale_color_distiller(palette = palette, direction = dir) +
            labs(color = var) +
            theme_bw() +
            theme(aspect.ratio = 1)
        
    } else {
        
        plotPCA(vst(dds, blind = FALSE), intgroup = var) +
            scale_color_brewer(palette = palette) +
            labs(color = var) +
            theme_bw() +
            theme(aspect.ratio = 1)
    }
}

deseq_ma <- function(res, p_co, lfc_co, lfc_plot_lim = 5) {
    
    res %<>%
        deseq_transform(p_co, lfc_co)
    
    res %>%
        mutate(shape = ifelse(log2FoldChange > lfc_plot_lim | log2FoldChange < -lfc_plot_lim, TRUE, FALSE),
               log2FoldChange = replace(log2FoldChange, log2FoldChange > lfc_plot_lim, lfc_plot_lim),
               log2FoldChange = replace(log2FoldChange, log2FoldChange < -lfc_plot_lim, -lfc_plot_lim)) %>%
        arrange(factor(significant, levels = c("Not Sig","Down","Up"))) %>%
        ggplot() +
        geom_point(aes(x = log10(baseMean), 
                       y = log2FoldChange, 
                       color = significant, 
                       shape = shape),
                   size = 2) +
        scale_color_manual(values = c("#1f78b4","#d9d9d9", "#e31a1c")) +
        scale_shape_manual(values = c(16, 17)) +
        theme_bw() +
        labs(y = expression(Log[2]~Fold~Change), x = expression(Log[10]~Mean~Normalized~Count)) +
        theme(legend.position = "none",
              plot.margin = unit(rep(1,4), "cm"),
              axis.text = element_text(size = 12),
              axis.title = element_text(size = 14)) +
        geom_hline(yintercept = 0, color = "#984ea3", size = 1.5, alpha = 0.5)
}



deseq_volcano <- function(res, p_co, lfc_co, 
                          p_plot_lim = 5, 
                          lfc_plot_lim = 5) {
    
    res %<>%
        deseq_transform(p_co, lfc_co)
    
    res %>%
        mutate(shape1 = ifelse(log2FoldChange > lfc_plot_lim | log2FoldChange < -lfc_plot_lim, TRUE, FALSE),
               log2FoldChange = replace(log2FoldChange, log2FoldChange > lfc_plot_lim, lfc_plot_lim),
               log2FoldChange = replace(log2FoldChange, log2FoldChange < -lfc_plot_lim, -lfc_plot_lim)) %>%
        mutate(shape2 = ifelse(-log10(padj) > p_plot_lim, TRUE, FALSE),
               padj = ifelse(-log10(padj) > p_plot_lim, 10^(-p_plot_lim), padj),
               shape = shape1 | shape2) %>%
        arrange(factor(significant, levels = c("Not Sig","Down","Up"))) %>%
        ggplot() +
        geom_point(aes(x = log2FoldChange,
                       y = -log10(padj),
                       color = significant,
                       shape = factor(shape)),
                   size = 2) +
        scale_color_manual(values = c("#1f78b4","#d9d9d9", "#e31a1c")) +
        scale_shape_manual(values = c(16, 17), guide = FALSE) +
        theme_bw() +
        labs(y = expression(-Log[10]~Adjusted~p-value), x = expression(Log[2]~Fold~Change)) +
        theme(aspect.ratio = 1,
              legend.position = "none",
              plot.margin = unit(rep(1,4), "cm"),
              axis.text = element_text(size = 12),
              axis.title = element_text(size = 14))
    
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
}


deseq_box <- function(dds, genes, var, palette) {
    
    df <- get_nm_count_dds(dds, genes, var)
    
    ggplot(df) +
        geom_boxplot(aes(x = !!sym(var),
                         y = log10(count),
                         fill = !!sym(var))) +
        geom_point(aes(x = !!sym(var),
                       y = log10(count))) +
        facet_wrap(~symbol) +
        scale_fill_brewer(palette = palette)
}


deseq_gsea <- function(res, pathways) {
    
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
    
}
