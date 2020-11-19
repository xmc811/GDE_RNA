
source("loading.R")
source("helpers.R")
source("ui_design.R")
source("tab_design.R")
source("visualization.R")

options(shiny.maxRequestSize = 500 * 1024 ^ 2)
options(spinner.color = "#3182bd")
options(spinner.type = 8)

# user interface
ui <- navbarPage(
    
    theme = shinytheme("yeti"),
    title = "Genomic Data Explorer",
    id = "tabs",
    
    # tab_design.R
    tab_file,
    tab_rna,
    tab_help,
    tab_about,
    
    tags$head(tags$link(rel = "stylesheet", 
                        type = "text/css", 
                        href = "style.css"))
)

# server function
server <- function(input, output, session) {
    
    # ----------
    # Data Input Tab
    # ----------
    
    # Raw Data Input - UI
    output$upload_panel <- renderUI({
        list(
            h3("RNA-Seq Data Input"),
            br(),
            h4("1. File Upload"),
            upload_widgets
        )
    })
    
    # Raw Data Input - Generation of Raw Counts and Metadata 
    mt_raw <- eventReactive(input$cts_upload_click, {
        if (input$cts_source == "Example") {
            read_csv("./data/metadata.csv")
        } else if (input$cts_source == "Upload") {
            read_csv(input$meta_file$datapath)
        } else {}
    })
    
    cts_raw <- eventReactive(input$cts_upload_click, {
        if (input$cts_source == "Example") {
            readRDS("./data/example_mtx.rds")
        } else if (input$cts_source == "Upload"){
            withProgress(message = "Loading Data..", value = 0.5, {
            htseq_to_mtx(input$cts_files)
            })
        } else {}
    })
    
    # Raw Data Input - Messages
    output$mt_message <- renderText({
        validate(need(input$cts_source, ""))
        if (input$cts_source == "Example") {} 
        else if (input$cts_source == "Upload") {
            validate(need(input$meta_file, "Please Upload Metadata"))
        } else {}
        paste0("Metadata uploaded: ", nrow(mt_raw()), " rows")
    })
    
    output$count_message <- renderText({
        validate(need(input$cts_source,""))
        if (input$cts_source == "Example") {} 
        else if (input$cts_source == "Upload") {
            validate(need(input$cts_files, "Please Upload Count Data"))
        } else {}
        get_count_message(cts_raw())
    })

    
    # Pre Processing of Raw Data - UI
    select_meta <- function(id, text, choices = colnames(mt_raw())) {
        selectInput(inputId = id, 
                    label = text, 
                    choices = choices, 
                    selected = choices[1])
    }
    
    
    output$cts_proc <- renderUI({
        validate(need(input$cts_source, ""))
        if (input$cts_source == "Example") {
        } else if (input$cts_source == "Upload") {
            validate(need(input$meta_file, "Please Upload Metadata"),
                     need(input$cts_files, "Please Upload Count Data"))
        } else {}
        list(
            br(),
            h4("2. Pre-Process Count and Metadata"),
            select_meta(id = "meta_sample_col", 
                        text = "Column with Sample Names"),
            select_meta(id = "meta_file_col", 
                        text = "Column with File Names"),
            number_cts_cutoff(),
            button_cts_process()
        )
    })
    
    
    # Pre Processing of Raw Data - Generation of Counts and Metadata
    cts <- eventReactive(input$cts_process_click, {
        mtx_name_match(cts_raw(), 
                       mt_raw(), 
                       input$meta_sample_col,
                       input$meta_file_col,
                       input$cts_cutoff)
    })

    mt <- reactiveVal()
    
    observeEvent(input$cts_process_click,{
        mt(filter_mt(cts(), mt_raw()))
    })
    
    # Pre Processing of Raw Data - Messages
    output$cts_summary <- renderPrint({
        validate(need(try(cts()), ""),
                 need(ncol(cts()) >= 1, "Name matching returns count matrix with 0 samples.\nPlease make sure the columns of sample names and files names are chosen correctly."))
        trubble(cts())
    })
    
    
    # Computation - UI
    output$compute_button <- renderUI({
        validate(need(try(cts()),""),
                 need(try(mt()),""))
        list(h4("3. Generate DESeq2 Data"),
             button_cts_compute())
    })
    
    
    # Computation - Generation of DDS and VSD
    dds <- reactiveVal()
    vsd <- reactiveVal()
    
    observeEvent(input$cts_compute_click, {
        dds(cts_to_dds(cts(), mt()))
        vsd(DESeq2::vst(dds(), blind = FALSE))
    })

    # Computation - Messages
    output$compute_message <- renderText({
        validate(need(!is.null(dds()), ""),
                 need(!is.null(vsd()), ""))
        paste0("DESeq2 data object generated.\nYou can explore your data in other panels.")
    })
    
    
    # PCA Sample Distance - UI
    output$pca_var_ui <- renderUI({
        validate(need(!is.null(vsd()), ""))
        selectInput(inputId = "pca_var", 
                    label = "Variable for PCA Plot and Sample Distance Heatmap",
                    choices = colnames(vsd()@colData),
                    selected = colnames(vsd()@colData)[1])
    })
    
    output$color_ui <- renderUI({
        list(palette_widgets,
             switch_palette_dir(""))
    })
    
    output$cluster_switch <- renderUI({
        switch_cluster_mode()
    })
    
    output$cluster_ui <- renderUI({
        validate(need(input$cluster_sw == TRUE, ""))
        cluster_widgets
    })
    
    # Plot Size
    output$plot_size <- renderUI({
        plot_size_widgets
    })
    
    plot_height <- reactive({
        validate(need(input$plot_height < 4000, message_plot_size("height")))
        return(input$plot_height)
    })
    
    plot_width <- reactive({
        validate(need(input$plot_width < 4000, message_plot_size("width")))
        return(input$plot_width)
    })
    
    # Plot Download Component
    button_download <- function(id) {
        renderUI({
            validate(need(try(vsd()), ""))
            downloadButton(id, label = "Download PDF")
        })
    }
    
    download_pdf <- function(file, type = NULL) {
        dpi <- 72
        if (is.null(type)) {
            pdf(file, height = viz_plot_height()/dpi, viz_plot_width()/dpi)
        } else {
            pdf(file, height = plot_height()/dpi, plot_width()/dpi)
        }
    }
    
    # PCA - Clustering
    km_res <- reactive({
        vsd_km(vsd(), input$n_cluster)
    })
    
    observeEvent(input$assign_cluster_click, {
        dds(assign_km_clu(dds(), km_res()))
        vsd(assign_km_clu(vsd(), km_res()))
        mt(assign_km_clu_col(mt(), km_res()))
    })
    
    # PCA - Plotting
    
    plots <- reactiveValues()
    
    validate_vsd <- function() {
        validate(need(try(vsd()), "VSD object not computed. Visualization not available."))
    }
    
    output$pca <- renderPlot({
        withProgress(message = "Plotting...", value = 0.5, {
        validate_vsd()
        plots$pca <- if (input$cluster_sw == TRUE) {
            plot_pca_vsd_km(vsd = vsd(), 
                            km_res = km_res(), 
                            pal = input$pal_categorical)
        } else if (is.numeric(dds()@colData[[input$pca_var]])) {
            plot_pca_vsd(vsd = vsd(), 
                         var = input$pca_var, 
                         pal = input$pal_continuous,
                         dir = ifelse(input$pal_dir, 1, -1))
        } else {
            plot_pca_vsd(vsd = vsd(), 
                         var = input$pca_var, 
                         pal = input$pal_categorical)
        }
        plots$pca
        })
    }, height = plot_height, width = plot_width)
    
    output$dl_pca_button <- button_download("dl_pca")
    
    output$dl_pca <- downloadHandler(filename = function() {"test.pdf"},
                                     content = function(file) {
                                         download_pdf(file, 1); print(plots$pca); dev.off()
                                         })
    
    # Sample Distance - Plotting
    output$hm <- renderPlot({
        withProgress(message = "Plotting...", value = 0.5, {
        validate_vsd()
        plots$hm <- plot_heatmap_vsd(vsd = vsd(), 
                         var = input$pca_var, 
                         pal = input$pal_continuous, 
                         dir = input$pal_dir)
        plots$hm
        })
    }, height = plot_height, width = plot_width)
    
    output$dl_hm_button <- button_download("dl_hm")
    
    output$dl_hm <- downloadHandler(filename = function() {"test.pdf"},
                                     content = function(file) {
                                         download_pdf(file, 1); print(plots$hm); dev.off()
                                     })
    
    
    # DGE
    output$dge_var_ui <- renderUI({
        validate(need(try(mt()),""))
            selectInput(inputId = "dge_var", 
                        label = "Variable for DGE",
                        choices = colnames(mt()))
    })
    
    output$dge_group1 <- renderUI({
        validate(need(try(mt()),""))
        selectInput(inputId = "dge_g1", 
                    label = "Group 1",
                    choices = unique(mt()[[input$dge_var]]))
    })
    
    output$dge_group2 <- renderUI({
        validate(need(try(mt()),""))
        selectInput(inputId = "dge_g2", 
                    label = "Group 2",
                    choices = setdiff(unique(mt()[[input$dge_var]]),
                                      input$dge_g1))
    })
    
    output$dge_button <- renderUI({
        validate(need(try(cts()),""),
                 need(try(mt()),""))
        button_dge()
    })
    
    dds_run <- reactiveVal()
    res <- reactiveVal()
    
    observeEvent(input$dge_click,{
        dds_run(cts_to_dds(cts(), mt(), input$dge_var))
        res(DESeq2::results(dds_run(), 
                            contrast = c(input$dge_var, 
                                         input$dge_g1, 
                                         input$dge_g2)))
    })
    
    output$dge_message <- renderText({
        validate(need(try(res()), ""))
        paste0("Differential gene expression (DGE) analysis done\n\n",
               res()@elementMetadata[2,2] %>%
                   str_remove("^.*:") %>%
                   str_trim(), "\n\n",
               "You can visualize your data in the 'DGE Visualization' tab.")
    })
    
    
    # DGE Visualization UI
    output$viz_plot_size <- renderUI({
        validate(need(try(res()), ""))
        viz_plot_size_widgets
    })
    
    output$dge_params <- renderUI({
        validate(need(try(res()), ""))
        list(
            h4("Differential Gene Expression Parameters"),
            dge_cutoff_widgets
            )
    })
    
    output$squash_params <- renderUI({
        validate(need(try(res()), ""))
        list(
            h4("Plotting Parameters"),
            dge_plot_limit_widgets
            )
    })
    
    output$rna_var <- renderUI({
        validate(need(try(dds_run()), ""))
        list(h4("Variable"),
             selectInput(inputId = "box_var", 
                         label = "Categorical Variable for Gene Boxplot",
                         choices = colnames(dds_run()@colData))
            )
    })
    
    output$viz_color_ui <- renderUI({
        validate(need(try(res()), ""))
        list(viz_palette_widgets,
             switch_palette_dir("viz_")
        )
    })
    
    output$viz_gene_ui <- renderUI({
        validate(need(try(res()), ""))
        gene_selection_widgets
        })
    
    output$pathway_ui <- renderUI({
        validate(need(try(res()), ""))
        pathway_selection_widgets
    })
    
    
    # DGE Visualization Parameters
    viz_plot_height <- reactive({
        validate(need(input$viz_plot_height < 4000, message_plot_size("height")))
        return(input$viz_plot_height)
    })
    
    viz_plot_width <- reactive({
        validate(need(input$viz_plot_width < 4000, message_plot_size("width")))
        return(input$viz_plot_width)
    })
    
    
    rna_genes <- eventReactive(
        c(input$rna_gene_read1,
          input$rna_gene_read2,
          input$rna_genes_file),
        {
            if(input$gene_source == 'Use Top Genes') {
                get_top_genes(res())[1:input$rna_gene_num]
            } else if (input$gene_source == 'Manual Input') {
                validate(need(input$rna_genes_man, "Please Input Gene List"))
                parse_rna_genes(input$rna_genes_man)
            } else {
                validate(need(input$rna_genes_file, "Please Upload Gene List"))
                readLines(input$rna_genes_file$datapath)
            }}
        )
    
    rna_pathways <- reactive({
        if(input$pathway_source == 'Use Hallmarks') {
            hmks_hs
        } else {
            validate(need(input$rna_pathway_file, "Please Upload Pathway File"))
            read_csv(input$rna_pathway_file$datapath, 
                     col_names = FALSE) %>% df_to_signature()
        }})
    
    
    # DGE Visualization Results
    validate_res <- function() {
        validate(need(try(res()), "No DGE results. Visualization not available."))
    }
    
    output$res_table <- DT::renderDataTable({
        validate_res()
        deseq_table(res(), 
                    input$p_co, 
                    input$lfc_co)
    })
    
    output$res_ma <- renderPlot({
        withProgress(message = "Plotting...", value = 0.5, {
        validate_res()
        plots$ma <- plot_deseq_ma(res(),
                                  input$p_co, 
                                  input$lfc_co,
                                  input$lfc_plot_lim)
        plots$ma
        })
    }, height = viz_plot_height, width = viz_plot_width)
    
    output$dl_ma_button <- button_download("dl_ma")
    
    output$dl_ma <- downloadHandler(filename = function() {"test.pdf"},
                                    content = function(file) {
                                        download_pdf(file); print(plots$ma); dev.off()
                                    })
    
    output$res_volcano <- renderPlot({
        withProgress(message = "Plotting...", value = 0.5, {
        validate_res()
        plots$volcano <- plot_deseq_volcano(res(), 
                                            input$p_co, 
                                            input$lfc_co,
                                            input$p_plot_lim,
                                            input$lfc_plot_lim)
        plots$volcano
        })
    }, height = viz_plot_height, width = viz_plot_width)
    
    output$dl_volcano_button <- button_download("dl_volcano")
    
    output$dl_volcano <- downloadHandler(filename = function() {"test.pdf"},
                                    content = function(file) {
                                        download_pdf(file); print(plots$volcano); dev.off()
                                    })
    
    output$res_box <- renderPlot({
        withProgress(message = "Plotting...", value = 0.5, {
        validate_res()
        plots$box <- plot_gene_boxplot(dds_run(), 
                                       rna_genes(),
                                       input$box_var, 
                                       input$viz_pal_categorical)
        plots$box
        })
    }, height = viz_plot_height, width = viz_plot_width)
    
    output$dl_box_button <- button_download("dl_box")
    
    output$dl_box <- downloadHandler(filename = function() {"test.pdf"},
                                         content = function(file) {
                                             download_pdf(file); print(plots$box); dev.off()
                                         })
    
    output$res_cluster <- renderPlot({
        withProgress(message = "Plotting...", value = 0.5, {
        validate_res()
        plots$cluster <- plot_sample_gene_mtx(dds_run(), 
                                              rna_genes(),
                                              input$viz_pal_continuous, 
                                              input$viz_pal_dir)
        plots$cluster
        })
    }, height = viz_plot_height, width = viz_plot_width)
    
    output$dl_cluster_button <- button_download("dl_cluster")
    
    output$dl_cluster <- downloadHandler(filename = function() {"test.pdf"},
                                     content = function(file) {
                                         download_pdf(file); print(plots$cluster); dev.off()
                                     })
    
    output$res_gsea <- renderPlot({
        withProgress(message = "Plotting...", value = 0.5, {
        validate_res()
        plots$gsea <- plot_deseq_gsea(res_to_gsea(res(), 
                                                  rna_pathways()))
        plots$gsea
        })
    }, height = viz_plot_height, width = viz_plot_width)
    
    output$dl_gsea_button <- button_download("dl_gsea")
    
    output$dl_gsea <- downloadHandler(filename = function() {"test.pdf"},
                                     content = function(file) {
                                         download_pdf(file); print(plots$gsea); dev.off()
                                     })
    
    
    # Help

    output$pdfview <- renderUI({
        tags$iframe(style = "height:700px; width:100%", 
                    src = "GDE_User_Guide_07102020.pdf")
    })
    
}

shinyApp(ui = ui, server = server)
