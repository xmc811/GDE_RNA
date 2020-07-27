
source("loading.R")
source("helpers.R")
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
    
    tab_file,
    tab_rna,
    tab_help,
    tab_about,
    
    tags$head(tags$link(rel="stylesheet", 
                        type="text/css", 
                        href="style.css"))
)

# server function
server <- function(input, output, session) {
    
    output$upload_panel <- renderUI({
        
        list(
            h3("RNA-seq Raw Data Input"),
            br(),
            h4("Data Source"),
            splitLayout(radioGroupButtons(inputId = "count_source",
                                          label = NULL,
                                          choices = c("Example","Upload","Select"),
                                          justified = TRUE),
                        actionButton(
                            inputId = "count_start",
                            label = "Upload",
                            icon = icon("bar-chart"),
                            style = "color: white; background-color: #0570b0;
                            float:right; margin-right: 5px;"),
                        
                        cellWidths = c("67%", "33%"))
        )
    })
    
    shinyFileChoose(input, 
                    id = 'count_upload', 
                    roots = c(home = "~/"),
                    filetypes = "tsv")
    
    mt_raw <- eventReactive(input$count_start, {
        
        if (input$count_source == "Example") {
            
            read_csv("./data/metadata.csv")
            
        } else if (input$count_source == "Upload") {
            
            read_csv(input$meta_input$datapath)
            
        } else {}
    })
    
    cts_raw <- eventReactive(input$count_start,{
        
        if (input$count_source == "Example") {
            
            readRDS("./data/example_mtx.rds")
            
        } else if (input$count_source == "Upload"){
            
            roots <- c(home = "~/")
            files_path <- parseFilePaths(roots, 
                                         input$count_upload)$datapath
            htseq_to_mtx(files_path)
            
        } else {}
    })
    
    output$mt_message <- renderText({
        validate(
            need(input$count_source,"")
        )
        if (input$count_source == "Example") {
            
        } else if (input$count_source == "Upload") {
            validate(
                need(input$meta_input, 
                     "Please Upload Metadata")
            )
        } else {}
        paste0("Metadata uploaded: ", nrow(mt_raw()), " rows")
    })
    
    output$count_message <- renderText({
        validate(
            need(input$count_source,"")
        )
        if (input$count_source == "Example") {
            
        } else if (input$count_source == "Upload") {
            validate(
                need(input$count_upload, 
                     "Please Upload Count Data")
            )
        } else {}
        get_count_message(cts_raw())
    })
    
    output$cts_proc <- renderUI({
        validate(
            need(input$count_source,"")
        )
        if (input$count_source == "Example") {
            
        } else if (input$count_source == "Upload") {
            validate(
                need(input$meta_input, "Please Upload Metadata"),
                need(input$count_upload, "Please Upload Count Data")
            )
        } else {}
        list(
            h4("File-Sample Name Matching"),
            selectInput(inputId = "meta_sample_col", 
                        label = "Sample Column",
                        choices = colnames(mt_raw())),
            selectInput(inputId = "meta_file_col", 
                        label = "File Name Column",
                        choices = colnames(mt_raw())),
            h4("Count cutoff (count data row sums < n)"),
            numericInput(inputId = "count_co",
                         label = NULL,
                         value = 10),
            actionButton(
                inputId = "cts_start",
                label = "Pre-Process",
                icon = icon("bar-chart"),
                style = "color: white; background-color: #2ca25f")
        )
    })
    
    cts <- eventReactive(input$cts_start, {
        mtx_name_match(cts_raw(), 
                       mt_raw(), 
                       input$meta_sample_col,
                       input$meta_file_col,
                       input$count_co)
    })
    
    mt <- reactiveVal()
    
    observeEvent(input$cts_start,{
        
        mt(filter_mt(cts(), mt_raw()))
    })
    
    
    output$cts_summary <- renderPrint({
        validate(
            need(try(cts()), ""),
            need(ncol(cts()) >= 1, "Name matching returns count matrix with 0 samples.\nPlease make sure the columns of sample names and files names are chosen correctly.")
        )
        trubble(cts())
    })
    
    output$compute_button <- renderUI({
        validate(
            need(try(cts()),""),
            need(try(mt()),"")
        )
        actionButton(
            inputId = "cts_compute",
            label = "Compute",
            icon = icon("bar-chart"),
            style = "color: white; background-color: #0570b0;")
        
    })
    
    dds <- reactiveVal()
    vsd <- reactiveVal()
    
    observeEvent(input$cts_compute, {
        
        dds(cts_to_dds(cts(), mt()))
        vsd(DESeq2::vst(dds(), blind = FALSE))

    })

    
    output$compute_message <- renderText({
        validate(
            need(!is.null(dds()), ""),
            need(!is.null(vsd()), "")
        )
        paste0("Initial computation done.\nYou can explore your data in other panels.")
    })
    
    # PCA Plot
    
    output$pca_var_ui <- renderUI({
        validate(
            need(!is.null(vsd()), "")
        )
        list(
            selectInput(inputId = "pca_var_ch", 
                        label = "Variable for PCA Plot and Heatmap",
                        choices = colnames(vsd()@colData))
        )
    })
    
    output$color_ui <- renderUI({
        list(
            splitLayout(selectInput(inputId = "pal_cat", 
                                    label = "Categorical Palette",
                                    choices = rownames(brewer.pal.info[brewer.pal.info$category == "qual",]),
                                    selected = "Set2"),
                        selectInput(inputId = "pal_con", 
                                    label = "Continuous Palette",
                                    choices = rownames(brewer.pal.info[brewer.pal.info$category != "qual",]),
                                    selected = "Spectral")),
            materialSwitch(
                inputId = "pal_dir",
                label = "Reverse Scale Color Direction",
                value = FALSE,
                right = TRUE)
        )
    })
    
    output$cluster_switch <- renderUI({
            materialSwitch(
                inputId = "cluster_sw",
                label = "K-means Clustering Mode",
                value = FALSE,
                right = TRUE)
    })
    
    output$cluster_ui <- renderUI({
        validate(
            need(input$cluster_sw == TRUE, "")
        )
        list(
            splitLayout(numericInput(inputId = "n_cluster",
                                     label = "No. of K-means Clusters",
                                     value = 3),
                        actionButton(
                            inputId = "assign_clu",
                            label = "Assign K-means Cluster",
                            icon = icon("bar-chart"),
                            style = "color: white; background-color: #2ca25f;margin-top: 25px;
                                    float:right;
                                    margin-right: 5px;"),
                        cellWidths = c("33%", "67%")
                        )
        )
    })
    
    
    km_res <- reactive({
        
        vsd_km(vsd(), input$n_cluster)
        
    })
    
    observeEvent(input$assign_clu, {
        
        dds(assign_km_clu(dds(), km_res()))
        vsd(assign_km_clu(vsd(), km_res()))
        mt(assign_km_clu_col(mt(), km_res()))
        
    })
    
    # Plot Size
    
    output$plot_size <- renderUI({
            splitLayout(numericInput("plot_height", 
                                 "Plot Height (px)", 
                                 value = 600),
                    numericInput("plot_width", 
                                 "Plot Width (px)", 
                                 value = 800),
                    cellWidths = c("50%", "50%"))
    })
    
    plot_height <- reactive({
        validate(
            need(input$plot_height < 4000, 
                 "Plot height shouldn't exceed 4000px.")
        )
        return(input$plot_height)
    })
    
    plot_width <- reactive({
        validate(
            need(input$plot_width < 4000, 
                 "Plot width shouldn't exceed 4000px.")
        )
        return(input$plot_width)
    })

    
    output$pca <- renderPlot({
        validate(
            need(vsd(), "VSD object not computed. PCA not available.")
        )
        if (input$cluster_sw == TRUE) {
            
            vsd_pca_km(vsd(), km_res(), input$pal_cat)
            
        } else if (is.numeric(dds()@colData[[input$pca_var_ch]])) {
            vsd_pca(dds(), 
                    vsd(), 
                    input$pca_var_ch, 
                    input$pal_con,
                    ifelse(input$pal_dir, 1, -1))
        } else {
            vsd_pca(dds(), 
                    vsd(), 
                    input$pca_var_ch, 
                    input$pal_cat)
        }
    },
    height = plot_height,
    width = plot_width
    )
    
    # DGE
    
    output$dge_var <- renderUI({
        validate(
            need(try(mt()),"")
        )
            selectInput(inputId = "dge_var_ch", 
                        label = "Variable for DGE",
                        choices = colnames(mt()))
    })
    
    output$dge_group1 <- renderUI({
        validate(
            need(try(mt()),"")
        )
        selectInput(inputId = "dge_g1", 
                    label = "Group 1",
                    choices = unique(mt()[[input$dge_var_ch]]))
        
    })
    
    output$dge_group2 <- renderUI({
        validate(
            need(try(mt()),"")
        )
        selectInput(inputId = "dge_g2", 
                    label = "Group 2",
                    choices = setdiff(unique(mt()[[input$dge_var_ch]]),
                                      input$dge_g1))
    })
    
    output$dge_button <- renderUI({
        validate(
            need(try(cts()),""),
            need(try(mt()),"")
        )
        actionButton(
            inputId = "dge_start",
            label = "Start DGE",
            icon = icon("bar-chart"),
            style = "color: white; background-color: #0570b0;")
    })
    
    dds_run <- reactiveVal()
    res <- reactiveVal()
    
    observeEvent(input$dge_start,{
        
        dds_run(cts_to_dds(cts(), mt(), input$dge_var_ch))
        res(DESeq2::results(dds_run(), 
                            contrast = c(input$dge_var_ch, 
                                         input$dge_g1, 
                                         input$dge_g2)))
    })
    
    
    output$dge_message <- renderText({
        validate(
            need(try(res()), "")
        )
        paste0("Differential gene expression (DGE) analysis done\n\n",
               res()@elementMetadata[2,2] %>%
                   str_remove("^.*:") %>%
                   str_trim(), "\n\n",
               "You can visualize your data in the 'DGE Visualization' tab.")
    })
    
    # Sample Distances
    
    output$hm <- renderPlot({
        validate(
            need(vsd(), "VSD object not computed. Heatmap not available.")
        )
        vsd_hm(vsd(), 
               input$pca_var_ch, 
               input$pal_con, 
               input$pal_dir)
    },
    height = plot_height,
    width = plot_width
    )
    
    
    # DGE Visualization UIs
    
    output$viz_plot_size <- renderUI({
        validate(
            need(try(res()), "")
        )
        splitLayout(numericInput("viz_plot_height", 
                                 "Plot Height (px)", 
                                 value = 600),
                    numericInput("viz_plot_width", 
                                 "Plot Width (px)", 
                                 value = 800),
                    cellWidths = c("50%", "50%"))
    })
    
    
    output$dge_params <- renderUI({
        validate(
            need(try(res()), "")
        )
        list(
            h4("Differential Gene Expression Parameters"),
            splitLayout(numericInput("p_co", 
                                     label = "Adjusted P-value Cutoff", 
                                     value = 0.05),
                        numericInput("lfc_co", 
                                     label = "Log2 Fold Change Cutoff", 
                                     value = 1))
        )
    })
    
    
    output$squash_params <- renderUI({
        validate(
            need(try(res()), "")
        )
        list(
            h4("Plotting Parameters"),
            splitLayout(numericInput("p_plot_lim", 
                                     label = "Adjusted P-value Squash", 
                                     value = 5),
                        numericInput("lfc_plot_lim", 
                                     label = "Log2 Fold Change Squash", 
                                     value = 5))
        )
    })
    
    
    output$rna_var <- renderUI({
        validate(
            need(try(dds_run()), "")
        )
        list(
            h4("Variable"),
            selectInput(inputId = "box_var", 
                        label = "Categorical Variable for Gene Boxplot",
                        choices = colnames(dds_run()@colData))
        )
    })
    
    
    output$viz_color_ui <- renderUI({
        validate(
            need(try(res()), "")
        )
        list(
            splitLayout(selectInput(inputId = "viz_pal_cat", 
                                    label = "Categorical Palette",
                                    choices = rownames(brewer.pal.info[brewer.pal.info$category == "qual",]),
                                    selected = "Set2"),
                        selectInput(inputId = "viz_pal_con", 
                                    label = "Continuous Palette",
                                    choices = rownames(brewer.pal.info[brewer.pal.info$category != "qual",]),
                                    selected = "Spectral")),
            materialSwitch(
                inputId = "viz_pal_dir",
                label = "Reverse Scale Color Direction",
                value = FALSE,
                right = TRUE)
        )
    })
    
    output$viz_gene_ui <- renderUI({
        validate(
            need(try(res()), "")
        )
        list(
            h4("Gene List"),
            radioGroupButtons(inputId = "viz_gene_src",
                              label = NULL,
                              choices = c("Use Top Genes",
                                          "Manual Input",
                                          "Upload File"),
                              justified = TRUE),
            conditionalPanel(
                condition = "input.viz_gene_src == 'Use Top Genes'",
                splitLayout(sliderInput(inputId = "rna_gene_num",
                                        label = "Number of Genes", 
                                        min = 1, max = 50, value = 6),
                            actionButton(
                                inputId = "rna_gene_read1",
                                label = "Plot",
                                icon = icon("check"),
                                style = "color: white; background-color: #737373; margin-top: 25px; float:right; margin-right: 5px;"),
                            cellWidths = c("75%", "25%")
                            )
                ),
            conditionalPanel(
                condition = "input.viz_gene_src == 'Manual Input'",
                splitLayout(textInput("rna_genes_man", 
                                      label = NULL, 
                                      value = ""),
                            actionButton(
                                inputId = "rna_gene_read2",
                                label = "Plot",
                                icon = icon("check"),
                                style = "color: white; background-color: #737373; float:right; margin-right: 5px;"),
                            cellWidths = c("75%", "25%")
                            )
                ),
            conditionalPanel(
                condition = "input.viz_gene_src == 'Upload File'",
                fileInput(inputId = "rna_genes_file",
                          label = NULL,
                          buttonLabel = "Browse..")
                ),
            br()
            )
        })
    
    
    # DGE Visualization Parameters
    
    viz_plot_height <- reactive({
        validate(
            need(input$plot_height < 4000, 
                 "Plot height shouldn't exceed 4000px.")
        )
        return(input$viz_plot_height)
    })
    
    viz_plot_width <- reactive({
        validate(
            need(input$plot_width < 4000, 
                 "Plot width shouldn't exceed 4000px.")
        )
        return(input$viz_plot_width)
    })
    
    
    rna_genes <- eventReactive(
        c(input$rna_gene_read1,
          input$rna_gene_read2,
          input$rna_genes_file),
        {
            if(input$viz_gene_src == 'Use Top Genes') {
                get_rna_genes(res())[1:input$rna_gene_num]
            } else if (input$viz_gene_src == 'Manual Input') {
                validate(
                    need(input$rna_genes_man, "Please Input Gene List")
                )
                parse_rna_genes(input$rna_genes_man)
            } else {
                validate(
                    need(input$rna_genes_file, "Please Upload Gene List")
                )
                readLines(input$rna_genes_file$datapath)
            }}
        )
    
    # DGE Visualization Panels
    
    output$res_table <- DT::renderDataTable({
        validate(
            need(try(res()), "No DGE results. Table not available.")
        )
        deseq_table(res(), 
                    input$p_co, 
                    input$lfc_co)
    })
    
    
    
    output$res_ma <- renderPlot({
        validate(
            need(try(res()), "No DGE results. Plot not available.")
        )
        deseq_ma(res(),
                 input$p_co, 
                 input$lfc_co,
                 input$lfc_plot_lim)
    }, 
    height = viz_plot_height, 
    width = viz_plot_width)
    
    
    
    output$res_volcano <- renderPlot({
        validate(
            need(try(res()), "No DGE results. Plot not available.")
        )
        deseq_volcano(res(), 
                      input$p_co, 
                      input$lfc_co,
                      input$p_plot_lim,
                      input$lfc_plot_lim)
    }, 
    height = viz_plot_height, 
    width = viz_plot_width)
    
    
    output$res_box <- renderPlot({
        validate(
            need(try(res()), "No DGE results. Plot not available.")
        )
        deseq_box(dds_run(), 
                  rna_genes(),
                  input$box_var, 
                  input$viz_pal_cat)
    }, 
    height = viz_plot_height, 
    width = viz_plot_width)
    
    
    # RNA
    
    observeEvent(input$rna_start, {
        
        rna_input <- if(input$data_source == "Example") {
            reactive({readRDS("./large_data/rnaseq.rds")})
        } else if (input$data_source == "Upload"){
            reactive({
                validate(
                    need(input$rna_input, 
                         "Please Upload Data")
                )
                readRDS(input$rna_input$datapath)})
        } else {
            reactive({
                readRDS(paste0("./large_data/",input$rna_select))})
        }
        
        rna_plot_height <- reactive({
            validate(
                need(input$rna_plot_height < 4000, "Plot height shouldn't exceed 4000px.")
            )
            return(input$rna_plot_height)
        })
        
        rna_plot_width <- reactive({
            validate(
                need(input$rna_plot_width < 4000, "Plot width shouldn't exceed 4000px.")
            )
            return(input$rna_plot_width)
        })
        
        
        
        
        
        
        output$deseq_cluster <- renderPlot({
            deseq_cluster(rna_input()[[1]], 
                          rna_genes(),
                          input$palette_con, 
                          input$palette_dir)
        }, 
        height = rna_plot_height, 
        width = rna_plot_width)
        
        
        rna_pathways <- reactive({
            
            if(input$rna_pathway_src == 'Use Hallmarks') {
                
                hmks_hs
                
            } else {
                
                validate(
                    need(input$rna_pathway_file, "Please Upload Pathway File")
                )
                read_csv(input$rna_pathway_file$datapath, 
                         col_names = FALSE) %>%
                    df_to_signature()
                
            }})
        
        output$deseq_gsea <- renderPlot({
            deseq_gsea(rna_input()[[2]],
                       rna_pathways())
        }, 
        height = rna_plot_height, 
        width = rna_plot_width)
        
    })
    
    output$pdfview <- renderUI({
        tags$iframe(style = "height:700px; width:100%", 
                    src = "GDE_User_Guide_07102020.pdf")
    })
    
}

shinyApp(ui = ui, server = server)
