
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
    
    metadata <- eventReactive(input$count_start, {
        
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
        if (input$count_source == "Example") {
            
        } else if (input$count_source == "Upload") {
            validate(
                need(input$meta_input, 
                     "Please Upload Metadata")
            )
        } else {}
        paste0("Metadata uploaded: ", nrow(metadata()), " rows")
    })
    
    output$count_message <- renderText({
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
                        choices = colnames(metadata())),
            selectInput(inputId = "meta_file_col", 
                        label = "File Name Column",
                        choices = colnames(metadata())),
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
                       metadata(), 
                       input$meta_sample_col,
                       input$meta_file_col,
                       input$count_co)
    })
    
    output$cts_summary <- renderPrint({
        validate(
            need(try(cts()), ""),
            need(ncol(cts()) >= 1, "Name matching returns count matrix with 0 samples.\nPlease make sure the columns of sample names and files names are chosen correctly.")
        )
        trubble(cts())
    })

    # RNA
    
    observeEvent(input$rna_start, {
        
        library(DESeq2)
        
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
        
        output$rna_var <- renderUI({
            validate(
                need(try(rna_input()), "")
            )
            list(
                h4("Variable"),
                selectInput(inputId = "pca_var", 
                            label = "Categorical Variable for Heatmap, PCA, and Gene Boxplot",
                            choices = colnames(rna_input()[[1]]@colData))
            )
        })
        
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
        
        output$deseq_hm <- renderPlot({
            validate(
                need(input$pca_var, "Please Upload Data")
            )
            deseq_heatmap(rna_input()[[1]], 
                          input$pca_var,
                          input$palette_con,
                          input$palette_dir)
        }, 
        height = rna_plot_height, 
        width = rna_plot_width)
        
        output$deseq_pca <- renderPlot({
            validate(
                need(input$pca_var, "Please Upload Data")
            )
            if (is.numeric(rna_input()[[1]]@colData[[input$pca_var]])) {
                deseq_pca(rna_input()[[1]], 
                          input$pca_var, 
                          input$palette_con,
                          ifelse(input$palette_dir, 1, -1))
            } else {
                deseq_pca(rna_input()[[1]], 
                          input$pca_var, 
                          input$palette_cat)
            }
        }, 
        height = rna_plot_height, 
        width = rna_plot_width)
        
        output$deseq_ma <- renderPlot({
            deseq_ma(rna_input()[[2]],
                     input$p_co, 
                     input$lfc_co,
                     input$lfc_plot_lim)
        }, 
        height = rna_plot_height, 
        width = rna_plot_width)
        
        output$deseq_volcano <- renderPlot({
            deseq_volcano(rna_input()[[2]], 
                          input$p_co, 
                          input$lfc_co,
                          input$p_plot_lim,
                          input$lfc_plot_lim)
        }, 
        height = rna_plot_height, 
        width = rna_plot_width)
        
        output$deseq_table <- DT::renderDataTable({
            deseq_table(rna_input()[[2]], 
                        input$p_co, 
                        input$lfc_co)
        })
        
        rna_genes <- eventReactive(
            
            c(input$rna_gene_read1,
              input$rna_gene_read2,
              input$rna_genes_file),
            
            {
                
                if(input$rna_gene_ls_src == 'Use Top Genes') {
                    
                    get_rna_genes(rna_input()[[2]])[1:input$rna_gene_num]
                    
                } else if (input$rna_gene_ls_src == 'Manual Input') {
                    
                    validate(
                        need(input$rna_genes_man, "Please Input Gene List")
                    )
                    parse_rna_genes(input$rna_genes_man)
                    
                } else {
                    
                    validate(
                        need(input$rna_genes_file, "Please Upload Gene List")
                    )
                    readLines(input$rna_genes_file$datapath)
                }})
        
        
        output$deseq_box <- renderPlot({
            validate(
                need(input$pca_var, "Please Upload Data")
            )
            deseq_box(rna_input()[[1]], 
                      rna_genes(),
                      input$pca_var, 
                      input$palette_cat)
        }, 
        height = rna_plot_height, 
        width = rna_plot_width)
        
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

