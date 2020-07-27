
hmks_hs <- fgsea::gmtPathways("./data/h.all.v7.0.symbols.gmt")


tab_file <- tabPanel(
    
    title = "Data Input",
    fluid = TRUE,
    value = "v_file",
    
    sidebarLayout(
        sidebarPanel = sidebarPanel(
            conditionalPanel(
                condition = "input.count_panel == 1",
                uiOutput("upload_panel")
            ),
            conditionalPanel(
                condition = "input.count_source == 'Upload' &&
                             input.count_panel == 1",
                shinyFilesButton('count_upload', 
                                 label='Select RNA Count Files', 
                                 title='Please select HTSeq RNA count files', 
                                 multiple = TRUE),
                br(),
                fileInput(inputId = "meta_input",
                          label = " ",
                          buttonLabel = "Upload Metadata..")
            ),
            
            conditionalPanel(
                condition = "input.count_source == 'Select' &&
                             input.count_panel == 1",
                fileInput(inputId = "meta_input",
                          label = NULL,
                          buttonLabel = "Upload Metadata.."),
                br()
            ),
            conditionalPanel(
                condition = "input.count_panel == 1",
                uiOutput("cts_proc"),
                br()
            ),
            conditionalPanel(
                condition = "input.count_panel == 1",
                uiOutput("compute_button")
            ),
            conditionalPanel(
                condition = "input.count_panel == 2 ",
                h3("PCA Plot")
            ),
            conditionalPanel(
                condition = "input.count_panel == 3 ",
                h3("Sample Distance Heatmap")
            ),
            conditionalPanel(
                condition = "input.count_panel == 4 ",
                h3("DGE Analysis")
            ),
            conditionalPanel(
                condition = "input.count_panel == 2 ||
                input.count_panel == 3",
                uiOutput("pca_var_ui")
            ),
            conditionalPanel(
                condition = "input.count_panel == 2 ||
                input.count_panel == 3",
                uiOutput("color_ui")
            ),
            conditionalPanel(
                condition = "input.count_panel == 2",
                uiOutput("cluster_switch")
            ),
            conditionalPanel(
                condition = "input.count_panel == 2",
                uiOutput("cluster_ui")
            ),
            conditionalPanel(
                condition = "input.count_panel == 2 ||
                input.count_panel == 3",
                uiOutput("plot_size")
            ),
            conditionalPanel(
                condition = "input.count_panel == 4",
                uiOutput("dge_var"),
                uiOutput("dge_group1"),
                uiOutput("dge_group2"),
                uiOutput("dge_button")
            ),
            tags$head(tags$style(HTML("
                              .shiny-split-layout > div {
                                overflow: visible;
                              }
                              .progress-bar {
                              height: 20px;
                              }
                              .shiny-notification {
                              width: 200px;
                              top: 50%;
                              right: 30%;
                              position: fixed;
                              }
                              ")))
        ),
        
        mainPanel(
            tabsetPanel(
                id = "count_panel",
                tabPanel(
                    value = 1,
                    title = "Count Data Summary",
                    br(),
                    verbatimTextOutput("mt_message"),
                    verbatimTextOutput("count_message"),
                    verbatimTextOutput("cts_summary"),
                    verbatimTextOutput("compute_message")
                ),
                tabPanel(
                    value = 2,
                    title = "PCA",
                    br(),
                    plotOutput("pca")
                ),
                tabPanel(
                    value = 3,
                    title = "Heatmap",
                    br(),
                    plotOutput("hm")
                ),
                tabPanel(
                    value = 4,
                    title = "DGE Run",
                    br(),
                    verbatimTextOutput("dge_message")
                )
            )
        )
    )
)


tab_rna <- tabPanel(
    
    title = "DGE Visualization",
    fluid = TRUE,
    value = "v_rna",
    
    sidebarLayout(
        
        sidebarPanel = sidebarPanel(
            h3("RNA DGE Visualization"),
            br(),
            
            conditionalPanel(
                condition = "input.data_source == 'Upload'",
                fileInput(inputId = "rna_input",
                          label = NULL,
                          buttonLabel = "Browse..")
            ),
            
            conditionalPanel(
                condition = "input.data_source == 'Select'",
                selectInput(inputId = "rna_select", 
                            label = "RNA-seq results",
                            choices = list.files("./large_data/"))
            ),
            br(),
            
            conditionalPanel(
                condition = "input.rna_panel == 1 || 
                            input.rna_panel == 2 || 
                            input.rna_panel == 6",
                
                uiOutput("rna_var"),
            ),
            
            
            conditionalPanel(
                condition = "input.rna_panel == 3 || 
                            input.rna_panel == 4 || 
                            input.rna_panel == 5",
                
                h4("Differential Gene Expression Parameters"),
                splitLayout(numericInput("p_co", 
                                         label = "Adjusted P-value Cutoff", 
                                         value = 0.05),
                            numericInput("lfc_co", 
                                         label = "Log2 Fold Change Cutoff", 
                                         value = 1)),
                br()
            ),
            
            conditionalPanel(
                condition = "input.rna_panel == 6 || 
                            input.rna_panel == 7",
                
                h4("Gene List"),
                radioGroupButtons(inputId = "rna_gene_ls_src",
                                  label = NULL,
                                  choices = c("Use Top Genes",
                                              "Manual Input",
                                              "Upload File"),
                                  justified = TRUE),
                conditionalPanel(
                    condition = "input.rna_gene_ls_src == 'Use Top Genes'",
                    splitLayout(sliderInput(inputId = "rna_gene_num",
                                            label = "Number of Genes", 
                                            min = 1, max = 50, value = 6),
                                actionButton(
                                    inputId = "rna_gene_read1",
                                    label = "Plot",
                                    icon = icon("check"),
                                    style = "color: white; 
                                    background-color: #737373;
                                    margin-top: 25px;
                                    float:right;
                                    margin-right: 5px;"),
                                cellWidths = c("75%", "25%")
                    )
                ),
                conditionalPanel(
                    condition = "input.rna_gene_ls_src == 'Manual Input'",
                    splitLayout(textInput("rna_genes_man", 
                                          label = NULL, 
                                          value = ""),
                                
                                actionButton(
                                    inputId = "rna_gene_read2",
                                    label = "Plot",
                                    icon = icon("check"),
                                    style = "color: white; background-color: #737373;
                            float:right; margin-right: 5px;"),
                                cellWidths = c("75%", "25%")
                    )
                    
                ),
                conditionalPanel(
                    condition = "input.rna_gene_ls_src == 'Upload File'",
                    fileInput(inputId = "rna_genes_file",
                              label = NULL,
                              buttonLabel = "Browse..")
                ),
                br()
            ),
            
            conditionalPanel(
                condition = "input.rna_panel == 8",
                
                h4("Pathway List"),
                splitLayout(radioGroupButtons(inputId = "rna_pathway_src",
                                              label = NULL,
                                              choices = c("Use Hallmarks",
                                                          "Upload File"),
                                              justified = TRUE),
                            cellWidths = "67%"),
                
                conditionalPanel(
                    condition = "input.rna_pathway_src == 'Upload File'",
                    fileInput(inputId = "rna_pathway_file",
                              label = NULL,
                              buttonLabel = "Browse..")
                )
            ),
            
            conditionalPanel(
                condition = "input.rna_panel == 1 || 
                            input.rna_panel == 2 ||
                            input.rna_panel == 6 || 
                            input.rna_panel == 7",
                
                h4("Plotting Parameters"),
                splitLayout(selectInput(inputId = "palette_cat", 
                                        label = "Categorical Palette",
                                        choices = rownames(brewer.pal.info[brewer.pal.info$category == "qual",]),
                                        selected = "Set2"),
                            selectInput(inputId = "palette_con", 
                                        label = "Continuous Palette",
                                        choices = rownames(brewer.pal.info[brewer.pal.info$category != "qual",]),
                                        selected = "Spectral")),
                materialSwitch(
                    inputId = "palette_dir",
                    label = "Reverse Scale Color Direction",
                    value = FALSE,
                    right = TRUE
                )
            ),
            
            conditionalPanel(
                condition = "input.rna_panel == 3 || 
                            input.rna_panel == 4",
                
                h4("Plotting Parameters"),
                splitLayout(numericInput("p_plot_lim", 
                                         label = "Adjusted P-value Squash", 
                                         value = 5),
                            numericInput("lfc_plot_lim", 
                                         label = "Log2 Fold Change Squash", 
                                         value = 5)),
            ),
            
            conditionalPanel(
                condition = "input.rna_panel != 5",
                
                splitLayout(numericInput("rna_plot_height", 
                                         "Plot Height (px)", 
                                         value = 600),
                            numericInput("rna_plot_width", 
                                         "Plot Width (px)", 
                                         value = 800),
                            cellWidths = c("50%", "50%"))
            ),
            
            tags$head(tags$style(HTML("
                              .shiny-split-layout > div {
                                overflow: visible;
                              }
                              .progress-bar {
                              height: 20px;
                              }
                              .shiny-notification {
                              width: 200px;
                              top: 50%;
                              right: 30%;
                              position: fixed;
                              }
                              ")))
        ),
        
        mainPanel(
            tabsetPanel(
                id = "rna_panel",
                tabPanel(
                    value = 5,
                    title = "Table",
                    br(),
                    DT::dataTableOutput("res_table")
                ),
                tabPanel(
                    value = 3,
                    title = "MA Plot",
                    br(),
                    plotOutput("deseq_ma", width = "100%") 
                ),
                tabPanel(
                    value = 4,
                    title = "Volcano Plot",
                    br(),
                    plotOutput("deseq_volcano")
                ),
                tabPanel(
                    value = 6,
                    title = "Gene Boxplot",
                    br(),
                    plotOutput("deseq_box", width = "100%")
                ),
                tabPanel(
                    value = 7,
                    title = "Gene Clustering",
                    br(),
                    plotOutput("deseq_cluster", width = "100%")
                ),
                tabPanel(
                    value = 8,
                    title = "GSEA",
                    br(),
                    plotOutput("deseq_gsea", width = "100%")
                )
            )
        )
    )
)

tab_help <- tabPanel(
    
    title = "Help",
    fluid = TRUE,
    value = "v_help",
    
    sidebarLayout(
        
        sidebarPanel = sidebarPanel(
            
        ),
        
        mainPanel(
            uiOutput("pdfview")
        )
    )
)

tab_about <- tabPanel(
    
    title = "About",
    fluid = TRUE,
    value = "v_about",
    
    sidebarLayout(
        
        sidebarPanel = sidebarPanel(
            
            
            
        ),
        
        mainPanel(
            br(),
            h4("Authors:"),
            br(),
            h4("Mingchu Xu"),
            br(),
            h4("Xiaogang (Sean) Wu"),
            br(),
            h4("Jianhua (John) Zhang")
        )
    )
    
)

