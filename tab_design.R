
tab_file <- tabPanel(
    
    title = "Data Input & Exploration",
    fluid = TRUE,
    value = "v_file",
    
    sidebarLayout(
        sidebarPanel = sidebarPanel(
            conditionalPanel(
                condition = "input.count_panel == 1",
                uiOutput("upload_panel")
            ),
            conditionalPanel(
                condition = "input.cts_source == 'Select' &&
                             input.count_panel == 1",
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
                h3("Data Exploration")
            ),
            conditionalPanel(
                condition = "input.count_panel == 3 ",
                h3("Data Exploration")
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
                uiOutput("dge_var_ui"),
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
                    title = "Principle Component Analysis",
                    br(),
                    uiOutput("dl_pca_button"),
                    plotOutput("pca"),
                ),
                tabPanel(
                    value = 3,
                    title = "Sample Distance Heatmap",
                    br(),
                    uiOutput("dl_hm_button"),
                    plotOutput("hm")
                ),
                tabPanel(
                    value = 4,
                    title = "Differential Gene Expression Analysis",
                    br(),
                    verbatimTextOutput("dge_message"),
                    verbatimTextOutput("dge_error")
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
                condition = "input.rna_panel == 6",
                uiOutput("rna_var"),
            ),
            conditionalPanel(
                condition = "input.rna_panel == 3 || 
                            input.rna_panel == 4 || 
                            input.rna_panel == 5",
                uiOutput("dge_params")
            ),
            conditionalPanel(
                condition = "input.rna_panel == 6 || 
                            input.rna_panel == 7",
                uiOutput("viz_gene_ui")
            ),
            conditionalPanel(
                condition = "input.rna_panel == 6 || 
                            input.rna_panel == 7",
                uiOutput("viz_color_ui")
            ),
            conditionalPanel(
                condition = "input.rna_panel == 3 || 
                            input.rna_panel == 4",
                uiOutput("squash_params")
            ),
            conditionalPanel(
                condition = "input.rna_panel == 8",
                uiOutput("pathway_ui")
            ),
            conditionalPanel(
                condition = "input.rna_panel != 5",
                uiOutput("viz_plot_size")
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
                    uiOutput("dl_ma_button"),
                    plotOutput("res_ma") 
                ),
                tabPanel(
                    value = 4,
                    title = "Volcano Plot",
                    br(),
                    uiOutput("dl_volcano_button"),
                    plotOutput("res_volcano")
                ),
                tabPanel(
                    value = 6,
                    title = "Gene Boxplot",
                    br(),
                    uiOutput("dl_box_button"),
                    plotOutput("res_box")
                ),
                tabPanel(
                    value = 7,
                    title = "Gene Clustering",
                    br(),
                    uiOutput("dl_cluster_button"),
                    plotOutput("res_cluster")
                ),
                tabPanel(
                    value = 8,
                    title = "GSEA",
                    br(),
                    uiOutput("dl_gsea_button"),
                    plotOutput("res_gsea")
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

