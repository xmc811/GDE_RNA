# ----------
# Basic shiny modules
# ----------

# Buttons
button_cts_upload <- function() {
    actionButton(inputId = "cts_upload_click", 
                 label = "Upload", 
                 icon = icon("upload"), 
                 style = "color:white; background-color:#0570b0; float:right; margin-right:5px; width:150px; max-width:90%;")
}

button_cts_process <- function() {
    actionButton(inputId = "cts_process_click", 
                 label = "Pre-Process", 
                 icon = icon("files-o"), 
                 style = "color:white; background-color:#0570b0; width:140px;")
}

button_cts_compute <- function() {
    actionButton(inputId = "cts_compute_click", 
                 label = "Compute", 
                 icon = icon("calculator"), 
                 style = "color:white; background-color:#2ca25f; width:140px;")
}

button_assign_cluster <- function() {
    actionButton(inputId = "assign_cluster_click", 
                 label = "Assign K-means Cluster", 
                 icon = icon("calculator"), 
                 style = "color:white; background-color:#2ca25f; margin-top:25px; float:right; margin-right:5px;")
}

button_dge <- function() {
    actionButton(inputId = "dge_click", 
                 label = "Start DGE",
                 icon = icon("bar-chart"),
                 style = "color:white; background-color:#0570b0;")
}

button_read_gene <- function(value) {
    actionButton(inputId = paste0("rna_gene_read", as.character(value)),
                 label = "Plot",
                 icon = icon("picture-o"),
                 style = "color:white; background-color:#737373; margin-top:25px; float:right; margin-right:5px;")
}

# Radio Selections
radio_source <- function() {
    radioGroupButtons(inputId = "cts_source",
                      label  = NULL,
                      choices = c("Example","Upload","Select"),
                      justified = TRUE)
}

radio_gene_source <- function() {
    radioGroupButtons(inputId = "gene_source",
                      label  = NULL,
                      choices = c("Use Top Genes","Manual Input","Upload File"),
                      justified = TRUE)
}

radio_pathway_source <- function() {
    radioGroupButtons(inputId = "pathway_source",
                      label = NULL,
                      choices = c("Use Hallmarks","Upload File"),
                      justified = TRUE)
}

# File Uploading
file_cts_upload <- function() {
    shinyFilesButton(id = 'cts_files', 
                     label = 'Select RNA Count Files', 
                     title = 'Please select HTSeq RNA count files', 
                     multiple = TRUE)
}

file_meta_upload <- function() {
    fileInput(inputId = "meta_file",
              label = " ",
              buttonLabel = "Upload Metadata..")
}

file_gene_pathway_upload <- function(id) {
    fileInput(inputId = id,
              label = NULL,
              buttonLabel = "Browse..")
}


# Switches
switch_palette_dir <- function(tab) {
    materialSwitch(inputId = paste0(tab, "pal_dir"), 
                   label = "Reverse Scale Color Direction",
                   right = TRUE)
}

switch_cluster_mode <- function() {
    materialSwitch(inputId = "cluster_sw", 
                   label = "K-means Clustering Mode",
                   right = TRUE)
}

# Selections
select_palette <- function(tab, type) {
    all_pals <- rownames(brewer.pal.info)
    if (type == "categorical") {
        pals <- all_pals[brewer.pal.info$category == "qual"]
        selected <- "Set2"
    } else {
        pals <- all_pals[brewer.pal.info$category != "qual"]
        selected <- "Spectral"
    }
    selectInput(inputId = paste0(tab, "pal_", type), 
                label = paste(str_to_title(type), "Palette"),
                choices = pals,
                selected = selected)
}

# Numbers
number_cts_cutoff <- function(){
    numericInput(inputId = "cts_cutoff", label = NULL, value = 10)
}

number_plot_size <- function(tab, type) {
    value <- ifelse(type == "height", 600, 800)
    numericInput(inputId = paste0(tab, "plot_", type), 
                 label = paste("Plot", str_to_title(type), "(px)"), 
                 value = value)
}

number_n_cluster <- function() {
    numericInput(inputId = "n_cluster",
                 label = "No. of K-means Clusters",
                 value = 3)
}

number_p_cutoff <- function() {
    numericInput("p_co", 
                 label = "Adjusted P-value Cutoff", 
                 value = 0.05)
}

number_lfc_cutoff <- function() {
    numericInput("lfc_co", 
                 label = "Log2 Fold Change Cutoff", 
                 value = 1)
}

number_p_limit <- function() {
    numericInput("p_plot_lim", 
                 label = "Adjusted P-value Squash", 
                 value = 5)
}

number_lfc_limit <- function() {
    numericInput("lfc_plot_lim", 
                 label = "Log2 Fold Change Squash", 
                 value = 5)
}

# Sliders
slider_gene_num <- function() {
    sliderInput(inputId = "rna_gene_num",
                label = "Number of Genes", 
                min = 0, max = 50, value = 6)
}

# Texts
text_gene_names <- function() {
    textInput("rna_genes_man", 
              label = NULL, 
              value = "")
}

# Messages
message_plot_size <- function(type) {
    paste("Plot", type, "shouldn't exceed 4000px.")
}

# ----------
# UI components
# ----------

upload_widgets <- splitLayout(radio_source(),
                              button_cts_upload(),
                              cellWidths = c("67%", "33%"))

palette_widgets <- splitLayout(select_palette("", "categorical"),
                               select_palette("", "continuous"))

viz_palette_widgets <- splitLayout(select_palette("viz_", "categorical"),
                                   select_palette("viz_", "continuous"))

cluster_widgets <- splitLayout(number_n_cluster(),
                               button_assign_cluster(),
                               cellWidths = c("33%", "67%"))

plot_size_widgets <- splitLayout(number_plot_size("", "height"),
                                 number_plot_size("", "width"),
                                 cellWidths = c("50%", "50%"))

viz_plot_size_widgets <- splitLayout(number_plot_size("viz_", "height"),
                                     number_plot_size("viz_", "width"),
                                     cellWidths = c("50%", "50%"))

dge_cutoff_widgets <- splitLayout(number_p_cutoff(),
                                  number_lfc_cutoff())

dge_plot_limit_widgets <- splitLayout(number_p_limit(),
                                      number_lfc_limit())

top_gene_widgets <- splitLayout(slider_gene_num(),
                                button_read_gene(1),
                                cellWidths = c("75%", "25%"))

manual_gene_widgets <- splitLayout(text_gene_names(),
                                   button_read_gene(2),
                                   cellWidths = c("75%", "25%"))


gene_selection_widgets <- list(h4("Gene List"),
                               radio_gene_source(),
                               conditionalPanel(
                                   condition = "input.gene_source == 'Use Top Genes'",
                                   top_gene_widgets),
                               conditionalPanel(
                                   condition = "input.gene_source == 'Manual Input'",
                                   manual_gene_widgets),
                               conditionalPanel(
                                   condition = "input.gene_source == 'Upload File'",
                                   file_gene_pathway_upload("rna_genes_file")),
                               br()
                               )

pathway_selection_widgets <- list(h4("Pathway List"),
                                  splitLayout(radio_pathway_source(), 
                                              cellWidths = "67%"),
                                  conditionalPanel(
                                      condition = "input.pathway_source == 'Upload File'",
                                      file_gene_pathway_upload("rna_pathway_file"))
                                  )





