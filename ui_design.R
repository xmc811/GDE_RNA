# ----------
# Basic shiny modules
# ----------

# Buttons
button_cts_upload <- function() {
    actionButton(inputId = "cts_upload_click", 
                 label = "Upload", 
                 icon = icon("upload"), 
                 style = "color: white; background-color: #0570b0; float:right; margin-right: 5px;")
}

button_cts_process <- function() {
    actionButton(inputId = "cts_process_click", 
                 label = "Pre-Process", 
                 icon = icon("files-o"), 
                 style = "color: white; background-color: #2ca25f")
}

button_cts_compute <- function() {
    actionButton(inputId = "cts_compute_click", 
                 label = "Compute", 
                 icon = icon("calculator"), 
                 style = "color: white; background-color: #0570b0;")
}

button_assign_cluster <- function() {
    actionButton(inputId = "assign_cluster_click", 
                 label = "Assign K-means Cluster", 
                 icon = icon("calculator"), 
                 style = "color: white; background-color: #2ca25f;margin-top: 25px; float:right; margin-right: 5px;")
}

button_dge <- function() {
    actionButton(inputId = "dge_click", 
                 label = "Start DGE",
                 icon = icon("bar-chart"),
                 style = "color: white; background-color: #0570b0;")
}

# Radio Selections
radio_source <- function() {
    radioGroupButtons(inputId = "cts_source",
                      label  = NULL,
                      choices = c("Example","Upload","Select"),
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


# Switches
switch_palette_dir <- function() {
    materialSwitch(inputId = "pal_dir", 
                   label = "Reverse Scale Color Direction",
                   right = TRUE)
}

switch_cluster_mode <- function() {
    materialSwitch(inputId = "cluster_sw", 
                   label = "K-means Clustering Mode",
                   right = TRUE)
}

#Selections
select_palette <- function(type) {
    all_pals <- rownames(brewer.pal.info)
    if (type == "categorical") {
        pals <- all_pals[brewer.pal.info$category == "qual"]
        selected <- "Set2"
    } else {
        pals <- all_pals[brewer.pal.info$category != "qual"]
        selected <- "Spectral"
    }
    selectInput(inputId = paste0("pal_", type), 
                label = paste(str_to_title(type), "Palette"),
                choices = pals,
                selected = selected)
}

# Numbers
number_cts_cutoff <- function(){
    numericInput(inputId = "cts_cutoff", label = NULL, value = 10)
}

number_plot_size <- function(type) {
    value <- ifelse(type == "height", 600, 800)
    numericInput(inputId = paste0("plot_", type), 
                 label = paste("Plot", str_to_title(type), "(px)"), 
                 value = value)
}

number_n_cluster <- function() {
    numericInput(inputId = "n_cluster",
                 label = "No. of K-means Clusters",
                 value = 3)
}


# ----------
# UI components
# ----------

upload_widgets <- splitLayout(radio_source(),
                              button_cts_upload(),
                              cellWidths = c("67%", "33%"))

palette_widgets <- splitLayout(select_palette(type = "categorical"),
                               select_palette(type = "continuous"))

cluster_widgets <- splitLayout(number_n_cluster(),
                               button_assign_cluster(),
                               cellWidths = c("33%", "67%"))

plot_size_widgets <- splitLayout(number_plot_size("height"),
                                 number_plot_size("width"),
                                 cellWidths = c("50%", "50%"))





