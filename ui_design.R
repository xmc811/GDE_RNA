
# Basic shiny modules

radio01 <- function(id, text, choices) {
    radioGroupButtons(inputId = id, label = text, choices = choices, justified = TRUE)
}

button01 <- function(id, text, icon, style) {
    actionButton(inputId = id, label = text, icon = icon(icon), style = style)
}

select01 <- function(id, text, choices, selected = NULL) {
    if (is.null(selected)) selected <- choices[1]
    selectInput(inputId = id, label = text, choices = choices, selected = selected)
}

switch01 <- function(id, text) {
    materialSwitch(inputId = id, label = text, value = FALSE, right = TRUE)
}

number01 <- function(id, text = NULL, value = NULL){
    numericInput(inputId = id, label = text, value = value)
}




number_cts_cutoff <- function(){
    numericInput(inputId = "cts_cutoff", label = NULL, value = 10)
}

button_cts_process <- function() {
    actionButton(inputId = "cts_process_click", 
                 label = "Pre-Process", 
                 icon = icon("files-o"), 
                 style = "color: white; background-color: #2ca25f")
}

radio_source <- function() {
    radioGroupButtons(inputId = "cts_source",
                      label  = NULL,
                      choices = c("Example","Upload","Select"),
                      justified = TRUE)
}

button_cts_upload <- function() {
    actionButton(inputId = "cts_upload_click", 
                 label = "Upload", 
                 icon = icon("upload"), 
                 style = "color: white; background-color: #0570b0; float:right; margin-right: 5px;")
}

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

button_cts_compute <- function() {
    actionButton(inputId = "cts_compute_click", 
                 label = "Compute", 
                 icon = icon("calculator"), 
                 style = "color: white; background-color: #0570b0;")
}


# UI components

upload_widgets <- splitLayout(radio_source(),
                              button_cts_upload(),
                              cellWidths = c("67%", "33%"))


palette_widgets <- splitLayout(select01(id = "pal_cat", 
                                        text = "Categorical Palette",
                                        choices = rownames(brewer.pal.info[brewer.pal.info$category == "qual",]),
                                        selected = "Set2"),
                               select01(id = "pal_con", 
                                        text = "Continuous Palette",
                                        choices = rownames(brewer.pal.info[brewer.pal.info$category != "qual",]),
                                        selected = "Spectral"))

cluster_widgets <- splitLayout(number01(id = "n_cluster",
                                        text = "No. of K-means Clusters",
                                        value = 3),
                               button01(id = "assign_clu",
                                        text = "Assign K-means Cluster",
                                        icon = "bar-chart",
                                        style = "color: white; background-color: #2ca25f;margin-top: 25px; float:right; margin-right: 5px;"),
                               cellWidths = c("33%", "67%"))

plot_size_widgets <- splitLayout(number01(id = "plot_height", 
                                          text = "Plot Height (px)", 
                                          value = 600),
                                 number01(id = "plot_width", 
                                          text = "Plot Width (px)", 
                                          value = 800),
                                 cellWidths = c("50%", "50%"))





