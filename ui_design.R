
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

# UI components

upload_widgets <- splitLayout(radio01(id = "count_source",
                                      text = NULL,
                                      choices = c("Example","Upload","Select")),
                              button01(id = "count_start",
                                       text = "Upload",
                                       icon = "upload",
                                       style = "color: white; background-color: #0570b0; float:right; margin-right: 5px;"),
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





