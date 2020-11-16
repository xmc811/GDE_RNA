
# This file contains the UI components

upload_widgets <- splitLayout(radioGroupButtons(inputId = "count_source",
                                            label = NULL,
                                            choices = c("Example","Upload","Select"),
                                            justified = TRUE),
                          actionButton(
                              inputId = "count_start",
                              label = "Upload",
                              icon = icon("upload"),
                              style = "color: white; background-color: #0570b0;
                              float:right; margin-right: 5px;"),
                          
                          cellWidths = c("67%", "33%"))


palette_widgets <- splitLayout(selectInput(inputId = "pal_cat", 
                                           label = "Categorical Palette",
                                           choices = rownames(brewer.pal.info[brewer.pal.info$category == "qual",]),
                                           selected = "Set2"),
                               selectInput(inputId = "pal_con", 
                                           label = "Continuous Palette",
                                           choices = rownames(brewer.pal.info[brewer.pal.info$category != "qual",]),
                                           selected = "Spectral"))


direction_widgets <- materialSwitch(inputId = "pal_dir",
                                    label = "Reverse Scale Color Direction",
                                    value = FALSE,
                                    right = TRUE)

