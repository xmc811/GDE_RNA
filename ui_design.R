
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
