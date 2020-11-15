
test_panel <- splitLayout(radioGroupButtons(inputId = "count_source",
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
