source("../app_core/source_app_core.R")

# Load packages ----
library(tibble)
library(dplyr)
library(shiny)
library(bslib)
library(DT)
library(patchwork)
library(gridExtra)
library(grid)
library(gt)
library(readr)
library(writexl)
library(ggplot2)

# Source helpers ----
source("helpers.R")

# ============================================================
# Plot typography settings for manuscript screenshots
# ============================================================

PLOT_TITLE_SIZE <- 26
AXIS_TITLE_SIZE <- 28
AXIS_TEXT_SIZE <- 24
LEGEND_TEXT_SIZE <- 24
PLOT_TAG_SIZE <- 32
RMST_ANNOTATION_SIZE <- 10
SUBGROUP_ANNOTATION_SIZE <- 11
RMST_TABLE_SIZE <- 28

large_plot_theme <- theme(
  plot.title = element_text(hjust = 0.5, size = PLOT_TITLE_SIZE, face = "bold"),
  legend.position = "bottom",
  legend.title = element_blank(),
  legend.text = element_text(size = LEGEND_TEXT_SIZE),
  legend.key.size = unit(1.2, "cm"),
  legend.spacing.x = unit(0.4, "cm"),
  axis.title = element_text(size = AXIS_TITLE_SIZE, face = "bold"),
  axis.text = element_text(size = AXIS_TEXT_SIZE)
)

large_rmst_theme <- theme(
  axis.title = element_text(size = AXIS_TITLE_SIZE, face = "bold"),
  axis.text = element_text(size = AXIS_TEXT_SIZE)
)

# user interface ----
ui <- navset_bar(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
  ),
  title = "Bayesian-ECD demo",
  navbar_options = navbar_options(
    bg = "#0062cc",
    underline = TRUE
  ),
  
  # Panel with uploaded data
  nav_panel(
    title = "Uploaded Data",
    layout_sidebar(
      sidebar = sidebar(
        width = 320,
        helpText("Read the variable descriptions carefully. 
                 Don't switch tabs until you fill out this page."),
        
        helpText("Name of the column containing the patient ID"),
        helpText(tags$ul(
          tags$li("Can just be numbers/letters")
        )),
        textInput("ID", "Patient ID: ", value = NULL),
        
        helpText("Name of the column containing overall survival in days"),
        helpText(tags$ul(
          tags$li("Calculated from Date of Randomization to Date of Death/Last Follow-up or 
                  calculated from Date of First Radiation Therapy to Date of Death/Last Follow-up")
        )),
        textInput("response", "Overall survival: ", value = NULL),
        
        helpText("Name of the column containing the censoring indicator"),
        helpText(tags$ul(
          tags$li("Use TRUE if the event/death was observed and FALSE if the patient was censored/alive at last follow-up.")
        )),
        textInput("censor_ind", "Censoring Indicator: ", value = NULL),
        
        helpText("Name of the column containing the patient treatment type"),
        helpText(tags$ul(
          tags$li("Make sure that your reference treatment is labeled as Control")
        )),
        textInput("trt_type", "Treatment type: ", value = NULL),
        
        helpText("Name of the column containing the cohort name"),
        helpText(tags$ul(
          tags$li("We only support one trial dataset at a time, so make sure this column only has one value")
        )),
        textInput("dat_type", "Cohort name: ", value = NULL),
        
        helpText("Name of the column containing baseline age in years"),
        helpText(tags$ul(
          tags$li("Preferably between 21.00 and 81.22 years")
        )),
        textInput("cont_vars", "Age: ", value = NULL),
        
        helpText("Name of the column containing sex"),
        helpText(tags$ul(
          tags$li("Expected factor levels after loading: Female, Male"),
          tags$li("Original coding assumed: 0 = Female, 1 = Male")
        )),
        textInput("sex", "Sex: ", value = NULL),
        
        helpText("Name of the column containing baseline KPS category"),
        helpText(tags$ul(
          tags$li("Expected factor levels after loading: > 80, (60, 80], <= 60"),
          tags$li("Original coding assumed: 0 = > 80, 1 = (60, 80], 2 = <= 60")
        )),
        textInput("kps", "KPS: ", value = NULL),
        
        helpText("Name of the column containing extent of resection"),
        helpText(tags$ul(
          tags$li("Expected factor levels after loading: GTR, STR, biopsy"),
          tags$li("Original coding assumed: 0 = GTR, 1 = STR, 2 = biopsy")
        )),
        textInput("eor", "EOR: ", value = NULL),
        
        helpText(
          "Please upload your data (.csv format), 
          make sure you pay attention to how things should be formatted, 
          and double check the relevant column names. To continue, click Submit."
        ),
        fluidPage(
          fileInput(
            "file1",
            "Choose CSV File",
            accept = c(
              "text/csv",
              "text/comma-separated-values",
              ".csv"
            )
          )
        ),
        submitButton("Submit")
      ),
      page_fillable(
        div(
          style = "width: 100%; overflow-x: auto;",
          DTOutput("fileinput")
        )
      )
    )
  ),
  
  # Optional concurrent-control compatibility diagnostic
  nav_panel(
    title = "Control Compatibility",
    layout_sidebar(
      helpText(
        "Optional concurrent-control compatibility diagnostic. This check is available when the uploaded dataset contains a concurrent control arm. If your uploaded dataset does not include a concurrent control arm, you can still continue with the Bayesian-ECD analysis using the other tabs."
      ),
      sidebar = sidebar(
        width = 320,
        helpText(
          "The concurrent control arm is assumed to be labeled exactly as: Control"
        ),
        numericInput(
          "compat_min_time",
          "Minimum time for plot:",
          value = 150,
          min = 1,
          max = 2000
        ),
        numericInput(
          "compat_max_time",
          "Maximum time for plot:",
          value = 1200,
          min = 10,
          max = 2500
        ),
        helpText("Set a specific width and/or height when downloading the plot. If not specified, the download button will use the default plot size."),
        numericInput(
          "compat_width",
          "Width (inches):",
          value = NA,
          min = 1,
          max = 50
        ),
        numericInput(
          "compat_height",
          "Height (inches):",
          value = NA,
          min = 1,
          max = 50
        ),
        submitButton("Run compatibility check")
      ),
      div(
        style = "width: 100%;",
        h4("Optional concurrent-control compatibility diagnostic"),
        uiOutput("control_compatibility_status"),
        br(),
        plotOutput("control_compatibility_plot", height = "760px"),
        br(),
        downloadButton(
          "downloadControlCompatibilityPlot",
          "Download Compatibility Plot (.png)",
          style = "width: 100%;"
        )
      )
    )
  ),
  
  # Panel with plots
  nav_panel(
    title = "Plot Output",
    layout_sidebar(
      helpText("Make your selections and click the submit button."),
      sidebar = sidebar(
        width = 320,
        selectInput(
          "ref_treatment_selection",
          label = "Select a reference treatment:",
          choices = NULL
        ),
        selectInput(
          "treatment_selection",
          label = "Select a comparison treatment:",
          choices = NULL
        ),
        selectInput(
          "plot_type",
          label = "Select a plot type",
          choices = c(
            "Survival probabilities vs. Time (Days)",
            "Hazard Ratio vs. Time (Days)"
          ),
          selected = "Survival probabilities vs. Time (Days)"
        ),
        numericInput("min_time", "Minimum time (Days):", 150, min = 10, max = 2000),
        numericInput("max_time", "Maximum time (Days):", 1200, min = 10, max = 2000),
        helpText("Set a specific width and/or height when downloading the plot.
                 If not specified, the download button will download the displayed plot."),
        numericInput("width", "Width (inches):", value = NA, min = 1, max = 50),
        numericInput("height", "Height (inches):", value = NA, min = 1, max = 50),
        submitButton("Submit")
      ),
      page_fillable(plotOutput("plot", height = "720px")),
      downloadButton("downloadresults", "Download Bayesian-ECD Results (.zip)", style = "width: 100%;"),
      downloadButton("downloadplot", "Download Plot (.png)", style = "width: 100%;")
    )
  ),
  
  # Panel with RMST Plots
  nav_panel(
    title = "RMST Output",
    layout_sidebar(
      helpText("Make your selections and click the submit button."),
      sidebar = sidebar(
        width = 320,
        selectInput(
          "ref_treatment_selection1",
          label = "Select a reference treatment:",
          choices = NULL
        ),
        selectInput(
          "treatment_selection1",
          label = "Select Treatment:",
          choices = NULL
        ),
        numericInput(
          "horizon_time",
          "Horizon time in days for RMST:",
          value = 730.5,
          min = 10,
          max = 2000
        ),
        helpText("Set a specific width and/or height when downloading the plot.
                 If not specified, the download button will download the displayed plot."),
        numericInput("width1", "Width (inches):", value = NA, min = 1, max = 50),
        numericInput("height1", "Height (inches):", value = NA, min = 1, max = 50),
        submitButton("Submit")
      ),
      page_fillable(plotOutput("rmst_plot", height = "760px")),
      downloadButton("downloadrmstplot", "Download RMST Plot (.png)", style = "width: 100%;")
    )
  ),
  
  # Panel with subgroup analysis output
  nav_panel(
    title = "Subgroup Analysis",
    layout_sidebar(
      helpText("Make your selections, and click the submit button."),
      sidebar = sidebar(
        width = 320,
        selectInput(
          "ref_treatment_selection2",
          label = "Select a reference treatment:",
          choices = NULL
        ),
        selectInput(
          "treatment_selection2",
          label = "Select Treatment:",
          choices = NULL
        ),
        selectInput(
          "plot_type1",
          label = "Select a plot type",
          choices = c(
            "Survival probabilities vs. Time (Days)",
            "Hazard Ratio vs. Time (Days)"
          ),
          selected = "Survival probabilities vs. Time (Days)"
        ),
        numericInput("min_time1", "Minimum time (Days):", 150, min = 10, max = 2000),
        numericInput("max_time1", "Maximum time (Days):", 1200, min = 10, max = 2000),
        
        selectInput(
          "sex_selection",
          label = "Select the patient sex you want to examine:",
          choices = c("Female", "Male"),
          selected = "Female",
          multiple = TRUE
        ),
        selectInput(
          "eor_selection",
          label = "Select which EOR categories you want to examine:",
          choices = c("GTR", "STR", "biopsy"),
          selected = NULL,
          multiple = TRUE
        ),
        selectInput(
          "kps_selection",
          label = "Select which KPS categories you want to examine:",
          choices = c("> 80", "(60, 80]", "<= 60"),
          selected = NULL,
          multiple = TRUE
        ),
        
        helpText("For age, we recommend not going below 21 and above 81.22 years"),
        helpText(tags$ul(
          tags$li("The minimum age will be included"),
          tags$li("The maximum age will not be included")
        )),
        numericInput("min_age", "Minimum age (Years):", value = 48, min = 0, max = 150),
        numericInput("max_age", "Maximum age (Years):", value = 81, min = 0, max = 150),
        helpText("For example, the selection above will include patients that are 48 years old, 
                 but exclude patients that are 81 years old."),
        
        helpText("Set a specific width and/or height when downloading the plot.
                 If not specified, the download button will download the displayed plot."),
        numericInput("width2", "Width (inches):", value = NA, min = 1, max = 50),
        numericInput("height2", "Height (inches):", value = NA, min = 1, max = 50),
        submitButton("Submit")
      ),
      page_fillable(plotOutput("subgroup_plot", height = "720px")),
      downloadButton("downloadsubgroupplot", "Download Plot (.png)", style = "width: 100%;")
    )
  ),
  
  # Panel with subgroup RMST output
  nav_panel(
    title = "Subgroup RMST Output",
    layout_sidebar(
      helpText("Make your selections in the Subgroup Analysis tab, then click Submit here."),
      sidebar = sidebar(
        width = 320,
        numericInput(
          "horizon_time1",
          "Horizon time in days for RMST:",
          value = 730.5,
          min = 10,
          max = 2000
        ),
        submitButton("Submit")
      ),
      page_fillable(DTOutput("subgroup_rmst_table")),
      downloadButton("downloadRMSTtable", "Download RMST Table (.csv)", style = "width: 100%;")
    )
  ),
  
  # Panel with prediction data
  nav_panel(
    title = "Prediction Data",
    layout_sidebar(
      sidebar = sidebar(
        width = 320,
        helpText("Upload a new prediction dataset containing the same covariates used in model fitting."),
        helpText("The prediction dataset only needs patient ID and covariates. It does not need survival outcome, censoring indicator, treatment, or cohort columns."),
        
        helpText("Name of the column containing the prediction patient ID"),
        textInput("pred_ID", "Prediction Patient ID: ", value = NULL),
        
        helpText("Please upload your prediction data (.csv format)."),
        fileInput(
          "pred_file",
          "Choose Prediction CSV File",
          accept = c(
            "text/csv",
            "text/comma-separated-values",
            ".csv"
          )
        ),
        submitButton("Submit")
      ),
      page_fillable(
        div(
          style = "width: 100%; overflow-x: auto;",
          DTOutput("pred_fileinput")
        )
      )
    )
  ),
  
  # Panel with prediction output
  nav_panel(
    title = "Prediction Output",
    layout_sidebar(
      helpText("Prediction output is generated from the uploaded prediction dataset. By default, the full prediction dataset is used unless subgroup filters are enabled."),
      sidebar = sidebar(
        width = 320,
        
        selectInput(
          "pred_ref_treatment_selection",
          label = "Select a prediction reference treatment:",
          choices = NULL
        ),
        selectInput(
          "pred_treatment_selection",
          label = "Select a prediction comparison treatment:",
          choices = NULL
        ),
        selectInput(
          "pred_plot_type",
          label = "Select a plot type",
          choices = c(
            "Survival probabilities vs. Time (Days)",
            "Hazard Ratio vs. Time (Days)"
          ),
          selected = "Survival probabilities vs. Time (Days)"
        ),
        numericInput("pred_min_time", "Minimum time (Days):", 150, min = 10, max = 2000),
        numericInput("pred_max_time", "Maximum time (Days):", 1200, min = 10, max = 2000),
        
        checkboxInput(
          "pred_use_subgroup",
          "Apply subgroup filters to prediction dataset",
          value = FALSE
        ),
        
        conditionalPanel(
          condition = "input.pred_use_subgroup == true",
          
          selectInput(
            "pred_sex_selection",
            label = "Select the patient sex you want to examine:",
            choices = c("Female", "Male"),
            selected = NULL,
            multiple = TRUE
          ),
          selectInput(
            "pred_eor_selection",
            label = "Select which EOR categories you want to examine:",
            choices = c("GTR", "STR", "biopsy"),
            selected = NULL,
            multiple = TRUE
          ),
          selectInput(
            "pred_kps_selection",
            label = "Select which KPS categories you want to examine:",
            choices = c("> 80", "(60, 80]", "<= 60"),
            selected = NULL,
            multiple = TRUE
          ),
          
          helpText("For age, the minimum age is included and the maximum age is excluded."),
          numericInput("pred_min_age", "Minimum age (Years):", value = NA, min = 0, max = 150),
          numericInput("pred_max_age", "Maximum age (Years):", value = NA, min = 0, max = 150)
        ),
        
        helpText("Set a specific width and/or height when downloading the plot.
                 If not specified, the download button will download the displayed plot."),
        numericInput("pred_width", "Width (inches):", value = NA, min = 1, max = 50),
        numericInput("pred_height", "Height (inches):", value = NA, min = 1, max = 50),
        submitButton("Submit")
      ),
      page_fillable(plotOutput("pred_plot", height = "720px")),
      downloadButton("downloadpredplot", "Download Prediction Plot (.png)", style = "width: 100%;")
    )
  ),
  
  # Panel with prediction RMST output
  nav_panel(
    title = "Prediction RMST Output",
    layout_sidebar(
      helpText("This RMST table uses the same prediction cohort/subgroup selected in the Prediction Output tab."),
      sidebar = sidebar(
        width = 320,
        numericInput(
          "pred_horizon_time",
          "Horizon time in days for RMST:",
          value = 730.5,
          min = 10,
          max = 2000
        ),
        submitButton("Submit")
      ),
      page_fillable(DTOutput("pred_rmst_table")),
      downloadButton("downloadPredRMSTtable", "Download Prediction RMST Table (.csv)", style = "width: 100%;")
    )
  )
)

# Server logic ----
server <- function(input, output, session) {
  
  session$onSessionEnded(function() {
    try(foreach::registerDoSEQ(), silent = TRUE)
    gc()
    stopApp()
  })
  
  make_rmst_combined_plot <- function(
      RMST_store,
      cntrl,
      trt,
      horizon_time
  ) {
    tmp1 <- RMST_store$RMST
    tmp2 <- RMST_store$Difference
    tmp2$`Difference (Days)` <- round(tmp2$`Difference (Days)`, 2)
    tmp3 <- RMST_store$Ratio
    tmp3$Ratio <- round(tmp3$Ratio, 2)
    
    p.trt <- RMST_store$p.surv.trt +
      annotate(
        "text",
        x = horizon_time / 2,
        y = 0.5,
        label = paste("RMST:", round(tmp1[1, 2]), "Days"),
        size = RMST_ANNOTATION_SIZE
      ) +
      labs(x = "Time") +
      theme_bw() +
      large_rmst_theme
    
    p.cntrl <- RMST_store$p.surv.cntrl +
      annotate(
        "text",
        x = horizon_time / 2,
        y = 0.5,
        label = paste("RMST:", round(tmp1[2, 2]), "Days"),
        size = RMST_ANNOTATION_SIZE
      ) +
      labs(x = "Time") +
      theme_bw() +
      large_rmst_theme
    
    ((p.trt + p.cntrl + plot_layout(axis_titles = "collect") &
        plot_annotation(tag_levels = list(c(trt, cntrl), "1")) &
        theme(plot.tag = element_text(size = PLOT_TAG_SIZE, face = "bold"))) +
        tableGrob(
          t(tmp2),
          cols = NULL,
          theme = ttheme_default(base_size = RMST_TABLE_SIZE)
        ) +
        tableGrob(
          t(tmp3),
          cols = NULL,
          theme = ttheme_default(base_size = RMST_TABLE_SIZE)
        )) +
      plot_layout(heights = c(7, 1))
  }
  
  input.specs <- reactive({
    req(input$response)
    req(input$censor_ind)
    req(input$trt_type)
    req(input$dat_type)
    req(input$cont_vars)
    req(input$sex)
    req(input$kps)
    req(input$eor)
    
    list(
      response = input$response,
      censor_ind = input$censor_ind,
      trt_type = input$trt_type,
      dat_type = input$dat_type,
      cont_vars = input$cont_vars,
      cat_vars = c(input$sex, input$kps, input$eor)
    )
  })
  
  # ------------------------------------------------------------
  # Subgroup list for original uploaded modeling data
  # ------------------------------------------------------------
  subgroup_list <- reactive({
    req(input$sex)
    req(input$eor)
    req(input$kps)
    req(input$cont_vars)
    
    subgrp_list <- list()
    subgrp_list[[input$sex]] <- NA
    subgrp_list[[input$eor]] <- NA
    subgrp_list[[input$kps]] <- NA
    subgrp_list[[input$cont_vars]] <- NA
    
    sex_map <- c("Female" = "0", "Male" = "1")
    if (!is.null(input$sex_selection) && length(input$sex_selection) > 0) {
      selected <- input$sex_selection[input$sex_selection %in% names(sex_map)]
      if (length(selected) > 0 && length(selected) < length(sex_map)) {
        subgrp_list[[input$sex]] <- paste0(sex_map[selected], collapse = "")
      }
    }
    
    eor_map <- c("GTR" = "0", "STR" = "1", "biopsy" = "2")
    if (!is.null(input$eor_selection) && length(input$eor_selection) > 0) {
      selected <- input$eor_selection[input$eor_selection %in% names(eor_map)]
      if (length(selected) > 0 && length(selected) < length(eor_map)) {
        subgrp_list[[input$eor]] <- paste0(eor_map[selected], collapse = "")
      }
    }
    
    kps_map <- c("> 80" = "0", "(60, 80]" = "1", "<= 60" = "2")
    if (!is.null(input$kps_selection) && length(input$kps_selection) > 0) {
      selected <- input$kps_selection[input$kps_selection %in% names(kps_map)]
      if (length(selected) > 0 && length(selected) < length(kps_map)) {
        subgrp_list[[input$kps]] <- paste0(kps_map[selected], collapse = "")
      }
    }
    
    if (!is.na(input$min_age) && is.na(input$max_age)) {
      subgrp_list[[input$cont_vars]] <- as.numeric(c(
        input$min_age,
        max(data()[[input$cont_vars]], na.rm = TRUE)
      ))
    } else if (is.na(input$min_age) && !is.na(input$max_age)) {
      subgrp_list[[input$cont_vars]] <- as.numeric(c(
        min(data()[[input$cont_vars]], na.rm = TRUE),
        input$max_age
      ))
    } else if (!is.na(input$min_age) && !is.na(input$max_age)) {
      subgrp_list[[input$cont_vars]] <- as.numeric(c(input$min_age, input$max_age))
    }
    
    if (all(is.na(subgrp_list[[input$cont_vars]]))) {
      subgrp_list[[input$cont_vars]] <- NA
    }
    
    subgrp_list
  })
  
  # ------------------------------------------------------------
  # Uploaded model-fitting data
  # ------------------------------------------------------------
  data <- reactive({
    req(input$file1)
    req(input.specs())
    
    df <- read_csv(
      input$file1$datapath,
      col_select = as.character(c(input$ID, unlist(input.specs()))),
      show_col_types = FALSE
    )
    
    # Ensure categorical variables are stored as factors for model/subgroup functions.
    cat_vars <- input.specs()$cat_vars
    
    if (!is.null(cat_vars) && length(cat_vars) > 0) {
      df[cat_vars] <- lapply(df[cat_vars], as.factor)
    }
    
    # Also treat ID, treatment, and cohort/source columns as factors.
    df[[input$ID]] <- as.factor(df[[input$ID]])
    df[[input.specs()$trt_type]] <- as.factor(df[[input.specs()$trt_type]])
    df[[input.specs()$dat_type]] <- as.factor(df[[input.specs()$dat_type]])
    
    df
  })
  
  data.mod <- reactive({
    cbind(
      lapply(data()[, input$ID], as.factor),
      data()[, input.specs()$cont_vars],
      lapply(data()[, input.specs()$cat_vars], as.factor),
      lapply(data()[, input.specs()$trt_type], as.factor),
      lapply(data()[, input.specs()$dat_type], as.factor),
      data()[, input.specs()$response],
      data()[, input.specs()$censor_ind]
    )
  })
  
  output$fileinput <- renderDT(
    { data.mod() },
    filter = "top",
    options = list(
      pageLength = 5,
      autoWidth = TRUE,
      scrollX = TRUE,
      dom = "ltp"
    ),
    rownames = FALSE
  )
  
  # ------------------------------------------------------------
  # Prediction data: ID + covariates only
  # ------------------------------------------------------------
  pred_data <- reactive({
    req(input$pred_file)
    req(input$pred_ID)
    req(input.specs())
    
    required_cols <- as.character(c(
      input$pred_ID,
      input.specs()$cont_vars,
      input.specs()$cat_vars
    ))
    
    pred_raw <- read_csv(
      input$pred_file$datapath,
      show_col_types = FALSE
    )
    
    missing_cols <- setdiff(required_cols, names(pred_raw))
    
    if (length(missing_cols) > 0) {
      stop(
        "Prediction dataset is missing required column(s): ",
        paste(missing_cols, collapse = ", ")
      )
    }
    
    df <- pred_raw[, required_cols, drop = FALSE]
    
    cat_vars <- input.specs()$cat_vars
    
    if (!is.null(cat_vars) && length(cat_vars) > 0) {
      df[cat_vars] <- lapply(df[cat_vars], as.factor)
    }
    
    df[[input$pred_ID]] <- as.factor(df[[input$pred_ID]])
    
    df
  })
  
  pred_data.mod <- reactive({
    req(pred_data())
    
    cbind(
      lapply(pred_data()[, input$pred_ID], as.factor),
      pred_data()[, input.specs()$cont_vars],
      lapply(pred_data()[, input.specs()$cat_vars], as.factor)
    )
  })
  
  output$pred_fileinput <- renderDT(
    { pred_data.mod() },
    filter = "top",
    options = list(
      pageLength = 5,
      autoWidth = TRUE,
      scrollX = TRUE,
      dom = "ltp"
    ),
    rownames = FALSE
  )
  
  pred_subgroup_list <- reactive({
    req(input.specs())
    req(pred_data())
    
    subgrp_list <- list()
    subgrp_list[[input$sex]] <- NA
    subgrp_list[[input$eor]] <- NA
    subgrp_list[[input$kps]] <- NA
    subgrp_list[[input$cont_vars]] <- NA
    
    if (!isTRUE(input$pred_use_subgroup)) {
      return(subgrp_list)
    }
    
    sex_map <- c("Female" = "0", "Male" = "1")
    if (!is.null(input$pred_sex_selection) && length(input$pred_sex_selection) > 0) {
      selected <- input$pred_sex_selection[input$pred_sex_selection %in% names(sex_map)]
      if (length(selected) > 0 && length(selected) < length(sex_map)) {
        subgrp_list[[input$sex]] <- paste0(sex_map[selected], collapse = "")
      }
    }
    
    eor_map <- c("GTR" = "0", "STR" = "1", "biopsy" = "2")
    if (!is.null(input$pred_eor_selection) && length(input$pred_eor_selection) > 0) {
      selected <- input$pred_eor_selection[input$pred_eor_selection %in% names(eor_map)]
      if (length(selected) > 0 && length(selected) < length(eor_map)) {
        subgrp_list[[input$eor]] <- paste0(eor_map[selected], collapse = "")
      }
    }
    
    kps_map <- c("> 80" = "0", "(60, 80]" = "1", "<= 60" = "2")
    if (!is.null(input$pred_kps_selection) && length(input$pred_kps_selection) > 0) {
      selected <- input$pred_kps_selection[input$pred_kps_selection %in% names(kps_map)]
      if (length(selected) > 0 && length(selected) < length(kps_map)) {
        subgrp_list[[input$kps]] <- paste0(kps_map[selected], collapse = "")
      }
    }
    
    if (!is.na(input$pred_min_age) && is.na(input$pred_max_age)) {
      subgrp_list[[input$cont_vars]] <- as.numeric(c(
        input$pred_min_age,
        max(pred_data()[[input$cont_vars]], na.rm = TRUE)
      ))
    } else if (is.na(input$pred_min_age) && !is.na(input$pred_max_age)) {
      subgrp_list[[input$cont_vars]] <- as.numeric(c(
        min(pred_data()[[input$cont_vars]], na.rm = TRUE),
        input$pred_max_age
      ))
    } else if (!is.na(input$pred_min_age) && !is.na(input$pred_max_age)) {
      subgrp_list[[input$cont_vars]] <- as.numeric(c(
        input$pred_min_age,
        input$pred_max_age
      ))
    }
    
    if (all(is.na(subgrp_list[[input$cont_vars]]))) {
      subgrp_list[[input$cont_vars]] <- NA
    }
    
    subgrp_list
  })
  
  observeEvent(input.specs(), {
    choices <- names(table(data()[, input.specs()$trt_type]))
    choices <- unique(c(choices, result.CAPPMx.base$Treatment_Levels))
    
    updateSelectInput(session, "ref_treatment_selection",
                      choices = choices,
                      selected = result.CAPPMx.base$Treatment_Levels[1])
    
    updateSelectInput(session, "ref_treatment_selection1",
                      choices = choices,
                      selected = result.CAPPMx.base$Treatment_Levels[1])
    
    updateSelectInput(session, "ref_treatment_selection2",
                      choices = choices,
                      selected = result.CAPPMx.base$Treatment_Levels[1])
    
    updateSelectInput(session, "treatment_selection", choices = choices)
    updateSelectInput(session, "treatment_selection1", choices = choices)
    updateSelectInput(session, "treatment_selection2", choices = choices)
    
    # Prediction treatment selections are intentionally separate from the main tabs.
    updateSelectInput(session, "pred_ref_treatment_selection",
                      choices = choices,
                      selected = result.CAPPMx.base$Treatment_Levels[1])
    updateSelectInput(session, "pred_treatment_selection", choices = choices)
  })
  
  result <- reactive({
    req(input$ref_treatment_selection)
    
    withProgress(
      message = "Building the Bayesian-ECD model",
      detail = "This may take a moment.",
      value = 0, {
        cappmx_extend_approx_fit(
          result_CAPPMx = result.CAPPMx.base,
          input_df = data(),
          input_specs = input.specs(),
          ref_trt = input$ref_treatment_selection,
          del_range_response_1 = c(0.005, 0.02) * 8,
          del_range_response_2 = c(0.005, 0.02) * 9,
          del_range_alp1 = c(0.1, 0.3) * 2.8
        )
      }
    )
  })
  
  pred_result <- reactive({
    req(input$pred_ref_treatment_selection)
    req(data())
    req(pred_data())
    req(input.specs())
    
    withProgress(
      message = "Building the Bayesian-ECD prediction model",
      detail = "This may take a moment.",
      value = 0, {
        cappmx_extend_approx_fit(
          result_CAPPMx = result.CAPPMx.base,
          input_df = data(),
          input_specs = input.specs(),
          ref_trt = input$pred_ref_treatment_selection,
          input_df_pred = pred_data(),
          del_range_response_1 = c(0.005, 0.02) * 8,
          del_range_response_2 = c(0.005, 0.02) * 9,
          del_range_alp1 = c(0.1, 0.3) * 2.8
        )
      }
    )
  })
  
  picube <- reactive({
    res <- result()
    if (!is.list(res) || is.null(res$picube)) return(NULL)
    
    num.cohort <- dim(res$picube)[3]
    
    tmp <- lapply(seq_len(num.cohort), function(i) {
      data.frame(res$picube[, , i])
    })
    
    names(tmp) <- res$Data_Levels
    tmp
  })
  
  Lognormal_Mu_Cube <- reactive({
    res <- result()
    if (!is.list(res) || is.null(res$Lognormal_Mu_Cube)) return(NULL)
    
    num.trt <- dim(res$Lognormal_Mu_Cube)[3]
    
    tmp <- lapply(seq_len(num.trt), function(i) {
      data.frame(res$Lognormal_Mu_Cube[, , i])
    })
    
    names(tmp) <- res$Treatment_Levels
    tmp
  })
  
  Lognormal_Sig_Cube <- reactive({
    res <- result()
    if (!is.list(res) || is.null(res$Lognormal_Sig_Cube)) return(NULL)
    
    num.trt <- dim(res$Lognormal_Sig_Cube)[3]
    
    tmp <- lapply(seq_len(num.trt), function(i) {
      data.frame(res$Lognormal_Sig_Cube[, , i])
    })
    
    names(tmp) <- res$Treatment_Levels
    tmp
  })
  
  Allocation_variables <- reactive({
    res <- result()
    if (!is.list(res) || is.null(res$Allocation_variables)) return(NULL)
    data.frame(res$Allocation_variables)
  })
  
  Variable_specifications <- reactive({
    res <- result()
    if (!is.list(res) || is.null(res$Variable_specifications)) return(NULL)
    res$Variable_specifications
  })
  
  output$downloadresults <- downloadHandler(
    filename = function() {
      paste0("results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip")
    },
    content = function(fname) {
      temp_dir <- tempdir()
      
      if (grepl("\\\\", temp_dir)) {
        temp_dir <- gsub("\\\\", "/", temp_dir)
      }
      
      files_to_zip <- c(
        file.path(temp_dir, "picube.xlsx"),
        file.path(temp_dir, "Lognormal_Mu_Cube.xlsx"),
        file.path(temp_dir, "Lognormal_Sig_Cube.xlsx"),
        file.path(temp_dir, "Allocation_variables.xlsx"),
        file.path(temp_dir, "Variable_specifications.txt")
      )
      
      write_xlsx(picube(), path = files_to_zip[1])
      write_xlsx(Lognormal_Mu_Cube(), path = files_to_zip[2])
      write_xlsx(Lognormal_Sig_Cube(), path = files_to_zip[3])
      write_xlsx(Allocation_variables(), path = files_to_zip[4])
      capture.output(print(Variable_specifications()), file = files_to_zip[5])
      
      zip::zip(zipfile = fname, files = files_to_zip, mode = "cherry-pick")
    }
  )
  
  plot_store <- reactive({
    req(input$ref_treatment_selection)
    req(input$treatment_selection)
    
    withProgress(
      message = "Building the plotting results.",
      detail = "This may take a moment.",
      value = 0, {
        plots_lognormal(
          result = result(),
          cntrl = input$ref_treatment_selection,
          trt = input$treatment_selection,
          timepoints = seq(input$min_time, input$max_time, 3),
          burnin = 200
        )
      }
    )
  })
  
  plot_store1 <- reactive({
    req(input$ref_treatment_selection2)
    req(input$treatment_selection2)
    
    withProgress(
      message = "Building the subgroup plotting results.",
      detail = "This may take a moment.",
      value = 0, {
        subgroup_data(
          result = result(),
          cntrl = input$ref_treatment_selection2,
          trt = input$treatment_selection2,
          input_df = data(),
          timepoints = seq(input$min_time1, input$max_time1, 3),
          burnin = 200,
          subgroup_list = subgroup_list(),
          use_pred = FALSE
        )
      }
    )
  })
  
  pred_plot_store <- reactive({
    req(input$pred_ref_treatment_selection)
    req(input$pred_treatment_selection)
    req(pred_data())
    
    withProgress(
      message = "Building the prediction plotting results.",
      detail = "This may take a moment.",
      value = 0, {
        subgroup_data(
          result = pred_result(),
          cntrl = input$pred_ref_treatment_selection,
          trt = input$pred_treatment_selection,
          input_df = pred_data(),
          timepoints = seq(input$pred_min_time, input$pred_max_time, 3),
          burnin = 200,
          subgroup_list = pred_subgroup_list(),
          use_pred = TRUE
        )
      }
    )
  })
  
  data_table_store <- reactive({
    req(input$horizon_time1)
    
    withProgress(
      message = "Building the subgroup RMST table.",
      detail = "This may take a moment.",
      value = 0, {
        RMST_result(
          result_subgroup_data = plot_store1(),
          time_horizons = input$horizon_time1
        )
      }
    )
  })
  
  pred_data_table_store <- reactive({
    req(input$pred_horizon_time)
    
    withProgress(
      message = "Building the prediction RMST table.",
      detail = "This may take a moment.",
      value = 0, {
        RMST_result(
          result_subgroup_data = pred_plot_store(),
          time_horizons = input$pred_horizon_time
        )
      }
    )
  })
  
  control_compatibility_has_control <- reactive({
    if (is.null(input$file1)) {
      return(FALSE)
    }
    
    if (
      is.null(input$trt_type) ||
      input$trt_type == ""
    ) {
      return(FALSE)
    }
    
    uploaded_df <- tryCatch(data(), error = function(e) NULL)
    
    if (is.null(uploaded_df)) {
      return(FALSE)
    }
    
    if (!input$trt_type %in% names(uploaded_df)) {
      return(FALSE)
    }
    
    sum(as.character(uploaded_df[[input$trt_type]]) == "Control") > 0
  })
  
  output$control_compatibility_status <- renderUI({
    compatibility_error <- control_compatibility_error_store()
    
    if (!is.null(compatibility_error)) {
      return(tags$div(
        class = "alert alert-danger",
        paste("Compatibility diagnostic failed:", compatibility_error)
      ))
    }
    
    if (is.null(input$file1)) {
      return(tags$div(
        class = "alert alert-info",
        "Please upload a dataset in the Uploaded Data tab before running this check. This diagnostic is optional and is only available when the uploaded dataset includes concurrent-control patients."
      ))
    }
    
    if (is.null(input$trt_type) || input$trt_type == "") {
      return(tags$div(
        class = "alert alert-warning",
        "Please specify the treatment column in the Uploaded Data tab. You can still continue with the Bayesian-ECD analysis using the other tabs."
      ))
    }
    
    uploaded_df <- tryCatch(data(), error = function(e) NULL)
    
    if (is.null(uploaded_df)) {
      return(tags$div(
        class = "alert alert-warning",
        "The uploaded dataset could not be read yet. Please check the Uploaded Data tab. You can still continue with other analyses after the uploaded data are valid."
      ))
    }
    
    if (!input$trt_type %in% names(uploaded_df)) {
      return(tags$div(
        class = "alert alert-warning",
        "The specified treatment column was not found in the uploaded dataset. You can still continue with the Bayesian-ECD analysis once the uploaded data fields are corrected."
      ))
    }
    
    control_n <- sum(
      as.character(uploaded_df[[input$trt_type]]) == "Control"
    )
    
    if (control_n == 0) {
      return(tags$div(
        class = "alert alert-warning",
        paste0(
          "No concurrent-control patients were found for control label '",
          "Control",
          "'. This diagnostic requires a concurrent control arm, but you may still continue with the Bayesian-ECD analysis using the other tabs."
        )
      ))
    }
    
    tags$div(
      class = "alert alert-success",
      paste0(
        "Found ",
        control_n,
        " uploaded concurrent-control patients. You can run the optional posterior predictive compatibility check."
      )
    )
  })
  
  control_compatibility_error_store <- reactiveVal(NULL)
  
  control_compatibility_result <- reactive({
    req(data())
    req(input.specs())
    req(input$compat_min_time)
    req(input$compat_max_time)
    
    control_compatibility_error_store(NULL)
    
    validate(
      need(
        control_compatibility_has_control(),
        "No concurrent-control patients labeled 'Control' were found. This optional diagnostic is skipped, but you can continue with the Bayesian-ECD analysis using the other tabs."
      ),
      need(
        input$compat_min_time < input$compat_max_time,
        "Minimum time must be smaller than maximum time."
      )
    )
    
    uploaded_all_df <- data()
    
    uploaded_control_df <- uploaded_all_df[
      as.character(uploaded_all_df[[input.specs()$trt_type]]) == "Control",
      ,
      drop = FALSE
    ]
    
    
    validate(
      need(
        nrow(uploaded_control_df) > 0,
        "No concurrent-control patients labeled 'Control' were found. This optional diagnostic is skipped, but you can continue with the Bayesian-ECD analysis using the other tabs."
      )
    )
    
    time_grid <- seq(
      from = input$compat_min_time,
      to = input$compat_max_time,
      length.out = 201
    )
    
    showNotification(
      "Starting optional concurrent-control compatibility diagnostic...",
      type = "message",
      duration = 5
    )
    
    tryCatch({
      withProgress(
        message = "Running optional concurrent-control compatibility diagnostic",
        detail = "Step 1 of 3: Building the Control-reference Bayesian-ECD model.",
        value = 0.1, {
          
          compatibility_model <- cappmx_extend_approx_fit(
            result_CAPPMx = result.CAPPMx.base,
            input_df = uploaded_control_df,
            input_specs = input.specs(),
            ref_trt = "Control",
            del_range_response_1 = c(0.005, 0.02) * 8,
            del_range_response_2 = c(0.005, 0.02) * 9,
            del_range_alp1 = c(0.1, 0.3) * 2.8
          )
          
          incProgress(
            amount = 0.35,
            detail = "Step 2 of 3: Computing the uploaded-control KM curve."
          )
          
          incProgress(
            amount = 0.25,
            detail = "Step 3 of 3: Computing historical posterior predictive survival bands and median interval."
          )
          
          out <- run_control_compatibility_check(
            result = compatibility_model,
            uploaded_data = uploaded_control_df,
            time_col = input.specs()$response,
            censor_col = input.specs()$censor_ind,
            trt_col = input.specs()$trt_type,
            control_label = "Control",
            time_grid = time_grid,
            burnin = 200L,
            n_posterior_draws = 1000L,
            conf_int = 0.95,
            seed = 20260706L
          )
          
          incProgress(
            amount = 0.30,
            detail = "Compatibility diagnostic complete."
          )
          
          showNotification(
            "Control compatibility diagnostic is complete.",
            type = "message",
            duration = 5
          )
          
          out
        }
      )
    }, error = function(e) {
      control_compatibility_error_store(conditionMessage(e))
      
      showNotification(
        paste("Compatibility diagnostic failed:", conditionMessage(e)),
        type = "error",
        duration = 10
      )
      
      stop(e)
    })
  })
  
  output$control_compatibility_plot <- renderPlot({
    validate(
      need(
        !is.null(control_compatibility_result()),
        "Click 'Run compatibility check' in the sidebar to generate the compatibility plot."
      )
    )
    
    compatibility_plot_result <- control_compatibility_result()
    
    compatibility_plot_result$plot_data <- compatibility_plot_result$plot_data[
      compatibility_plot_result$plot_data$time >= input$compat_min_time &
        compatibility_plot_result$plot_data$time <= input$compat_max_time,
      ,
      drop = FALSE
    ]
    
    combined_grob <- make_control_compatibility_combined_grob(
      compatibility_result = compatibility_plot_result,
      x_axis_min = 0,
      x_axis_max = input$compat_max_time + 100,
      plot_theme = large_plot_theme,
      table_base_size = 18
    )
    
    grid::grid.draw(combined_grob)
  })
  
  output$downloadControlCompatibilityPlot <- downloadHandler(
    filename = function() {
      paste0(
        "control_compatibility_plot_",
        format(Sys.time(), "%Y%m%d_%H%M%S"),
        ".png"
      )
    },
    content = function(file) {
      req(control_compatibility_result())
      
      compatibility_plot_result <- control_compatibility_result()
      
      compatibility_plot_result$plot_data <- compatibility_plot_result$plot_data[
        compatibility_plot_result$plot_data$time >= input$compat_min_time &
          compatibility_plot_result$plot_data$time <= input$compat_max_time,
        ,
        drop = FALSE
      ]
      
      combined_grob <- make_control_compatibility_combined_grob(
        compatibility_result = compatibility_plot_result,
        x_axis_min = 0,
        x_axis_max = input$compat_max_time + 100,
        plot_theme = large_plot_theme,
        table_base_size = 18
      )
      
      plot_width <- ifelse(
        is.na(input$compat_width),
        14,
        input$compat_width
      )
      
      plot_height <- ifelse(
        is.na(input$compat_height),
        8,
        input$compat_height
      )
      
      ggplot2::ggsave(
        filename = file,
        plot = combined_grob,
        width = plot_width,
        height = plot_height,
        dpi = 300
      )
    }
  )
  
  make_main_survival_or_hr_plot <- function(
      plot_store_result,
      plot_type,
      cntrl,
      trt,
      min_time,
      max_time
  ) {
    plt <- switch(
      plot_type,
      "Survival probabilities vs. Time (Days)" = "Survival",
      "Hazard Ratio vs. Time (Days)" = "HR"
    )
    
    if (plt == "Survival") {
      plot_store_result$p.surv +
        scale_fill_discrete(labels = c(trt, cntrl)) +
        scale_color_discrete(labels = c(trt, cntrl)) +
        labs(x = "Time") +
        theme_bw() +
        large_plot_theme +
        scale_x_continuous(limits = c(min_time - 10, max_time + 100))
    } else {
      plot_store_result$p.hazarad.ratio.withCI +
        labs(x = "Time") +
        theme_bw() +
        large_plot_theme +
        scale_x_continuous(limits = c(min_time - 10, max_time + 100))
    }
  }
  
  
  make_subgroup_survival_or_hr_plot <- function(
      subgroup_plot_result,
      plot_type,
      cntrl,
      trt,
      min_time,
      max_time,
      cohort_label = "Subgroup sample size:"
  ) {
    plt <- switch(
      plot_type,
      "Survival probabilities vs. Time (Days)" = "Survival",
      "Hazard Ratio vs. Time (Days)" = "HR"
    )
    
    n_subgroup <- length(subgroup_plot_result$subgroup_index)
    
    if (plt == "Survival") {
      subgroup_plot_result$p.surv +
        scale_fill_discrete(labels = c(trt, cntrl)) +
        scale_color_discrete(labels = c(trt, cntrl)) +
        annotate(
          "text",
          x = max_time / 1.5,
          y = max(subgroup_plot_result$surv.data$surv.pred.97.5, na.rm = TRUE) / 1.5,
          label = paste(cohort_label, n_subgroup),
          size = SUBGROUP_ANNOTATION_SIZE
        ) +
        labs(x = "Time") +
        theme_bw() +
        large_plot_theme +
        scale_x_continuous(limits = c(min_time - 10, max_time + 100))
    } else {
      subgroup_plot_result$p.hazarad.ratio.withCI +
        annotate(
          "text",
          x = max_time / 1.5,
          y = max(subgroup_plot_result$HR.data$HR.pred.97.5, na.rm = TRUE) / 1.5,
          label = paste(cohort_label, n_subgroup),
          size = SUBGROUP_ANNOTATION_SIZE
        ) +
        labs(x = "Time") +
        theme_bw() +
        large_plot_theme +
        scale_x_continuous(limits = c(min_time - 10, max_time + 100))
    }
  }
  
  
  save_plot_with_optional_size <- function(
      file,
      plot,
      width_input,
      height_input,
      output_width_px,
      output_height_px
  ) {
    if (is.na(width_input) && is.na(height_input)) {
      ggsave(
        file,
        plot = plot,
        device = "png",
        width = as.numeric(output_width_px),
        height = as.numeric(output_height_px),
        units = "px",
        dpi = 96
      )
    } else if (!is.na(width_input) && is.na(height_input)) {
      ggsave(
        file,
        plot = plot,
        device = "png",
        width = width_input,
        units = "in"
      )
    } else if (is.na(width_input) && !is.na(height_input)) {
      ggsave(
        file,
        plot = plot,
        device = "png",
        height = height_input,
        units = "in"
      )
    } else {
      ggsave(
        file,
        plot = plot,
        device = "png",
        width = width_input,
        height = height_input,
        units = "in"
      )
    }
  }
  
  
  output$plot <- renderPlot({
    req(input$ref_treatment_selection)
    req(input$treatment_selection)
    
    withProgress(
      message = "Building the survival/hazard ratio plot",
      detail = "This may take a moment.",
      value = 0, {
        main_plot_result <- plot_store()
        
        make_main_survival_or_hr_plot(
          plot_store_result = main_plot_result,
          plot_type = input$plot_type,
          cntrl = input$ref_treatment_selection,
          trt = input$treatment_selection,
          min_time = input$min_time,
          max_time = input$max_time
        )
      }
    )
  })
  
  rmst_store <- reactive({
    req(input$ref_treatment_selection1)
    req(input$treatment_selection1)
    req(input$horizon_time)
    
    withProgress(
      message = "Building the RMST results",
      detail = "This may take a moment.",
      value = 0, {
        rmst_lognormal(
          result(),
          input$ref_treatment_selection1,
          input$treatment_selection1,
          input$horizon_time,
          burnin = 200
        )
      }
    )
  })
  
  output$rmst_plot <- renderPlot({
    req(input$ref_treatment_selection1)
    req(input$treatment_selection1)
    
    withProgress(
      message = "Building the RMST plot",
      detail = "This may take a moment.",
      value = 0, {
        make_rmst_combined_plot(
          RMST_store = rmst_store(),
          cntrl = input$ref_treatment_selection1,
          trt = input$treatment_selection1,
          horizon_time = input$horizon_time
        )
      }
    )
  })
  
  output$subgroup_plot <- renderPlot({
    req(input$ref_treatment_selection2)
    req(input$treatment_selection2)
    
    withProgress(
      message = "Building the subgroup survival/hazard ratio plot",
      detail = "This may take a moment.",
      value = 0, {
        subgroup_plot_result <- plot_store1()
        
        make_subgroup_survival_or_hr_plot(
          subgroup_plot_result = subgroup_plot_result,
          plot_type = input$plot_type1,
          cntrl = input$ref_treatment_selection2,
          trt = input$treatment_selection2,
          min_time = input$min_time1,
          max_time = input$max_time1,
          cohort_label = "Subgroup sample size:"
        )
      }
    )
  })
  
  output$pred_plot <- renderPlot({
    req(input$pred_ref_treatment_selection)
    req(input$pred_treatment_selection)
    
    withProgress(
      message = "Building the prediction survival/hazard ratio plot",
      detail = "This may take a moment.",
      value = 0, {
        prediction_plot_result <- pred_plot_store()
        
        cohort_label <- if (isTRUE(input$pred_use_subgroup)) {
          "Prediction subgroup sample size:"
        } else {
          "Prediction sample size:"
        }
        
        make_subgroup_survival_or_hr_plot(
          subgroup_plot_result = prediction_plot_result,
          plot_type = input$pred_plot_type,
          cntrl = input$pred_ref_treatment_selection,
          trt = input$pred_treatment_selection,
          min_time = input$pred_min_time,
          max_time = input$pred_max_time,
          cohort_label = cohort_label
        )
      }
    )
  })
  
  output$downloadplot <- downloadHandler(
    filename = function() {
      paste0("plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      main_plot_result <- plot_store()
      
      plt_str <- make_main_survival_or_hr_plot(
        plot_store_result = main_plot_result,
        plot_type = input$plot_type,
        cntrl = input$ref_treatment_selection,
        trt = input$treatment_selection,
        min_time = input$min_time,
        max_time = input$max_time
      )
      
      save_plot_with_optional_size(
        file = file,
        plot = plt_str,
        width_input = input$width,
        height_input = input$height,
        output_width_px = session$clientData$output_plot_width,
        output_height_px = session$clientData$output_plot_height
      )
    }
  )
  
  output$downloadrmstplot <- downloadHandler(
    filename = function() {
      paste0("rmst_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      plt_str2 <- make_rmst_combined_plot(
        RMST_store = rmst_store(),
        cntrl = input$ref_treatment_selection1,
        trt = input$treatment_selection1,
        horizon_time = input$horizon_time
      )
      
      save_plot_with_optional_size(
        file = file,
        plot = plt_str2,
        width_input = input$width1,
        height_input = input$height1,
        output_width_px = session$clientData$output_rmst_plot_width,
        output_height_px = session$clientData$output_rmst_plot_height
      )
    }
  )
  
  output$downloadsubgroupplot <- downloadHandler(
    filename = function() {
      paste0("subgroup_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      subgroup_plot_result <- plot_store1()
      
      plt_str1 <- make_subgroup_survival_or_hr_plot(
        subgroup_plot_result = subgroup_plot_result,
        plot_type = input$plot_type1,
        cntrl = input$ref_treatment_selection2,
        trt = input$treatment_selection2,
        min_time = input$min_time1,
        max_time = input$max_time1,
        cohort_label = "Subgroup sample size:"
      )
      
      save_plot_with_optional_size(
        file = file,
        plot = plt_str1,
        width_input = input$width2,
        height_input = input$height2,
        output_width_px = session$clientData$output_subgroup_plot_width,
        output_height_px = session$clientData$output_subgroup_plot_height
      )
    }
  )
  
  output$downloadpredplot <- downloadHandler(
    filename = function() {
      paste0("prediction_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      prediction_plot_result <- pred_plot_store()
      
      cohort_label <- if (isTRUE(input$pred_use_subgroup)) {
        "Prediction subgroup sample size:"
      } else {
        "Prediction sample size:"
      }
      
      plt_str <- make_subgroup_survival_or_hr_plot(
        subgroup_plot_result = prediction_plot_result,
        plot_type = input$pred_plot_type,
        cntrl = input$pred_ref_treatment_selection,
        trt = input$pred_treatment_selection,
        min_time = input$pred_min_time,
        max_time = input$pred_max_time,
        cohort_label = cohort_label
      )
      
      save_plot_with_optional_size(
        file = file,
        plot = plt_str,
        width_input = input$pred_width,
        height_input = input$pred_height,
        output_width_px = session$clientData$output_pred_plot_width,
        output_height_px = session$clientData$output_pred_plot_height
      )
    }
  )
  
  output$subgroup_rmst_table <- renderDT({
    trt_names <- rev(gsub("mu_", "", names(plot_store1()$thetas[1:2])))
    
    as.data.frame(
      round(data_table_store()[[1]], 2),
      row.names = c(trt_names[2], trt_names[1], "Difference", "Ratio")
    )
  },
  options = list(
    pageLength = 5,
    autoWidth = TRUE,
    dom = "t"
  ))
  
  output$pred_rmst_table <- renderDT({
    trt_names <- rev(gsub("mu_", "", names(pred_plot_store()$thetas[1:2])))
    
    as.data.frame(
      round(pred_data_table_store()[[1]], 2),
      row.names = c(trt_names[2], trt_names[1], "Difference", "Ratio")
    )
  },
  options = list(
    pageLength = 5,
    autoWidth = TRUE,
    dom = "t"
  ))
  
  output$downloadRMSTtable <- downloadHandler(
    filename = function() {
      paste0("data_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      trt_names <- rev(gsub("mu_", "", names(plot_store1()$thetas[1:2])))
      
      write.csv(
        as.data.frame(
          data_table_store()[[1]],
          row.names = c(trt_names[2], trt_names[1], "Difference", "Ratio")
        ),
        file
      )
    }
  )
  
  output$downloadPredRMSTtable <- downloadHandler(
    filename = function() {
      paste0("prediction_rmst_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      trt_names <- rev(gsub("mu_", "", names(pred_plot_store()$thetas[1:2])))
      
      write.csv(
        as.data.frame(
          pred_data_table_store()[[1]],
          row.names = c(trt_names[2], trt_names[1], "Difference", "Ratio")
        ),
        file
      )
    }
  )
}

# Run the app ----
shinyApp(ui, server)