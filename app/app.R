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
          tags$li("If the patient is alive at the end of follow-up, this value is TRUE; otherwise FALSE.")
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
      downloadButton("downloadresults", "Download Bayesian-ECD Results (.zip)"),
      downloadButton("downloadplot", "Download Plot (.png)")
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
      downloadButton("downloadrmstplot", "Download RMST Plot (.png)")
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
      downloadButton("downloadsubgroupplot", "Download Plot (.png)")
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
      downloadButton("downloadRMSTtable", "Download RMST Table (.csv)")
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
      downloadButton("downloadpredplot", "Download Prediction Plot (.png)")
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
      downloadButton("downloadPredRMSTtable", "Download Prediction RMST Table (.csv)")
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
  
  output$plot <- renderPlot({
    req(input$ref_treatment_selection)
    req(input$treatment_selection)
    
    withProgress(
      message = "Building the survival/hazard ratio plot",
      detail = "This may take a moment.",
      value = 0, {
        plt <- switch(
          input$plot_type,
          "Survival probabilities vs. Time (Days)" = "Survival",
          "Hazard Ratio vs. Time (Days)" = "HR"
        )
        
        cntrl <- input$ref_treatment_selection
        trt <- input$treatment_selection
        
        if (plt == "Survival") {
          plot_store()$p.surv +
            scale_fill_discrete(labels = c(trt, cntrl)) +
            scale_color_discrete(labels = c(trt, cntrl)) +
            labs(x = "Time") +
            theme_bw() +
            large_plot_theme +
            scale_x_continuous(limits = c(input$min_time - 10, input$max_time + 100))
        } else {
          plot_store()$p.hazarad.ratio.withCI +
            labs(x = "Time") +
            theme_bw() +
            large_plot_theme +
            scale_x_continuous(limits = c(input$min_time - 10, input$max_time + 100))
        }
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
        cntrl <- input$ref_treatment_selection1
        trt <- input$treatment_selection1
        
        RMST_store <- rmst_lognormal(
          result(),
          cntrl,
          trt,
          input$horizon_time,
          burnin = 200
        )
        
        tmp1 <- RMST_store$RMST
        tmp2 <- RMST_store$Difference
        tmp2$`Difference (Days)` <- round(tmp2$`Difference (Days)`, 2)
        tmp3 <- RMST_store$Ratio
        tmp3$Ratio <- round(tmp3$Ratio, 2)
        
        p.trt <- RMST_store$p.surv.trt +
          annotate(
            "text",
            x = input$horizon_time / 2,
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
            x = input$horizon_time / 2,
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
            tableGrob(t(tmp2), cols = NULL, theme = ttheme_default(base_size = RMST_TABLE_SIZE)) +
            tableGrob(t(tmp3), cols = NULL, theme = ttheme_default(base_size = RMST_TABLE_SIZE))) +
          plot_layout(heights = c(7, 1))
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
        plt <- switch(
          input$plot_type1,
          "Survival probabilities vs. Time (Days)" = "Survival",
          "Hazard Ratio vs. Time (Days)" = "HR"
        )
        
        cntrl <- input$ref_treatment_selection2
        trt <- input$treatment_selection2
        
        if (plt == "Survival") {
          plot_store1()$p.surv +
            scale_fill_discrete(labels = c(trt, cntrl)) +
            scale_color_discrete(labels = c(trt, cntrl)) +
            annotate(
              "text",
              x = input$max_time1 / 1.5,
              y = max(plot_store1()$surv.data$surv.pred.97.5, na.rm = TRUE) / 1.5,
              label = paste("Subgroup sample size:", length(plot_store1()$subgroup_index)),
              size = SUBGROUP_ANNOTATION_SIZE
            ) +
            labs(x = "Time") +
            theme_bw() +
            large_plot_theme +
            scale_x_continuous(limits = c(input$min_time1 - 10, input$max_time1 + 100))
        } else {
          plot_store1()$p.hazarad.ratio.withCI +
            annotate(
              "text",
              x = input$max_time1 / 1.5,
              y = max(plot_store1()$HR.data$HR.pred.97.5, na.rm = TRUE) / 1.5,
              label = paste("Subgroup sample size:", length(plot_store1()$subgroup_index)),
              size = SUBGROUP_ANNOTATION_SIZE
            ) +
            labs(x = "Time") +
            theme_bw() +
            large_plot_theme +
            scale_x_continuous(limits = c(input$min_time1 - 10, input$max_time1 + 100))
        }
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
        plt <- switch(
          input$pred_plot_type,
          "Survival probabilities vs. Time (Days)" = "Survival",
          "Hazard Ratio vs. Time (Days)" = "HR"
        )
        
        cntrl <- input$pred_ref_treatment_selection
        trt <- input$pred_treatment_selection
        
        cohort_label <- if (isTRUE(input$pred_use_subgroup)) {
          "Prediction subgroup sample size:"
        } else {
          "Prediction sample size:"
        }
        
        if (plt == "Survival") {
          pred_plot_store()$p.surv +
            scale_fill_discrete(labels = c(trt, cntrl)) +
            scale_color_discrete(labels = c(trt, cntrl)) +
            annotate(
              "text",
              x = input$pred_max_time / 1.5,
              y = max(pred_plot_store()$surv.data$surv.pred.97.5, na.rm = TRUE) / 1.5,
              label = paste(cohort_label, length(pred_plot_store()$subgroup_index)),
              size = SUBGROUP_ANNOTATION_SIZE
            ) +
            labs(x = "Time") +
            theme_bw() +
            large_plot_theme +
            scale_x_continuous(limits = c(input$pred_min_time - 10, input$pred_max_time + 100))
        } else {
          pred_plot_store()$p.hazarad.ratio.withCI +
            annotate(
              "text",
              x = input$pred_max_time / 1.5,
              y = max(pred_plot_store()$HR.data$HR.pred.97.5, na.rm = TRUE) / 1.5,
              label = paste(cohort_label, length(pred_plot_store()$subgroup_index)),
              size = SUBGROUP_ANNOTATION_SIZE
            ) +
            labs(x = "Time") +
            theme_bw() +
            large_plot_theme +
            scale_x_continuous(limits = c(input$pred_min_time - 10, input$pred_max_time + 100))
        }
      }
    )
  })
  
  output$downloadplot <- downloadHandler(
    filename = function() {
      paste0("plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      plt <- switch(
        input$plot_type,
        "Survival probabilities vs. Time (Days)" = "Survival",
        "Hazard Ratio vs. Time (Days)" = "HR"
      )
      
      cntrl <- input$ref_treatment_selection
      trt <- input$treatment_selection
      
      if (plt == "Survival") {
        plt_str <- plot_store()$p.surv +
          scale_fill_discrete(labels = c(trt, cntrl)) +
          scale_color_discrete(labels = c(trt, cntrl)) +
          labs(x = "Time") +
          theme_bw() +
          large_plot_theme +
          scale_x_continuous(limits = c(input$min_time - 10, input$max_time + 100))
      } else {
        plt_str <- plot_store()$p.hazarad.ratio.withCI +
          labs(x = "Time") +
          theme_bw() +
          large_plot_theme +
          scale_x_continuous(limits = c(input$min_time - 10, input$max_time + 100))
      }
      
      if (is.na(input$width) && is.na(input$height)) {
        plot_width <- as.numeric(session$clientData$output_plot_width)
        plot_height <- as.numeric(session$clientData$output_plot_height)
        ggsave(file, plot = plt_str, device = "png", width = plot_width,
               height = plot_height, units = "px", dpi = 96)
      } else if (!is.na(input$width) && is.na(input$height)) {
        ggsave(file, plot = plt_str, device = "png", width = input$width, units = "in")
      } else if (is.na(input$width) && !is.na(input$height)) {
        ggsave(file, plot = plt_str, device = "png", height = input$height, units = "in")
      } else {
        ggsave(file, plot = plt_str, device = "png",
               width = input$width, height = input$height, units = "in")
      }
    }
  )
  
  output$downloadrmstplot <- downloadHandler(
    filename = function() {
      paste0("plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      cntrl <- input$ref_treatment_selection1
      trt <- input$treatment_selection1
      
      RMST_store <- rmst_lognormal(
        result(),
        cntrl,
        trt,
        input$horizon_time,
        burnin = 200
      )
      
      tmp1 <- RMST_store$RMST
      tmp2 <- RMST_store$Difference
      tmp2$`Difference (Days)` <- round(tmp2$`Difference (Days)`, 2)
      tmp3 <- RMST_store$Ratio
      tmp3$Ratio <- round(tmp3$Ratio, 2)
      
      p.trt <- RMST_store$p.surv.trt +
        annotate(
          "text",
          x = input$horizon_time / 2,
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
          x = input$horizon_time / 2,
          y = 0.5,
          label = paste("RMST:", round(tmp1[2, 2]), "Days"),
          size = RMST_ANNOTATION_SIZE
        ) +
        labs(x = "Time") +
        theme_bw() +
        large_rmst_theme
      
      plt_str2 <- ((p.trt + p.cntrl + plot_layout(axis_titles = "collect") &
                      plot_annotation(tag_levels = list(c(trt, cntrl), "1")) &
                      theme(plot.tag = element_text(size = PLOT_TAG_SIZE, face = "bold"))) +
                     tableGrob(t(tmp2), cols = NULL, theme = ttheme_default(base_size = RMST_TABLE_SIZE)) +
                     tableGrob(t(tmp3), cols = NULL, theme = ttheme_default(base_size = RMST_TABLE_SIZE))) +
        plot_layout(heights = c(7, 1))
      
      if (is.na(input$width1) && is.na(input$height1)) {
        plot_width <- as.numeric(session$clientData$output_rmst_plot_width)
        plot_height <- as.numeric(session$clientData$output_rmst_plot_height)
        ggsave(file, plot = plt_str2, device = "png", width = plot_width,
               height = plot_height, units = "px", dpi = 96)
      } else if (!is.na(input$width1) && is.na(input$height1)) {
        ggsave(file, plot = plt_str2, device = "png", width = input$width1, units = "in")
      } else if (is.na(input$width1) && !is.na(input$height1)) {
        ggsave(file, plot = plt_str2, device = "png", height = input$height1, units = "in")
      } else {
        ggsave(file, plot = plt_str2, device = "png",
               width = input$width1, height = input$height1, units = "in")
      }
    }
  )
  
  output$downloadsubgroupplot <- downloadHandler(
    filename = function() {
      paste0("plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      plt <- switch(
        input$plot_type1,
        "Survival probabilities vs. Time (Days)" = "Survival",
        "Hazard Ratio vs. Time (Days)" = "HR"
      )
      
      cntrl <- input$ref_treatment_selection2
      trt <- input$treatment_selection2
      
      if (plt == "Survival") {
        plt_str1 <- plot_store1()$p.surv +
          scale_fill_discrete(labels = c(trt, cntrl)) +
          scale_color_discrete(labels = c(trt, cntrl)) +
          labs(x = "Time") +
          theme_bw() +
          large_plot_theme +
          scale_x_continuous(limits = c(input$min_time1 - 10, input$max_time1 + 100))
      } else {
        plt_str1 <- plot_store1()$p.hazarad.ratio.withCI +
          labs(x = "Time") +
          theme_bw() +
          large_plot_theme +
          scale_x_continuous(limits = c(input$min_time1 - 10, input$max_time1 + 100))
      }
      
      if (is.na(input$width2) && is.na(input$height2)) {
        plot_width <- as.numeric(session$clientData$output_subgroup_plot_width)
        plot_height <- as.numeric(session$clientData$output_subgroup_plot_height)
        ggsave(file, plot = plt_str1, device = "png", width = plot_width,
               height = plot_height, units = "px", dpi = 96)
      } else if (!is.na(input$width2) && is.na(input$height2)) {
        ggsave(file, plot = plt_str1, device = "png", width = input$width2, units = "in")
      } else if (is.na(input$width2) && !is.na(input$height2)) {
        ggsave(file, plot = plt_str1, device = "png", height = input$height2, units = "in")
      } else {
        ggsave(file, plot = plt_str1, device = "png",
               width = input$width2, height = input$height2, units = "in")
      }
    }
  )
  
  output$downloadpredplot <- downloadHandler(
    filename = function() {
      paste0("prediction_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      plt <- switch(
        input$pred_plot_type,
        "Survival probabilities vs. Time (Days)" = "Survival",
        "Hazard Ratio vs. Time (Days)" = "HR"
      )
      
      cntrl <- input$pred_ref_treatment_selection
      trt <- input$pred_treatment_selection
      
      cohort_label <- if (isTRUE(input$pred_use_subgroup)) {
        "Prediction subgroup sample size:"
      } else {
        "Prediction sample size:"
      }
      
      if (plt == "Survival") {
        plt_str <- pred_plot_store()$p.surv +
          scale_fill_discrete(labels = c(trt, cntrl)) +
          scale_color_discrete(labels = c(trt, cntrl)) +
          annotate(
            "text",
            x = input$pred_max_time / 1.5,
            y = max(pred_plot_store()$surv.data$surv.pred.97.5, na.rm = TRUE) / 1.5,
            label = paste(cohort_label, length(pred_plot_store()$subgroup_index)),
            size = SUBGROUP_ANNOTATION_SIZE
          ) +
          labs(x = "Time") +
          theme_bw() +
          large_plot_theme +
          scale_x_continuous(limits = c(input$pred_min_time - 10, input$pred_max_time + 100))
      } else {
        plt_str <- pred_plot_store()$p.hazarad.ratio.withCI +
          annotate(
            "text",
            x = input$pred_max_time / 1.5,
            y = max(pred_plot_store()$HR.data$HR.pred.97.5, na.rm = TRUE) / 1.5,
            label = paste(cohort_label, length(pred_plot_store()$subgroup_index)),
            size = SUBGROUP_ANNOTATION_SIZE
          ) +
          labs(x = "Time") +
          theme_bw() +
          large_plot_theme +
          scale_x_continuous(limits = c(input$pred_min_time - 10, input$pred_max_time + 100))
      }
      
      if (is.na(input$pred_width) && is.na(input$pred_height)) {
        plot_width <- as.numeric(session$clientData$output_pred_plot_width)
        plot_height <- as.numeric(session$clientData$output_pred_plot_height)
        ggsave(file, plot = plt_str, device = "png", width = plot_width,
               height = plot_height, units = "px", dpi = 96)
      } else if (!is.na(input$pred_width) && is.na(input$pred_height)) {
        ggsave(file, plot = plt_str, device = "png", width = input$pred_width, units = "in")
      } else if (is.na(input$pred_width) && !is.na(input$pred_height)) {
        ggsave(file, plot = plt_str, device = "png", height = input$pred_height, units = "in")
      } else {
        ggsave(file, plot = plt_str, device = "png",
               width = input$pred_width, height = input$pred_height, units = "in")
      }
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