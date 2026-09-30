#260915

library(shiny)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ggplot2)
library(lubridate)
library(iglu)
library(FNN)

# Shiny's default upload limit is 5 MB. Allow CGM CSV files up to 30 MB.
options(shiny.maxRequestSize = 30 * 1024^2)

key_metrics <- c(
  "TITR", "70_to_180", "above_180", "below_70",
  "hbgi", "lbgi", "cv_glu", "mage"
)

# =====================================================
# Load server-side reference CSV
# =====================================================

REF_PATH <- "ref_umap.csv"

if (!file.exists(REF_PATH)) {
  stop(
    "ref_umap.csv was not found. Place it in the same folder as app.R before deployment."
  )
}

ref_df_server <- read.csv(
  REF_PATH,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_ref_columns <- c(
  "State", "DDR_dim1", "DDR_dim2", "Pseudotime", key_metrics
)

missing_ref_columns <- setdiff(required_ref_columns, colnames(ref_df_server))

if (length(missing_ref_columns) > 0) {
  stop(
    paste(
      "ref_umap.csv is missing required columns:",
      paste(missing_ref_columns, collapse = ", ")
    )
  )
}
# -----------------------------
# State -> Clinical cluster mapping
# -----------------------------
state_to_cluster2 <- c(
  "1" = "EHHR",
  "2" = "SE",
  "3" = "MH",
  "4" = "HHV",
  "5" = "HLV"
)

ref_df_server$State <- as.character(ref_df_server$State)

ref_df_server$clinical_cluster2 <- unname(
  state_to_cluster2[ref_df_server$State]
)

# -----------------------------
# Helper functions
# -----------------------------

cluster_descriptions <- tibble(
  assigned_cluster = c(
    "SE",
    "EHHR",
    "MH",
    "HHV",
    "HLV"
  ),
  
  cluster_label = c(
    "SE",
    "EHHR",
    "MH",
    "HHV",
    "HLV"
  ),
  
  cluster_name = c(
    "Stable Euglycemic",
    "Euglycemic high-hypoglycemia risk",
    "Mild hyperglycemic",
    "Hyperglycemia high Variability",
    "Hyperglycemia low Variability"
  ),
  
  cluster_description = c(
    "Glucose levels remain consistently within the target range. Preserved glucose homeostasis with maintained β-cell function.",
    "Overall glucose levels remain within the target range. Increased hypoglycemia risk and glycemic variability. Associated with lower muscle mass and body water.",
    "Mild elevation of glucose levels. Early metabolic dysfunction with declining β-cell function and reduced insulin sensitivity.",
    "Frequent and large glucose fluctuations. β-cell dysfunction with relatively preserved insulin sensitivity.",
    "Persistent hyperglycemia with limited glucose variability. β-cell dysfunction accompanied by insulin resistance. Associated with higher adiposity and visceral fat accumulation."
  )
)
read_raw_cgm <- function(path, patient_id = "uploaded_patient") {
  
  # -----------------------------
  # First, read first line to detect file type
  # -----------------------------
  first_line <- readLines(path, n = 1, encoding = "UTF-8")
  
  # -----------------------------
  # Type 2: header contains Korean column names
  # time = column 8
  # glucose = column 9
  # -----------------------------
  if (grepl("측정 및 기록 시간", first_line) || grepl("혈당 값", first_line)) {
    
    df <- read.csv(
      path,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      fileEncoding = "UTF-8-BOM"
    )
    
    if (ncol(df) < 9) {
      stop("Type 2 CGM file must have at least 9 columns.")
    }
    
    out <- df %>%
      transmute(
        id = patient_id,
        time = as.POSIXct(
          .[[8]],
          format = "%Y.%m.%d %H:%M:%S",
          tz = "UTC"
        ),
        gl = suppressWarnings(as.numeric(.[[9]]))
      ) %>%
      drop_na(time, gl) %>%
      arrange(time) %>%
      distinct(time, .keep_all = TRUE)
    
    return(out)
  }
  
  # -----------------------------
  # Type 1: Libre format
  # skip first 3 rows
  # time = column 3
  # glucose = column 5 or 6
  # -----------------------------
  df <- read.csv(
    path,
    header = FALSE,
    skip = 3,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    fileEncoding = "UTF-8-BOM"
  )
  
  if (ncol(df) < 6) {
    stop("Type 1 CGM file must have at least 6 columns after skipping first 3 rows.")
  }
  
  out <- df %>%
    transmute(
      id = patient_id,
      time = as.POSIXct(
        V3,
        format = "%Y-%m-%d %H:%M",
        tz = "UTC"
      ),
      gl = dplyr::coalesce(
        suppressWarnings(as.numeric(V5)),
        suppressWarnings(as.numeric(V6))
      )
    ) %>%
    drop_na(time, gl) %>%
    arrange(time) %>%
    distinct(time, .keep_all = TRUE)
  
  return(out)
}

calculate_iglu_8metrics <- function(raw_df) {
  raw_df <- raw_df %>%
    mutate(
      id = as.character(id),
      time = as.POSIXct(time),
      gl = as.numeric(gl)
    ) %>%
    drop_na(id, time, gl)
  
  titr_value <- mean(raw_df$gl >= 70 & raw_df$gl <= 140) * 100
  tir_value <- mean(raw_df$gl >= 70 & raw_df$gl <= 180) * 100
  above_180_value <- mean(raw_df$gl > 180) * 100
  below_70_value <- mean(raw_df$gl < 70) * 100
  
  hbgi_value <- iglu::hbgi(raw_df)$HBGI[1]
  lbgi_value <- iglu::lbgi(raw_df)$LBGI[1]
  cv_value   <- iglu::cv_glu(raw_df)$CV[1]
  mage_value <- iglu::mage(raw_df)$MAGE[1]
  
  above_180_value <- log1p(above_180_value)
  below_70_value  <- log1p(below_70_value)
  hbgi_value      <- log1p(hbgi_value)
  lbgi_value      <- log1p(lbgi_value)
  mage_value      <- log1p(mage_value)
  
  tibble(
    id = raw_df$id[1],
    TITR = titr_value,
    `70_to_180` = tir_value,
    above_180 = above_180_value,
    below_70 = below_70_value,
    hbgi = hbgi_value,
    lbgi = lbgi_value,
    cv_glu = cv_value,
    mage = mage_value
  )
}
detect_cgm_followup_periods <- function(raw_df,
                                        gap_days = 60,
                                        window_days = 14,
                                        max_periods = 2,
                                        visit_gap_days = 7) {
  
  raw_df <- raw_df %>%
    mutate(
      time = as.POSIXct(time, tz = "UTC"),
      gl = as.numeric(gl)
    ) %>%
    drop_na(time, gl) %>%
    arrange(time)
  
  overall_start <- min(raw_df$time, na.rm = TRUE)
  overall_end   <- max(raw_df$time, na.rm = TRUE)
  overall_span_days <- as.numeric(difftime(overall_end, overall_start, units = "days"))
  
  # Detect smaller CGM blocks/visits using 7-day gap by default
  raw_df <- raw_df %>%
    mutate(
      gap_from_prev_days = as.numeric(difftime(time, lag(time), units = "days")),
      new_visit = ifelse(
        is.na(gap_from_prev_days),
        0,
        ifelse(gap_from_prev_days >= visit_gap_days, 1, 0)
      ),
      visit_id = cumsum(new_visit) + 1
    )
  
  visit_summary <- raw_df %>%
    group_by(visit_id) %>%
    summarise(
      visit_start = min(time, na.rm = TRUE),
      visit_end = max(time, na.rm = TRUE),
      n_points = n(),
      duration_days = as.numeric(difftime(max(time), min(time), units = "days")) + 1,
      .groups = "drop"
    ) %>%
    arrange(visit_start)
  
  # --------------------------------------------------
  # Case 1: overall span is short -> use only latest 2 weeks
  # --------------------------------------------------
  if (overall_span_days < gap_days) {
    
    selected_windows <- tibble(
      period_order = 1,
      period_label = "Uploaded CGM",
      period_start = overall_end - lubridate::days(window_days),
      period_end = overall_end
    )
    
  } else {
    
    # --------------------------------------------------
    # Case 2: enough total follow-up span
    # Select earliest CGM block and latest CGM block
    # --------------------------------------------------
    first_visit <- visit_summary %>%
      arrange(visit_start) %>%
      slice(1)
    
    last_visit <- visit_summary %>%
      arrange(desc(visit_end)) %>%
      slice(1)
    
    months_gap <- round(
      as.numeric(difftime(last_visit$visit_end, first_visit$visit_end, units = "days")) / 30.44
    )
    
    baseline_end <- first_visit$visit_end
    baseline_start <- max(
      first_visit$visit_start,
      baseline_end - lubridate::days(window_days)
    )
    
    followup_end <- last_visit$visit_end
    followup_start <- max(
      last_visit$visit_start,
      followup_end - lubridate::days(window_days)
    )
    
    selected_windows <- tibble(
      period_order = c(1, 2),
      period_label = c(
        "Baseline 2 weeks",
        paste0("Follow-up ~", months_gap, "M")
      ),
      period_start = c(baseline_start, followup_start),
      period_end = c(baseline_end, followup_end)
    )
  }
  
  period_data <- purrr::map_dfr(seq_len(nrow(selected_windows)), function(i) {
    
    w <- selected_windows[i, ]
    
    raw_df %>%
      filter(
        time > w$period_start,
        time <= w$period_end
      ) %>%
      mutate(
        period_order = w$period_order,
        period_label = w$period_label,
        period_start = w$period_start,
        period_end = w$period_end
      )
  })
  
  list(
    visit_summary = visit_summary,
    selected_visits = selected_windows,
    period_data = period_data
  )
}


preprocess_cgm_metrics <- function(df, ref_mean = NULL, ref_sd = NULL) {
  x <- df[, key_metrics]
  
  x <- x %>%
    mutate(across(everything(), as.numeric))
  
  #  for (col in log_cols) {
  #    x[[col]] <- log1p(x[[col]])
  #  }
  
  if (is.null(ref_mean)) {
    ref_mean <- sapply(x, mean, na.rm = TRUE)
  }
  
  if (is.null(ref_sd)) {
    ref_sd <- sapply(x, sd, na.rm = TRUE)
    ref_sd[ref_sd == 0] <- 1
  }
  
  x_scaled <- sweep(x, 2, ref_mean, "-")
  x_scaled <- sweep(x_scaled, 2, ref_sd, "/")
  
  list(
    x_scaled = as.data.frame(x_scaled),
    ref_mean = ref_mean,
    ref_sd = ref_sd
  )
}

project_to_ddrtree <- function(new_metrics_df, ref_df, k = 10) {
  
  required_ref <- c(
    key_metrics,
    "DDR_dim1",
    "DDR_dim2",
    "clinical_cluster2"
  )
  
  missing_ref <- setdiff(required_ref, colnames(ref_df))
  if (length(missing_ref) > 0) {
    stop(paste("Reference file missing:", paste(missing_ref, collapse = ", ")))
  }
  
  missing_new <- setdiff(key_metrics, colnames(new_metrics_df))
  if (length(missing_new) > 0) {
    stop(paste("New metrics missing:", paste(missing_new, collapse = ", ")))
  }
  
  ref_df <- ref_df %>%
    drop_na(all_of(required_ref))
  
  ref_processed <- preprocess_cgm_metrics(ref_df)
  
  new_processed <- preprocess_cgm_metrics(
    new_metrics_df,
    ref_mean = ref_processed$ref_mean,
    ref_sd = ref_processed$ref_sd
  )
  
  k_use <- min(k, nrow(ref_df))
  
  nn <- FNN::get.knnx(
    data = ref_processed$x_scaled,
    query = new_processed$x_scaled,
    k = k_use
  )
  
  idx <- nn$nn.index[1, ]
  dist <- nn$nn.dist[1, ]
  
  neighbors <- ref_df[idx, ]
  
  projected <- tibble(
    id = new_metrics_df$id[1],
    assigned_cluster = neighbors$clinical_cluster2[1],
    DDR_dim1_projected = neighbors$DDR_dim1[1],
    DDR_dim2_projected = neighbors$DDR_dim2[1],
    nearest_distance_mean = mean(dist),
    nearest_distance_min = min(dist)
  )
  
  projected %>%
    left_join(cluster_descriptions, by = "assigned_cluster")
}

# -----------------------------
# UI
# -----------------------------

ui <- fluidPage(
  titlePanel("CGM Heterogeneity DDRTree Projection"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      fileInput(
        "raw_cgm_file",
        "Upload raw CGM CSV",
        accept = c(".csv")
      ),
      
      helpText(
        "Supported formats: Libre CSV (three metadata rows, time in column 3, glucose in columns 5/6) or the Korean-header CGM CSV (time in column 8, glucose in column 9)."
      ),
      
      tags$div(
        style = "margin: 12px 0; padding: 10px; border-left: 4px solid #d9534f; background: #fff7f7; font-size: 13px;",
        tags$strong("Privacy notice"),
        tags$br(),
        "Upload de-identified CGM files only. This app does not intentionally save uploaded files after the active session."
      ),
      
      helpText(
        "Research-use only. The projected phenotype is not a medical diagnosis or treatment recommendation."
      ),
      
      
      selectInput(
        "color_var",
        "Color by",
        choices = c(
          "CGM Glycemic Phenotype" = "clinical_cluster2",
          "Pseudotime" = "Pseudotime"
        ),
        selected = "clinical_cluster2"
      ),
      
      actionButton(
        "run_btn",
        "Run IGLU + Projection"
      ),
      
      checkboxInput(
        "followup_mode",
        "Auto-detect follow-up periods",
        value = TRUE
      ),
      
      numericInput(
        "gap_days",
        "Minimum gap between periods (days)",
        value = 60,
        min = 30,
        max = 365,
        step = 15
      ),
      
      numericInput(
        "window_days",
        "CGM window per period (days)",
        value = 14,
        min = 7,
        max = 21,
        step = 1
      ),
      
      
      
    ),
    
    mainPanel(
      width = 9,
      
      tabsetPanel(
        tabPanel(
          "DDRTree Plot",
          plotOutput(
            "ddrtree_plot",
            height = "720px"
          ),
          br(),
          uiOutput("cluster_interpretation_box")
        ),
        
        tabPanel(
          "AGP Plot",
          plotOutput(
            "agp_plot",
            height = "650px"
          )
        )
      )
    )
  )
)

# -----------------------------
# Server
# -----------------------------
server <- function(input, output, session) {
  
  ref_df <- reactive({
    ref_df_server
  })
  
  raw_cgm_df <- eventReactive(input$run_btn, {
    req(input$raw_cgm_file)
    
    validate(
      need(
        identical(tolower(tools::file_ext(input$raw_cgm_file$name)), "csv"),
        "Please upload a CSV file."
      )
    )
    
    out <- tryCatch(
      read_raw_cgm(
        path = input$raw_cgm_file$datapath,
        patient_id = "uploaded_patient"
      ),
      error = function(e) {
        validate(
          need(FALSE, paste("Could not read the CGM file:", conditionMessage(e)))
        )
      }
    )
    
    validate(
      need(nrow(out) >= 2, "No usable glucose observations were found in the uploaded file.")
    )
    
    out
  })
  
  followup_periods <- eventReactive(input$run_btn, {
    
    if (isTRUE(input$followup_mode)) {
      detect_cgm_followup_periods(
        raw_df = raw_cgm_df(),
        gap_days = input$gap_days,
        window_days = input$window_days,
        max_periods = 2
      )
    } else {
      list(
        visit_summary = NULL,
        selected_visits = NULL,
        period_data = raw_cgm_df() %>%
          mutate(
            period_order = 1,
            period_label = "Uploaded CGM",
            period_start = min(time, na.rm = TRUE),
            period_end = max(time, na.rm = TRUE)
          )
      )
    }
  })
  metrics_df <- eventReactive(input$run_btn, {
    
    period_data <- followup_periods()$period_data
    
    period_data %>%
      group_by(period_order, period_label, period_start, period_end) %>%
      group_modify(~ {
        out <- calculate_iglu_8metrics(.x)
        out$id <- paste0("uploaded_patient_", .y$period_label)
        out
      }) %>%
      ungroup()
  })
  
  projected_df <- eventReactive(input$run_btn, {
    
    m <- metrics_df()
    
    purrr::map_dfr(seq_len(nrow(m)), function(i) {
      
      proj <- project_to_ddrtree(
        new_metrics_df = m[i, ],
        ref_df = ref_df(),
        k = 10
      )
      
      proj %>%
        mutate(
          period_order = m$period_order[i],
          period_label = m$period_label[i],
          period_start = m$period_start[i],
          period_end = m$period_end[i]
        )
    })
  })
  output$agp_plot <- renderPlot({
    req(raw_cgm_df())
    
    agp_df <- raw_cgm_df() %>%
      mutate(
        time = as.POSIXct(time, tz = "UTC"),
        tod_min = hour(time) * 60 + minute(time),
        tod_bin = floor(tod_min / 15) * 15   # 15분 bin
      ) %>%
      group_by(tod_bin) %>%
      summarise(
        p05 = quantile(gl, 0.05, na.rm = TRUE),
        p25 = quantile(gl, 0.25, na.rm = TRUE),
        p50 = quantile(gl, 0.50, na.rm = TRUE),
        p75 = quantile(gl, 0.75, na.rm = TRUE),
        p95 = quantile(gl, 0.95, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        clock_time = as.POSIXct("2000-01-01 00:00:00", tz = "UTC") + tod_bin * 60
      )
    x_min <- min(agp_df$clock_time, na.rm = TRUE)
    p <- ggplot(agp_df, aes(x = clock_time)) +
      annotate(
        "rect",
        xmin = min(agp_df$clock_time, na.rm = TRUE),
        xmax = max(agp_df$clock_time, na.rm = TRUE),
        ymin = 70,
        ymax = 180,
        fill = "gray85",
        alpha = 0.25
      ) +
      
      geom_ribbon(aes(ymin = p05, ymax = p95),
                  fill = "#dbeeff", alpha = 0.8) +
      geom_ribbon(aes(ymin = p25, ymax = p75),
                  fill = "#7ec9f5", alpha = 0.9) +
      geom_line(aes(y = p50), color = "black", linewidth = 1.4) +
      geom_hline(yintercept = 180, color = "#D62728", linewidth = 0.8) +
      geom_hline(yintercept = 70, color = "#D62728", linewidth = 0.8) +
      geom_hline(yintercept = 54, color = "#4db8ff", linetype = "dashed", linewidth = 0.8) +
      geom_hline(yintercept = 250, color = "#4db8ff", linetype = "dashed", linewidth = 0.8) +
      scale_x_datetime(
        date_labels = "%I %p",
        date_breaks = "3 hours",
        expand = c(0.01, 0.01)
      ) +
      coord_cartesian(
        ylim = c(0, 350),
        clip = "off"
      ) +
      labs(
        title = "Ambulatory Glucose Profile (AGP)",
        x = "Time of day",
        y = NULL
      ) +
      annotate(
        "text",
        x = x_min,
        y = 125,
        label = "TIR",
        angle = 90,
        hjust = 0.5,
        vjust = 2.5,
        size = 5,
        fontface = "bold",
        color = "#D62728"
      ) +
      theme_bw(base_size = 14) +
      theme(
        axis.title.y = element_blank(),
        plot.title = element_text(face = "bold", size = 18),
        axis.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 13, face = "bold"),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 95)
      )
    
    print(p)
    
    
    grid::grid.text(
      "Glucose (mg/dL)",
      x = grid::unit(0.065, "npc"),
      y = grid::unit(0.52, "npc"),
      rot = 90,
      gp = grid::gpar(
        fontsize = 14,
        fontface = "bold",
        col = "black"
      )
    )
  })
  
  
  output$cluster_interpretation_box <- renderUI({
    req(projected_df())
    
    x <- projected_df() %>% arrange(period_order)
    
    tagList(
      lapply(seq_len(nrow(x)), function(i) {
        div(
          style = "
            margin-top: 20px;
            padding: 18px 22px;
            border: 1px solid #b8d4ff;
            border-radius: 12px;
            background-color: #f5f9ff;
            box-shadow: 0 2px 6px rgba(0,0,0,0.08);
            font-size: 16px;
            line-height: 1.6;
          ",
          h4(
            style = "margin-top: 0; color: #0b3d91; font-weight: 700;",
            paste0("Projected Clinical Cluster — ", x$period_label[i])
          ),
          div(
            style = "font-size: 18px; font-weight: 700; margin-bottom: 8px;",
            paste0(x$cluster_label[i], " — ", x$cluster_name[i])
          ),
          div(
            style = "font-size: 16px;",
            x$cluster_description[i]
          )
        )
      })
    )
  })
  
  output$ddrtree_plot <- renderPlot({
    df <- ref_df()
    
    validate(
      need("DDR_dim1" %in% colnames(df), "DDR_dim1 column not found"),
      need("DDR_dim2" %in% colnames(df), "DDR_dim2 column not found")
    )
    
    color_var <- input$color_var
    
    LEGEND_TITLE <- "CGM Glycemic Phenotype"
    LEGEND_TITLE_SIZE <- 16
    LEGEND_TEXT_SIZE <- 14
    
    p <- ggplot(df, aes(x = DDR_dim1, y = DDR_dim2))
    
    if (color_var == "clinical_cluster2") {
      
      p <- p +
        geom_point(aes(color = clinical_cluster2), alpha = 0.7, size = 1.8) +
        scale_color_discrete(
          name = "CGM Glycemic Phenotype",
          labels = c(
            "EHHR" = "Euglycemic high-hypoglycemia risk",
            "SE"   = "Stable Euglycemic",
            "MH"   = "Mild hyperglycemic",
            "HHV"  = "Hyperglycemia high Variability",
            "HLV"  = "Hyperglycemia low Variability"
          )
        )
      
    } else if (color_var == "Pseudotime") {
      
      p <- p +
        geom_point(aes(color = Pseudotime), alpha = 0.75, size = 1.8) +
        scale_color_viridis_c(
          name = "Pseudotime",
          option = "plasma",
          direction = -1
        )
    }
    
    if (input$run_btn > 0) {
      
      proj <- projected_df() %>%
        arrange(period_order)
      
      # trajectory line
      if (nrow(proj) >= 2) {
        p <- p +
          geom_path(
            data = proj,
            aes(
              x = DDR_dim1_projected,
              y = DDR_dim2_projected,
              group = 1
            ),
            inherit.aes = FALSE,
            color = "black",
            linewidth = 1.1,
            arrow = grid::arrow(length = grid::unit(0.25, "cm"))
          )
      }
      
      p <- p +
        geom_point(
          data = proj,
          aes(
            x = DDR_dim1_projected,
            y = DDR_dim2_projected
          ),
          inherit.aes = FALSE,
          shape = 21,
          size = 6,
          stroke = 1.4,
          color = "black",
          fill = "yellow"
        ) +
        geom_text(
          data = proj,
          aes(
            x = DDR_dim1_projected,
            y = DDR_dim2_projected,
            label = period_label
          ),
          inherit.aes = FALSE,
          hjust = -0.1,
          vjust = -0.7,
          size = 5,
          fontface = "bold"
        )
    }
    
    p +
      theme_bw() +
      coord_equal() +
      theme(
        legend.title = element_text(
          size = LEGEND_TITLE_SIZE,
          face = "bold"
        ),
        legend.text = element_text(
          size = LEGEND_TEXT_SIZE
        )
      ) +
      labs(
        title = paste("DDRTree Projection colored by CGM Glycemic Phenotype"),
        x = "DDR dim1",
        y = "DDR dim2",
        color = LEGEND_TITLE
        
      )
  })
}
shinyApp(ui = ui, server = server)
