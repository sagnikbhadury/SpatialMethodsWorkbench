library(shiny)
library(bslib)
library(ggplot2)
library(DT)

invisible(lapply(list.files("R", pattern = "\\.R$", full.names = TRUE), source, local = FALSE))

app_theme <- bs_theme(version = 5, bg = "#f6f3ec", fg = "#112d35", primary = "#176b68",
                      secondary = "#d46245", base_font = font_collection("Inter", "system-ui"),
                      heading_font = font_collection("Georgia", "serif"))

ui <- fluidPage(
  theme = app_theme,
  tags$head(
    tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
    tags$link(rel = "stylesheet", href = "styles.css"),
    tags$title("Spatial Methods Workbench")
  ),
  div(class = "app-shell",
    tags$header(class = "app-header",
      div(class = "brand", span(class = "brand-mark", "SB"), div(strong("Spatial Methods Workbench"), tags$small("Public-method analysis portal"))),
      div(class = "header-links", a(href = "https://sagnikbhadury.github.io", "Research site ↗"), a(href = "https://github.com/sagnikbhadury", "GitHub ↗"))
    ),
    tags$section(class = "app-intro",
      div(p(class = "kicker", "UPLOAD · VALIDATE · CHOOSE · RUN"), h1("One careful path through spatial and structured data."),
          p("Upload a table, map its columns, and run only analyses whose input requirements are satisfied. Results stay tied to their assumptions and reproducible settings.")),
      tags$aside(strong("Research use only"), p("This application does not provide clinical decisions. Uploaded data are held only in the active Shiny session; do not upload protected health information to a public deployment."))
    ),
    div(class = "workbench",
      tags$aside(class = "control-panel",
        div(class = "step-label", "01 · DATA"),
        checkboxInput("use_demo", "Use synthetic spatial demonstration", TRUE),
        fileInput("file", "Or upload CSV", accept = c(".csv", "text/csv"), buttonLabel = "Choose CSV"),
        helpText("One row per observation. Coordinates, annotations, outcomes, and numeric marker/feature columns may be combined."),
        hr(),
        div(class = "step-label", "02 · COLUMN MAP"),
        selectInput("x_col", "X coordinate", choices = character()),
        selectInput("y_col", "Y coordinate", choices = character()),
        selectInput("z_col", "Z / section coordinate", choices = character()),
        selectInput("zone_col", "Region / zone", choices = character()),
        selectInput("id_col", "Sample / contour / volume ID", choices = character()),
        selectInput("outcome_col", "Outcome", choices = character()),
        selectInput("exposure_col", "Exposure", choices = character()),
        selectInput("mediator_col", "Mediator", choices = character()),
        selectizeInput("feature_cols", "Numeric markers / features", choices = character(), multiple = TRUE,
                       options = list(plugins = list("remove_button"), placeholder = "Select at least one feature")),
        hr(),
        div(class = "step-label", "03 · ANALYSIS"),
        selectInput("analysis", "Method", choices = character()),
        uiOutput("method_status"),
        uiOutput("method_params"),
        checkboxInput("citation_ack", "I will cite the workbench and applicable method publications in outputs that use these results.", value = FALSE),
        p(class = "citation-note", "Citation metadata is written into every downloaded result bundle."),
        actionButton("run", "Run selected analysis", class = "run-button", width = "100%"),
        downloadButton("download_bundle", "Download reproducibility bundle", class = "download-button")
      ),
      tags$main(class = "result-panel",
        uiOutput("dataset_summary"),
        navset_card_tab(
          nav_panel("Readiness", uiOutput("readiness")),
          nav_panel("Figure", div(class = "plot-wrap", plotOutput("result_plot", height = "610px"))),
          nav_panel("Table", DTOutput("result_table")),
          nav_panel("Discuss results", div(class = "collaboration-panel",
            h3("Interpretation requires scientific context."),
            p("This application does not explain, narrate, or draw scientific conclusions from results."),
            p("For study-specific interpretation, methodological guidance, or collaboration, contact Sagnik Bhadury."),
            a(class = "btn btn-primary", href = "mailto:bhadury@umich.edu?subject=Spatial%20Methods%20Workbench%20collaboration", "Email about collaboration"),
            a(class = "collaboration-link", href = "https://sagnikbhadury.github.io/work-with-me/", target = "_blank", rel = "noopener", "Collaboration information ↗")
          )),
          nav_panel("Reproducibility", verbatimTextOutput("manifest")),
          id = "results_tabs"
        )
      )
    ),
    tags$footer(class = "app-footer",
      p("Static guidance and fast exploratory modules run in-session. Intensive Bayesian engines appear only when installed on the deployment."),
      p(a(href = "https://doi.org/10.1038/s41598-026-35341-8", "ISPat publication ↗"), " · ", a(href = "https://github.com/sagnikbhadury/SpatialMethodsWorkbench", "Application source ↗"))
    )
  )
)

server <- function(input, output, session) {
  uploaded <- reactive({
    if (isTRUE(input$use_demo) || is.null(input$file)) demo_spatial_data() else read_uploaded_csv(input$file$datapath)
  })

  observeEvent(input$file, { if (!is.null(input$file)) updateCheckboxInput(session, "use_demo", value = FALSE) })

  observeEvent(uploaded(), {
    data <- uploaded(); all_cols <- names(data); nums <- numeric_columns(data); cats <- categorical_columns(data)
    none <- c("Not selected" = "")
    guess <- function(pattern, candidates = all_cols) { hit <- grep(pattern, candidates, ignore.case = TRUE, value = TRUE); if (length(hit)) hit[1] else "" }
    updateSelectInput(session, "x_col", choices = c(none, nums), selected = guess("^(x|coord_x|x_coord)$", nums))
    updateSelectInput(session, "y_col", choices = c(none, nums), selected = guess("^(y|coord_y|y_coord)$", nums))
    updateSelectInput(session, "z_col", choices = c(none, nums), selected = guess("^(z|section|depth|coord_z)$", nums))
    updateSelectInput(session, "zone_col", choices = c(none, all_cols), selected = guess("zone|region|cluster", all_cols))
    updateSelectInput(session, "id_col", choices = c(none, all_cols), selected = guess("sample.*id|patient.*id|volume.*id|^id$", all_cols))
    updateSelectInput(session, "outcome_col", choices = c(none, all_cols), selected = guess("outcome|response|disease|class|label", all_cols))
    updateSelectInput(session, "exposure_col", choices = c(none, all_cols), selected = guess("exposure|treatment|group", all_cols))
    updateSelectInput(session, "mediator_col", choices = c(none, nums), selected = guess("mediator", nums))
    excluded <- unique(c(guess("^(x|coord_x|x_coord)$", nums), guess("^(y|coord_y|y_coord)$", nums), guess("^(z|section|depth|coord_z)$", nums)))
    defaults <- head(setdiff(nums, excluded), 8)
    updateSelectizeInput(session, "feature_cols", choices = nums, selected = defaults, server = TRUE)
  }, ignoreInit = FALSE)

  mapping <- reactive(coerce_mapping(uploaded(), list(x = input$x_col, y = input$y_col, z = input$z_col,
    zone = input$zone_col, id = input$id_col, outcome = input$outcome_col,
    exposure = input$exposure_col, mediator = input$mediator_col)))
  features <- reactive(intersect(input$feature_cols %||% character(), numeric_columns(uploaded())))
  validation <- reactive(validate_dataset(uploaded(), mapping()))
  recommendations <- reactive(recommend_methods(uploaded(), mapping(), features()))

  observe({
    registry <- analysis_registry(); rec <- recommendations()
    labels <- vapply(names(registry), function(id) {
      status <- if (rec[[id]]$compatible) "Ready" else "Needs input/engine"
      paste0(registry[[id]]$label, " — ", status)
    }, character(1))
    choices <- stats::setNames(names(registry), labels)
    selected <- input$analysis
    if (is.null(selected) || !selected %in% names(registry)) selected <- names(registry)[which(vapply(rec, `[[`, logical(1), "compatible"))[1]] %||% names(registry)[1]
    updateSelectInput(session, "analysis", choices = choices, selected = selected)
  })

  output$dataset_summary <- renderUI({
    v <- validation()
    div(class = "metric-row",
      div(class = "metric", strong(format(v$n, big.mark = ",")), span("observations")),
      div(class = "metric", strong(v$p), span("columns")),
      div(class = "metric", strong(length(v$numeric)), span("numeric")),
      div(class = "metric", strong(sprintf("%.1f%%", 100 * v$missing_rate)), span("missing"))
    )
  })

  output$method_status <- renderUI({
    req(input$analysis); method <- analysis_registry()[[input$analysis]]; status <- recommendations()[[input$analysis]]
    div(class = paste("method-note", if (status$compatible) "method-note--ready" else "method-note--blocked"),
        strong(if (status$compatible) "Ready to run" else "Not ready"),
        p(method$short), tags$small(paste(method$runtime, "runtime ·", method$engine)),
        if (!status$compatible) tags$small(status$reason))
  })

  output$method_params <- renderUI({
    req(input$analysis)
    switch(input$analysis,
      spatial_qc = tagList(selectInput("param_feature", "Feature to map", choices = features()), sliderInput("param_neighbors", "Nearest neighbors", 3, 20, 8), sliderInput("param_permutations", "Permutations", 99, 999, 199, step = 100)),
      conditional_network = tagList(sliderInput("param_threshold", "|Partial correlation| threshold", 0, .8, .15, step = .05), sliderInput("param_ridge", "Ridge regularization", .01, .5, .05, step = .01), checkboxInput("param_by_zone", "Estimate one network per selected zone", !is.null(mapping()$zone))),
      ispat = tagList(selectInput("param_kernel", "Spatial kernel", c("Matern", "RBF")), selectInput("param_spatial_fit", "Spatial fit", c("VB", "MLE")), selectInput("param_factor_fit", "Factor fit", c("CAVI", "SVI")), sliderInput("param_threshold", "Display threshold", 0, .8, .15, .05), numericInput("param_cores", "Cores", 1, 1, 8)),
      gpghs = tagList(sliderInput("param_basis", "Basis functions per dimension", 2, 8, 4), selectInput("param_nu", "Matérn smoothness", c(.5, 1.5, 2.5), 1.5), numericInput("param_nmc", "MCMC iterations", 1200, 500, 10000, 100), numericInput("param_burn", "Burn-in", 400, 100, 5000, 100), selectInput("param_symmetry", "Symmetry rule", c("AND", "OR")), numericInput("param_cores", "Cores", 1, 1, 8)),
      ispat3d = tagList(selectInput("param_kernel", "Volumetric kernel", c("Matern", "RBF")), selectInput("param_factor_fit", "Factor fit", c("CAVI", "SVI")), sliderInput("param_threshold", "Display threshold", 0, .8, .15, .05), numericInput("param_cores", "Cores", 1, 1, 8)),
      spatial_ml = numericInput("param_seed", "Random seed", 2026, 1, 999999),
      mediation = tagList(sliderInput("param_bootstrap", "Bootstrap replicates", 100, 1000, 300, 100), numericInput("param_seed", "Random seed", 2026, 1, 999999)),
      shape_pca = sliderInput("param_landmarks", "Resampled landmarks", 20, 100, 40, 5)
    )
  })

  params <- reactive({
    list(feature = input$param_feature, neighbors = input$param_neighbors, permutations = input$param_permutations,
         threshold = input$param_threshold, ridge = input$param_ridge, by_zone = input$param_by_zone,
         kernel = input$param_kernel, spatial_fit = input$param_spatial_fit, factor_fit = input$param_factor_fit,
         cores = input$param_cores, basis = input$param_basis, nu = as.numeric(input$param_nu),
         nmc = input$param_nmc, burn = input$param_burn, thin = 4L, symmetry = input$param_symmetry,
         seed = input$param_seed, bootstrap = input$param_bootstrap, landmarks = input$param_landmarks)
  })

  output$readiness <- renderUI({
    v <- validation(); registry <- analysis_registry(); rec <- recommendations()
    div(class = "readiness",
      if (!v$ok) div(class = "alert alert-danger", lapply(v$issues, tags$p)) else div(class = "alert alert-success", "The file passed structural validation."),
      if (length(v$warnings)) div(class = "alert alert-warning", lapply(v$warnings, tags$p)),
      h3("Compatible analysis map"),
      div(class = "readiness-grid", lapply(names(registry), function(id) {
        item <- registry[[id]]; status <- rec[[id]]
        tags$article(class = paste("readiness-card", if (status$compatible) "is-ready" else "is-blocked"),
                span(class = "status-dot"), h4(item$label), p(item$short), tags$small(if (status$compatible) paste("Ready ·", item$runtime) else status$reason))
      }))
    )
  })

  result <- eventReactive(input$run, {
    req(validation()$ok, input$analysis)
    validate(need(isTRUE(input$citation_ack), "Please acknowledge the citation requirement before running an analysis."))
    status <- recommendations()[[input$analysis]]
    validate(need(status$compatible, status$reason))
    withProgress(message = "Running analysis", value = 0, {
      incProgress(.15, detail = "Preparing validated input")
      tryCatch({
        out <- run_analysis(input$analysis, uploaded(), mapping(), features(), params())
        incProgress(.85, detail = "Preparing results")
        out
      }, error = function(e) list(error = conditionMessage(e), method = analysis_registry()[[input$analysis]]$label))
    })
  }, ignoreInit = TRUE)

  output$result_plot <- renderPlot({ req(result()); validate(need(is.null(result()$error), result()$error), need(!is.null(result()$plot), "This result has no figure.")); result()$plot })
  output$result_table <- renderDT({ req(result()); validate(need(is.null(result()$error), result()$error)); datatable(result()$table, filter = "top", options = list(pageLength = 12, scrollX = TRUE)) })
  output$manifest <- renderText({ req(result()); if (!is.null(result()$error)) return(result()$error); jsonlite::toJSON(result_manifest(result()), auto_unbox = TRUE, pretty = TRUE, null = "null") })
  output$download_bundle <- downloadHandler(filename = function() paste0("spatial-workbench-", Sys.Date(), ".zip"),
    content = function(file) { req(result()); validate(need(is.null(result()$error), result()$error)); write_result_bundle(result(), file) }, contentType = "application/zip")
}

shinyApp(ui, server)
