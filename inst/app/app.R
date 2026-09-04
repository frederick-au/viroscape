library(shiny)
library(bslib)
library(ggplot2)
library(viroscape)

vs_metrics <- viroscape::distance_metrics()
vs_wb <- viroscape::wb_indicators()

`%||%` <- function(x, y) if (is.null(x)) y else x

## An empty selectInput reports "" rather than NULL, and the choices are
## refreshed only after the client answers, so the metric input can be blank
## or left over from a differently-typed dataset when divergence is first
## computed. Resolve it against the alignment instead of trusting the input.
resolve_metric <- function(metric, molecule) {
  allowed <- vs_metrics$metric[vs_metrics$molecule %in% c(molecule, "both")]
  if (length(metric) == 1L && !is.na(metric) && nzchar(metric) &&
      metric %in% allowed) return(metric)
  if (identical(molecule, "nucleotide")) "ml_nt" else "ml_aa"
}

dt <- function(x) {
  DT::datatable(as.data.frame(x), rownames = FALSE,
                options = list(dom = "tip", pageLength = 10, scrollX = TRUE)) |>
    DT::formatRound(columns = which(vapply(as.data.frame(x), is.numeric, logical(1))),
                    digits = 4)
}

ui <- bslib::page_navbar(
  title = "viroscape",
  theme = bslib::bs_theme(version = 5, primary = "#2c6e8f"),
  ## Not fillable. Every plot here has an explicit pixel height, so filling buys
  ## nothing, and on a crowded panel it costs: bslib shrinks fill-enabled
  ## children to share the viewport, the plot container reaches zero height, and
  ## the graphics device is asked for a zero-size canvas - "invalid quartz()
  ## size" on macOS, an equivalent error elsewhere. The panels scroll instead.
  fillable = FALSE,

  sidebar = bslib::sidebar(
    width = 290,
    title = "Data",

    radioButtons("source", NULL,
      choices = c("Bundled example" = "example",
                  "Local files" = "local",
                  "Live retrieval" = "live"),
      selected = "example"),

    conditionalPanel(
      "input.source == 'local'",
      fileInput("fasta", "Sequence FASTA", accept = c(".fasta", ".fa", ".fas", ".txt")),
      fileInput("envcsv", "Environmental CSV", accept = ".csv"),
      radioButtons("molecule", "Molecule", c("protein", "nucleotide"), inline = TRUE)
    ),

    conditionalPanel(
      "input.source == 'live'",
      selectizeInput("live_countries", "Countries",
        choices = c(names(viroscape::sea_countries()), "Vietnam", "China",
                    "Bangladesh", "India", "Egypt", "Nigeria", "Japan"),
        selected = names(viroscape::sea_countries()), multiple = TRUE),
      sliderInput("live_years", "Years", min = 1996, max = 2026,
                  value = c(2003, 2022), sep = ""),
      textInput("subtype", "HA subtype", "H5"),
      radioButtons("live_molecule", "Molecule", c("protein", "nucleotide"), inline = TRUE),
      numericInput("max_records", "Maximum records", 800, min = 50, max = 8000, step = 50),
      selectizeInput("wb_vars", "Environmental indicators",
        choices = stats::setNames(vs_wb$variable, vs_wb$label),
        selected = vs_wb$variable[vs_wb$default], multiple = TRUE,
        options = list(plugins = list("remove_button"))),
      helpText("Rows missing any selected indicator are dropped before modelling, so adding a sparsely reported indicator can cost you country-years."),
      checkboxInput("refresh", "Ignore cache", FALSE),
      helpText("Requires the rentrez package and network access. Results are cached on disk.")
    ),

    actionButton("load", "Load data", class = "btn-primary w-100"),
    hr(),
    uiOutput("filters"),
    hr(),
    div(class = "small text-muted", textOutput("status"))
  ),

  bslib::nav_panel(
    "Sampling",
    bslib::layout_columns(
      ## One card per row. Two sidebars already take ~550px, so a two-column
      ## split leaves each plot about 350px wide - narrow enough that axis
      ## labels are clipped and a twelve-country boxplot is unreadable. Full
      ## width and scroll beats side by side and squinting.
      col_widths = c(12, 12),
      bslib::card(bslib::card_header("Sampling effort"),
                  plotOutput("p_sampling", height = 460),
                  bslib::card_footer(class = "small text-muted",
                    "Uneven sampling is the main threat to interpreting country coefficients. A country represented by a single year, or a handful of sequences, cannot support a reliable estimate.")),
      bslib::card(bslib::card_header("By country"), DT::DTOutput("t_sampling"))
    ),
    bslib::card(bslib::card_header("Sequence metadata"), DT::DTOutput("t_meta"))
  ),

  bslib::nav_panel(
    "Divergence",
    bslib::layout_sidebar(
      sidebar = bslib::sidebar(
        width = 260, position = "right",
        selectInput("aligner", "Alignment backend", choices = NULL),
        selectInput("model_criterion", "Model selection criterion",
                    c("BIC", "AIC"), selected = "BIC"),
        checkboxInput("gamma", "Gamma rate heterogeneity", TRUE),
        numericInput("gamma_k", "Rate categories", 4, min = 1, max = 8),
        selectInput("metric", "Divergence metric",
                    choices = vs_metrics$metric, selected = "ml_aa"),
        actionButton("run_seq", "Recompute divergence", class = "btn-primary w-100"),
        hr(),
        checkboxInput("do_rate", "Estimate substitution rate (root-to-tip)", FALSE),
        helpText("Divergence from a consensus is not a rate. This roots a tree and regresses root-to-tip distance on collection year.")
      ),
      bslib::layout_columns(
        col_widths = c(12, 12),
        bslib::card(bslib::card_header("Substitution model ranking"),
                    plotOutput("p_models", height = 420)),
        bslib::card(bslib::card_header("Divergence from consensus"),
                    plotOutput("p_dist", height = 420))
      ),
      bslib::card(bslib::card_header("Model comparison"), DT::DTOutput("t_models")),
      bslib::card(bslib::card_header("Root-to-tip rate by country"),
                  DT::DTOutput("t_rate"),
                  bslib::card_footer(class = "small text-muted",
                    "Slope is substitutions per site per year. Tips are not independent, so these standard errors are optimistic."))
    )
  ),

  bslib::nav_panel(
    "Selection",
    bslib::layout_sidebar(
      sidebar = bslib::sidebar(
        width = 260, position = "right",
        radioButtons("level", "Model at",
          c("Country-year (recommended)" = "cluster", "Sequence" = "sequence"),
          selected = "cluster"),
        helpText("Environmental predictors vary by country-year. Selecting per sequence treats near-identical sequences as independent and retains predictors that explain nothing."),
        numericInput("perm", "Calibration permutations", 199, min = 0, max = 2000, step = 50),
        radioButtons("structure", "Country term",
          c("Choose by criterion" = "auto", "Fixed effect" = "fixed",
            "Random intercept" = "random", "Omit" = "none"), selected = "auto"),
        checkboxInput("use_year", "Include year covariate", TRUE),
        selectizeInput("candidates", "Candidate predictors", choices = NULL, multiple = TRUE),
        selectInput("criterion", "Criterion", c("AIC", "BIC")),
        numericInput("delta", "Improvement required (delta)", -2, max = 0, step = 0.5),
        numericInput("gvif", "Scaled GVIF threshold", 10, min = 1, step = 1),
        numericInput("stability", "Max proportional shift", 0.5, min = 0.05, max = 5, step = 0.05),
        checkboxInput("allow_base", "Exempt base terms from the collinearity rule", TRUE),
        actionButton("run_model", "Run forward selection", class = "btn-primary w-100")
      ),
      bslib::layout_columns(
        col_widths = c(4, 4, 4),
        bslib::value_box("Rows modelled", textOutput("v_rows")),
        bslib::value_box("Design effect", textOutput("v_deff")),
        bslib::value_box("Calibrated p", textOutput("v_calp"))
      ),
      bslib::layout_columns(
        col_widths = c(12, 12),
        bslib::card(bslib::card_header("Criterion path"), plotOutput("p_path", height = 420)),
        bslib::card(bslib::card_header("Structural comparison"), DT::DTOutput("t_structure"))
      ),
      bslib::card(bslib::card_header("Selection calibrated against permuted predictors"),
                  plotOutput("p_calib", height = 340),
                  bslib::card_footer(class = "small text-muted",
                    "Whole country-year predictor vectors are reassigned to other country-years, so the clustering and collinearity survive and only the association with divergence is destroyed. If the observed improvement sits inside this distribution, the search found nothing.")),
      bslib::card(bslib::card_header("Step log"),
                  DT::DTOutput("t_steps"),
                  bslib::card_footer(class = "small text-muted",
                    "A candidate is accepted only if it improves the criterion and leaves the existing coefficients stable. Rejections record which rule stopped them."))
    )
  ),

  bslib::nav_panel(
    "Diagnostics",
    bslib::layout_columns(
      col_widths = c(12, 12),
      bslib::card(bslib::card_header("Collinearity"), plotOutput("p_gvif", height = 420)),
      bslib::card(bslib::card_header("Sequential likelihood ratio tests"), DT::DTOutput("t_lrt"))
    ),
    bslib::card(bslib::card_header("Coefficient stability at the final step"),
                DT::DTOutput("t_stability"),
                bslib::card_footer(class = "small text-muted",
                  "Estimates before and after the last accepted predictor. Sign flips or large proportional shifts indicate the predictors are not separately identified."))
  ),

  bslib::nav_panel(
    "Results",
    bslib::layout_columns(
      col_widths = c(3, 3, 3, 3),
      bslib::value_box("Adjusted R-squared", textOutput("v_r2"), theme = "primary"),
      bslib::value_box("Observations", textOutput("v_n")),
      bslib::value_box("Predictors retained", textOutput("v_k")),
      bslib::value_box("Substitution model", textOutput("v_model"))
    ),
    bslib::layout_columns(
      col_widths = c(12, 12),
      bslib::card(bslib::card_header("Coefficients"), plotOutput("p_coef", height = 460)),
      bslib::card(bslib::card_header("Predictor relationship"),
                  selectInput("focus", NULL, choices = NULL),
                  plotOutput("p_focus", height = 400))
    ),
    bslib::card(
      bslib::card_header("Downloads"),
      div(class = "d-flex gap-2",
        downloadButton("dl_report", "Report (txt)", class = "btn-outline-primary"),
        downloadButton("dl_data", "Modelling dataset (csv)", class = "btn-outline-primary"),
        downloadButton("dl_steps", "Step log (csv)", class = "btn-outline-primary"))
    )
  ),

  bslib::nav_spacer(),
  bslib::nav_item(tags$a("Source", href = "https://github.com/frederickau/viroscape", target = "_blank"))
)

server <- function(input, output, session) {

  state <- reactiveValues(seqs = NULL, aln = NULL, msel = NULL, dists = NULL,
                          env = NULL, ds = NULL, mds = NULL, rate = NULL,
                          struct = NULL, sel = NULL, calib = NULL,
                          status = "Press Load data to begin.")

  observe({
    av <- viroscape::available_aligners()
    ok <- av$backend[av$available]
    updateSelectInput(session, "aligner", choices = c("auto", ok), selected = "auto")
  })

  output$status <- renderText(state$status)

  ## ---- load ----------------------------------------------------------------
  observeEvent(input$load, {
    withProgress(message = "Loading data", value = 0, {
      tryCatch({
        if (input$source == "example") {
          incProgress(0.3, detail = "bundled example")
          state$seqs <- viroscape::example_sequences()
          state$env <- viroscape::example_environment()
        } else if (input$source == "local") {
          req(input$fasta, input$envcsv)
          incProgress(0.3, detail = "reading files")
          state$seqs <- viroscape::read_sequences(input$fasta$datapath,
                                                  molecule = input$molecule)
          state$env <- viroscape::read_environment(input$envcsv$datapath)
        } else {
          req(input$live_countries)
          incProgress(0.2, detail = "querying NCBI")
          state$seqs <- viroscape::fetch_sequences(
            countries = input$live_countries, years = input$live_years,
            subtype = input$subtype, molecule = input$live_molecule,
            max_records = input$max_records, refresh = isTRUE(input$refresh))
          incProgress(0.4, detail = "querying the World Bank")
          state$env <- viroscape::fetch_environment(
            countries = input$live_countries, years = input$live_years,
            variables = if (length(input$wb_vars)) input$wb_vars else NULL,
            refresh = isTRUE(input$refresh))
        }
        incProgress(0.4, detail = "done")
        state$aln <- NULL; state$dists <- NULL; state$ds <- NULL; state$sel <- NULL
        state$status <- sprintf("%d sequences loaded.", length(state$seqs$seqs))
        mol <- state$seqs$molecule
        allowed <- vs_metrics$metric[vs_metrics$molecule %in% c(mol, "both")]
        updateSelectInput(session, "metric", choices = allowed,
                          selected = if (mol == "protein") "ml_aa" else "ml_nt")
      }, error = function(e) {
        state$status <- paste("Error:", conditionMessage(e))
        showNotification(conditionMessage(e), type = "error", duration = 10)
      })
    })
  })

  output$filters <- renderUI({
    req(state$seqs)
    meta <- state$seqs$meta
    ctry <- sort(unique(stats::na.omit(meta$country)))
    yr <- range(meta$year, na.rm = TRUE)
    hosts <- sort(unique(as.character(meta$host_class)))
    tagList(
      selectizeInput("f_countries", "Countries", choices = ctry, selected = ctry, multiple = TRUE),
      sliderInput("f_years", "Years", min = yr[1], max = yr[2], value = yr, sep = ""),
      selectizeInput("f_hosts", "Host classes", choices = hosts, selected = hosts, multiple = TRUE)
    )
  })

  filtered <- reactive({
    req(state$seqs)
    s <- state$seqs
    if (!is.null(input$f_countries)) {
      s <- try(viroscape::filter_sequences(s, countries = input$f_countries,
                                           years = input$f_years,
                                           hosts = input$f_hosts), silent = TRUE)
      if (inherits(s, "try-error")) return(NULL)
    }
    s
  })

  ## ---- divergence ----------------------------------------------------------
  recompute <- function() {
    s <- filtered(); req(s)
    withProgress(message = "Computing divergence", value = 0, {
      incProgress(0.15, detail = "aligning")
      state$aln <- viroscape::align_sequences(s, method = input$aligner %||% "auto")
      incProgress(0.35, detail = "ranking substitution models")
      state$msel <- viroscape::select_substitution_model(
        state$aln, criterion = input$model_criterion %||% "BIC",
        gamma = isTRUE(input$gamma), k = input$gamma_k %||% 4)
      incProgress(0.3, detail = "distances")
      cons <- viroscape::build_consensus(state$aln)
      mol <- state$aln$molecule
      metric <- resolve_metric(input$metric, mol)
      updateSelectInput(session, "metric",
                        choices = vs_metrics$metric[vs_metrics$molecule %in% c(mol, "both")],
                        selected = metric)
      state$dists <- viroscape::distance_from_consensus(
        state$aln, cons, metric = metric, model = state$msel)
      incProgress(0.15, detail = "joining environment")
      state$ds <- viroscape::join_environment(state$dists, state$env)
      state$rate <- NULL
      if (isTRUE(input$do_rate)) {
        state$rate <- try(viroscape::divergence_rate(state$aln, by = "country",
                                                     model = state$msel), silent = TRUE)
        if (inherits(state$rate, "try-error")) state$rate <- NULL
      }
      updateSelectizeInput(session, "candidates", choices = state$ds$predictors,
                           selected = state$ds$predictors)
      updateSelectInput(session, "focus", choices = state$ds$predictors)
      state$status <- sprintf("%d observations ready for modelling.", nrow(state$ds$data))
    })
  }

  observeEvent(input$load, {
    req(state$seqs, state$env)
    tryCatch(recompute(), error = function(e) {
      state$status <- paste("Error:", conditionMessage(e))
      showNotification(conditionMessage(e), type = "error", duration = 10)
    })
  }, priority = -10)

  observeEvent(input$run_seq, {
    tryCatch(recompute(), error = function(e) {
      showNotification(conditionMessage(e), type = "error", duration = 10)
    })
  })

  ## ---- modelling -----------------------------------------------------------
  observeEvent(input$run_model, {
    req(state$ds)
    withProgress(message = "Forward selection", value = 0.2, {
      tryCatch({
        ds <- if (identical(input$level, "sequence")) state$ds else
          viroscape::aggregate_to_cluster(state$ds)
        state$mds <- ds
        if (input$structure == "auto") {
          st <- viroscape::select_structure(ds, criterion = input$criterion,
                                            time = if (isTRUE(input$use_year)) "year" else NULL)
          state$struct <- st
          base <- st$base_terms
          rand <- st$random_term
        } else {
          state$struct <- NULL
          base <- switch(input$structure, fixed = "country", random = character(), none = character())
          if (isTRUE(input$use_year)) base <- c(base, "year")
          rand <- if (input$structure == "random") "(1|country)" else NULL
        }
        incProgress(0.4)
        state$calib <- NULL
        if (isTRUE(input$perm > 0)) {
          setProgress(message = "Calibrating against permuted predictors")
          state$calib <- try(viroscape::calibrate_selection(
            ds, base_terms = base, candidates = input$candidates,
            random_term = rand, criterion = input$criterion,
            B = input$perm), silent = TRUE)
          if (inherits(state$calib, "try-error")) state$calib <- NULL
        }
        incProgress(0.4)
        state$sel <- viroscape::forward_select(
          ds, base_terms = base, candidates = input$candidates,
          random_term = rand, criterion = input$criterion,
          delta_threshold = input$delta, gvif_threshold = input$gvif,
          stability_threshold = input$stability,
          allow_collinear_base = isTRUE(input$allow_base),
          calibrate = FALSE, calibration = state$calib, verbose = FALSE)
        state$status <- sprintf("Selected: %s",
          if (length(state$sel$selected)) paste(state$sel$selected, collapse = ", ") else "no predictors")
      }, error = function(e) {
        showNotification(conditionMessage(e), type = "error", duration = 10)
      })
    })
  })

  ## ---- outputs -------------------------------------------------------------
  output$p_sampling <- renderPlot({ req(state$ds); viroscape::plot_sampling(state$ds) })
  output$t_sampling <- DT::renderDT({
    req(state$ds); dt(viroscape::sampling_summary(state$ds)$by_country)
  })
  output$t_meta <- DT::renderDT({
    req(filtered())
    m <- filtered()$meta
    dt(m[c("id", "strain", "country", "year", "host", "host_class", "subtype")])
  })

  output$t_rate <- DT::renderDT({ req(state$rate); dt(as.data.frame(state$rate)) })
  output$v_rows <- renderText({
    req(state$mds); sprintf("%d", nrow(state$mds$data))
  })
  output$v_deff <- renderText({
    req(state$ds)
    de <- viroscape::cluster_summary(state$ds)$design_effect
    if (is.null(de) || !is.finite(de)) "-" else sprintf("%.1f", de)
  })
  output$v_calp <- renderText({
    if (is.null(state$calib)) "not run" else signif(state$calib$p_value, 3)
  })
  output$p_calib <- renderPlot({
    req(state$calib)
    ## a zero-width container means the panel has not been laid out yet; drawing
    ## into it is what produces the invalid-device error
    req(isTRUE(session$clientData$output_p_calib_width > 0))
    cal <- state$calib
    ggplot(data.frame(delta = cal$null), aes(x = delta)) +
      geom_histogram(bins = 40, fill = "grey80", colour = "white") +
      geom_vline(xintercept = cal$threshold_05, linetype = "dashed") +
      geom_vline(xintercept = cal$observed, colour = "#b2182b", linewidth = 1) +
      labs(x = sprintf("best first-round delta %s under permutation", cal$criterion),
           y = "permutations",
           subtitle = sprintf("observed %.1f (red), 5%% threshold %.1f (dashed), empirical p = %.3g",
                              cal$observed, cal$threshold_05, cal$p_value)) +
      theme_minimal(base_size = 12)
  })

  output$p_models <- renderPlot({ req(state$msel); viroscape::plot_model_selection(state$msel) })
  output$t_models <- DT::renderDT({ req(state$msel); dt(state$msel$table) })
  output$p_dist <- renderPlot({ req(state$ds); viroscape::plot_distances(state$ds) })

  output$p_path <- renderPlot({ req(state$sel); viroscape::plot_selection_path(state$sel) })
  output$t_structure <- DT::renderDT({ req(state$struct); dt(state$struct$table) })
  output$t_steps <- DT::renderDT({ req(state$sel); dt(state$sel$steps) })

  output$p_gvif <- renderPlot({ req(state$sel); viroscape::plot_gvif(state$sel, threshold = input$gvif) })
  output$t_lrt <- DT::renderDT({ req(state$sel); dt(viroscape::lrt_chain(state$sel)) })
  output$t_stability <- DT::renderDT({
    req(state$sel)
    ch <- state$sel$chain
    if (length(ch) < 2) return(dt(data.frame(note = "No predictor was added.")))
    dt(viroscape::coef_stability(ch[[length(ch) - 1]], ch[[length(ch)]],
                                 threshold = input$stability))
  })

  output$p_coef <- renderPlot({ req(state$sel); viroscape::plot_coefficients(state$sel) })
  output$p_focus <- renderPlot({
    req(state$ds, input$focus)
    viroscape::plot_predictor(state$ds, input$focus)
  })

  fit <- reactive({ req(state$sel); viroscape::model_fit_summary(state$sel$model) })
  output$v_r2 <- renderText(sprintf("%.3f", fit()$adj_r_squared))
  output$v_n <- renderText(format(fit()$n, big.mark = ","))
  output$v_k <- renderText(as.character(length(state$sel$selected)))
  output$v_model <- renderText({ req(state$msel); state$msel$best })

  analysis <- reactive({
    req(state$sel)
    structure(list(
      sequences = state$seqs, alignment = state$aln,
      model_selection = state$msel, consensus = NULL,
      distances = state$dists, environment = state$env,
      dataset = state$ds, model_dataset = state$mds,
      level = input$level, calibration = state$calib,
      structure = state$struct, selection = state$sel,
      lrt = viroscape::lrt_chain(state$sel),
      coefficients = viroscape::tidy_coefficients(state$sel),
      gvif = viroscape::gvif_table(state$sel$model),
      fit = fit(), sampling = viroscape::sampling_summary(state$ds)
    ), class = "vs_analysis")
  })

  output$dl_report <- downloadHandler(
    filename = function() sprintf("viroscape-report-%s.txt", Sys.Date()),
    content = function(file) {
      a <- analysis()
      if (is.null(a$structure)) a$structure <- list(table = data.frame())
      writeLines(viroscape::report_analysis(a), file)
    })
  output$dl_data <- downloadHandler(
    filename = function() sprintf("viroscape-dataset-%s.csv", Sys.Date()),
    content = function(file) utils::write.csv(as.data.frame(state$ds$data), file, row.names = FALSE))
  output$dl_steps <- downloadHandler(
    filename = function() sprintf("viroscape-steps-%s.csv", Sys.Date()),
    content = function(file) utils::write.csv(state$sel$steps, file, row.names = FALSE))
}

shinyApp(ui, server)
