server <- function(input, output, session) {

  # load config data (default and optionally custom ones)
  CONFIG <- gDRcomponents::get_config(dump = TRUE)

  #-##### spinner
  gDRcomponents::show_spinner()
  #-#####

  rv <- reactiveValues()
  sv <- list() # static values

  output$restart <- renderUI({
    a("RESTART", href = paste0("javascript:history.go(0)"))
  })

  # currently only shinyproxy deployments provide logout functionality
  output$logout <- renderUI({
    hostname <- session$clientData$url_hostname
    if (gDRcomponents::is_deployed_on(session, c("shinyproxy"))) {
      a("LOGOUT", href = paste0("https://", hostname, "/logout"))
    }
  })

  #-#### front page
  frontPageSERVER("frontpage", parent_session = session)

  #-####

  #-#### logs
  callModule(moduleDisplayLogs, "displayLogs", filepath = LOG_FILE)
  #-####

  #-##### admin mode
  isAdminMode <- isAdminUser(session)
  if (isAdminMode) {
    initAdminMode(input, CONFIG$plugins)
  }
  #-#####

  #-#### datasets
  # load projects and metaprojects data
  project_table <- get_project_table()

  metametaproject <- data.table::data.table(
    `project_id` = 0,
    `Project name` = "All projects",
    `#` = "?",
    `Description` = glue::as_glue("All projects"),
    `Projects` = ""
  )
  project_table$meta_project <-
    rbind(
      project_table$meta_project,
      metametaproject
    )
  # load filtered projects and metaprojects data
  project_table_filtered <-
    gDRcomponents::filter_project_table(project_table = project_table, metaprojects = NEW_TAB_FILTRATION)

  #-##### additional dataset (present/absent) with logic based on envs defined in global.R
  # TODO create a module that support multiple tabs/project-views 
  #   and/or upgrade while switching to elasticsearch-based solution
  observe({
    if (any(c(NEW_TAB_DESC, NEW_TAB_VALUE, NEW_TAB_FILTRATION) == "")) {
      shinyjs::addClass(selector = sprintf("#tabs > li > a[data-value='%s']", NEW_TAB_VALUE), class = "hidden")
    } else {
      shinyjs::removeClass(selector = sprintf("#tabs > li > a[data-value='%s']", NEW_TAB_VALUE), class = "hidden")
    }
  })
  #-#####

  #-##### example datasets
  # TODO: wrap this chunk into a helper function
  # TODO: switch to MAEs
  # list files that are available for viewing
  futile.logger::flog.trace("Main App: getting list of example files", name = "trace.logger")
  example_files <- USED_DATA
  previous_tab <- reactiveVal("")

  # For testing/dev purposes, you can define GDRVIZ_TEST_DATA envvar to quickly load a dataset
  test_dataset <- Sys.getenv("GDRVIZ_TEST_DATA")
  if (test_dataset != "") {
    if (test_dataset == "1") {
      test_dataset <- "medium"
    }
    shinyjs::delay(150, updateTabsetPanel(session, "tabs", selected = "example"))
    test_dataset_idx <- which(names(example_files) == test_dataset)
    if (length(test_dataset_idx) == 1) {
      example_files <- c(example_files[test_dataset_idx], example_files[-test_dataset_idx])
    }
  }
  #-#####

  #-##### global reactiveValues
  # Define reactiveValues to communicate across the modules
  mae_reactive <- reactiveValues(
      mae = NULL,
      meta_params = list(),
      out = NULL,
      reset = NULL,
      metaproject = NULL,
      project = NULL
  )

  # Define reactiveValues to communicate across the modules
  # - exps: experiments found in given mae/dataset
  # - sel_exps: experiments selected in given dataset (for visualisations/etc)
  .mae_filtering <- reactiveValues(
      exps = NULL,
      sel_exps = NULL
  )

  rective_simple_search_filters <- reactiveValues(
    mae_list = NULL,
    mae = NULL,
    mae_out = NULL,
    projects_DT = NULL,
    cell_lines = NULL,
    drugs = NULL
  )
  #-#####

  #-##### reset when dataset changed
  # check if data were selected in another module, if so - trigger reset in another one
  observe({
    if (!is.null(mae_reactive$mae)) {
      rective_simple_search_filters$mae <- NULL
    } else if (!is.null(rective_simple_search_filters$mae)) {
      mae_reactive$mae <- NULL
    }
  })
  #-#####

  #-##### too-big dataset
  # TODO: rewrite the logic for this part
  # the current logic for enabling/disabling some modules based on the dataset size
  # is overly complex and outdated (prepared for single SE, not valid for MAE with potentially three experiments)
  issue_with_mae <- reactiveValues(
    too_big_mae = FALSE
  )
  observeEvent(c(rective_simple_search_filters$error, mae_reactive$error), {
    issue_with_mae$too_big_mae <- if (input$tabs == "simple" && !is.null(rective_simple_search_filters$error)) {
      TRUE
    } else if (isTRUE(input$tabs %in% c("dash_mae", NEW_TAB_VALUE)) && !is.null(mae_reactive$error)) {
      TRUE
    } else {
      FALSE
    }

    updateTabItems(session, "main_tabs", "gDRsearch")
  })

  ######################
  #TODO: fixme - GDR-1642
  sv$reactive_filters < rective_simple_search_filters
  sv$search_simple_projects <- SEARCH_SIMPLE_PROJECTS
  sv$search_simple_celllines_cols <- SEARCH_SIMPLE_CELLLINES_COLS
  sv$search_simple_drug_cols <- SEARCH_SIMPLE_DRUG_COLS
  sv$example_files <- example_files
  sv$big_matrix_list <- mae_reactive
  sv$project_table <- project_table
  sv$project_table_filtered <- project_table_filtered
  sv$base_url <- BASE_URL
  sv$data_source <- DATA_SOURCE
  sv$log_file <- LOG_FILE
  sv$save_folder <- SAVE_FOLDER
  sv$useDSDB <- gDRcomponents::getBoolEnv("GDR_USE_DSDB", TRUE)
  callModule(gDRsearch::gDRsearch_SERVER, "search", rv, sv)
  #######################

  observeEvent(input$main_tabs, {
    if (input$main_tabs == "gDRin") {
      gDRutils::reset_env_identifiers()
    }
  }, ignoreInit = TRUE)

  observeEvent(input$`search_mae-project_DT_row_last_clicked`, {
    selected_project_id <- unlist(project_table$meta_project[input$`search_mae-meta_project_DT_row_last_clicked`]$Projects)[[ # nolint
      input$`search_mae-project_DT_row_last_clicked`]]
    combo_project_ids <- unlist(project_table_filtered$meta_project$Projects)
    if (selected_project_id %in% combo_project_ids) {
      shinyjs::hide(selector = '[data-value="heatmaps"]')
      shinyjs::hide(selector = '[data-value="boxes"]')
      shinyjs::hide(selector = '[data-value="bars"]')
      shinyjs::hide(selector = '[data-value="scatters"]')
      shinyjs::hide(selector = '[data-value="grid"]')
      shinyjs::hide(selector = '[data-value="curve"]')
    } else {
      shinyjs::hide(selector = '[data-value="combo"]')
    }
  })

  observeEvent(input$`search_mae2-project_DT_rows_selected`, {
    selected_project_id <- unlist(project_table_filtered$meta_project[input$`search_mae2-meta_project_DT_rows_selected`]$Projects)[[ # nolint
      input$`search_mae2-project_DT_rows_selected`]]
    combo_project_ids <- unlist(project_table_filtered$meta_project$Projects)
    if (selected_project_id %in% combo_project_ids) {
      shinyjs::hide(selector = '[data-value="heatmaps"]')
      shinyjs::hide(selector = '[data-value="boxes"]')
      shinyjs::hide(selector = '[data-value="bars"]')
      shinyjs::hide(selector = '[data-value="scatters"]')
      shinyjs::hide(selector = '[data-value="grid"]')
      shinyjs::hide(selector = '[data-value="curve"]')
    }
  })
  #-#####

  #-##### load dataset from any data source/(aka tab )

  observeEvent(input$main_tabs, {
    if (input$main_tabs == "frontpage") {
      shinyjs::hide(id = "miniSummary")
    } else {
      shinyjs::show(id = "miniSummary")
    }
    #TODO: we can add modal to tell user about possibility to lost all progress after switch input plugin
    if (input$main_tabs == "gDRsearch") {
      rv$mae_source <- "gDRsearch"
    } else if (input$main_tabs == "gDRin") {
      rv$mae_source <- "gDRin"
    }
  })

  # TODO: update the logic
  # as there are 'cotreatment' assays possible, next to 'single-agent' and 'matrix' assays
  .ds_raw <- reactive({
    s_dtype <- c("MultiAssayExperiment")
    .mae_filtering$sel_exps <- NULL

    req(rv$mae_source)
    my_ds <- if (rv$mae_source == "gDRsearch") {
      req(rv$mae_search())
    } else if (rv$mae_source == "gDRin") {
      rv$mae
    }
    req(my_ds)
    if (inherits(my_ds, "MultiAssayExperiment")) {
      my_ds
    } else {
      stop(sprintf("Unsupported dataset type. Currently supported: '%s'", toString(s_dtype)))
    }
  })
  #-#####

  #-##### filter experiments inside dataset
  filter_mae_mod <- gDRcomponents::filter_mae(
    id = "filter_mae",
    .exps = reactive(names(.ds_raw()))
  )

  # when filter_mae$trigger changes, update related values in mae_reactive
  observeEvent(filter_mae_mod$trigger, priority = -20, {

    req(.ds_raw)
    req(filter_mae_mod$trigger > 0)
    .mae_filtering$sel_exps  <- filter_mae_mod$sel_exps
  })



  #-##### load summary module
  # TODO: refactor for MAE dataset (show data for all experiments not only the first one)
  callModule(miniSummary, id = "miniSummary", ds = .ds)

  #TODO: fixme - there is duplicated values, it will be fix in GDR-1639
  rv$pidfs <- pidfs <- reactiveVal(gDRutils::get_prettified_identifiers(simplify = TRUE))

  #-##### validate and transform metadata of the dataset
  # TODO: get rid of 'drug_name' temporary fix (once dataset are updated)
  # TODO: wrap the code into functions/reactives
  # - validate MAE
  # - disambiguate cell-line/drugs data
  # - substitute 'drugname' to 'drug_name' column in metadataa
  # - check if 'too-big-mae' condition applies
  .ds <- reactive({
    futile.logger::flog.trace("Main App: \t disambiguating", name = "trace.logger")
    # why do we start with req for the first element? performnce-wise?
    req(any(vapply(.mae_filtering$sel_exps, function(x) {
      !is.null(x)
    }, logical(1))))
    req(.ds_raw()[[1]])
    sel_exps <- as.character(unlist(.mae_filtering$sel_exps))
    mae <- .ds_raw()[k = sel_exps]
    # checking if MAE identifiers are the same env
    mae_idfs <- gDRutils::get_MAE_identifiers(mae)
    gDRutils::update_env_idfs_from_mae(mae_idfs)
    pidfs(gDRutils::get_prettified_identifiers(simplify = TRUE))

    # validate MAE
    v <- tryCatch(gDRutils::validate_MAE(mae),
                  error = function(e) {
                    e
                  })
    # show popup if validation failed
    if (!is.null(unlist(unname(v)))) {
        futile.logger::flog.trace(sprintf("Main App: \t 'validate_mae' failed with '%s'", v$message),
                                  name = "trace.logger")
        # inform the user that given dataset in not gDRviz compatible and can't be used
        shinyalert::shinyalert(
            sprintf(
                "Selected dataset is not compatible with gDRviz. Please select another one."
            ),
            confirmButtonCol = "#3c8dbc"
        )
        req(FALSE)
    }


    for (i in seq_along(mae)) {
      SummarizedExperiment::colData(mae[[i]]) <-
        disambiguate_df_identifiers(SummarizedExperiment::colData(mae[[i]]), "cellline_name")
      SummarizedExperiment::rowData(mae[[i]]) <-
        disambiguate_df_identifiers(SummarizedExperiment::rowData(mae[[i]]), "drug_name")
      
    # renaming drugnameX identifiers to drug_nameX
    # this is a temporary fix until the identifier in the incoming data and testData is updated
      if (any(grepl("drugname", names(metadata(mae[[i]])$identifiers)))) {
        names(metadata(mae[[i]])$identifiers) <-
          gsub("drugname", "drug_name", names(metadata(mae[[i]])$identifiers))
      }
      
      # refine colData
      SummarizedExperiment::colData(mae[[i]]) <-
        gDRutils::refine_coldata(SummarizedExperiment::colData(mae[[i]]), mae[[i]])
      # refine rowData
      SummarizedExperiment::rowData(mae[[i]]) <-
        gDRutils::refine_rowdata(SummarizedExperiment::rowData(mae[[i]]), mae[[i]])
    }

    search_simple_max_cell_lines <- gDRcomponents::check_env_for_value(NULL, "search_simple_max_cell_lines")
    search_simple_max_drugs <- gDRcomponents::check_env_for_value(NULL, "search_simple_max_drugs")

    for (i in seq_along(mae)) {
      if (nrow(SummarizedExperiment::colData(mae[[i]])) < as.numeric(search_simple_max_cell_lines) /
          20 &&
          nrow(SummarizedExperiment::rowData(mae[[i]])) < as.numeric(search_simple_max_drugs) /
          20) {
        shinyjs::js$enableTab("shiny-tab-grid")
      } else {
        shinyjs::js$disableTab("shiny-tab-grid")
      }
    }

    issue_with_mae$too_big_mae  <- FALSE
    return(mae)
  })
  #-#####

  #-##### fit source
  rv$fit_source <- reactive({
    futile.logger::flog.trace("Main App: extracting fit sources", name = "trace.logger")
    data <- req(assay_metrics_initial())
    req(.mae_filtering$sel_exps)
    fit_sources <-
      lapply(data, function(x) {
        unique(na.omit(x$fit_source))
      })
    lapply(fit_sources, function(x) {
      checkmate::assert_true(length(x) > 0)
    })
    
    if (length(.mae_filtering$sel_exps) == 1L) {
      if ("single-agent" %in% names(.mae_filtering$sel_exps)) {
        list(`single-agent` = fit_sources[[.mae_filtering$sel_exps[["single-agent"]]]])
      } else if ("matrix" %in% names(.mae_filtering$sel_exps)) {
        list(matrix = fit_sources[[.mae_filtering$sel_exps[["matrix"]]]])
      } else {
        NULL
      }
    } else if (length(.mae_filtering$sel_exps > 1L)) {
      list(
        `single-agent` = fit_sources[[.mae_filtering$sel_exps[["single-agent"]]]],
        `matrix` = fit_sources[[.mae_filtering$sel_exps[["matrix"]]]]
      )
    }

  })

  #-##### transform dataset assays
  # -- for all assays (including combo-matrix):
  # - convert BumpyMatrice(s) to data.table(s)
  # -- for non-combo-matrix assays:
  # - dirty hack for ec50/c50 assay data
  # - flatten with normalization type/fit source
  # - prettify metrics
  # - merge assay-data with drug/cell_line data
  # - drop excessive columns (drug2-related if only untreated data found)
  # - run capVals (on metrics only)
  # TODO: make the logic valid for 'cotreatmnent' assay as well
  # TODO: split into helper function (see especialy repetitive code for non-matrix assays)
  # TODO: wrap the code into functions/reactives/shiny modules

  assay_metrics_initial <- reactive({
    futile.logger::flog.trace("Main App: extracting initial Metrics data", name = "trace.logger")
    ds <- req(.ds())
    shinybusy::show_modal_spinner(spin = "orbit", text = "Loading the project")
    gDRutils::MAEpply(
      ds,
      gDRcomponents::convert_se_assay_to_custom_dt,
      assay_name = "Metrics",
      output_table = "Metrics_initial"
    )
  })

  rv$metrics_wish_list <- reactive({
    ds <- req(.ds_raw())
    gDRcomponents::extract_additional_perturbations(ds)
  })

  assay_metrics_raw <- reactive({
    futile.logger::flog.trace("Main App: extracting Metrics data", name = "trace.logger")
    req(assay_metrics_initial())
    ds <- req(.ds())

    data <-
      gDRutils::MAEpply(
        ds,
        gDRcomponents::convert_se_assay_to_custom_dt,
        assay_name = "Metrics",
        output_table = "Metrics_raw"
      )
    gDRcomponents::remove_spinner()
    return(data)
  })

  assay_metrics <- reactive({
    futile.logger::flog.trace("Main App: \t capping values", name = "trace.logger")
    assay_metrics <- assay_metrics_raw()
    lapply(assay_metrics, capVals)
  })

  assay_normalized <- reactive({
    futile.logger::flog.trace("Main App: extracting Normalized data", name = "trace.logger")
    ds <- req(.ds())

    gDRutils::MAEpply(
      ds,
      gDRcomponents::convert_se_assay_to_custom_dt,
      assay_name = "Normalized"
    )
  })
  assay_averaged <- reactive({
    futile.logger::flog.trace("Main App: extracting Averaged data", name = "trace.logger")
    ds <- req(.ds())

    futile.logger::flog.trace("Main App: \t extracting", name = "trace.logger")

    gDRutils::MAEpply(
      ds,
      gDRcomponents::convert_se_assay_to_custom_dt,
      assay_name = "Averaged"
    )
  })

  combo_dts <- reactive({
    ds <- req(.ds())
    for (i in seq_along(ds)) {
      if (isTRUE(gDRcomponents::is_combo_data(ds[[i]]))) {
        combo_dt <- gDRutils::convert_combo_data_to_dt(ds[[i]])
      } else {
        combo_dt <- NULL
      }
    }
    combo_dt
  })

  # create object that stores all assays
  assay_objects <- reactiveVal()
  observe({
    req(.ds())
    futile.logger::flog.trace("Main App: building list of assays", name = "trace.logger")

    newObjects <- list()
    for (i in seq_along(.ds())) {
      newObjects[[i]] <- list(
        metrics_raw = assay_metrics_raw()[[i]],
        metrics = assay_metrics()[[i]],
        normalized = assay_normalized()[[i]],
        averaged = assay_averaged()[[i]]
      )
      if (!is.null(combo_dts()) && names(.ds())[i] %in% gDRutils::get_experiment_groups("matrix")) {
        futile.logger::flog.trace("Main App: appending combo data", name = "trace.logger")
        newObjects[[i]] <- c(newObjects[[i]], combo_dts())
      }
    }
    names(newObjects) <- names(.ds())

    assay_objects(newObjects)
  })

  rv$assay_object_manage_data <- reactive({
    req(assay_objects())
    req(.mae_filtering$sel_exps)
    rv$vars_wish_list <- c("drug_name2", "duration", "data_source", rv$metrics_wish_list())
    if (length(.mae_filtering$sel_exps) == 1L) {
      if ("single-agent" %in% names(.mae_filtering$sel_exps)) {
        list(`single-agent` = assay_objects()[[.mae_filtering$sel_exps[["single-agent"]]]])
      } else if ("matrix" %in% names(.mae_filtering$sel_exps)) {
        list(`matrix` = assay_objects()[[.mae_filtering$sel_exps[["matrix"]]]])
      } else {
        NULL
      }
    } else if (length(.mae_filtering$sel_exps > 1L)) {
      list(
        `single-agent` = assay_objects()[[.mae_filtering$sel_exps[["single-agent"]]]],
        `matrix` = assay_objects()[[.mae_filtering$sel_exps[["matrix"]]]]
      )
    }
  })

  #-##### disable grid module if too many drug vs cell-line combinations
  observeEvent(input$elementExist, {

    assay_metrics <- assay_objects()[["metrics"]]
    n_cell_lines <- length(unique(assay_metrics$clid))
    n_drugs <- length(unique(assay_metrics$Gnumber))

    search_simple_max_cell_lines <- as.numeric(gDRcomponents::check_env_for_value(NULL, "search_simple_max_cell_lines"))
    search_simple_max_drugs <- as.numeric(gDRcomponents::check_env_for_value(NULL, "search_simple_max_drugs"))


    if (n_cell_lines < search_simple_max_cell_lines / 20 &&
        n_drugs < search_simple_max_drugs / 20) {
      shinyjs::js$enableTab("shiny-tab-grid")
    } else {
      shinyjs::js$disableTab("shiny-tab-grid")
    }
    
    # handle tooltips below
    toltip_text <-
      sprintf(
        "Dose Response Overview is available for max %i drugs and %i cell lines",
        as.numeric(search_simple_max_drugs) / 20,
        as.numeric(search_simple_max_cell_lines) / 20
      )
    shinyjs::runjs(sprintf('$("#gDRvizTabs > ul > li:nth-child(6) > a").attr("title_dis", "%s");', toltip_text))

    shinyjs::runjs(sprintf('$("#gDRvizTabs > ul > li:nth-child(8) > a").attr("title_dis", "%s");',
                           "Drug Combo Treatments is available for combo data only"))
  })
  #-#####

  observeEvent(.mae_filtering$sel_exps, ignoreInit = FALSE, ignoreNULL = FALSE, {

    plugins_table <- data.table::rbindlist(CONFIG$plugins, fill = TRUE, idcol = TRUE)
    plugins_table$.id <- gDRcomponents:::plugin_make_id(parent_id = NULL, plugins_table$.id)

    combo_plugins <- plugins_table[plugins_table$type == "output - combo"]$.id
    sa_plugins <- plugins_table[plugins_table$type == "output - sa"]$.id

    if (is.null(.mae_filtering$sel_exps)) {
      # hiding/showing is analogous in every condition
      # collapse tabs
      shinyjs::hide("out_plugins")
      shinyjs::hide("manage_plugins")
      # hide tabs
      shinyjs::js$hideSidebarTab("shiny-tab-out_plugins")
      shinyjs::js$hideSidebarTab("shiny-tab-manage_plugins")
      # hide subtab for each single-agent and combo plugin
      lapply(combo_plugins, shinyjs::hide)
      lapply(sa_plugins, shinyjs::hide)
    } else if (!is.null(.mae_filtering$sel_exps[["matrix"]]) && !is.null(.mae_filtering$sel_exps[["single-agent"]])) {
      shinyjs::js$showSidebarTab("shiny-tab-out_plugins")
      shinyjs::show("out_plugins")
      shinyjs::js$showSidebarTab("shiny-tab-manage_plugins")
      shinyjs::show("manage_plugins")
      lapply(combo_plugins, shinyjs::show)
      lapply(sa_plugins, shinyjs::show)
    } else {
      shinyjs::js$showSidebarTab("shiny-tab-out_plugins")
      shinyjs::show("out_plugins")
      shinyjs::js$showSidebarTab("shiny-tab-manage_plugins")
      shinyjs::show("manage_plugins")
      if (!is.null(.mae_filtering$sel_exps[["matrix"]])) {
        lapply(combo_plugins, shinyjs::show)
        lapply(sa_plugins, shinyjs::hide)
      } else {
        lapply(combo_plugins, shinyjs::hide)
        lapply(sa_plugins, shinyjs::show)
      }
    }
  })

  rv$response_metrics <- reactive({
    req(assay_objects())
    req(rv$assay_object_manage_data())
    rv$assay_object_manage_data()[["single-agent"]][["metrics"]]
  })

  rv$response_data <- reactive({
    req(assay_objects())
    req(rv$assay_object_manage_data())
    rv$assay_object_manage_data()[["single-agent"]][["normalized"]]
  })

  observe({
    req(.mae_filtering$sel_exps[["matrix"]])
    req(.ds()[[.mae_filtering$sel_exps[["matrix"]]]])
    se_as_dt <- gDRutils::convert_combo_data_to_dt(.ds()[[.mae_filtering$sel_exps[["matrix"]]]])
    rv$combo_object <- se_as_dt
  })

  #-##### render visualisation tabs


  #-##### plugins
  gDRcomponents::init_plugins(
    plugins_config = CONFIG$plugins,
    session = session,
    rv = rv,
    sv = sv
  )

  #-##### spinner removal
  shinyjs::delay(1500, gDRcomponents::remove_spinner()) # remove it when done
  #-#####

}
