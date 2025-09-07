rm(list = ls())

packages <- c(
  "this.path",
  "comprehenr",
  "dplyr",
  "readr",
  "stringr",
  "tidyverse",
  "tidyr",
  "EnhancedVolcano",
  "ggplot2",
  "plotly",
  "org.Hs.eg.db",
  "clusterProfiler",
  "enrichR",
  "enrichplot",
  "shiny",
  "shinyBS",
  "bslib",
  "DT"
)

for(pkg in packages){
  library(pkg, character.only = TRUE)
}

source("R/functions.R")



options(shiny.maxRequestSize=160*1024^2)
ui <- fluidPage(
  uiOutput("ui_body"),
  tags$footer(
    div(
      style = "position: fixed; bottom: 0; width: 100%; background-color: #f8f9fa; padding: 10px; text-align: center; border-top: 1px solid #ddd;",
      "© 2025 Debojyoti Das Bioinformactics Unit, Core Facility & Clinical Genomics Linköping, Linköping University"
    )
  )
)


server <- function(input, output, session) {

  output$ui_body <- renderUI({
    navbarPage("Mass spectrometry data",
               # Page 1: Volcano plot and Data tables
               tabPanel("About", fluidPage(
                 page_sidebar(
                  includeMarkdown("README.md")
                 )
                 )
                 ),
               tabPanel("DE visualization", fluidPage(
                 page_sidebar(
                   "",
                   sidebar = sidebar(
                     bsCollapse(id = "sidebar_collapse_tab1", open = c("Data Uploading"), multiple = TRUE,
                                bsCollapsePanel("Data Uploading", style = "primary",
                                                div(style = "display: flex; align-items: center;",
                                                    fileInput("file", "Upload DE Table", accept = ".tsv")
                                                ),
                                                uiOutput("log2fc_col_ui"),
                                                uiOutput("pvalue_col_ui"),
                                                uiOutput("protein_name_col_ui"),
                                                uiOutput("gene_name_col_ui"),
                                                uiOutput("uniprotid_col_ui"),
                                                uiOutput("protein_description_col_ui"),
                                                uiOutput("log2FC_cutoff_ui"),
                                                uiOutput("pval_cutoff_ui"),
                                                uiOutput("process_data_ui")
                                ),
                                bsCollapsePanel("Volcano", style = "success",
                                                uiOutput("plot_volcano_ui"),
                                                uiOutput("volcano_title_ui"),
                                                uiOutput("genes_to_label_ui"),
                                                uiOutput("volcano_filetype_ui"),
                                                uiOutput("download_volcano_ui")

                                )
                     )
                   ),
                   tabsetPanel(
                     tabPanel("Data table: reformatted", DT::dataTableOutput("data_reformatted")),
                     tabPanel("Data table: raw", DT::dataTableOutput("data_uploaded")),
                     tabPanel("Data table: clean", DT::dataTableOutput("cleaned_data")),
                     tabPanel("Data table: DE genes", DT::dataTableOutput("de_data")),
                     tabPanel("Volcano Plot", plotOutput("volcanoPlot", height = "800px", width = "100%")),
                     tabPanel("Data Tables: duplicate genes", DT::dataTableOutput("duplicated_rows"))
                   )
                 )
               )
               ),
               # Page 2: Enrichment results
               tabPanel("GO over-representation analysis", fluidPage(
                 page_sidebar(
                   "",
                   sidebar = sidebar(
                     bsCollapse(id = "sidebar_collapse_tab2", open = c("Parameters"), multiple = TRUE,
                                bsCollapsePanel("Parameters", style = "primary",
                                                uiOutput("pval_threshold_ego_ui"),
                                                uiOutput("ont_ui"),
                                                uiOutput("gene_set_ego_ui"),
                                                uiOutput("pval_ego_cutoff_ui"),
                                                uiOutput("enrich_ui")
                                ),
                                bsCollapsePanel("Download", style = "success",
                                                uiOutput("download_options_ui")
                                )
                     )
                   ),
                   tabsetPanel(id = "go_plots",
                               tabPanel("Enrich GO results", DT::dataTableOutput("enrich_go_table")),
                               tabPanel("Heat Plot", plotOutput("heatmapPlot", height = "800px", width = "100%")),
                               tabPanel("Barplot", plotOutput("barplot", height = "800px", width = "100%")),
                               tabPanel("Dotplot", plotOutput("dotplot", height = "800px", width = "100%")),
                               tabPanel("Cnetplot", plotOutput("cnetplot", height = "800px", width = "100%")),
                               tabPanel("Treeplot", plotOutput("treeplot", height = "800px", width = "100%"))
                   )
                 )
               )
               ),
               # Page 3: Enrichment results
               tabPanel("KEGG pathway over-representation analysis", fluidPage(
                 page_sidebar(
                   "",
                   sidebar = sidebar(
                     bsCollapse(id = "sidebar_collapse_tab3", open = c("Parameters"), multiple = TRUE,
                                bsCollapsePanel("Parameters", style = "primary",
                                                uiOutput("pval_threshold_ekegg_ui"),
                                                uiOutput("gene_set_ekegg_ui"),
                                                uiOutput("pval_ekegg_cutoff_ui"),
                                                uiOutput("enrich_kegg_ui")
                                ),
                                bsCollapsePanel("Download", style = "success",
                                                uiOutput("download_kegg_options_ui")
                                )
                     )
                   ),
                   tabsetPanel(id = "kegg_plots",
                               tabPanel("Enrich KEGG results", DT::dataTableOutput("enrich_ekegg_table")),
                               tabPanel("Heat Plot", plotOutput("heatmapPlot_kegg", height = "800px", width = "100%")),
                               tabPanel("Barplot", plotOutput("barplot_kegg", height = "800px", width = "100%")),
                               tabPanel("Dotplot", plotOutput("dotplot_kegg", height = "800px", width = "100%")),
                               tabPanel("Cnetplot", plotOutput("cnetplot_kegg", height = "800px", width = "100%")),
                               tabPanel("Treeplot", plotOutput("treeplot_kegg", height = "800px", width = "100%"))
                   )
                 )
               )
               ),
               # Page 4: Enrichment results
               tabPanel("REACTOME pathway over-representation analysis", fluidPage(
                 page_sidebar(
                   "",
                   sidebar = sidebar(
                     bsCollapse(id = "sidebar_collapse_tab4", open = c("Parameters"), multiple = TRUE,
                                bsCollapsePanel("Parameters", style = "primary",
                                                uiOutput("pval_threshold_ereactome_ui"),
                                                uiOutput("gene_set_ereactome_ui"),
                                                uiOutput("pval_ereactome_cutoff_ui"),
                                                uiOutput("enrich_reactome_ui")
                                ),
                                bsCollapsePanel("Download", style = "success",
                                                uiOutput("download_reactome_options_ui")
                                )
                     )
                   ),
                   tabsetPanel(id = "reactome_plots",
                               tabPanel("Enrich REACTOME results", DT::dataTableOutput("enrich_ereactome_table")),
                               tabPanel("Heat Plot", plotOutput("heatmapPlot_reactome", height = "800px", width = "100%")),
                               tabPanel("Barplot", plotOutput("barplot_reactome", height = "800px", width = "100%")),
                               tabPanel("Dotplot", plotOutput("dotplot_reactome", height = "800px", width = "100%")),
                               tabPanel("Cnetplot", plotOutput("cnetplot_reactome", height = "800px", width = "100%")),
                               tabPanel("Treeplot", plotOutput("treeplot_reactome", height = "800px", width = "100%"))
                           )
                         )
                       )
                    ),
               # Page 5: Enrichment results
               tabPanel("Panther/ KEGG: enrichR over-representation analysis", fluidPage(
                 page_sidebar(
                   "",
                   sidebar = sidebar(
                     bsCollapse(id = "sidebar_collapse_tab5", open = c("Parameters"), multiple = TRUE,
                                bsCollapsePanel("Parameters", style = "primary",
                                                uiOutput("pval_threshold_enrichr_ui"),
                                                uiOutput("gene_set_enrichr_ui"),
                                                uiOutput("pval_enrichr_cutoff_ui"),
                                                uiOutput("database_enrichR_ui"),
                                                uiOutput("enrich_enrichr_ui")
                                ),
                                bsCollapsePanel("Download", style = "success",
                                                uiOutput("download_enrichr_options_ui")
                                )
                     )
                   ),
                   tabsetPanel(id = "enrichr_plots",
                               tabPanel("enrichR results", DT::dataTableOutput("enrich_enrichr_table")),
                               tabPanel("Barplot", plotOutput("barplot_enrichr", height = "800px", width = "100%")),
                               tabPanel("Dotplot", plotOutput("dotplot_enrichr", height = "800px", width = "100%"))
                          )
                        )
                      )
                    )
                  )
                })


  # Reactive expression to read the uploaded file
  raw_data_upload  <- reactive({
    req(input$file)
    df_raw <- read_delim(input$file$datapath) %>% as.data.frame()
    make_numeric_if_possible(df_raw)
  })


  column_types <- reactive({
    req(raw_data_upload())

    column_classes <- sapply(raw_data_upload(), class)
    numeric_columns <- names(column_classes[column_classes %in% c("numeric", "integer")])
    non_numeric_columns <- names(column_classes[!(column_classes %in% c("numeric", "integer"))])

    list(numeric_columns, non_numeric_columns)
  })
  # Generate UI for selecting log2FC and p-value columns
  output$log2fc_col_ui <- renderUI({
    req(column_types())
    selectInput("log2fc_col", "Select log2FC Column", choices = column_types()[[1]])
  })

  output$pvalue_col_ui <- renderUI({
    req(column_types())
    selectInput("pvalue_col", "Select p-value Column", choices = column_types()[[1]])
  })

  output$protein_name_col_ui <- renderUI({
    req(column_types())
    selectInput("protein_name_col", "Select protein column", choices = column_types()[[2]])
  })

  output$gene_name_col_ui <- renderUI({
    req(column_types())
    selectInput("gene_name_col", "Select gene column", choices = column_types()[[2]])
  })

  output$uniprotid_col_ui <- renderUI({
    req(column_types())
    selectInput("uniprotid_col", "Select uniprotid column", choices = column_types()[[2]])
  })

  output$protein_description_col_ui <- renderUI({
    req(column_types())
    selectInput("protein_description_col", "Select protein description column", choices = column_types()[[2]])
  })

  output$log2FC_cutoff_ui <- renderUI({
    req(raw_data_upload(), input$log2fc_col, input$pvalue_col)
    numericInput("log2FC_cutoff", "log2FC Cutoff", min = 0, value = 1.0)
  })

  output$pval_cutoff_ui <- renderUI({
    req(raw_data_upload(), input$log2fc_col, input$pvalue_col)
    numericInput("pval_cutoff", "P-value Cutoff", min = 0, value = 0.05)
  })

  output$process_data_ui <- renderUI({
    req(
        input$file,
        input$log2FC_cutoff,
        input$pval_cutoff,
        input$log2fc_col,
        input$pvalue_col,
        input$protein_name_col,
        input$gene_name_col,
        input$uniprotid_col,
        input$protein_description_col
        )
    div(style = "display: flex; align-items: center;",
        actionButton("process_data", "Process data",
                     style = "padding: 10px 10px; font-size: 16px; line-height: 1; border: 2px solid #ccc; border-radius: 4px; background-color: #f9f9f9; color: #333; margin-left: 3px;  margin-top: 3px; cursor: pointer;")
    )
  })

  output$plot_volcano_ui <- renderUI({
    req(data())
    div(style = "display: flex; align-items: center;",
        actionButton("plot_volcano", "Generate Volcano", style = "padding: 12px 12px; font-size: 16px; line-height: 1; border: 1px solid #ccc; border-radius: 4px; background-color: #f9f9f9; color: #333; margin-left: 3px; cursor: pointer;")
    )
  })

  output$volcano_title_ui <- renderUI({
    req(plot_volcano())
    textInput("volcano_title", "title for volcano plot", value = "")
  })

  output$genes_to_label_ui <- renderUI({
    req(plot_volcano(), clean_data())
    selectizeInput("genes_to_label", "Select Genes to Label", choices = clean_data()[["GeneName"]],
                   multiple = TRUE, selected = NULL,
                   options = list(placeholder = "Select from significant genes"))
  })

  output$volcano_filetype_ui <- renderUI({
    req(plot_volcano())
    selectInput("volcano_filetype", "Select file type for Volcano Plot", choices = c("pdf", "png", "jpeg"))
  })

  output$download_volcano_ui <- renderUI({
    req(plot_volcano())
    downloadButton("download_volcano", "Download Volcano Plot")})

  data <- reactive({
    req(input$file,
        input$log2FC_cutoff,
        input$pval_cutoff,
        input$log2fc_col,
        input$pvalue_col,
        input$protein_name_col,
        input$gene_name_col,
        input$uniprotid_col,
        input$protein_description_col,
        input$process_data
        )
    read_mass_spec_tsv(input$file$datapath,
                       log2FC_cutoff = input$log2FC_cutoff,
                       pval_cutoff = input$pval_cutoff,
                       logfc_col2use = input$log2fc_col,
                       pval_col2use = input$pvalue_col,
                       protein_name_col = input$protein_name_col,
                       gene_name_col = input$gene_name_col,
                       uniprotid_col = input$uniprotid_col,
                       protein_description_col = input$protein_description_col)
  })

  # Generate and display data tables
  output$data_reformatted <- DT::renderDataTable({
    req(data(), input$process_data)
    datatable(
      data(),
      options = list(
        pageLength = 10,
        autoWidth = TRUE,
        filter = 'top' # Enable column filters
      ),
      filter = 'top'
    )
  })

  output$data_uploaded <- DT::renderDataTable({
    req(raw_data_upload())
    datatable(
      raw_data_upload(),
      options = list(
        pageLength = 10,
        autoWidth = TRUE,
        filter = 'top' # Enable column filters
      ),
      filter = 'top'
    )
  })

  clean_data <- reactive({
    req(data(), input$log2fc_col, input$log2FC_cutoff)
    filter_duplicate_genes(data(),
                           x = "log2FC",
                           y = "Pvalue",
                           log2FC_cutoff = input$log2FC_cutoff)
  })

  output$cleaned_data <- DT::renderDataTable({
    req(clean_data(), input$process_data)
    datatable(clean_data() %>% as.data.frame(),
              filter = "top",
              extensions = 'Buttons',
              options = list(
                dom = 'Bflrtip',
                buttons = list('copy', 'csv', 'excel', list(
                  extend = "collection",
                  text = 'Show All',
                  action = DT::JS("function ( e, dt, node, config ) {
                                    dt.page.len(-1);
                                    dt.ajax.reload();
                                }")
                )
                )
              )
    )
  })

  duplicated_data <- reactive({
    req(data())
    duplicated_gene_rows(data(), x = "log2FC")
  })

  output$duplicated_rows <- DT::renderDataTable({
    req(duplicated_data(), input$process_data)
    datatable(
      duplicated_data() %>%
        as.data.frame(),
      filter = "top",
      extensions = 'Buttons',
      options = list(
        dom = 'Bflrtip',
        buttons = list('copy', 'csv', 'excel', list(
          extend = "collection",
          text = 'Show All',
          action = DT::JS("function ( e, dt, node, config ) {
                                    dt.page.len(-1);
                                    dt.ajax.reload();
                                }")
        )
        )
      )
    )
  })


  # Reactive to hold the filtered data
  significant_data <- eventReactive(input$process_data, {
    req(
      clean_data(),
      input$log2FC_cutoff,
      input$pval_cutoff
    )
    data_to_use <- clean_data()

    sig_data <- data_to_use %>%
      dplyr::filter(abs(log2FC) > input$log2FC_cutoff) %>%
      dplyr::filter(Pvalue < input$pval_cutoff)

    sig_data
  })

  output$de_data <- DT::renderDataTable({
    req(significant_data(), input$process_data)
    datatable(
      significant_data() %>%
        dplyr::select(all_of(c("log2FC", "Pvalue", "GeneName", "ProteinName")))  %>%
        dplyr::mutate(across(where(is.numeric), ~ round(., 3))),
      filter = "top",
      extensions = 'Buttons',
      options = list(
        dom = 'Bflrtip',
        buttons = list('copy', 'csv', 'excel', list(
          extend = "collection",
          text = 'Show All',
          action = DT::JS("function ( e, dt, node, config ) {
                                    dt.page.len(-1);
                                    dt.ajax.reload();
                                }")
        )
        )
      )
    )
  })


  # Update the choices for genes_to_label based on significant data
  observe({
    req(
      significant_data(),
      input$log2FC_cutoff,
      input$pval_cutoff
    )

    sig_data <- significant_data()

    if (nrow(sig_data) > 0) {
      updateSelectizeInput(
        session,
        "genes_to_label",
        choices = sig_data$GeneName,
        selected = NULL,
        server = TRUE
      )
    } else {
      updateSelectizeInput(
        session,
        "genes_to_label",
        choices = NULL,
        selected = NULL,
        server = TRUE
      )
    }
  })

  # Generate the volcano plot
  plot_volcano <- eventReactive(input$plot_volcano, {
    req(clean_data())
    # Modify volcano_plot to allow choosing genes significant for both p-value and log2FC
    volcano_plot(
      clean_data(),
      x = "log2FC",
      y = "Pvalue",
      log2FC_cutoff = input$log2FC_cutoff,
      pval_cutoff = input$pval_cutoff,
      label_genes = input$genes_to_label,
      title = input$volcano_title
    )  # You may add this argument for conditional labeling
  })

  # Code for generating volcano plot based on selected genes
  output$volcanoPlot <- renderPlot({
    req(plot_volcano())
    plot_volcano()
  })

  output$download_volcano <- downloadHandler(
    filename = function() {
      paste("volcano_plot.", input$volcano_filetype, sep = "")
    },
    content = function(file) {
      filetype <- input$volcano_filetype
      if (!is.null(plot_volcano())) {
        filetype <- input$heatmap_filetype
        ggsave(
          file,
          plot = plot_volcano(),
          device = filetype,
          width = 8,
          height = 7,
          dpi = 300
        )
      }
    }
  )

  # Perform enrichment analysis and generate plots
  output$pval_threshold_ego_ui <- renderUI({
    req(clean_data())
    numericInput("pval_threshold_ego", "P-value Threshold for Enrichment", min = 0, value = 0.05)
  })
  output$ont_ui <- renderUI({
    req(clean_data())
    selectInput("ont", "Ontology", choices = ont_list)
  })
  output$gene_set_ego_ui <- renderUI({
    req(clean_data())
    selectInput("gene_set_ego", "gene set", choices = c("all", "up", "down"))
  })
  output$enrich_ui <- renderUI({
    req(data(), input$log2fc_col)
    actionButton("enrich", "Perform Enrichment Analysis")
  })

  enrich_results <- eventReactive(input$enrich, {
    req(clean_data(), input$ont)
    enrichGO_pathway(
      data(),
      ont = input$ont,
      x = "log2FC",
      y_pval = "Pvalue",
      pval_cutoff = input$pval_cutoff,
      ego_pval_threshold = input$pval_threshold_ego,
      log2FC_cutoff = input$log2FC_cutoff,
      gene_set = input$gene_set_ego
    )
  })

  ###########
  # UI for download options (conditional on active tab)
  output$download_options_ui <- renderUI({
    req(input$go_plots, enrich_results())  # Ensure a tab is selected
    tagList(
      h4(paste("Download", input$go_plots, "Plot")),
      numericInput("download_plot_width", "Width (inches)", value = 12, min = 1),
      numericInput("download_plot_height", "Height (inches)", value = 8, min = 1),
      selectInput("download_plot_type", "File Type", choices = c("pdf", "png", "jpg")),
      downloadButton("download_plot", "Download Plot")
    )
  })

  output$enrich_go_table <- DT::renderDataTable({
      req(enrich_results())
      datatable(enrich_results()[[1]] %>% as.data.frame(),
                caption = "enrich GO results",
                filter = "top",
                extensions = 'Buttons',
                options = list(
                  dom = 'Bflrtip',
                  buttons = list('copy', 'csv', 'excel', list(
                    extend = "collection",
                    text = 'Show All',
                    action = DT::JS("function ( e, dt, node, config ) {
                                    dt.page.len(-1);
                                    dt.ajax.reload();
                                }")
                      )
                    )
                  )
        )
    })


  # Download handler
  output$download_plot <- downloadHandler(
    filename = function() {
      paste0("GO_plot_", input$go_plots, "_", Sys.Date(), ".", input$download_plot_type)
    },
    content = function(file) {
      plot_obj <- switch(input$go_plots,
                         "Heat Plot" = enrich_results()[[2]],
                         "Barplot" = enrich_results()[[3]],
                         "Dotplot" = enrich_results()[[4]],
                         "Cnetplot" = enrich_results()[[5]],
                         "Treeplot" = enrich_results()[[7]],
                         NULL)

      req(plot_obj)  # Ensure the plot exists

      ggsave(
        filename = file,
        plot = plot_obj,
        device = input$download_plot_type,
        width = input$download_plot_width,
        height = input$download_plot_height
      )
    }
  )
  ###########

  output$heatmapPlot <- renderPlot({
    req(enrich_results())
    if (is.null(enrich_results()[[2]])) {
      # Fallback plot with a message
      ggplot() +
        annotate(
          "text",
          x = 0.5,
          y = 0.5,
          label = "No significantly enriched terms found using enrichGO",
          size = 12,
          hjust = 0.5,
          vjust = 0.5
        ) +
        theme_void()
    } else {
      enrich_results()[[2]]
    }
  })

  output$barplot <- renderPlot({
    req(enrich_results())
    enrich_results()[[3]]
  })

  output$dotplot <- renderPlot({
    req(enrich_results())
    enrich_results()[[4]]
  })

  output$cnetplot <- renderPlot({
    req(enrich_results())
    enrich_results()[[5]]
  })

  output$treeplot <- renderPlot({
    req(enrich_results())
    enrich_results()[[7]]
  })

  # enrichKEGG
  # Perform enrichment analysis and generate plots
  output$pval_threshold_ekegg_ui <- renderUI({
    req(clean_data())
    numericInput("pval_threshold_ekegg", "P-value Threshold for Enrichment", min = 0, value = 0.05)
  })

    output$gene_set_ekegg_ui <- renderUI({
    req(clean_data())
    selectInput("gene_set_ekegg", "gene set", choices = c("all", "up", "down"))
  })
  output$enrich_kegg_ui <- renderUI({
    req(clean_data())
    actionButton("enrich_kegg", "Perform Enrichment Analysis")
  })

  enrich_results_kegg <- eventReactive(input$enrich_kegg, {
    req(clean_data(), input$log2fc_col)
    enrichKEGG_pathway(
      clean_data(),
      x = "log2FC",
      y_pval = "Pvalue",
      pval_cutoff = input$pval_cutoff,
      ekegg_pval_threshold = input$pval_threshold_ekegg,
      log2FC_cutoff = input$log2FC_cutoff,
      gene_set = input$gene_set_ekegg
    )
  })

  # UI for download options (conditional on active tab)
  output$download_kegg_options_ui <- renderUI({
    req(input$kegg_plots)  # Ensure a tab is selected
    tagList(
      h4(paste("Download", input$kegg_plots, "Plot")),
      numericInput("download_ekegg_plot_width", "Width (inches)", value = 12, min = 1),
      numericInput("download_ekegg_plot_height", "Height (inches)", value = 8, min = 1),
      selectInput("download_ekegg_plot_type", "File Type", choices = c("pdf", "png", "jpg")),
      downloadButton("download_ekegg_plot", "Download Plot")
    )
  })

  output$enrich_ekegg_table <- DT::renderDataTable({
    req(enrich_results_kegg())
    datatable(enrich_results_kegg()[[1]] %>% as.data.frame(),
              caption = "enrich KEGG results",
              filter = "top",
              extensions = 'Buttons',
              options = list(
                dom = 'Bflrtip',
                buttons = list('copy', 'csv', 'excel', list(
                  extend = "collection",
                  text = 'Show All',
                  action = DT::JS("function ( e, dt, node, config ) {
                                    dt.page.len(-1);
                                    dt.ajax.reload();
                                }")
                )
                )
              )
    )
  })


  # Download handler
  output$download_ekegg_plot <- downloadHandler(
    filename = function() {
      paste0("KEGG_plot_", input$kegg_plots, "_", Sys.Date(), ".", input$download_ekegg_plot_type)
    },
    content = function(file) {
      plot_obj_ekegg <- switch(input$kegg_plots,
                         "Heat Plot" = enrich_results_kegg()[[2]],
                         "Barplot" = enrich_results_kegg()[[3]],
                         "Dotplot" = enrich_results_kegg()[[4]],
                         "Cnetplot" = enrich_results_kegg()[[5]],
                         "Treeplot" = enrich_results_kegg()[[7]],
                         NULL)

      req(plot_obj_ekegg)  # Ensure the plot exists

      ggsave(
        filename = file,
        plot = plot_obj_ekegg,
        device = input$download_ekegg_plot_type,
        width = input$download_ekegg_plot_width,
        height = input$download_ekegg_plot_height
      )
    }
  )
  ###########
  output$heatmapPlot_kegg <- renderPlot({
    req(enrich_results_kegg())
    if (is.null(enrich_results_kegg()[[2]])) {
      # Fallback plot with a message
      ggplot() +
        annotate(
          "text",
          x = 0.5,
          y = 0.5,
          label = "No significantly enriched terms found using KEGG",
          size = 12,
          hjust = 0.5,
          vjust = 0.5
        ) +
        theme_void()
    } else {
      enrich_results_kegg()[[2]]
    }
  })

  output$barplot_kegg <- renderPlot({
    req(enrich_results_kegg())
    enrich_results_kegg()[[3]]
  })

  output$dotplot_kegg <- renderPlot({
    req(enrich_results_kegg())
    enrich_results_kegg()[[4]]
  })

  output$cnetplot_kegg <- renderPlot({
    req(enrich_results_kegg())
    enrich_results_kegg()[[5]]
  })

  output$treeplot_kegg <- renderPlot({
    req(enrich_results_kegg())
    enrich_results_kegg()[[7]]
  })

  # UI for download options (conditional on active tab)
  output$download_reactome_options_ui <- renderUI({
    req(input$reactome_plots)  # Ensure a tab is selected
    tagList(
      h4(paste("Download", input$reactome_plots, "Plot")),
      numericInput("download_ereactome_plot_width", "Width (inches)", value = 12, min = 1),
      numericInput("download_ereactome_plot_height", "Height (inches)", value = 8, min = 1),
      selectInput("download_ereactome_plot_type", "File Type", choices = c("pdf", "png", "jpg")),
      downloadButton("download_ereactome_plot", "Download Plot")
    )
  })

  output$enrich_ereactome_table <- DT::renderDataTable({
    req(enrich_results_reactome())
    datatable(enrich_results_reactome()[[1]] %>% as.data.frame(),
              caption = "enrich REACTOME results",
              filter = "top",
              extensions = 'Buttons',
              options = list(
                dom = 'Bflrtip',
                buttons = list('copy', 'csv', 'excel', list(
                  extend = "collection",
                  text = 'Show All',
                  action = DT::JS("function ( e, dt, node, config ) {
                                    dt.page.len(-1);
                                    dt.ajax.reload();
                                }")
                )
                )
              )
    )
  })


  # Download handler
  output$download_ereactome_plot <- downloadHandler(
    filename = function() {
      paste0("REACTOME_plot_", input$reactome_plots, "_", Sys.Date(), ".", input$download_ereactome_plot_type)
    },
    content = function(file) {
      plot_obj_ereactome <- switch(input$reactome_plots,
                                   "Heat Plot" = enrich_results_reactome()[[2]],
                                   "Barplot" = enrich_results_reactome()[[3]],
                                   "Dotplot" = enrich_results_reactome()[[4]],
                                   "Cnetplot" = enrich_results_reactome()[[5]],
                                   "Treeplot" = enrich_results_reactome()[[7]],
                                   NULL)

      req(plot_obj_ereactome)  # Ensure the plot exists

      ggsave(
        filename = file,
        plot = plot_obj_ereactome,
        device = input$download_ereactome_plot_type,
        width = input$download_ereactome_plot_width,
        height = input$download_ereactome_plot_height
      )
    }
  )



  # Reactome Pathway Analysis
  # Perform enrichment analysis and generate plots
  output$pval_threshold_ereactome_ui <- renderUI({
    req(clean_data())
    numericInput("pval_threshold_ereactome", "P-value Threshold for Enrichment", min = 0, value = 0.05)
  })

  output$gene_set_ereactome_ui <- renderUI({
    req(clean_data())
    selectInput("gene_set_ereactome", "gene set", choices = c("all", "up", "down"))
  })
  output$enrich_reactome_ui <- renderUI({
    req(clean_data())
    actionButton("enrich_reactome", "Perform Enrichment Analysis")
  })

  enrich_results_reactome <- eventReactive(input$enrich_reactome, {
    req(clean_data())
    enrichReactome_pathway(
      data(),
      x = "log2FC",
      y_pval = "Pvalue",
      pval_cutoff = input$pval_cutoff,
      ereactome_pval_threshold = input$pval_threshold_ereactome,
      log2FC_cutoff = input$log2FC_cutoff,
      gene_set = input$gene_set_ereactome
    )
  })

  output$heatmapPlot_reactome <- renderPlot({
    req(enrich_results_reactome())
    if (is.null(enrich_results_reactome()[[2]])) {
      # Fallback plot with a message
      ggplot() +
        annotate(
          "text",
          x = 0.5,
          y = 0.5,
          label = "No significantly enriched terms found using Reactome",
          size = 12,
          hjust = 0.5,
          vjust = 0.5
        ) +
        theme_void()
    } else {
      enrich_results_reactome()[[2]]
    }
  })

  output$barplot_reactome <- renderPlot({
    req(enrich_results_reactome())
    enrich_results_reactome()[[3]]
  })

  output$dotplot_reactome <- renderPlot({
    req(enrich_results_reactome())
    enrich_results_reactome()[[4]]
  })

  output$cnetplot_reactome <- renderPlot({
    req(enrich_results_reactome())
    enrich_results_reactome()[[5]]
  })

  output$treeplot_reactome <- renderPlot({
    req(enrich_results_reactome())
    enrich_results_reactome()[[7]]
  })




  # enrichPathway
  # Perform enrichment analysis and generate plots
  output$pval_threshold_enrichr_ui <- renderUI({
    req(clean_data())
    numericInput("pval_threshold_enrichR", "P-value Threshold for Enrichment", min = 0, value = 0.05)
  })

  output$gene_set_enrichr_ui <- renderUI({
    req(clean_data())
    selectInput("gene_set_enrichr", "gene set", choices = c("all", "up", "down"))
  })

  output$database_enrichR_ui <- renderUI({
    req(clean_data())
    selectInput("database_enrichR", "database", choices = c("Panther_2016", "Panther_2015", "KEGG_2021_Human", "KEGG_2019_Human"))
  })

  output$enrich_enrichr_ui <- renderUI({
    req(clean_data())
    actionButton("enrich_enrichR", "Perform Enrichment Analysis")
  })

  enrich_results_enrichR <- eventReactive(input$enrich_enrichR, {
    req(
      clean_data(),
      input$database_enrichR,
      input$pval_threshold_enrichR,
      input$log2FC_cutoff
    )
    enrichR_pathway(
      clean_data(),
      database = input$database_enrichR,
      x = "log2FC",
      y_pval = "Pvalue",
      pval_cutoff = input$pval_cutoff,
      enrichr_pval_threshold = input$pval_threshold_enrichr,
      log2FC_cutoff = input$log2FC_cutoff,
      gene_set = input$gene_set_enrichr
    )
  })


  output$barplot_enrichR <- renderPlot({
    req(enrich_results_enrichR())
    if (is.null(enrich_results_enrichR()[[2]])) {
      # Fallback plot with a message
      ggplot() +
        annotate(
          "text",
          x = 0.5,
          y = 0.5,
          label = "No significantly enriched terms found using enrichR",
          size = 12,
          hjust = 0.5,
          vjust = 0.5
        ) +
        theme_void()
    } else {
      enrich_results_enrichR()[[2]]
    }
  })

  output$dotplot_enrichR <- renderPlot({
    req(enrich_results_enrichR())
    enrich_results_enrichR()[[3]]
  })


  # UI for download options (conditional on active tab)
  output$download_enrichr_options_ui <- renderUI({
    req(input$enrichr_plots)  # Ensure a tab is selected
    tagList(
      h4(paste("Download", input$enrichr_plots, "Plot")),
      numericInput("download_enrichr_plot_width", "Width (inches)", value = 12, min = 1),
      numericInput("download_enrichr_plot_height", "Height (inches)", value = 8, min = 1),
      selectInput("download_enrichr_plot_type", "File Type", choices = c("pdf", "png", "jpg")),
      downloadButton("download_enrichr_plot", "Download Plot")
    )
  })

  output$enrich_enrichr_table <- DT::renderDataTable({
    req(enrich_results_enrichR())
    datatable(enrich_results_enrichR()[[1]] %>% as.data.frame(),
              caption = "enrichR results",
              filter = "top",
              extensions = 'Buttons',
              options = list(
                dom = 'Bflrtip',
                buttons = list('copy', 'csv', 'excel', list(
                  extend = "collection",
                  text = 'Show All',
                  action = DT::JS("function ( e, dt, node, config ) {
                                    dt.page.len(-1);
                                    dt.ajax.reload();
                                }")
                )
                )
              )
    )
  })


  # Download handler
  output$download_enrichr_plot <- downloadHandler(
    filename = function() {
      paste0("enrichr_plot_", input$enrichr_plots, "_", Sys.Date(), ".", input$download_enrichr_plot_type)
    },
    content = function(file) {
      plot_obj_enrichr <- switch(input$enrichr_plots,
                                   "Heat Plot" = enrich_results_enrichR()[[2]],
                                   "Barplot" = enrich_results_enrichR()[[3]],
                                   "Dotplot" = enrich_results_enrichR()[[4]],
                                   "Cnetplot" = enrich_results_enrichR()[[5]],
                                   "Treeplot" = enrich_results_enrichR()[[7]],
                                   NULL)

      req(plot_obj_enrichr)  # Ensure the plot exists

      ggsave(
        filename = file,
        plot = plot_obj_enrichr,
        device = input$download_enrichr_plot_type,
        width = input$download_enrichr_plot_width,
        height = input$download_enrichr_plot_height
      )
    }
  )

}


shiny::shinyApp(ui = ui, server = server)


