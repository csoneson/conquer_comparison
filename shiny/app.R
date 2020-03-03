suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(pheatmap))

trsf <- function(x) {
  if (x == "none") "."
  else x
}

topdir <- ".."

## Load summarized data
res <- readRDS(paste0(topdir, "/export_results/shiny_results.rds"))
aspects <- c("fracna", "nbrgenes", "type1error", "fdp", "tpr", "auroc", "timing", "origvsmock")
aspects2 <- c("decharac")

## Define app
summary_app <- function(res, aspects, aspects2) {

  p_layout <- function(request) {
    shinydashboard::dashboardPage(
      skin = "purple", 
      
      shinydashboard::dashboardHeader(title = "Exploration of results", 
                                      titleWidth = 300),
      
      shinydashboard::dashboardSidebar(
        selectInput(inputId = "keep_methods", label = "Methods to include",
                    choices = res$allmethods,
                    selected = res$allmethods,
                    multiple = TRUE, selectize = TRUE),
        
        selectInput(inputId = "highlight_method", label = "Method to highlight",
                    choices = c("none", res$allmethods),
                    selected = "none",
                    multiple = FALSE, selectize = TRUE)
      ),
      
      shinydashboard::dashboardBody(fluidRow(
        do.call(shinydashboard::tabBox,
          c(
            width = 12, 
            
            list(shiny::tabPanel("About",
                                 includeMarkdown("about.md"),
                                 value = "about")),
            
            list(shiny::tabPanel(
              "Select inputs",
              fluidRow(
                column(width = 4, 
                       selectInput(inputId = "keep_ncells",
                                   label = "Keep instances with the following number of cells per group",
                                   choices = res$allncells,
                                   selected = res$allncells, 
                                   multiple = TRUE, selectize = TRUE)),
                column(width = 4,
                       selectInput(inputId = "keep_datasets",
                                   label = "Keep instances from the following data sets",
                                   choices = res$alldatasets,
                                   selected = res$alldatasets,
                                   multiple = TRUE, selectize = TRUE)),
                column(width = 4,
                       selectInput(inputId = "keep_filterings",
                                   label = "Keep instances with the following filter settings",
                                   choices = res$allfilterings,
                                   selected = res$allfilterings,
                                   multiple = TRUE, selectize = TRUE))
              )
            )), 

            list(shiny::tabPanel(
              "Performance summary",
              fluidRow(
                column(width = 10,
                       selectInput(inputId = "performance_summary_criteria",
                                   label = "Evaluation criteria to include",
                                   choices = colnames(res$perfsummary),
                                   selected = colnames(res$perfsummary),
                                   multiple = TRUE, selectize = TRUE)),
                column(width = 2,
                       numericInput(inputId = "performance_summary_plotheight",
                                    label = "Plot height (numeric, in pixels)",
                                    value = 600, min = 200, max = 20000, step = 10)
              )),
              fluidRow(
                column(width = 12,
                       uiOutput("performance_summary_plot_ui"))
              )
            )),

            lapply(aspects, function(w)
              shiny::tabPanel(
                res[[w]]$id,
                fluidRow(
                  column(width = 3,
                         selectInput(inputId = paste0(w, "_facet1"),
                                     label = "Facet by (rows)",
                                     choices = c("none", "dtype", "dataset", "ncells"),
                                     selected = "none", 
                                     multiple = FALSE, 
                                     selectize = FALSE),
                         selectInput(inputId = paste0(w, "_facet2"),
                                     label = "Facet by (columns)",
                                     choices = res[[w]]$facet2,
                                     selected = res[[w]]$facet2[1], 
                                     multiple = FALSE, 
                                     selectize = FALSE)),
                  column(width = 3,
                         checkboxInput(inputId = paste0(w, "_dosqrt"),
                                       label = "Square root transform y-axis",
                                       value = FALSE)),
                  column(width = 2,
                         selectInput(inputId = paste0(w, "_xaxis"),
                                     label = "x-axis",
                                     choices = c("ncells", "method"),
                                     selected = "method",
                                     multiple = FALSE,
                                     selectize = FALSE)),
                  column(width = 2, 
                         sliderInput(inputId = paste0(w, "_yrange"),
                                     label = "y-axis range",
                                     min = res[[w]]$ymin, max = res[[w]]$ymax, 
                                     value = c(res[[w]]$ymin, res[[w]]$ymax))),
                  column(width = 2,
                         numericInput(inputId = paste0(w, "_plotheight"),
                                      label = "Plot height (numeric, in pixels)",
                                      value = 600, min = 200, max = 20000, step = 10))
                ), 
                fluidRow(
                  column(width = 12, DT::dataTableOutput(paste0(w, "_hover_info")))
                ), 
                fluidRow(
                  column(width = 12, uiOutput(paste0(w, "_plot_ui")))
                )
              )),
            
            ## DE characteristics
            lapply(aspects2, function(w)
              shiny::tabPanel(
                res[[w]]$id,
                fluidRow(
                  column(width = 6,
                         selectInput(inputId = paste0(w, "_facet"),
                                     label = "Facet by",
                                     choices = c("characteristic", "method"),
                                     selected = "characteristic", 
                                     multiple = FALSE, 
                                     selectize = FALSE)),
                  column(width = 4,
                         selectInput(inputId = paste0(w, "_xaxis"),
                                     label = "x-axis",
                                     choices = c("characteristic", "method"),
                                     selected = "method",
                                     multiple = FALSE,
                                     selectize = FALSE)),
                  column(width = 2,
                         numericInput(inputId = paste0(w, "_plotheight"),
                                      label = "Plot height (numeric, in pixels)",
                                      value = 600, min = 200, max = 20000, step = 10))
                ), 
                fluidRow(
                  column(width = 12, DT::dataTableOutput(paste0(w, "_hover_info")))
                ), 
                fluidRow(
                  column(width = 12, uiOutput(paste0(w, "_plot_ui")))
                )
              )),
            
            list(shiny::tabPanel(
              "Cross-method similarity",
              fluidRow(
                column(width = 12,
                       uiOutput("crossmethodcons_plot_ui"))
              )
            ))
            
          )
        ))
      )
    )
  }
  
  server_function <- function(input, output, session) {
    ## ====================================================================== ##
    ## Filter results based on input selections
    ## ====================================================================== ##
    filtvals <- reactiveValues()
    observe({
      for (w in c(aspects, aspects2)) {
        filtvals[[w]] <- 
          res[[w]]$data %>% dplyr::filter(method %in% input$keep_methods) %>% 
          dplyr::filter(ncells %in% input$keep_ncells) %>%
          dplyr::filter(dataset %in% input$keep_datasets) %>%
          dplyr::filter(filtering %in% input$keep_filterings)
      }
    })
    
    ## ====================================================================== ##
    ## Cross-method consistency plot and UI object
    ## ====================================================================== ##
    output$crossmethodcons_plot_ui <- renderUI({
      plotOutput("crossmethodcons_plot", width = "100%")
    })
    
    output$crossmethodcons_plot <- renderPlot({
      cmcons <- res$crossmethodcons %>% 
        dplyr::filter(method1 %in% input$keep_methods & method2 %in% input$keep_methods) %>%
        dplyr::filter(dataset %in% input$keep_datasets) %>%
        dplyr::filter(filtering %in% input$keep_filterings) %>%
        dplyr::filter(ncells %in% input$keep_ncells) %>%
        dplyr::group_by(method1, method2) %>%
        dplyr::summarize(meanAUCs = mean(AUCs)) %>% as.data.frame()
      cmcons <- rbind(cmcons, data.frame(method1 = unique(c(cmcons$method1, cmcons$method2)),
                                         method2 = unique(c(cmcons$method1, cmcons$method2)),
                                         meanAUCs = 1, stringsAsFactors = FALSE))
      cmcons <- reshape2::dcast(cmcons, method1 ~ method2, value.var = "meanAUCs")
      rownames(cmcons) <- cmcons$method1
      cmcons$method1 <- NULL
      stopifnot(all(rownames(cmcons) == colnames(cmcons)))
      stopifnot(all((cmcons == t(cmcons))[!is.na(cmcons == t(cmcons))]))
      tmpdist <- 1 - cmcons
      tmpdist[is.na(tmpdist)] <- 1
      hcl_average <- hclust(as.dist(tmpdist))
      plot(hcl_average, main = "", hang = -1, sub = "", xlab = "(complete linkage)", ylab = "")
    })
    
    ## ====================================================================== ##
    ## Performance summary plot and UI object
    ## ====================================================================== ##
    output$performance_summary_plot_ui <- renderUI({
      plotOutput("performance_summary_plot", width = "100%",
                 height = paste0(input$performance_summary_plotheight, "px"))
    })
    
    output$performance_summary_plot <- renderPlot({
      allperf <- res$perfsummary[match(input$keep_methods, rownames(res$perfsummary)), 
                                 match(input$performance_summary_criteria, colnames(res$perfsummary)), 
                                 drop = FALSE]
      allperf <- allperf[order(rowMeans(allperf, na.rm = TRUE), decreasing = TRUE), , drop = FALSE]
      pheatmap(allperf, cluster_rows = FALSE, cluster_cols = FALSE,
               color = c("#E8601C", "#F6C141", "#90C987"), breaks = c(-0.5, 0.5, 1.5, 2.5),
               scale = "none", legend_breaks = c(0, 1, 2),
               legend_labels = c("poor", "intermediate", "good"), fontsize = 14,
               gaps_col = seq_len(ncol(allperf)),
               gaps_row = seq_len(nrow(allperf)))
      
    })
    
    ## ====================================================================== ##
    ## Generate data tables with hover information
    ## ====================================================================== ##
    Map(function(w) {
      output[[paste0(w, "_hover_info")]] <-
        DT::renderDataTable({
          nearPoints(filtvals[[w]], 
                     input[[paste0(w, "_hover")]],
                     threshold = 500, maxpoints = 1, addDist = FALSE, 
                     xvar = input[[paste0(w, "_xaxis")]], 
                     yvar = res[[w]]$yvar) %>%
            dplyr::select(-repl)
        }, rownames = FALSE)
    }, c(aspects, aspects2))
    
    ## ====================================================================== ##
    ## Generate UI objects for plots
    ## ====================================================================== ##
    Map(function(w) {
      output[[paste0(w, "_plot_ui")]] <-
        renderUI({
          plotOutput(paste0(w, "_plot"), width = "100%",
                     height = paste0(input[[paste0(w, "_plotheight")]], "px"),
                     hover = paste0(w, "_hover"))
        })
    }, c(aspects, aspects2))
    
    ## ====================================================================== ##
    ## Generate plots
    ## ====================================================================== ##
    Map(function(w) {
      output[[paste0(w, "_plot")]] <-
        renderPlot({
          p <- ggplot(filtvals[[w]],
                      aes_string(x = input[[paste0(w, "_xaxis")]], 
                                 y = res[[w]]$yvar, color = "method", 
                                 group = "method"))
          if (input[[paste0(w, "_xaxis")]] == "method")  ## boxplot + points
            p <- p + geom_boxplot(outlier.size = -1) + 
              geom_point(position = position_jitter(width = 0.2), size = 0.5) + 
              guides(color = "none")
          else  ## points + smooth
            p <- p + geom_point(size = 1, aes(alpha = method)) + 
              geom_line(stat = "smooth", size = 0.75, method = "loess", span = 1, 
                        na.rm = TRUE, aes(alpha = method)) + 
              theme(legend.position = "bottom") + 
              guides(color = guide_legend(title = "", nrow = 4), alpha = "none")
          p <- p + theme_bw() + xlab("") + 
            scale_color_manual(values = res$cols) + 
            scale_alpha_manual(values = structure(
              c(0.2, 1)[(res$allmethods == input$highlight_method) + 1] - 
                max(c(0.2, 1)[(res$allmethods == input$highlight_method) + 1]) + 1, 
              names = res$allmethods)) + 
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
                  axis.text.y = element_text(size = 12),
                  axis.title.y = element_text(size = 13),
                  strip.text = element_text(size = 12),
                  legend.position = "bottom") + ylab(res[[w]]$id)
          if (input[[paste0(w, "_dosqrt")]]) 
            p <- p + scale_y_sqrt(limits = c(input[[paste0(w, "_yrange")]][1], 
                                             input[[paste0(w, "_yrange")]][2]))
          else 
            p <- p + ylim(input[[paste0(w, "_yrange")]][1], input[[paste0(w, "_yrange")]][2])
          if (any(c(input[[paste0(w, "_facet1")]], input[[paste0(w, "_facet2")]]) != "none"))
            p <- p + facet_grid(as.formula(paste0(trsf(input[[paste0(w, "_facet1")]]), 
                                                  " ~ ", trsf(input[[paste0(w, "_facet2")]]))))
          p
        })
    }, aspects)
    
    ## DE characteristics
    Map(function(w) {
      output[[paste0(w, "_plot")]] <-
        renderPlot({
          p <- ggplot(filtvals[[w]], 
                      aes_string(x = input[[paste0(w, "_xaxis")]], y = res[[w]]$yvar, 
                                 color = "method")) + 
            geom_hline(yintercept = 0) + theme_bw() + xlab("") + 
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
                  axis.text.y = element_text(size = 12),
                  axis.title.y = element_text(size = 13),
                  legend.position = "none") + 
            scale_color_manual(values = res$cols) + 
            geom_boxplot(outlier.size = -1) + 
            geom_point(position = position_jitter(width = 0.2), size = 0.5) + 
            facet_wrap(as.formula(paste0("~ ", input[[paste0(w, "_facet")]])),
                       scales = "free_y") + ylab("SNR")
          p
        })
    }, aspects2)
    
  }
  
  shinyApp(ui = p_layout, server = server_function)
}

#print(summary_app(res, aspects, aspects2))
