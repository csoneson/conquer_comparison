suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

topdir <- ".."

## Load help functions and summarized data
res <- readRDS(paste0(topdir, "/figures/multi_dataset/fracNA/summary_fracNA_real.rds"))

## Define app
summary_app <- function(res) {

  p_layout <- function(request) {
    shinydashboard::dashboardPage(
      skin = "blue", 
      
      shinydashboard::dashboardHeader(title = "Exploration of results", titleWidth = 300),
      
      shinydashboard::dashboardSidebar(
        checkboxGroupInput(inputId = "keep_methods", label = "Methods to include",
                           choices = unique(res[["fracna_comb_"]]$data$method),
                           selected = unique(res[["fracna_comb_"]]$data$method))
      ),
      
      shinydashboard::dashboardBody(fluidRow(
        shinydashboard::tabBox(
          width = 12, 
          
          # tabPanel("About",
          #          includeMarkdown(paste0(topdir, "shiny/about_app.md"))),
          # 
          tabPanel("Fraction NA adjusted p-values",
                   fluidRow(
                     column(width = 6, uiOutput("fracna.plot.1.ui")),
                     column(width = 6, uiOutput("fracna.plot.2.ui"))
                   ),
                   fluidRow(
                     column(width = 8,
                            checkboxGroupInput(inputId = "keep_ncells_fracna",
                                               label = "Number of cells",
                                               choices = as.character(sort(as.numeric(as.character(unique(res[["fracna_comb_"]]$data$ncells_fact))))),
                                               selected = as.character(sort(as.numeric(as.character(unique(res[["fracna_comb_"]]$data$ncells_fact))))),
                                               inline = TRUE)),
                     column(width = 2, checkboxInput(inputId = "dosqrt",
                                                     label = "Square root transformation",
                                                     value = FALSE)),
                     column(width = 2, sliderInput(inputId = "yrange",
                                                   label = "y axis range",
                                                   min = 0, max = 1,
                                                   value = c(0, 1)))
                   )
          )
          
        )
      ))
    )
  }
  
  server_function <- function(input, output, session) {
    # ===================== PCA plot =========================
    output$fracna.plot.1 <- renderPlot({
      p <- ggplot(res[["fracna_comb_"]]$data %>% 
                    dplyr::filter(method %in% input$keep_methods) %>%
                    dplyr::filter(ncells_fact %in% input$keep_ncells_fracna),
                  aes(x = method, y = fracNA, color = plot_color)) +
        geom_boxplot(outlier.size = -1) +
        geom_point(position = position_jitter(width = 0.2), size = 0.5, aes(shape = plot_char)) +
        theme_bw() + xlab("") + ylab("Fraction of NA adjusted p-values") +
        scale_color_identity() + 
        scale_shape_identity() + 
        guides(color = guide_legend(ncol = 2, title = ""),
               shape = guide_legend(ncol = 4, title = "Number of \ncells per group")) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
              axis.text.y = element_text(size = 12),
              axis.title.y = element_text(size = 13)) + 
        ggtitle("Without filtering")
      if (input$dosqrt) p <- p + scale_y_sqrt(limits = c(input$yrange[1], input$yrange[2]))
      else p <- p + ylim(input$yrange[1], input$yrange[2])
      p
    })
      
    output$fracna.plot.2 <- renderPlot({
      p <- ggplot(res[["fracna_comb_TPM_1_25p"]]$data %>% 
                    dplyr::filter(method %in% input$keep_methods) %>%
                    dplyr::filter(ncells_fact %in% input$keep_ncells_fracna),
                  aes(x = method, y = fracNA, color = plot_color)) +
        geom_boxplot(outlier.size = -1) +
        geom_point(position = position_jitter(width = 0.2), size = 0.5, aes(shape = plot_char)) +
        theme_bw() + xlab("") + ylab("Fraction of NA adjusted p-values") +
        scale_color_identity() + 
        scale_shape_identity() + 
        guides(color = guide_legend(ncol = 2, title = ""),
               shape = guide_legend(ncol = 4, title = "Number of \ncells per group")) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
              axis.text.y = element_text(size = 12),
              axis.title.y = element_text(size = 13)) + 
        ggtitle("After filtering (TPM > 1 in > 25% of cells)")
      if (input$dosqrt) p <- p + scale_y_sqrt(limits = c(input$yrange[1], input$yrange[2]))
      else p <- p + ylim(input$yrange[1], input$yrange[2])
      p
    })
    
    output$fracna.plot.1.ui <- renderUI({
      plotOutput("fracna.plot.1", height = "600px")
    })
    output$fracna.plot.2.ui <- renderUI({
      plotOutput("fracna.plot.2", height = "600px")
    })
    
  }
  
  shinyApp(ui = p_layout, server = server_function)
}

#print(summary_app(res))
