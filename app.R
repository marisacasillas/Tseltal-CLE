library(shiny)
source("1-GUI-TseltalCLE-findings.R")
model.res.fig <-"c_o.tpm_random_log_gaus.res.plot.png"

# Define UI for data upload app ----
ui <- fluidPage(

  # App title ----
  titlePanel("Tseltal child language environment estimates"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      # Input: Sample ----
      selectizeInput("sample", "Which clips?",
                   choices = c("Random", "High turn-taking", "High vocal activity", "All"),
                   options = list(
                     placeholder = 'Select a clip type below',
                     onInitialize = I('function() { this.setValue(""); }'))),

      # Input: Data version ----
      selectizeInput("version", "Which version of the dataset?",
                   choices = c("Casillas, Brown, & Levinson (submitted March 2019)"),
                   options = list(
                     placeholder = 'Select a dataset version below',
                     onInitialize = I('function() { this.setValue(""); }'))),

      # Input: Measures ----
      selectizeInput("measures", "Which measure?",
                   choices = c("Target-child-directed speech (TCDS) min/hr",
                               "All child-directed speech (CDS) min/hr",
                               "All other-directed-speech (ODS) min/hr",
                               "All speech (XDS) min/hr",
                               "Other-to-Target Child turn transitions/min",
                               "Target Child-to-Other turn transitions/min",
                               "Interactional sequence duration (sec)"),
                   options = list(
                     placeholder = 'Select an measure below',
                     onInitialize = I('function() { this.setValue(""); }'))),

      # Input: Models ----
      selectizeInput("model", "Wanna see model output(s)?",
                   choices = c("Yes", "No"),
                   options = list(
                     placeholder = 'Select an option below',
                     onInitialize = I('function() { this.setValue(""); }'))),

      # Submit button:
      actionButton("submit", "Update")
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      uiOutput("report")
    )
  )
)

# Define server logic to read selected file ----
server <- function(input, output) {
  report <- eventReactive(input$submit, {
    req(input$sample, input$version, input$measures, input$model)
    retrieve.summary(input$sample,
                  input$version,
                  input$measures,
                  input$model)
  })

  output$report <- renderUI({
    req(report())
    
    tagList(
      tags$h1("Summary statistics"),
      renderTable(report()$sum.stat.tbl),
      tags$h1("Graphical summary"),
      renderPlot(report()$sum.stat.fig),
      if (model.res.fig != "No") {
        tags$h1("Graphical summary")
        renderPlot(report()$sum.stat.fig)
      },
      tags$br()
    )
  })

  }

# Create Shiny app ----
shinyApp(ui, server)