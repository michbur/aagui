library(shiny)

shinyUI(pageWithSidebar(
  headerPanel("Aminoacids clustering"),
  sidebarPanel(
    selectInput("x_ii", "Class of parameters (x-axis):", choices = names(ii), "phys"),
    selectInput("y_ii", "Class of parameters (y-axis):", choices = names(ii), "phys"),
    selectInput("z_ii", "Class of parameters (z-axis):", choices = names(ii), "phys"),
    selectInput("w_ii", "Class of parameters (w-axis):", choices = names(ii), "phys"),
    uiOutput("x_ii_ctrl"),
    uiOutput("y_ii_ctrl"),
    uiOutput("z_ii_ctrl"),
    uiOutput("w_ii_ctrl"),
    sliderInput("trim", "Trim:", min = 0, max = 1, value = 0, step= 0.05),
    checkboxInput(inputId = "spheres", label = "Spheres instead of names:", value = TRUE)
  ),
  mainPanel(
    tabsetPanel(
      tabPanel("Scatter plot", webGLOutput("sctPlot"), plotOutput("col_scale", height = 210),
               plotOutput("clust"), tableOutput("groups")),
      tabPanel("Data table", tableOutput("table"), tableOutput("table_n"), 
               verbatimTextOutput("test")),
      tabPanel("dput", verbatimTextOutput("text"))
    )
  )
))

