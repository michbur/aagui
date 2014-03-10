options(rgl.useNULL=TRUE)

library(shiny)
library(shinyRGL)
library(rgl)
library(tree)

data(ii)
data(aa_prop)
data(aa_nprop)
data(prop_names)
short_prop_names <- 1L:length(prop_names[["short"]])
names(short_prop_names) <- prop_names[["short"]]

shinyServer(function(input, output) {
  # Defines the input for the total number of positive molecules (k) and the total number of 
  # partitions (n).
  
  output$x_ii_ctrl <- renderUI({
    selectInput("x_axs", "X axis:", choices = short_prop_names[ii[[input$x_ii]]], 1)
  })
  output$y_ii_ctrl <- renderUI({
    selectInput("y_axs", "Y axis:", choices = short_prop_names[ii[[input$y_ii]]], 2)
  })
  output$z_ii_ctrl <- renderUI({
    selectInput("z_axs", "Z axis:", choices = short_prop_names[ii[[input$z_ii]]], 3)
  })
  
  output$w_ii_ctrl <- renderUI({
    selectInput("w_axs", "W axis:", choices = short_prop_names[ii[[input$w_ii]]], 3)
  })
  
  output$table_n <- renderTable({
    aa_nprop[as.numeric(c(input$x_axs, input$y_axs, input$z_axs)), ]
  })
  
  output$table <- renderTable({
    aa_prop[as.numeric(c(input$x_axs, input$y_axs, input$z_axs)), ]
  })
  
  output$test <- renderPrint({
    cat(paste0(c("Wybrane kolumny", input$x_axs, input$y_axs, input$z_axs), collapse = " "))
  })
  
  output$col_scale <- renderPlot({
    z <- 0L:10/10
    image(z = matrix(z, ncol = 1), col = sapply(z, function(i) adjustcolor("red", red.f = i)), 
          xlim = c(0.035, 1.035), yaxt = "n")
  })
  
  output$sctPlot <- renderWebGL({
    if (input$spheres == FALSE) {
      for (i in 1L:20) {
        spheres3d(aa_nprop[as.numeric(input$x_axs), i],
                  aa_nprop[as.numeric(input$y_axs), i],
                  aa_nprop[as.numeric(input$z_axs), i],
                  radius = 0.03, 
                  col = adjustcolor("red", red.f = aa_nprop[as.numeric(input$w_axs), i]))
      }
    } else {
      for (i in 1L:20) {
        text3d(aa_nprop[as.numeric(input$x_axs), i],
               aa_nprop[as.numeric(input$y_axs), i],
               aa_nprop[as.numeric(input$z_axs), i],
               text = colnames(aa_prop)[i], 
               col = adjustcolor("red", red.f = aa_nprop[as.numeric(input$w_axs), i]))
      }
    }
    
    
    axes3d()
    grid3d(c("x", "y", "z"))
    title3d('','', "X","Y","Z")
  })
  
  output$clust <- renderPlot({
    inds <- as.numeric(c(input$x_axs, input$y_axs, input$z_axs, input$w_axs))
    cl <- hclust(dist(t(aa_nprop[inds, ])))
    if (input$trim > 0) 
      cl <- contrim(cl, input$trim)
    plot_cluster(cl)
    
  })
  
  output$groups <- renderTable({
    inds <- as.numeric(c(input$x_axs, input$y_axs, input$z_axs, input$w_axs))
    cl <- hclust(dist(t(aa_nprop[inds, ])))
    gr <- cutree(cl, h = input$trim)
    data.frame(Name = unique(gr), Aas = sapply(unique(gr), function(i) 
      paste0(names(gr[gr == i]), collapse = " ")))
  })
  
  output$text <- renderPrint({
    inds <- as.numeric(c(input$x_axs, input$y_axs, input$z_axs, input$w_axs))
    cl <- hclust(dist(t(aa_nprop[inds, ])))
    gr <- cutree(cl, h = input$trim)
    gr <- lapply(unique(gr), function(i) 
      names(gr[gr == i]))
    names(gr) <- 1:length(gr)
    dput(gr)
  })
})