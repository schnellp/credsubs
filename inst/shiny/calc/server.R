library(shiny)

shinyServer(function(input, output) {
  
  load("config.RData")
  
  output$title <- renderUI({
    titlePanel(title)
  })
  
  output$instructions <- renderText({
    paste(instructions, "<br /><br />")
  })
  
  for (j in 1:length(cov.space)) {
    if (!is.factor(cov.space[, j]) && 
        (is.null(attr(cov.space[, j], "scaled:scale")) ||
        is.null(attr(cov.space[, j], "scaled:center")))) {
      attr(cov.space[, j], "scaled:scale") <- 1
      attr(cov.space[, j], "scaled:center") <- 0
    }
  }
  
  output$predictors <- renderUI({
    inputs <- list()
    for (predictor in colnames(cov.space)) {
      if (is.factor(cov.space[, predictor])) {
        inputs[[predictor]] <-
          selectInput(predictor,
                      predictor,
                      unique(cov.space[, predictor])
          )
      } else {
        inputs[[predictor]] <-
          selectInput(predictor,
                      predictor,
                      unique(cov.space[, predictor]) *
                        attr(cov.space[, predictor], "scaled:scale") +
                        attr(cov.space[, predictor], "scaled:center")
          )
      }
      
    }
    inputs
  })
  
  output$status <- renderText({
    matches <- 1:nrow(cov.space)
    
    for (predictor in names(cov.space)) {
      if (is.factor(cov.space[, predictor])) {
        matches <-
          intersect(matches,
                    which(cov.space[, predictor]  == input[[predictor]]))
      } else {
        matches <-
          intersect(matches,
                    which(cov.space[, predictor] *
                            attr(cov.space[, predictor], "scaled:scale") +
                            attr(cov.space[, predictor], "scaled:center") ==
                            input[[predictor]]))
      }
      
    }
    
    if (length(matches) == 0) {
      return("Predictor combination untested.")
    } else if (credsubs.level$sign[matches[1]] == 0) {
      return(paste("<b>No conclusion</b> may be drawn",
                    "for the covariate point selected",
                    "<b>at any credible level</b>."))
    } else {
      return(paste0("The covariate point selected is <b>",
               ifelse(credsubs.level$sign[matches] == 1,
                      paste0("<font color=\"#00AA00\">",
                             "in the exclusive credible subset",
                             "</font>"),
                      paste0("<font color=\"#AA0000\">",
                             "not in the inclusive credible subset",
                             "</font>")),
               "</b> at a maximum credible level of <b>",
               round(credsubs.level$level[matches] * 100, digits=2),
               "%</b>. ",
               "At higher credible levels, no conclusion may be drawn."))
    }
  })
  
  observeEvent(input$closeApp, {
    stopApp()
  })
})
