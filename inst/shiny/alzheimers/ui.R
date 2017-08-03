shinyUI(fluidPage(
  
  uiOutput("title"),
  
  htmlOutput("instructions"),
  
  uiOutput("predictors"),
  
  h3("Result"),
  htmlOutput("status")
))
