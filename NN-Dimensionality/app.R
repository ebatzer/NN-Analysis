#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny); library(tidyverse); library(plot3D); 
library(ggrepel); library(ggpmisc); library(DT); library(shinythemes)

# Reading in datasets
SS_pvals <- read.csv("Data/tradeoffs_RRPP_Pvals_Summary.csv", row.names = 1) %>%
  
  filter(site != "saline.us")

specscores <- read.csv("Data/tradeoffs_specscores_long.csv", row.names = 1)   %>%
  
  filter(site != "saline.us")

sitenames <- as.character(unique(SS_pvals$site))
  
# Creating the projection matrix
dim1 <- c(1, 1, 1) / sqrt(sum(c(1, 1, 1)^2))
dim2 <- c(0, -1, 1) / sqrt(sum(c(0, -1, 1)^2))
dim3 <- c(-1, .5, .5) / sqrt(sum(c(-1, .5, .5)^2))

projmat <- matrix(c(dim1, dim2, dim3), nrow = 3)

# Define UI for application that draws a histogram
ui <- fixedPage(theme = shinytheme("readable"),

  fixedRow(
     column(12,
            h1("Nutrient Network Site Response Dimensionality"),
            br(),
            selectInput("selected_site", "Chosen Site", 
                        sitenames, selected = "hopl.us", multiple = FALSE,
                        selectize = TRUE, width = NULL, size = NULL),
            h2("Model Fitting and Sum of Squares Decomposition"),
            column(8, 
                   "To estimate species responses to different single nutrient addition treatments,
                   a multiple linear regression model was fit to log2-transformed species abundance data 
                   at all sites containing at least 5 years of observations.",
                   br(), br(),
                   "Model Formula:",
                   br(), br(),
                   "Y ~ Year_Fctr + Block + Treatment : Log(Year_Num)",
                   br(), br(), 
                   "Where the abundance of a species in a given observation, Y, is a function of site-level interannual variation (Year_Fctr), spatial heterogeneity (Block), and a time-saturating relationship with nutrient addition treatment (Treatment : Log(Year_Num)).",
                   br(),br(), 
                   "The significance of model terms was calculated using a randomized residual permutation procedure (RRPP). Resulting sums-of-squares decomposition and P-values for different coefficients within a site are shown to the right.",
                   br(),br(), 
                   "While sums-of-squares can be difficult to interpret in a multivariate dataset using transformed data, general sources of variation within a site can be checked here. Do sites with large year-to-year variation at the site scale show significance of the year_trt term? Do sites with strong responses to different nutrient treatments show significance of each of these terms? Are more impactful nutrient addition treatments associated with larger sums of squares?"),
            column(4, 
                   tableOutput('SStable')))),
  column(12,
         br(),
         h2("Pairwise Comparisons of Species Responses"),
         br(),
         "Following model fitting, treatment coefficients were extracted to estimate species responses to nutrient enrichment. Below shows a visualization of these fitted coefficients, filtered to just those species occurring in >20% of all community observations within a site. This figure can be manipulated to show pairwise contrasts (if two different variables are selected), or a dot plot displaying the magnitude of response organized alphabetically (if both selected variables are the same). Points are colored by functional group.",
         br(), br(),
         "Are treatments with correlated patterns of community change (positive or negative slopes) as expected?",
         br(), br()),
 
  fluidRow(column(6,
                  selectInput("X_chosen", "Variable 1 (X axis)", 
                              c("N", "P", "K"), selected = "N", multiple = FALSE,
                              selectize = TRUE, width = NULL, size = NULL)),
           column(6,
                  selectInput("Y_chosen", "Variable 2 (Y axis)", 
                              c("N", "P", "K"), selected = "P", multiple = FALSE,
                              selectize = TRUE, width = NULL, size = NULL))),
   fluidRow(
     column(12, 
            plotOutput('specviz'))
   ),
  column(12,
         br(), br(),
         h2("3D Visualization of Response Vectors"),
         br(),
         "To visualize differences in estimated response trajectories, 3D visualizations of filtered treatment response coefficients are presented below. In this figure, each point corresponds to the estimated change in cover (at the log2 scale) for each species within a site, per log(year) of treatment.",
         br(), br(),
         "The 1:1:1 line presented illustrates a one-dimensional (neutral) expectation of nutrient enrichment response, where species responses are consistent across all nutrient types. Differences between species observed responses and their projection on this line are highlighted by dashed lines.",
         br(), br(),
         "Total variance of species responses captured by this projection will range between 0-1. When this neutral expectation fails to characterized response patterns, observed responses will deviate strongly from this 1:1:1 line (variance captured = 0); consistently one-dimensional responses across species will fall largely long this 1:1:1 line (variance captured = 1)."
  ),
  fixedRow(column(12, align = "center",
                  plotOutput('threedviz', height = 500, width = 500),
                  br(), br(),
                  textOutput('formula'))
  ),
  
  column(12,
         br(), br(),
         h2("Table of Fitted Response Coefficients"),
         "A summary table of all (unfiltered) fitted responses and species functional group assignments is presented below. Much of this information can be gleaned from previous figures, but may benefit from a quick inspection of general patterns. Are the most extreme responders (on the log2 scale) as expected?",
         br(), br(),
         "For perspective, a value of '+1' indicates that a species doubles in cover for each log-year of treatment, meaning that this species is estimated to increase in cover by ~70% after 1 year of treatment, by ~100% after 2 years of treatment, and by ~150% after 5 years of treatment.",
         br(), br(),
         "Range of each set of values is as follows:",
         tableOutput('resprange'),
         br(), br()),
   fluidRow(column(12), 
            DTOutput('resptable')
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
   
   output$SStable <- renderTable({
     
      # subset rows
      ss_row <- SS_pvals[SS_pvals$site == input$selected_site,]
      
      ss_row %>% 
        mutate(P_value = format.pval(P_Value, digits = 3))

   })
   
   output$threedviz <- renderPlot({

     # Filter Site Dat
     sitedat <- specscores %>% 
       
       filter(site == input$selected_site) %>%
       
       filter(Occ_Total > 0.2) %>%
       
       # Rescaling responses
       mutate(trt_K = scale(trt_K, scale = TRUE),
              trt_N = scale(trt_N, scale = TRUE),
              trt_P = scale(trt_P, scale = TRUE)) 

     # 1:1:1 Projection
     site_proj <- data.frame(as.matrix(sitedat %>% select(trt_N, trt_P, trt_K)) %*% projmat)[,1] / sqrt(3)
     
     # Plotting parameters
     max_line <- (matrix(c(-3,3,-3,3,-3,3),
                        nrow = 2) %*% projmat)[,1] / sqrt(3)
     theta = 40
     phi = 30
     
     par(mar = c(2, 2, 2, 0))
     
     scatter3D(z = sitedat$trt_N, x = sitedat$trt_P, y = sitedat$trt_K, pch = 19, cex = 2,
               theta = theta, phi = phi, xlim = c(-3, 3), ylim = c(-3,3), zlim = c(-3,3),
               xlab = "Stdized P Response", ylab = "Stdized K Response", zlab = "Stdized N Response", type = "p",
               main = "Site-Level Treatment Responses", 
               colkey = FALSE, ticktype = "detailed",
               col = "black")
     
     for(i in 1:length(site_proj)){
       scatter3D(z = c(site_proj[i], sitedat$trt_N[i]), 
                 x = c(site_proj[i], sitedat$trt_P[i]), 
                 y = c(site_proj[i], sitedat$trt_K[i]), pch = 19, lwd = 2,
                 theta = theta, phi = phi, col = "grey40", add = TRUE, type = "l", lty = 2)
     }

     
     scatter3D(z = max_line, x = max_line, y = max_line, pch = 19, lwd = 2,
               theta = theta, phi = phi, col = "black", add = TRUE, type = "l", lty = 1)
     
   })
   
   output$formula <- renderText({
     
     sitedat <- specscores %>% 
       
       filter(site == input$selected_site) %>%
       
       filter(Occ_Total > 0.2) %>%
       
       # Rescaling responses
       mutate(trt_K = scale(trt_K, scale = TRUE),
              trt_N = scale(trt_N, scale = TRUE),
              trt_P = scale(trt_P, scale = TRUE)) 
     
     
     siteproj <- as.matrix(sitedat %>% select(trt_N, trt_P, trt_K)) %*% projmat     
      
     var_frac <- apply(siteproj, MARGIN = 2, function(x) sum(x^2))[1] / 
       sum(apply(siteproj, MARGIN = 2, function(x)sum(x^2)))
   
     paste("Fraction of variance in site responses captured by neutral (1:1:1) expectation:",
           round(var_frac,3))  
   })
   
   output$specviz <- renderPlot({
     
     sitedat <- specscores %>% 
       
       filter(site == input$selected_site) %>% 
       
       filter(Occ_Total > .2) %>%
       
       # Rescaling responses
       mutate(trt_K = scale(trt_K, scale = TRUE),
              trt_N = scale(trt_N, scale = TRUE),
              trt_P = scale(trt_P, scale = TRUE)) 
     
     X_chosen = paste("trt_", input$X_chosen, sep = "")
     Y_chosen = paste("trt_", input$Y_chosen, sep = "")
     
     chosen_dat <- sitedat %>%
       select(X_chosen, Y_chosen, Species, functional_group)
       
     chosen_dat %>% 
       ggplot(aes(x = !!sym(X_chosen),
                  y = !!sym(Y_chosen),
                  label = Species,
                  color = functional_group)) +
       theme_minimal() + 
       stat_smooth(method = "lm", se = FALSE, color = "black") +
       geom_point(size = 4) +
       stat_dens2d_filter(geom = "label_repel", keep.fraction = 0.2) +
       xlab(paste("Estimated", input$X_chosen, "Response (Change in Log2 Cover / Log Year of Treatment)")) +
       ylab(paste("Estimated", input$Y_chosen, "Response (Change in Log2 Cover / Log Year of Treatment)"))

   })
   
   output$resprange <-  renderTable({
     specscores %>% filter(site == input$selected_site) %>%
       select(trt_K, trt_P, trt_N) %>%
       mutate_if(is.numeric,  round, 4) %>%
       apply(., MARGIN = 2, range) %>%
       as.data.frame(., row.names = c("min", "max")) 
   })
   
   output$resptable <- renderDT({
     
     sitedat <- specscores %>% filter(site == as.character(input$selected_site)) %>%
       select(site, Species, trt_K, trt_P, trt_N, functional_group, Occ_Total) %>%
       mutate_if(is.numeric,  round, 4) %>%
       rename("Site Code" = "site",
              "K Response" = "trt_K",
              "N Response" = "trt_N",
              "P Response" = "trt_P",
              "Functional Group" = "functional_group",
              "Fraction Occupancy" = "Occ_Total")
     
     DT::datatable(sitedat, filter = "top",
                          options = list(pageLength = 20))
     })
}

# Run the application 
shinyApp(ui = ui, server = server)
