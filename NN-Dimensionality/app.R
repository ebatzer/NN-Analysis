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

SS_pvals <- read.csv("Data/tradeoffs_RRPP_Pvals_Summary.csv", row.names = 1)

specscores_full <- read.csv("Data/tradeoffs_specscores_long.csv", row.names = 1) 

specscores <- specscores_full %>%
  filter(Occ_Total > 0.33)

sitenames <- as.character(unique(SS_pvals$site))
  

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
         "Following model fitting, treatment coefficients were extracted to estimate species responses to nutrient enrichment. Below shows a visualization of these fitted coefficients, filtered to just those species occurring in >33% of all community observations within a site. This figure can be manipulated to show pairwise contrasts (if two different variables are selected), or a dot plot displaying the magnitude of response organized alphabetically (if both selected variables are the same). Points are colored by functional group.",
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
         "A best fit line generated through principal components decomposition is highlighted in black. The best fit line generated from data across all sites is presented as a dashed line. To better assess how a site's patterns deviate from this global average, the angle of difference between vectors is presented below.",
         br(), br(),
         "The right panel shows a screeplot of the proportion of variance captured by different principal components. Larger variance captured by PC1 indicates more 1-dimensional patterns than variance spread equally across all three PCs.",
         br(), br()          
  ),
  fixedRow(  withMathJax(),
             column(8, 
                    plotOutput('threedviz'),
                    textOutput('formula')),
             column(4, 
                    plotOutput('screeplot'))
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

     # All site PCA
     as_PCA <- prcomp(scale(specscores %>% dplyr::select(trt_N, trt_K, trt_P)))

     # Filter Site Dat
     sitedat <- specscores %>% filter(site == input$selected_site)

     # Site-level PCA
     site_PCA <- prcomp(scale(sitedat %>% dplyr::select(trt_N, trt_K, trt_P)))
     
     theta = 40
     phi = 30
     
     pc1_loadings <- site_PCA$rotation[,1]
     pc1_fit <- data.frame(t((pc1_loadings) %*% t(c(-4, 4))))
     names(pc1_fit) <- names(pc1_loadings)
     
     all_pc1_loadings <- as_PCA$rotation[,1]
     all_pc1_fit <- data.frame(t((all_pc1_loadings) %*% t(c(-4, 4))))
     names(all_pc1_fit) <- names(all_pc1_loadings)
     
     # par(mfrow = c(2,1))
     
     par(mar = c(2, 2, 2, 0))
     
     scatter3D(z = sitedat$trt_N, x = sitedat$trt_P, y = sitedat$trt_K, pch = 19, cex = 1,
               theta = theta, phi = phi, xlim = c(-3, 3), ylim = c(-3,3), zlim = c(-3,3),
               xlab = "P Response", ylab = "K Response", zlab = "N Response", type = "p",
               main = "Site-Level Treatment Responses", 
               colkey = FALSE, ticktype = "detailed",
               col = "black")
     
     scatter3D(z = pc1_fit$trt_N, x = pc1_fit$trt_P, y = pc1_fit$trt_K, pch = 19, lwd = 2,
               theta = theta, phi = phi, col = "black", add = TRUE, type = "l", lty = 1)
     
     scatter3D(z = all_pc1_fit$trt_N, x = all_pc1_fit$trt_P, y = all_pc1_fit$trt_K, pch = 19, lwd = 2,
               theta = theta, phi = phi, col = "grey40", add = TRUE, type = "l", lty = 2)
     
   })
   
   output$formula <- renderText({
     
     # All site PCA
     as_PCA <- prcomp(scale(specscores %>% dplyr::select(trt_N, trt_K, trt_P)))
     
     # Filter Site Dat
     sitedat <- specscores %>% filter(site == input$selected_site)
     
     # Site-level PCA
     site_PCA <- prcomp(scale(sitedat %>% dplyr::select(trt_N, trt_K, trt_P)))
     
     a <- as_PCA$rotation[,1] 
     b <- site_PCA$rotation[,1]
     my_calculated_value <-  acos( sum(a*b) / ( sqrt(sum(a * a)) * sqrt(sum(b * b)) ) ) * (180/pi)
     
     if(my_calculated_value > 90){
       my_calculated_value = 180 - my_calculated_value 
     }
     
     paste("Angle between site-level and global PC vectors is ", round(my_calculated_value, 3), "degrees")
   })
   
   output$screeplot <- renderPlot({
     
     # Filter Site Dat
     sitedat <- specscores %>% filter(site == input$selected_site)
     
     # Site-level PCA
     site_PCA <- prcomp(scale(sitedat %>% dplyr::select(trt_N, trt_K, trt_P)))
     
     data.frame(x = c(1, 2, 3), y = summary(site_PCA)$importance[2,]) %>%
       ggplot(aes(x = x,
                  y = y)) +
       geom_point(size = 2) +
       geom_line() +
       ylim(0, 1) +
       theme_bw() +
       ylab("Proportion of Variance Explained") +
       xlab("") +
       scale_x_continuous(breaks = c(1,2,3),
                          labels = c("PC1", "PC2", "PC3"))
   })
   
   output$specviz <- renderPlot({
     
     sitedat <- specscores %>% filter(site == input$selected_site)
     
     X_chosen = paste("trt_", input$X_chosen, sep = "")
     Y_chosen = paste("trt_", input$Y_chosen, sep = "")
     
     chosen_dat <- sitedat %>%
       select(X_chosen, Y_chosen, Species, functional_group)
       
     chosen_dat %>% 
       ggplot(aes(x = chosen_dat[,1],
                  y = chosen_dat[,2],
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
     
     sitedat <- specscores_full  %>% filter(site == input$selected_site) %>%
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
