#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny); library(tidyverse); library(shinythemes); library(DT)

# Reading in datasets
sb_features <- read.csv("data/sb_features.csv", row.names = 1)
spec_summary <- read.csv("data/spec_summary.csv", row.names = 1, stringsAsFactors = FALSE)
proposed_adjustments <- read.csv("data/proposed_adjustments.csv", row.names = 1)

# Setting initial parameters
sitenames <- unique(sb_features$site_code)

# Define UI for application that draws a histogram
ui <- fixedPage(theme = shinytheme("readable"),
                
                tags$style(type='text/css', "

          .selectize{closeAfterSelect: true;}
          .selectize-dropdown.form-control{height: 50vh !important; } 
          .selectize-dropdown-content{ max-height: 100% !important; height: 100% !important; }"),
   
   # Application title
   titlePanel("NutNet Taxonomic Adjustments"),
   
   # First data chunk                   
   fixedRow(
     column(12, 
            h3("Choose Site:"), 
            selectInput("selected_site", "Site Code", 
                        sitenames, selected = "hopl.us", multiple = FALSE,
                        selectize = TRUE, width = NULL, size = NULL)),
     column(12, h3("Site Overview Table:")),
     column(6, 
            "The following table summarizes general characteristics of each site, including
       the total years of observation, first and last recorded years, total number of blocks,
       and total number of plots."),
     column(6,
            tableOutput('site_tab')),
     column(12,
            h3("Visualization of Species Cover / Frequency"),
            "To help determine overlapping taxa, the following plot can be used to determine
            site-level mean cover (across all plots) or total # of occurences (across all plots)
     in each year of observation."),
     column(9,
            selectInput("plot_type", "Show me:", 
                        c("Cover", "Frequency"), selected = "cover", multiple = FALSE,
                        selectize = TRUE, width = NULL, size = NULL),
            plotOutput('summary_figure')
            ),
     column(3,
            selectizeInput(inputId = "Species", label = "Select Species", selected = "ALL",
                           choices = NA, multiple = TRUE,
                           width = 900, size = 25),
            selectizeInput(inputId = "Family", label = "Select Family", selected = "ALL",
                           choices = NA, multiple = TRUE,
                           width = 900, size = 25),
            selectizeInput(inputId = "Treatment", label = "Select Treatment", selected = "ALL",
                           choices = NA, multiple = TRUE,
                           width = 900, size = 25),
            "Choose which observations to plot. Note that this selection can take multiple 
            inputs. 'ALL' will plot all species/families/treatments present within a site. 
            Names can be added or deleted using this text box.",
            br(), br(),
            "By default, plot legend will not be included if # taxa > 10"
            ),
     column(12, 
            h3("Suggested Taxonomic Changes"),
            "Previously suggested taxonomic changes, grouping one or more raw taxon names (old_taxon) to
            a new identifier (new_taxon).",
            br(), br(), br(),
            DTOutput('suggested_changes'))
     
     )
)
# Define server logic required to draw a histogram
server <- function(input, output, session) {
   
  output$site_tab <- renderTable(
    sb_features %>% filter(site_code == input$selected_site)
  )
  
  # Update selectizer for figures
  observe({
    
    # Species
    updateSelectizeInput(session = session, inputId = "Species", selected = "ALL",
                         choices = c("ALL", as.character(unique(choices_species()$Taxon))))
    
    # Family
    updateSelectizeInput(session = session, inputId = "Family", selected = "ALL",
                         choices = c("ALL", as.character(unique(choices_species()$Family))))

    # Treatment
    updateSelectizeInput(session = session, inputId = "Treatment", selected = "ALL",
                         choices = c("ALL", as.character(unique(choices_species()$trt))))
  })
  
  # Generate figure
  plotY <-reactive({
    sym(if_else(input$plot_type == "Cover", "mean_cover", "obs_freq"))
    })
  
  # Family
  selectedfamily <- reactive({
    if("ALL" %in% as.character(input$Family)){
      as.character(unique(spec_summary$Family[spec_summary$site_code == input$selected_site]))
    }else{
      as.character(input$Family)
    }
  })  
  
  selectedspecs <- reactive({
    # Species
    if("ALL" %in% as.character(input$Species)){
      as.character(unique(spec_summary$Taxon[spec_summary$site_code == input$selected_site]))
    }else{
      as.character(input$Species)
    }
  })

    # Treatment
  selectedtrt <- reactive({
    if("ALL" %in% as.character(input$Treatment)){
      as.character("ALL")
    }else{
      as.character(input$Treatment)
    }
  })
  
  choices_species <- reactive({
    spec_summary %>%
      filter(site_code == input$selected_site)
  })
  
  plotparams <- reactive({
    sb_features %>% filter(site_code == input$selected_site)
  })
  
  output$summary_figure <- renderPlot(

    (newdf <- spec_summary %>% 
      filter(site_code == input$selected_site) %>%
      filter(Taxon %in% selectedspecs(),
             Family %in% selectedfamily(),
             trt %in% selectedtrt()) %>%
      group_by(Taxon, site_code, year) %>%
      summarise(mean_cover = mean(mean_cover),
                obs_freq = sum(obs_freq))) %>%
      ggplot(aes(x = year,
                 y = !!plotY(),
                 color = Taxon,
                 fill = Taxon)) +
      geom_line(size = 1) +
      geom_point(pch = 21, size = 5, color = 'black') + 
      theme_bw() + xlab("Year") + 
      ylab(input$plot_type) +
      theme(text = element_text(size=20)) + 
      scale_x_continuous(limits = c(plotparams()$first_year, plotparams()$last_year),
                         breaks = seq(plotparams()$first_year, plotparams()$last_year, by = 1)) +
      if(length(unique(newdf$Taxon)) > 10){guides(color = FALSE, fill = FALSE)}

    )
  
  
  output$suggested_changes <- renderDT(
    
    proposed_adjustments %>% filter(site_code == input$selected_site) %>%
      arrange(new_taxon)
  )
}

# Run the application 
shinyApp(ui = ui, server = server)

