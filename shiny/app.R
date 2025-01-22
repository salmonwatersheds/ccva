###############################################################################
# Packages
###############################################################################

library(shiny)
library(shinymanager)
library(shinyBS)
library(shinyWidgets)
library(leaflet)
library(colorRamps)
library(mapproj)
library(proj4)
library(sf)
library(dplyr)
library(wesanderson)
library(ncdf4) # package for netcdf manipulation
library(zoo) # package with rollmean function
library(shiny)
library(leaflet)

cols <- wes_palette("Darjeeling1")

# gpclibPermit()

###############################################################################
# Load data
###############################################################################

#------------------------------------------------------------------------------
# Conservation Unit Boundaries
#------------------------------------------------------------------------------
cu_list <- read.csv("data/cu_list.csv")
cu_boundary <- readRDS("data/cu_boundary_fraser.rds")

#------------------------------------------------------------------------------
# spatial distribution
#------------------------------------------------------------------------------

# Load grid cell ids for each CU and life stage
incl_FW <- readRDS("data/PCIC_incl.rds") # freshwater stages (and CU-specific)
incl_EM <- read.csv("data/incl_NOAA_EM.csv")[, 2] # early marine (same for all CUs)
incl_MR <- read.csv("data/incl_NOAA_marRear.csv")[, 2] # marine residence (same for all CUs)

# Load grid cells
grid_polys_fw <- readRDS("data/grid_polys_fw.rds")
grid_polys_mar <- readRDS("data/grid_polys_mar.rds")

# load spatial output
spat <- readRDS("data/spat_combined.rds")

#------------------------------------------------------------------------------
# Timing
#------------------------------------------------------------------------------

timing <- read.csv("data/timing-fraser_days.csv")
DOY <- c(1:365)
xDate <- as.Date(paste(2019, DOY), format = "%Y %j")
labDate <- as.Date(c("2019-01-01", "2019-04-01", "2019-07-01", "2019-10-01", "2019-12-31"))
labName <- c("Jan", "Apr", "Jul", "Oct", "Jan")

stages <- unique(timing$stage)

overall_exposure <- read.csv("data/overall_exposure_rcp45_mid.csv")
stage_exposure <- read.csv("data/stage_exposure_rcp45_mid.csv")
stagevar_exposure <- read.csv("data/stagevariable_exposure_rcp45_mid.csv")

#------------------------------------------------------------------------------
# set colour palette
#------------------------------------------------------------------------------

col_levels <- seq(0, 1, 0.05)
n <- length(col_levels)
col_palette <- paste0(colorRampPalette(colors = cols[c(5,3,1)])(n = length(col_levels) + 1))

col_palette_fun <- colorNumeric(palette = col_palette, domain = c(0,1))

###############################################################################
# RUn map app
###############################################################################
species_opt <- unique(cu_list$species_pooled)
stage_opt <- c("Incubation", "Freshwater rearing", "Early marine", "Marine rearing", "Adult freshwater migration", "Spawning")
stage_disp <- c("Incubation", "Freshwater\nrearing", "Early\nmarine", "Marine\nrearing", "Adult\nfreshwater\nmigration", "Spawning")
var_lookup <- cbind(c("Stream temperature", "Low flow", "SST", "Salinity"), c("stream_temp", "stream_flow", "SST", "SSS"))

shinyApp(
  
  ui = fluidPage(
    
    tags$head(
      includeCSS("www/bootstrap_custom.css")
    ),
    
    HTML("<h2 style=' background: #0B1C30; padding: 12px'><img src='SWP-logo-2023r12-horizontal-reverse-trans-2000w.png' alt='PSF logo' width='300'>  </h2>"),
    
    HTML("<h2 style='color: #3F3F46; padding: 6px'>CLIMATE CHANGE EXPOSURE <p><small>of Pacific Salmon and Steelhead Conservation Units </small></p> </h2>"),
    
    HTML("This interactive mapping tool provides information on the spatial distribution of exposure to climate changes for Pacific salmon and steelhead Conservation Units in the Fraser River basin. Exposure is quantified as the proportion of time that each life stage encounters potentially harmful high stream temperatures, low stream flows, high sea surface temperature (SST), and low salinity. Outputs shown here are based on Global Climate Model projections for the mid-century period (2040-2069) under the business-as-usual emissions scenario (RCP 4.5). See the ABOUT tab for further information.<br>"),
    
    br(),
    tabsetPanel(type = "tabs",
                selected = "DATA",
                tabPanel("DATA",
                         br(),
                         
                         fluidRow(
                           
                           column(width = 1),
                           column(width = 2,
                                  selectInput(inputId = "sp", label = "Species",
                                              choices = species_opt,
                                              selected = "Sockeye")
                           ),
                           column(width = 3, 
                                  selectInput(inputId = "cu", label =  "Conservation Unit",
                                              choices = cu_list$cu_name_pse)
                           ),
                           column(width = 3,
                                  selectInput(inputId = "stage", label =   "Life stage",
                                              choices = stage_opt,
                                              selected = "Adult freshwater migration")
                           ),
                           column(width = 3,
                                  selectInput("var", "Climate variable",
                                              choices = c("Stream temperature", "Low flow", "SST", "Salinity"),
                                              selected = "Stream temperature")
                           )
                         ), # end fluidRow
                         # input <- list(sp = "Chinook", cu = "Boundary Bay (Fall 4-1)", stage = "Adult freshwater migration", var = "Stream temperature")
                         # input <- list(sp = "Sockeye", cu = "Nahatlatch-Early Summer", stage = "Adult freshwater migration", var = "Stream temperature")
                         # input <- list(sp = "Coho", cu = "Boundary Bay", stage = "Adult freshwater migration", var = "Low flow")
                         
                         fluidRow( # Leaflet map
                           column(width = 7,
                                  leafletOutput("map", width = "100%", height = 600)
                           ),
                           column(width = 5,
                                  fluidRow(
                                    plotOutput(outputId = "overall", height = 200)),
                                  fluidRow(
                                    plotOutput(outputId = "timing", height = 200)),
                                  fluidRow(
                                    plotOutput(outputId = "stage_exposure", height = 200))
                           )
                         ) # end fluidRow
                ),# end tab DATA
                tabPanel("ABOUT",
                         br(),
                         HTML("<p>There is a pressing need to understand and account for the challenges (and opportunities) that climate change presents for salmon, and consider this information in the development of recovery and rebuilding plans. These quantitative assessments of future conditions that salmon will be exposed to can help us understand which Conservation Units will be most at risk and why.</p>"),
                         HTML("<h3>METHODS OVERVIEW</h3>"),
                         HTML("<p>This assessment considers exposure of 60 Conservation Units (CUs) in the Fraser river basin. Conservation Units are ecologically, geographically, and genetically unique groups of wild Pacific salmon defined under <a href = 'https://www.pac.dfo-mpo.gc.ca/fm-gp/salmon-saumon/wsp-pss/index-eng.html'>Canada's Wild Salmon Policy</a>, and are an appropriate unit for considering risks to biodiversity.</p><p> We considered climate change exposure of six life stages within each CU: incubation, freshwater rearing, early marine, marine rearing, adult freshwater migration, and spawning. Life stages were defined based on timing, spatial distribution, and sensitivity to climate variables. </p><p>Within each life stage, we assesed exposure to two climate variables: stream temperature and low stream flow for freshwater stages, and sea surface temperature and low salinity for marine stages. Exposure was quantified as the proportion of of the stage duration (for each CU) that was above or below a threshold value, based on the spatial distribution and timing of the life stage for that CU. Threshold values were unique for each climate variable:<ul>
                              <li>stream temperatures above species- and life-stage- specific thermal tolerance as defined by the Province of BC <a href = 'https://www2.gov.bc.ca/assets/gov/environment/air-land-water/water/waterquality/water-quality-guidelines/approved-wqgs/temperature-tech.pdf'>(Oliver and Fidler 2001)</a>; </li>
                              <li>stream flows below 20% historical mean annual discharge (MAD) for the location </li>
                              <li>sea surface temperatures (SSTs) above the current core thermal range of each species in the marine environment <a href = 'https://onlinelibrary.wiley.com/doi/abs/10.1111/faf.12825'>(Langan et al, 2024)</a> </li>
                              <li>sea surface salinity below the 2.5th percentile of historical values for the location.</li></ul></p><p>Exposure was calculated from Global Climate Model projections from the CMIP5 ensemble in four 30-year periods (historical, early century, mid century, and late century) under RCP 4.5 and RCP 8.5 emissions scenarios. <b>At this time, this app shows output for the mid-century period (2041-2070) under RCP 4.5 only.</b> Exposure was calculated independently for each of six GCMs, and the median value among GCMs is shown by the colour of the grid cell on the map.</p><p>Plots on the right show the median exposure for the CU (and, on the bottom, the life stage) among all grid cells across which the CU (or life stage) is distributed. Within a stage, exposure was averaged across the two climate variables (e.g., stream temperature and low flow for freshwater life stages) to yield the stage exposure (bottom-right plot). Exposure was further averaged across life stages to yield an overall exposure score for the CU (top-right plot).</p>"),
                         HTML("<h3>MORE INFORMATION</h3>"),
                         HTML("<p>Full details of the analysis are forthcoming in Peacock et al. (in. prep) and the associated supplemental information is available <a href='https://bookdown.org/salmonwatersheds/climate-change-exposure-supp/online-supplement.html'>online</a>. For more information, contact Stephanie Peacock (speacock@psf.ca).</p>"),
                         HTML("<h3>DATA SOURCES</h3>"),
                         HTML("<p>The climate exposure assessments rely on a number of data sources including:
                              <ul>
                              <li><a href='https://pacificclimate.org/data/gridded-hydrologic-model-output'>Gridded hydrologic model output</a> provided by the <a href='https://pacificclimate.org'>Pacific Climate Impacts Consortium</a> (Markus Schnorbus, pers. comm.)</li>
                              <li>Global Climate Model Projections from the CMIP5 ensemble, from <a href='https://psl.noaa.gov/ipcc/cmip5/'>NOAA's Climate Change Web Portal</a> (Jamie Scott, pers. comm).</li>
                              <li>Salmon distribution modelling, including spawning and rearing habitat for each species, from <a href='https://smnorris.github.io/bcfishpass/'> bcfishpass</a> (Simon Norris, pers. comm.)</li>
                              <li><a href='https://bookdown.org/salmonwatersheds/life-cycle-timing'>Life-history timing of Pacific salmon Conservation Units</a>, compiled by Wilson and Peacock (in press, CJFAS)</li>
                              <li>Spatial boundaries of Pacific salmon and steelhead Conservation Units, from PSF's <a href = 'https://data.salmonwatersheds.ca/result?datasetid=104'>Salmon Data Library</a></li>
                              </ul>
                              </p>"),
                HTML("<h3>ACKNOWLEDGEMENTS</h3>"),
                HTML("<p> This research was led by Stephanie J. Peacock (PSF), William Cheung (UBC), Brendan Connors (DFO), Lisa Crozier (NOAA), Sue Grant (DFO), Eric Hertz (PSF), Brian Hunt (UBC), Josephine Iacarella (DFO), Cory Lagasse (DFO), Dan Moore (UBC), Jon Moore (SFU), Francois Nicolas-Robinne (PSF), Marc Porter (PSF), Markus Schnorbus (PCIC), Samantha M. Wilson (SFU), and Katrina Connors (PSF).</p><p>We acknowledge the World Climate Research Programme’s Working Group on Coupled Modelling, which is responsible for CMIP, and we thank the climate modeling groups (listed <a href='https://bookdown.org/salmonwatersheds/climate-change-exposure-supp/online-supplement.html#2_Global_Climate_Models'>here</a>) for producing and making available their model output. For CMIP the U.S. Department of Energy’s Program for Climate Model Diagnosis and Intercomparison provides coordinating support and led development of software infrastructure in partnership with the Global Organization for Earth System Science Portals.</p><p>We thank Simon Norris for providing output from the bcfishpass model that was used to define spawning, rearing, and migration habitats for each CU. We thank Travis Tai and Diana Dobson for early input on this work through PSF’s Climate Science Advisory Committee. We also thank members of PSF’s Population Science Advisory Committee and Habitat Science Advisory Committee for their feedback.<p></p> This research was funded by the British Columbia Salmon Restoration and Innovation Fund, project #2020_279 led by the Pacific Salmon Foundation (PSF).</p>")
                ) #end tab ABOUT
    )# end tab panel
    
  ),
  
  ###############################################################################
  server = function(input, output) {
    
    #----------------------------------------------------------------------------
    # Limit list of CUs based on species selection
    #----------------------------------------------------------------------------
    limitCUs <- reactive({
      spp <- input$sp
      cus <- sort(cu_list$cu_name_pse[cu_list$species_pooled == spp])
      return(cus)
    })
    
    observe({
      updateSelectInput(inputId = "cu", label = "Conservation Unit:", choices = limitCUs() 
      )
    })
    
    #----------------------------------------------------------------------------
    # Limit list of climate variables based on life stage
    #----------------------------------------------------------------------------
    # c("Incubation", "Freshwater rearing", "Early Marine", "Marine Rearing", "Adult freshwater migration", "Spawning")
    # c("Stream temperature", "Low flow", "SST", "Salinity")
    limitVar <- reactive({
      stage <- input$stage
      varOpt <- list(c("Stream temperature", "Low flow"),
                     c("Stream temperature", "Low flow"),
                     c("SST", "Salinity"),
                     c("SST", "Salinity"),
                     c("Stream temperature", "Low flow"),
                     c("Stream temperature", "Low flow")
      )[[match(stage, stage_opt)]]
      return(varOpt)
    })
    
    observe({
      updateSelectInput(inputId = "var", label = "Climate variable", choices = limitVar())
    })
    
  
    #----------------------------------------------------------------------------
    # Plot map  
    #----------------------------------------------------------------------------
    
    output$map <- renderLeaflet({
      
      # Set selections
      i <- which(cu_list$cu_name_pse == input$cu & cu_list$species_pooled == input$sp)
      cu_boundary.i <- cu_boundary %>%
        filter(CUID == cu_list$cuid[i])

      j <- match(input$stage, stage_opt)
      
      k <- match(input$var, list(c("Stream temperature", "Low flow"),
                                c("Stream temperature", "Low flow"),
                                c("SST", "Salinity"),
                                c("SST", "Salinity"),
                                c("Stream temperature", "Low flow"),
                                c("Stream temperature", "Low flow")
      )[[match(input$stage, stage_opt)]])
       
     if(input$stage == "Marine rearing"){
       bounds1 <- as.numeric(st_bbox(grid_polys_mar))
     } else if(input$stage == "Early marine"){
       bounds1 <-  as.numeric(st_bbox(grid_polys_mar[grid_polys_mar$id %in% incl_EM,]))
     } else {
         # Mouth of Fraser is approx. 49.305188, -123.240498
      bounds1 <- as.numeric(st_bbox(cu_boundary.i))
      bounds1[1] <- min(bounds1[1], -123.240498)
      bounds1[2] <- min(bounds1[2], 49.305188)
     }
      
      # Select grid
      
      # Marine: highlight grid cells
      if(input$stage == "Early marine"){
        grid_polys.ij <- grid_polys_mar[grid_polys_mar$id %in% incl_EM, ]
      } else if(input$stage == "Marine rearing"){
        grid_polys.ij <- grid_polys_mar[grid_polys_mar$id %in% incl_MR, ]
      } else {
        j_FW <- match(input$stage, c("Adult freshwater migration", "Spawning", "Incubation", "Freshwater rearing")) # order of incl_FW
        grid_polys.ij <- grid_polys_fw[which(grid_polys_fw$id %in% incl_FW[[i, j_FW]]), ]
  
      }
      
      opacity.ijk <- ifelse(is.na(spat[[i,j]][k,]), 0, 0.8)
      
      leaflet() %>%
        addProviderTiles(providers$Esri.WorldTopoMap,
                         options = providerTileOptions(noWrap = FALSE)) %>%
        fitBounds(lng1 = bounds1[1], lat1 = bounds1[2], lng2 = bounds1[3], lat2 = bounds1[4])  %>%
        addLegend("topleft", pal = col_palette_fun, values = c(0,1),
                  title = "Exposure",
                  opacity = 1
        ) %>%
        addPolygons(data = grid_polys.ij, color = "#000000", fillColor = col_palette_fun(spat[[i,j]][k,]), weight = 1, opacity = 0.4, fillOpacity = opacity.ijk) %>%
        addPolygons(data = cu_boundary.i, color = "#000000", weight = 1.5, opacity = 0.8, fillOpacity = 0)
        
      
    })
    
    output$timing = renderPlot({
      # Set selections
      i <- which(cu_list$cu_name_pse == input$cu & cu_list$species_pooled == input$sp)
      j <- match(input$stage, stage_opt)
      
      timing.i <- timing[timing$cuid == cu_list$cuid[i], ]
      
      par(mar = c(2,0,4,5), family = "Trebuchet MS")
      plot(range(xDate), c(1,6), "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
      mtext(side = 3, "Life-history timing for selected CU", line = 1, cex = 1.3)
      axis(side = 1, at = labDate, labels = labName)
      abline(v = labDate, lty = 3, col = grey(0.8))
      for(jj in 1:6){
        col_stage <- ifelse(j == jj, "#000000", "#00000040")
        text(xDate[365] + 7, c(6:1)[jj], stage_disp[jj], xpd = NA, adj = 0, cex = 0.9, col = col_stage)
        
        timing.j <- timing.i %>% filter(stage == stages[[jj]])
        
        if(timing.j$duration > 0){
          
        segments(x0 = xDate[timing.j$start],
                 y0 = 6 - jj + 1,
                 x1 =xDate[min(365, timing.j$start + timing.j$duration)],
                 col = col_stage, lwd = 8)
        
        
        days_left <- timing.j$duration - (365 - timing.j$start)
        lower_line <- 0.2
        while(days_left > 0){
          segments(x0 = xDate[1],
                   y0 = 6 - jj + 1 - lower_line,
                   x1 =xDate[min(365, days_left)],
                   col = col_stage, lwd = 8)
          days_left <- days_left - 365
          lower_line <- lower_line + 0.2
        } # end while
        if(j == jj){
          
          points(xDate[timing.j$start], 6 - jj + 1, pch = 21, bg = "#83B687", cex = 1.2, xpd = NA)
          text(xDate[timing.j$start], 6 - jj + 1, strftime(xDate[timing.j$start], format = "%b %d"), pos = 2, cex = 0.9, xpd = NA, col = "#83B687")
          points(xDate[days_left + 365], 6 - jj + 1 - lower_line + 0.2, pch = 21, bg = "#C06363", cex = 1.2, xpd = NA)
          text(xDate[days_left + 365], 6 - jj + 1 - lower_line + 0.2, strftime(xDate[days_left + 365], format = "%b %d"), pos = 4, cex = 0.9, xpd = NA, col = "#C06363")
        }
        } # end if duration > 0
      } # end jj
      
    })
    
    output$overall = renderPlot({
      
      # Set selections
      i <- which(cu_list$cu_name_pse == input$cu & cu_list$species_pooled == input$sp)
      species_num <- as.numeric(as.factor(cu_list$species_pooled))
      par(mar = c(2,0,4,5), family = "Trebuchet MS")
      plot(c(0,0.5), c(-0.32,0.22), "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
      mtext(side = 3, "Overall exposure", line = 1, cex = 1.3)
      abline(v = seq(0, 1, 0.1), lty = 3, col = grey(0.8))
      abline(h = c(0.2, 0.1, 0, -0.1, -0.2, -0.3), lty = 3, col = grey(0.8))
      text(0.52, c(0.2, 0.1, 0, -0.1, -0.2, -0.3), c("Chinook", "Chum", "Coho", "Pink", "Sockeye", "Steelhead"), xpd = NA, adj = 0, cex = 0.9, col = ifelse(c("Chinook", "Chum", "Coho", "Pink", "Sockeye", "Steelhead") == input$sp, "#000000", "#00000040"))
      points(overall_exposure$x, c(0.2, 0.1, 0, -0.1, -0.2, -0.3)[species_num], pch = species_num, col = col_palette_fun(overall_exposure$x), lwd = 1.2, cex = 1.5)
      points(overall_exposure$x[i], c(0.2, 0.1, 0, -0.1, -0.2, -0.3)[species_num[i]], pch = species_num[i], lwd = 2, cex= 2)
      text(overall_exposure$x[i], c(0.2, 0.1, 0, -0.1, -0.2, -0.3)[species_num[i]], "Selected CU", pos = 3, cex = 0.9, xpd = NA)
    })
    
    
    output$stage_exposure = renderPlot({
      
      # Set selections
      i <- which(cu_list$cu_name_pse == input$cu & cu_list$species_pooled == input$sp)
      j <- match(input$stage, stage_opt)
      
      # stage_exposure_i <- stage_exposure %>% filter(cuid == cu_list$cuid[i])
      # 
      # par(mar = c(2,0,4,5), family = "Trebuchet MS")
      # plot(c(0,1), c(0.75,6.25), "n", bty = "n", yaxt = "n", xlab = "", ylab = "")
      # mtext(side = 3, "Stage exposure for selected CU", line = 1, cex = 1.3)
      # abline(v = seq(0, 1, 0.1), lty = 3, col = grey(0.8))
      # abline(h = c(1:6), lty = 3, col = grey(0.8))
      # 
      # text(1.02, c(6:1), stage_disp, xpd = NA, adj = 0, cex = 0.9, col = ifelse(stage_opt == input$stage, "#000000", "#00000040"))
      # points(stage_exposure_i$p_j_median[match(stages, stage_exposure_i$stage)], c(6:1), col = col_palette_fun(stage_exposure_i$p_j_median[match(stages, stage_exposure_i$stage)]), pch = 19, cex = 3)
      # points(stage_exposure_i$p_j_median[match(stages, stage_exposure_i$stage)][j], c(6:1)[j], cex = 3, lwd = 1.5)
      # text(stage_exposure_i$p_j_median[match(stages, stage_exposure_i$stage)][j], c(6:1)[j], pos = 3, "Selected stage", xpd = NA, cex =0.9)
      
      # Alternative
      k <- match(input$var, list(c("Stream temperature", "Low flow"),
                                 c("Stream temperature", "Low flow"),
                                 c("SST", "Salinity"),
                                 c("SST", "Salinity"),
                                 c("Stream temperature", "Low flow"),
                                 c("Stream temperature", "Low flow")
      )[[match(input$stage, stage_opt)]])
      
      V <- list(c("stream_temp", "stream_flow"),
                c("stream_temp", "stream_flow"),
                c("SST", "SSS"),
                c("SST", "SSS"),
                c("stream_temp", "stream_flow"),
                c("stream_temp", "stream_flow"))
      
      stage_exposure_i <- stagevar_exposure %>% filter(cuid == cu_list$cuid[i])
      
      par(mar = c(4,0,4,0), family = "Trebuchet MS")
      plot(c(0,6), c(0,1.2), "n", bty = "n", yaxt = "n", xlab = "", ylab = "", xaxt = "n")
      text(c(1:6)-0.5, rep(1.1, 6), stage_disp, xpd = NA, cex = 0.9, col = ifelse(stage_opt == input$stage, "#000000", "#00000040"))
      mtext(side = 3, "Stage exposure for selected CU", line = 1, cex = 1.3)
      mtext(side = 3, "median among grid cells", line = 0.4, cex = 0.9)
      
      for(jj in 1:6){
        stage_exposure_ijj <- stage_exposure_i %>% filter(stage == stages[jj])
        for(v in 1:2){
          polygon(x = jj - list(rep(c(1, 0.5), each = 2), rep(c(0.5, 0), each = 2))[[v]],
                  y = c(0.2, 0.8, 0.8, 0.2),
                  col = col_palette_fun(stage_exposure_ijj$exp_median[match(V[[jj]][v], stage_exposure_ijj$variable)]))
          text(jj - c(0.75, 0.25)[v], 0.5, sprintf(fmt = "%.2f", stage_exposure_ijj$exp_median[match(V[[jj]][v], stage_exposure_ijj$variable)]), cex = 0.9)
          } # end v 
      } # end jj
      
      # Highlight selected
      polygon(x = j - list(rep(c(1, 0.5), each = 2), rep(c(0.5, 0), each = 2))[[k]],
              y = c(-0.5, 0.8, 0.8, -0.5),
              col = NA, lwd = 5, xpd = NA)
      
      segments(x0 = c(0:6), y0 = 0.2, x1 = c(0:6), y1 = 1.2, lwd = 1.2, xpd = NA)
      segments(x0 = seq(0, 6, 0.5), y0 = 0.8, x1 = seq(0, 6, 0.5), y1 = -0.3, xpd = NA, lty = 3)
      text(seq(0, 5.5, 0.5) + 0.25, y = rep(0.1, 12), srt = 90, c(rep(c("stream\ntemp", "low\nflow"), 2), rep(c("SST", "salinity"), 2), rep(c("stream\ntemp", "low\nflow"), 2)), xpd = NA, adj = 1, cex = 0.9)
    })
    
  },
  options = list(height = 700)
)
