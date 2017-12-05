library(shiny)
library(shinyjs)
library(shinyBS)
library(shinythemes)

# Start of Shiny HTML5/JS page
shinyUI(navbarPage("UmbrellaEffect",
                   theme = shinytheme("flatly"),
position = "fixed-top",
collapsible = TRUE,
                   # theme = shinytheme("spacelab"),
                   # theme = shinytheme("sandstone"),
                   # theme = shinytheme("yeti"),
## model input data tab  
tabPanel("Model input data", align = "center", #style = "vertical-align: top;",
tags$style(type="text/css", "body {padding-top: 70px;}"),
    wellPanel(
    fluidRow(
    column(6,
    h3("Program settings"),
      # app directories
      fluidRow(
      column(6,
      textInput("inputDir", label = "Working directory input", value = "~/temp/UM/")),
      column(6,
      textInput("outputDir", label = "Working directory output", value = "~/temp/UM/"))
      ),
## tooltip
bsTooltip(id = 'inputDir', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'outputDir', title = "", placement = "top", trigger = "hover"),
## end tooltips
      fluidRow(
        column(6,
    h4("Columns with spatial data", align = "center")),
        column(4,
    h4("Data column", align = "center")),
        column(2,
    h4("Output", align = "center"))
      ),
      fluidRow(
      column(2,
      numericInput("spat_col_x", label = "x:", value = NA)),
      column(2,
      numericInput("spat_col_y", label = "y:", value = NA)),
      column(2,
      numericInput("spat_col_z", label = "z:", value = 2)),
      column(2,
      numericInput("data_col", label = "for csv:", value = 3)),
      column(2,
      numericInput("data_tsf", label = "for tsf:", value = 13)),
      column(2,
      radioButtons('outputType', 'filetype',
       c(csv='csv',
         Rdata='.rdata'),
       '.rdata'))
      ),
## tooltip
bsTooltip(id = 'spatial_col_x', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'spatial_col_y', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'spatial_col_z', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'data_col', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'data_tsf', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'outputType', title = "", placement = "top", trigger = "hover"),
## end tooltips
    h3("Gravimeter", align = "center"),
      fluidRow(
        column(3,
        numericInput("grav_x", label = "Location x:", value = 3)),
        column(3,
        numericInput("grav_y", label = "Location y:", value = 3)),
        column(3,
        numericInput("grav_z", label = "Location z:", value = 0)),
        column(3,
        numericInput("grav_height", label = "Sensor height:", value = 1.05))
      ),
## tooltip
bsTooltip(id = 'grav_x', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'grav_y', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'grav_z', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'grav_height', title = "", placement = "top", trigger = "hover"),
## end tooltips
    h3("Input files", align = "center"),
      fluidRow(
      column(4,
      textInput("inputFile_DEM", label = "DEM", value = "WE_UP_TO_300m_05m.asc")),
      column(4,
      textInput("inputFile_gObs", label = "Gravity observations", value = "SG030_TS_1month.tsf")),
      column(4,
      textInput("inputFile_SM", label = "Soil moisture time series", value = "SMdata_TS_1d.rData"))
      ),
## tooltip
bsTooltip(id = 'inputFile_DEM', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'inputFile_gObs', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'inputFile_SM', title = "", placement = "top", trigger = "hover"),
## end tooltips
      h3("Model domain", align = "center"),
      fluidRow(
        column(3,
    h4("")),
        column(3,
    h4("x coordinate", align = "center")),
        column(3,
    h4("y coordinate", align = "center")),
        column(3,
    h4("z coordinate", align = "center"))
      ),
      fluidRow(
        column(3,
    h4("minimal", align = "center")),
        column(3,
        numericInput("modelDomain_xmin", label = "", value = 0)),
        column(3,
        numericInput("modelDomain_ymin", label = "", value = 0)),
        column(3,
        numericInput("modelDepth_min", label = "", value = -3))
      ),
      fluidRow(
        column(3,
    h4("maximal", align = "center")),
        column(3,
        numericInput("modelDomain_xmax", label = "", value = 6)),
        column(3,
        numericInput("modelDomain_ymax", label = "", value = 6)),
        column(3,
        numericInput("modelDepth_max", label = "", value = 0))
      ),
      fluidRow(
        column(3,
    h4("discretization", align = "center")),
      column(3,
      numericInput("modelDiscr_x", label = "", value = 0.5)),
      column(3,
      numericInput("modelDiscr_y", label = "", value = 0.5)),
      column(3,
      numericInput("modelDiscr_z", label = "", value = 0.5))
      ),
## tooltip
bsTooltip(id = 'modelDomain_xmin', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'modelDomain_xmax', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'modelDomain_ymin', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'modelDomain_ymax', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'modelDiscr_x', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'modelDiscr_y', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'modelDiscr_z', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'modelDepth_min', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'modelDepth_max', title = "", placement = "top", trigger = "hover")
## end tooltips
    ), # end left input data column
    ## input data: specific umbrella settings
    column(6,
    h3("Umbrella settings", align = "center"),
      fluidRow(
      column(2,
        h4("")),
      column(6,
        h4("Building walls")),
      column(4,
        h4("Building baseplate"))
      ),
      fluidRow(
      column(2),
      column(2,
      numericInput("building_walls_x", label = "x", value = 0.6)),
      column(2,
      numericInput("building_walls_y", label = "y", value = 0.6)),
      column(2,
      numericInput("building_walls_z", label = "z", value = 1.5)),
      column(2,
      numericInput("building_baseplate_z_min", label = "z min", value = -0.5)),
      column(2,
      numericInput("building_baseplate_z_max", label = "z max", value = 0))
      ),
## tooltip
bsTooltip(id = 'pillar_r', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'pillar_d', title = "", placement = "top", trigger = "hover"),
## end tooltips
      fluidRow(
      column(7, align = "right",
        h4("Pillar correction with:")),
      column(5, align = "left",
      radioButtons('pillar_correction_type', '',
       c(rectangular='rectangular',
         circle='circle'),
       'rectangular'))
      ),
      fluidRow(
      column(2,
        h4("")),
      column(6,
        h4("rectangle")),
      column(4,
        h4("circle"))
      ),
      fluidRow(
      column(2),
      column(2,
      numericInput("pillar_x_min", label = "x min", value = 2)),
      column(2,
      numericInput("pillar_y_min", label = "y min", value = 2)),
      column(2,
      numericInput("pillar_z_min", label = "z min", value = -1)),
      column(2,
      numericInput("pillar_r", label = "Radius:", value = 1)),
      column(2,
      numericInput("pillar_d", label = "Depth below surf:", value = 1.2))
      ),
      fluidRow(
      column(2),
      column(2,
      numericInput("pillar_x_max", label = "x max", value = 4)),
      column(2,
      numericInput("pillar_y_max", label = "y max", value = 4)),
      column(2,
      numericInput("pillar_z_max", label = "z max", value = 0)),
      column(4)
      ),
    h3("Surounding conditions"),
      fluidRow(
        column(6,
          h4("SG position")),
        column(6,
          selectInput("SGposition", "",
                      c("Center" = "Center",
                        "Corner" = "Corner",
                        "Side" = "Side",
                        "Trans" = "Trans",
                        "Wettzell" = "Wettzell"),
                      selected = "Center"))
      ),
      fluidRow(
        column(6,
          h4("Hydrological site conditon")),
        column(6,
          selectInput("Hydro_condition", "",
                      c("Soil type sandy loam" = "Soil type sandy loam",
                        "Soil type silty loam" = "Soil type silty loam",
                        "Soil type clay loam" = "Soil type clay loam",
                        "Anisotropy sandy loam" = "Anisotropy sandy loam",
                        "Anisotropy silty loam" = "Anisotropy sandy loam",
                        "Anisotropy clay loam" = "Anisotropy sandy loam",
                        "Climate AI=0.35" = "Climate AI=0.35",
                        "Climate AI=0.575" = "Climate AI=0.575",
                        "Climate AI=1.0" = "Climate AI=1.0",
                        "Climate AI=2.0" = "Climate AI=2.0"),
                      selected = "Soil type sandy loam"))
      ),
      fluidRow(
        column(6,
          h4("Vertical extension")),
        column(6,
      numericInput("Vertical_extent", label = "", value = 5))
      )
    )
    ) # end model parameters
    ),#end wellPanel
#     fluidRow(
#         column(5),
#         column(2,
# ## run model button
# actionButton(
#         inputId = "run_model",
#         label = "Run infiltration model"
#       )),
#         column(5)
#     ),
    fluidRow()
# ) # end fluidRow around wellPanel
# ))
#
), #end tap:Model input
tabPanel("RUN model", align = "center",
    ## trigger model run
    fluidRow(
        column(3,
        ## run model button
        actionButton(
                inputId = "run_reduction",
                label = "Reduce time series for umbrella effect"
              )),
        column(9,
            shinyjs::useShinyjs(),
            div(style = "height:100px; overflow-y: scroll", verbatimTextOutput("text"))
               )
    ),
    fluidRow(
        column(4,
            h4("Surface grid")),
        column(4,
            h4("Gravity component grid (from above)")),
        column(4,
            h4("Gravity component grid (vertically)"))
      ),
    fluidRow(
        column(4,
            plotOutput(outputId = "surfaceGrid", width = "95%")),
        column(4,
            plotOutput(outputId = "gCompGrid", width = "95%")),
        column(4,
            plotOutput(outputId = "gCompGrid_PillarRemoved", width = "95%"))
            ),
    # fluidRow(
    #          h4("Optimized parameters: "),
    #         tableOutput("parameter_output")
    #         
    #     ),
fluidRow(
        column(12,
            h4("Gravity time series of input, calculations and reduction"))#,
        # column(4,
        #     h4("Soil moisture over depth and time")),
        # column(4,
        #     h4("Soil moisture transect at different times"))
        ),
fluidRow(
        column(12,
            plotOutput(outputId = "plot_all_gravity_data"))#,
        # column(4,
        #     plotOutput(outputId = "plot_soilMoisture_run")),
        # column(4,
        #     plotOutput(outputId = "plot_soilMoisture_transect_run"))
        )
),
tabPanel("Model results", align = "center",
    wellPanel(
    fluidRow(
    column(8,
      h3("Plot settings", align = "center"),
      # app directories
      fluidRow(
      column(3,
      textInput("plot_inputDir", label = "Working directory input", value = "~/temp/UM/")),
      column(3,
      textInput("plot_outputDir", label = "Working directory output", value = "~/temp/UM/")),
      column(4,
      textInput("plot_inputFile_gObs", label = "Gravity observations", value = "SG030_TS_1month.tsf")),
      column(2,
      numericInput("plot_data_tsf", label = "for tsf:", value = 13))
      )
      ),
      column(4, 
         actionButton(
             inputId = "plot_results",
             label = "Read output and plot it"))
      )
    ),
fluidRow(
        column(12,
            h4("Gravity time series of input, calculations and reduction"))
        ),
fluidRow(
        column(12,
            plotOutput(outputId = "plot_all_gravity_data_results"))
        )
) # end tab: Model results

)# end navbarPage
)# end shinyUI
