library(shiny)
library(shinyjs)
library(shinyBS)
library(shinythemes)

# Start of Shiny HTML5/JS page
shinyUI(navbarPage("GravityInfiltraitonProcessExplorer",
                   theme = shinytheme("flatly"),
                   # theme = shinytheme("spacelab"),
                   # theme = shinytheme("sandstone"),
                   # theme = shinytheme("yeti"),
  
  # titlePanel("gravity infiltration process detection"), # App title
tabPanel("Model input data",
## introduction, explication what to do here..
# fluidRow(
#     h3("Explanation")
#     ),
# fluidRow(
# column(12,
## input data
	wellPanel(
    fluidRow(
    column(6,
	h3("Program settings"),
      # app directories
      fluidRow(
      column(6,
      textInput("inputDir", label = "Working directory input", value = "~/temp/GI/")),
      column(6,
      textInput("outputDir", label = "Working directory output", value = "~/temp/GI/"))
      ),
## tooltip
bsTooltip(id = 'inputDir', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'outputDir', title = "", placement = "top", trigger = "hover"),
## end tooltips
      fluidRow(
        column(6,
    h4("Columns with spatial data")),
        column(4,
    h4("Data column")),
        column(2,
    h4("Output"))
      ),
      fluidRow(
      column(2,
      numericInput("spat_col_x", label = "x:", value = 1)),
      column(2,
      numericInput("spat_col_y", label = "y:", value = 2)),
      column(2,
      numericInput("spat_col_z", label = "z:", value = NA)),
      column(2,
      numericInput("data_col", label = "for csv:", value = 3)),
      column(2,
      numericInput("data_tsf", label = "for tsf:", value = 7)),
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
	h3("Site parameters"),
      fluidRow(
      column(8,
	h4("Gravimeter")),
      column(4,
	h4("Pillar"))
      ),
      fluidRow(
      column(2,
      numericInput("grav_x", label = "Location x:", value = 0)),
      column(2,
      numericInput("grav_y", label = "Location y:", value = 0)),
      column(2,
      numericInput("grav_z", label = "Location z:", value = 0)),
      column(2,
      numericInput("grav_height", label = "Sensor height:", value = 0)),
      column(2,
      numericInput("pillar_r", label = "Radius:", value = 0)),
      column(2,
      numericInput("pillar_d", label = "Depth below surf:", value = 0))
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
## tooltip
bsTooltip(id = 'pillar_r', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'pillar_d', title = "", placement = "top", trigger = "hover"),
## end tooltips
	h3("Input files"),
      fluidRow(
      column(6,
      # fileInput("inputFile_DEM", label = "DEM")),
      textInput("inputFile_DEM", label = "DEM", value = "")),
      column(6,
      # fileInput("inputFile_gObs", label = "Gravity observations"))
      textInput("inputFile_gObs", label = "Gravity observations", value = "iGrav006_obs_60sec.tsf"))
      ),
      fluidRow(
      column(6,
      # fileInput("inputFile_IntDistr", label = "Intensitry distribution")),
      textInput("inputFile_IntDistr", label = "Intensitry distribution", value = "waterIntensity_measured.rData")),
      column(3,
	  radioButtons('IntpMethod', 'Interpolation method',
	   c(iwd='IDW',
	     kriging='Kriging'),
	   'IDW')),
      column(3,
      numericInput("ZeroBorderDensity", label = "Zero fillup borders", value = 0.2))
      ),
## tooltip
bsTooltip(id = 'inputFile_DEM', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'inputFile_gObs', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'inputFile_IntDistr', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'IntpMethod', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'ZerosBorderDensity', title = "", placement = "top", trigger = "hover"),
## end tooltips
      h3("Model domain"),
      fluidRow(
        column(3,
    h4("")),
        column(3,
    h4("x coordinate")),
        column(3,
    h4("y coordinate")),
        column(3,
    h4("z coordinate"))
      ),
      fluidRow(
        column(3,
    h4("minimal")),
        column(3,
        numericInput("modelDomain_xmin", label = "", value = -7.5)),
        column(3,
        numericInput("modelDomain_ymin", label = "", value = -7.5)),
        column(3,
        numericInput("modelDepth_min", label = "", value = -3))
      ),
      fluidRow(
        column(3,
    h4("maximal")),
        column(3,
        numericInput("modelDomain_xmax", label = "", value = 7.5)),
        column(3,
        numericInput("modelDomain_ymax", label = "", value = 7.5)),
        column(3,
        numericInput("modelDepth_max", label = "", value = 0))
      ),
      fluidRow(
        column(3,
    h4("discretization")),
      column(3,
      numericInput("modelDiscr_x", label = "", value = 0.25)),
      column(3,
      numericInput("modelDiscr_y", label = "", value = 0.25)),
      column(3,
      numericInput("modelDiscr_z", label = "", value = 0.1))
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
    ## model input data and parametrization
    column(6,
    h3("Model settings"),
      fluidRow(
        column(3,
	    radioButtons('Modeling_mode', 'Model run method',
	     c(conversion='conversion',
	      inverse='inverse'),
	     'inverse')),
        column(3,
	    radioButtons('del_prevModRuns', 'Delete previous runs?',
	     c(yes='yes',
	      no='no'),
	     'yes')),
        column(3,
        numericInput("model_iterations", label = "Number of iterations:", value = 10)),
        column(3,
        numericInput("mb_permitted_error", label = "Max error for mass balance:", value = 0.05))
      ),
## tooltips
bsTooltip(id = 'Modeling_mode', title = "this is just a test [min]", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'del_prevModRuns', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'model_iterations', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'mb_permitted_error', title = "", placement = "top", trigger = "hover"),
## end tooltips
      fluidRow(
        column(3,
        numericInput("precip_time", label = "Sprinkling time:", value = 10)),
        column(3,
        numericInput("water_vol_min", label = "Water volume per timestep:", value = 0.035)),
        column(3,
        numericInput("plot_interval", label = "Plot interval:", value = 60)),
        column(3,
        numericInput("plot_transect_loc", label = "Plot transect location:", value = 0))
      ),
## tooltips
bsTooltip(id = 'precip_time', title = "this is just a test [min]", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'water_vol_min', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'plot_interval', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'plot_transect_loc', title = "", placement = "top", trigger = "hover"),
## end tooltips


      fluidRow(
        column(3,
        numericInput("macroporeLayers", label = "Number of macro pore layers", value = 2))
		# radioButtons('macroporeLayers', 'Number of macro pore layers',
		#  c(0='0',
        #    1='1',
		#    2='2'),
		#    '2'))
        ),



    # checkboxInput('macrolayers', 'Use macro pore layers?', TRUE)
    # checkboxInput('2macrolayers', '', TRUE)

    h3("Hydrological model parameters"),
      fluidRow(
        column(4,
    h4("")),
        column(2,
    h4("minimal")),
        column(2,
    h4("maximal")),
        column(2,
    h4("start value"))
      ),
      fluidRow(
        column(4,
        h5("Theta macro layer 1 [%]")),
        column(2,
        numericInput("theta_macro_min", label = "", value = 0.05)),
        column(2,
        numericInput("theta_macro_max", label = "", value = 0.15)),
        column(2,
        numericInput("theta_macro_start", label = "", value = 0.1))
      ),
## tooltip
bsTooltip(id = 'theta_macro_min', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'theta_macro_max', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'theta_macro_start', title = "", placement = "top", trigger = "hover"),
## end tooltips
      fluidRow(
        column(4,
        h5("Theta macro layer 2 [%]")),
        column(2,
        numericInput("theta_macro2_min", label = "", value = 0.05)),
        column(2,
        numericInput("theta_macro2_max", label = "", value = 0.15)),
        column(2,
        numericInput("theta_macro2_start", label = "", value = 0.1))
      ),
## tooltip
bsTooltip(id = 'theta_macro2_min', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'theta_macro2_max', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'theta_macro2_start', title = "", placement = "top", trigger = "hover"),
## end tooltips
      fluidRow(
        column(4,
        h5("Theta other process [%]")),
        column(2,
        numericInput("theta_other_min", label = "", value = 0.05)),
        column(2,
        numericInput("theta_other_max", label = "", value = 0.15)),
        column(2,
        numericInput("theta_other_start", label = "", value = 0.1))
      ),
## tooltip
bsTooltip(id = 'theta_other_min', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'theta_other_max', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'theta_other_start', title = "", placement = "top", trigger = "hover"),
## end tooltips
      fluidRow(
        column(4,
        h5("Depth macro layer 1 [m]")),
        column(2,
        numericInput("depth_macro_min", label = "", value = -0.5)),
        column(2,
        numericInput("depth_macro_max", label = "", value = -0.1)),
        column(2,
        numericInput("depth_macro_start", label = "", value = -0.25))
      ),
## tooltip
bsTooltip(id = 'depth_macro_min', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'depth_macro_max', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'depth_macro_start', title = "", placement = "top", trigger = "hover"),
## end tooltips
      fluidRow(
        column(4,
        h5("Depth macro layer 2 [m]")),
        column(2,
        numericInput("depth_macro2_min", label = "", value = -1.5)),
        column(2,
        numericInput("depth_macro2_max", label = "", value = -0.3)),
        column(2,
        numericInput("depth_macro2_start", label = "", value = -1.0))
      ),
## tooltip
bsTooltip(id = 'depth_macro2_min', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'depth_macro2_max', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'depth_macro2_start', title = "", placement = "top", trigger = "hover"),
## end tooltips
      fluidRow(
        column(4,
        h5("Depth other process [m]")),
        column(2,
        numericInput("depth_other_min", label = "", value = -2.5)),
        column(2,
        numericInput("depth_other_max", label = "", value = -0.1)),
        column(2,
        numericInput("depth_other_start", label = "", value = -1.5))
      ),
## tooltip
bsTooltip(id = 'depth_other_min', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'depth_other_max', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'depth_other_start', title = "", placement = "top", trigger = "hover"),
## end tooltips
      fluidRow(
        column(4,
        h5("Infiltration process")),
        column(2,
        numericInput("inf_dynamics_min", label = "", value = 1)),
        column(2,
        numericInput("inf_dynamics_max", label = "", value = 3)),
        column(2,
        numericInput("inf_dynamics_start", label = "", value = 1))
      ),
## tooltip
bsTooltip(id = 'inf_dynamic_min', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'inf_dynamic_max', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'inf_dynamic_start', title = "", placement = "top", trigger = "hover"),
## end tooltips
      fluidRow(
        column(4,
        h5("Lateral flow factor [-]")),
        column(2,
        numericInput("latflow_fac_min", label = "", value = 0)),
        column(2,
        numericInput("latflow_fac_max", label = "", value = 1)),
        column(2,
        numericInput("latflow_fac_start", label = "", value = 0.5))
      ),
## tooltip
bsTooltip(id = 'latflow_fac_min', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'latflow_fac_max', title = "", placement = "top", trigger = "hover"),
## end tooltips
## tooltip
bsTooltip(id = 'latflow_fac_start', title = "", placement = "top", trigger = "hover")
## end tooltips

    ) # end model parameters
    )
	),#end wellPanel
    fluidRow(
        column(5),
        column(2,
## run model button
actionButton(
        inputId = "run_model",
        label = "Run infiltration model"
      )),
        column(5)
    ),
    fluidRow()
# ) # end fluidRow around wellPanel
# ))
##
), #end tap:Model input
tabPanel("Model runtime output", 
    # run model button
    # fluidRow(
    #     column(2),
    #     column(2,
    #     actionButton(
    #         inputId = "run_model",
    #         label = "Run infiltration model"
    #     )),
    #     column(8)
    # ),
    # console output
    fluidRow(
        column(8,

            verbatimTextOutput("messages")
            # verbatimTextOutput("console_output")
            
        ),
        column(4,
            h3("Corresponding figures"),
            plotOutput(outputId = "test", width = "90%", height = "90%"),
            plotOutput(outputId = "testgrid", width = "90%", height = "90%"),
            plotOutput(outputId = "tt", width = "90%", height = "90%")
        )
    )
),
tabPanel("Model results", 
fluidRow(
h3("Result figures"),
plotOutput(outputId = "plot_gravity_responses", width = "90%", height = "90%")

# plotOutput(outputId, width = "100%", height = "400px", click = NULL,
#   dblclick = NULL, hover = NULL, hoverDelay = NULL,
#   hoverDelayType = NULL, brush = NULL, clickId = NULL, hoverId = NULL,
#   inline = FALSE)
)
) # end tab: Model results

)# end navbarPage
)# end shinyUI
