######################################
## Infiltration model: shiny server ##
######################################

###################
## Loading libraries
###################
library(shiny)
library(shinyjs)
#libraries
library(devtools)
load_all("/home/mreich/Dokumente/written/ResearchRepos/UmbrellaEffect")
load_all("/home/mreich/Dokumente/written/ResearchRepos/gravityInf")
# library(UmbrellaEffect)
# library(gravityInf)
# other Rpackages necessary
library(zoo); Sys.setenv(TZ = "GMT")
library(xts)
library(dplyr)
library(raster)
# library(UmbrellaEffect)
library(reshape2)
library(ggplot2)
library(viridis)
library(gstat)
library(ptinpoly)
#
library(plot3D)
library(grid)
library(gridExtra)
library(scales)
library(hydroGOF)
library(data.table)
library(ppso)
# for kriging
library(spacetime)
library(sp)
library(automap)
library(akima)

###################


## start shiny server
shinyServer(function(input, output) {
###################

# print_input = reactive({
# rbind(
#     dir_input = input$inputDir,
#     dir_output = input$outputDir,
#     output_type = input$outputType,
#     SG_x = input$grav_x,
#     SG_y = input$grav_y,
#     SG_Z = input$grav_z,
#     SG_SensorHeight = input$grav_height,
#     modelDomain_xmin = input$modelDomain_xmin,
#     modelDomain_xmax = input$modelDomain_xmax,
#     modelDomain_ymin = input$modelDomain_ymin,
#     modelDomain_ymax = input$modelDomain_ymax,
#     modelDiscr_x = input$modelDiscr_x,
#     modelDiscr_y = input$modelDiscr_y,
#     modelDiscr_z = input$modelDiscr_z,
#     modelDepth_min = input$modelDepth_min,
#     modelDepth_max = input$modelDepth_max,
#     thres_radius = input$pillar_r,
#     spatial_col_x = input$spat_col_x,
#     spatial_col_y = input$spat_col_y,
#     spatial_col_z = input$spat_col_z,
#     data_col = input$data_col,
#     data_tsf = input$data_tsf,
#     DEM_input_file = input$inputFile_DEM,
#     IntensityDistribution_file = input$inputFile_IntDistr,
#     interpolation_method = input$IntpMethod,
#     Zeros_border_density = input$ZeroBorderDensity,
#     gravityObservations_file = input$inputFile_gObs,
#     ModelingMode = input$Modeling_mode,
#     model_runs = input$model_iterations,
#     del_prevModRuns = input$del_prevModRuns,
#     number_macroLayers = input$macroporeLayers,
#     inf_dynamics_min = input$inf_dynamics_min,
#     precip_time = input$precip_time,
#     water_vol_min = input$water_vol_min,
#     mb_permitted_error = input$mb_permitted_error,
#     dtheta_macro_min = input$theta_macro_min,
#     dtheta_macro_max = input$theta_macro_max,
#     dtheta_macro2_min = input$theta_macro2_min,
#     dtheta_macro2_max = input$theta_macro2_max,
#     mdepth_min = input$depth_macro_min,
#     mdepth_max = input$depth_macro_max,
#     mdepth2_min = input$depth_macro2_min,
#     mdepth2_max = input$depth_macro2_max,
#     dtheta_other_min = input$theta_other_min,
#     dtheta_other_max = input$theta_other_max,
#     pdepth_min = input$depth_other_min,
#     pdepth_max = input$depth_other_max,
#     latflow_fac_min = input$latflow_fac_min,
#     latflow_fac_max = input$latflow_fac_max,
#     dtheta_macro_start = input$theta_macro_start,
#     dtheta_macro2_start = input$theta_macro2_start,
#     mdepth_start = input$depth_macro_start,
#     mdepth2_start = input$depth_macro2_start,
#     dtheta_other_start = input$theta_other_start,
#     latflow_fac_start = input$latflow_fac_start,
#     inf_dynamics_start = input$inf_dynamics_start,
#     pdepth_start = input$depth_other_start,
#     plot_interval = input$plot_interval,
#     plot_transect_loc = input$plot_transect_loc
# )
# })

####################
## model run EXECUTION CODE
####################
observeEvent(
      eventExpr = input[["run_model"]],
      handlerExpr = {
## start code, when "run model"-button is pressed
## tracking R console output and pass it to shiny text output
withCallingHandlers({
        shinyjs::html("text", "")

###################
## SETUP: assigning values
###################
## Output and input settings
# Directory
# path should be absolute
# (if not, it will be relative to the packages library path)
# use "test-data" for dir_input to use supplied example files within the package
# dir_input = "test-data"
# assign_inputDir = reactive({
#   input$inputDir
# })
dir_input = input$inputDir
dir_output = input$outputDir
# Output file type
# set to "csv", if output should also be saved as .csv (besides .rData)
if(input$outputType == 'csv'){ output_type = "csv"
}else{
output_type = ""
}
# Plotting option: should a plot of all time series be shown (and saved) in the end?
plot_data = TRUE
#
## Gravimeter location
# in [m]
# relativ, local coordinate sytem
SG_x = input$grav_x
SG_y = input$grav_y
SG_Z = input$grav_z
SG_SensorHeight = input$grav_height
# UTM coordinate system
# SG_x = 4564082.00
# SG_y = 5445669.70
# SG_Z = 609.755
# SG_SensorHeight = 1.5 
#
## Model domain
# in [m]
# local grid or UTM, depending on the coordinates of the SG !
sprinklingArea_x = c(input$modelDomain_xmin, input$modelDomain_xmax) # min, max
sprinklingArea_y = c(input$modelDomain_ymin, input$modelDomain_ymax) # min, max
# grid3d_depth = c(-3, 0) # min, max
# UTM
# sprinklingArea_x = c(SG_x - 7.5, SG_x + 7.5) # min, max
# sprinklingArea_y = c(SG_y - 7.5, SG_y + 7.5) # min, max
# grid3d_depth = c(SG_Z, SG_Z - 3) # min, max
#
## Model discretization
# in [m]
grid3d_discr = data.frame(x = input$modelDiscr_x, y = input$modelDiscr_y, z = input$modelDiscr_z) # 10 runs; 5.30 mins
# grid3d_discr = data.frame(x = .1, y = .1, z = .1) # 10 runs;  mins
# in REAL modeling so far: allDir = 0.1; depth up to 5 (?10) m;
grid3d_depth = c(input$modelDepth_min, input$modelDepth_max) # min, max
#
## Boundaries of SG pillar
# please use same units as in DEM and model domain
# if pillar has the structure of a rectangular pillar
# local grid
# Building_SGpillar_x = c(-1, 1) # min, max
# Building_SGpillar_y = c(-1, 1) # min, max
# Building_SGpillar_z = c(-1.2, 0) # min, max
# # UTM
# Building_SGpillar_x = c(SG_x - 1, SG_x + 1) # min, max
# Building_SGpillar_y = c(SG_y - 1, SG_y + 1) # min, max
# Building_SGpillar_z = c(SG_Z - 1.2, SG_Z) # min, max
# if pillar has the structure of a cylinder
# this is independent of local grid or UTM coordinates
# in [m]
thres_radius = input$pillar_r
thres_depth = input$pillar_d
#
## Input files
## general settings
# in case using .csv data, the special information has to supplied, in which columns the spatial information is stored
# the settings below are valid for 2d data files
# in the vector, the order is: x, y, z
spatial_col = c(input$spat_col_x, input$spat_col_y, input$spat_col_z)
# in all cases, a column has to be specified, containing the observation data
# columns of observation data (in .csv and .rData files)
data_col = input$data_col
# columns of observation data (in .tsf files)
data_tsf = input$data_tsf
# if the .csv has special characters for separating columns
# or decimal places, etc.
# the have to be EXPLICITLY specified in the read_data-function
# using sep = "??"
# using dec = "??"
# for further usage see ?read.csv
#
## DEM input file
# file name including its path
# should be absolute
# if left empty, a flat topographie will be assumed
# if(input$inputFile_DEM == ""){
# DEM_input_file = ""
# }else{
DEM_input_file = input$inputFile_DEM
# }
# DEM_input_file = "WE_UP_TO_300m_05m.asc"
#
## Water intensity distribution file
# IntensityDistribution_file = "waterIntensity_measured.csv"
IntensityDistribution_file = input$inputFile_IntDistr
# which interpolation algorithm should be used:
# inverse distance weight (IDW) or kriging (krige)
interpolation_method = input$IntpMethod
# set percentage of how many zeros should be added at border / side of grid
# if not desired, set to 0
Zeros_border_density = input$ZeroBorderDensity
#
## Observed gravity data time series
# this is optional and can be left empty if no automatized reduction is desired
gravityObservations_file = input$inputFile_gObs
# gravityObservations_input_file = "iGrav006_obs_60sec.rData"
## Information for optimization algorithm
# use inverse or conversion mode
# if set to FALSE, a single conversion run of the infiltration model is exectuted
# in this case, infiltration parameters supplied in 'starting values' are used as model input
if(input$Modeling_mode == "inverse"){
inverse = TRUE
}else{
inverse = FALSE
}
# number of iterations of the algorithm
model_runs = input$model_iterations
# delete results and log-file of previous inverse model runs?
# recommended to turn this on, experienced problems with the algoritum to 
# deal with previous input
# but it should be possible to even continue and extend previous runs
if(input$del_prevModRuns == "yes"){
delete_prevRuns = TRUE
}else{
delete_prevRuns = FALSE
}
#
## set infiltration model parameters
#
# set scenarios to include in optimization routine
# use macro pore flow on top?
# use_macro = TRUE
# 2 macro pore layers area only used if two_macro is set TRUE
# two_macro = TRUE 

if(input$macroporeLayers == 0){
        use_macro = FALSE
        two_macro = FALSE
}
if(input$macroporeLayers == 1){
        use_macro = TRUE
        two_macro = FALSE
}
if(input$macroporeLayers >= 2){
        use_macro = TRUE
        two_macro = TRUE 
}
# switch(input$macroporeLayers,
#        0 = {
#         use_macro = FALSE
#         two_macro = FALSE
#        },
#        1 = {
#         use_macro = TRUE
#         two_macro = FALSE
#        },
#        2 = {
#         use_macro = TRUE
#         two_macro = TRUE 
#        }
#        )

# further infiltration dynamics
# possible options are: wetting front advancement (wfa), by-pass flow and perched water table
# in the optimization, input is treated and dealt with in terms of real numbers
# consequentially, discrete values as infiltration process descriptions have to be translated accordingly
# including or exlucing a process is therefore realized via adjusting the boundary values of the following vector
# Wfa = 1
# perched water table = 2
# by-pass flow = 3
# the example setting therefore includes both Wfa and perched water table scenarios
inf_dynamics_min = input$inf_dynamics_min
inf_dynamics_max = input$inf_dynamics_max
#
## modeling time (duration of sprinkling experiment)
# [min]
# precip_time = 360 
precip_time = input$precip_time
## water input per timestep
# in [m³/min]
water_vol_min = input$water_vol_min
#
# set permitted error for mass balance
mb_permitted_error = input$mb_permitted_error
#
## Defines soil parameter boundaries
# min and max values, defining the search boundaries for the optimization algorithm
# Saturation deficit (dtheta)
# macropore layer
dtheta_macro_min = input$theta_macro_min #[VWC]
dtheta_macro_max = input$theta_macro_max #[VWC]
dtheta_macro2_min = input$theta_macro2_min #[VWC]
dtheta_macro2_max = input$theta_macro2_max #[VWC]
# Depth of processes
# in the case of no macro pore flow layers, the parameter 'mdepth' will be used for
# the depth of the single infiltration processes
mdepth_min = input$depth_macro_min #[m]
mdepth_max = input$depth_macro_max #[m]
mdepth2_min = input$depth_macro2_min #[m]
mdepth2_max = input$depth_macro2_max #[m]
# secondary infiltration process
# vertical bounaries
# if use_macro is set FALSE, this will be the only process used
dtheta_other_min = input$theta_other_min #[VWC]
dtheta_other_max = input$theta_other_max #[VWC]
# other infiltration processes ("pipe") are now spatially DIRECTLY connected below the macro pore layer
# if other is desired (e.g. gap between macro and pipe), this has to be activated again !
# with the following lines uncommented, a gap between macro and pipe layer is allowed
# break up criteria when pipedepth < mdepth is implemented in objective function (inf_model_3d)
# pipedepth_min = 0.2 #[m]
# pipedepth_max = 4.5 #[m]
pdepth_min = input$depth_other_min
pdepth_max = input$depth_other_max
#
# min max values for dividing water into horizontal / lateral parts (factor)
# in [%]
latflow_fac_min = input$latflow_fac_min #[1]
latflow_fac_max = input$latflow_fac_max #[1]
#
# Starting values of above set infiltration model parameters
dtheta_macro_start = input$theta_macro_start
dtheta_macro2_start = input$theta_macro2_start
mdepth_start = input$depth_macro_start
mdepth2_start = input$depth_macro2_start
dtheta_other_start = input$theta_other_start
# pipedepth_start = 0.4
latflow_fac_start = input$latflow_fac_start
# infiltration process
inf_dynamics_start = input$inf_dynamics_start
# other process starting depth
pdepth_start = input$depth_other_start
#
## plotting options
plot_interval = input$plot_interval
plot_transect_loc = input$plot_transect_loc
#

# ## test
# output$tt = renderPlot({
#   plot(pdepth_start)
# })
# ##

message("All input read..running infiltration model..")
## set working directory
setwd(dir_input)

## test
# plotting_output = reactive({
#   plot(input$plot_interval : 100)
# })

# output$plot_running_images = renderPlot({
#   plot(input$plot_interval : 100)
# })

#########################################
## Gravimeter location
#########################################
SG_z = SG_Z + SG_SensorHeight
SGloc = data.frame(x=SG_x, y=SG_y, z=SG_z)
#
#########################################
## Generate cropped DEM and surface grid
#########################################
message("Generate cropped DEM and surface grid..")
#
surface_grid = surface_grid(
            DEM = DEM_input_file,
            grid_domain_x = sprinklingArea_x,
            grid_domain_y = sprinklingArea_y,
            grid_discretization = grid3d_discr,
            input_dir = dir_input,
            output_dir = dir_output
            # , sep = "a", etc.
)
#
if(!is.null(surface_grid)){
  if(plot_data){
  surface.gg =  ggplot(surface_grid, 
           aes(x=x, y=y)) + 
           geom_tile(aes(fill = z))
 output$surfaceGrid = renderPlot({ 
    surface.gg
  })
  }
}
#
message("done.")
#########################################
# Generate 3d gravity component grid 
#########################################
message("Generate 3d gravity component grid..")
#
gravity_component_grid3d = gravity_comp_grid(
            surface = surface_grid,
            SG_coordinates = SGloc,
            grid_discretization = grid3d_discr,
            grid_depth = grid3d_depth,
            range_coords_x = sprinklingArea_x,
            range_coords_y = sprinklingArea_y
)
#
if(plot_data){
# message("plotting.")
    gravity_grid.gg = ggplot(gravity_component_grid3d, 
         aes(x=x, y=y)) + 
         geom_tile(aes(fill = z))
         # geom_point(aes(color = z))
 output$gCompGrid = renderPlot({ 
     gravity_grid.gg
  })
}
# 
message("done.")
#########################################
## Correct gravity component grid for SG pillar 
#########################################
message("Removing SG pillar from gravity component grid..")
#
gravity_component_grid3d = correct_SGpillar(
            gravity_comp3d = gravity_component_grid3d,
            # Pillar_x = NA,
            # Pillar_y = NA,
            # Pillar_z = NA,
            correct_radius = thres_radius,
            correct_depth = thres_depth,
            SG_X = SG_x,
            SG_Y = SG_y #,
            # grid_discretization = NA
)
#
# save grid
save(gravity_component_grid3d, file = paste0(dir_output, "gravity_component_grid3d.rData"))
#
if(plot_data){
  message("Plotting transect of gravity component grid and saving plot to output directory..")
  gravity_grid_mod.gg = plot_gcomp_grid(
                  grid_input = gravity_component_grid3d,
                  yloc = SG_y,
                  output_dir = dir_output,
                  grid_discretization = grid3d_discr
  )
 output$gCompGrid_PillarRemoved = renderPlot({ 
    gravity_grid_mod.gg
  })
}
#
message("done.")
#########################################
## Create water intensity distribution grid
#########################################
message("Creating grid for water intensity distribution..")
#
Intensity_distribution_interpolated = create_WaterIntensityGrid(
            input_file = IntensityDistribution_file,
            intp_method = interpolation_method,
            surface = gravity_component_grid3d,
            zerosBorder = Zeros_border_density,
            spat_col = spatial_col,
            dat_col = data_col,
            UTM_gridcenter_x = SG_x,
            UTM_gridcenter_y = SG_y,
            input_dir = dir_input, 
            output_dir = dir_output
)
#
if(plot_data){
  waterIntGird.gg = ggplot(Intensity_distribution_interpolated, 
         aes(x=x, y=y)) + 
         geom_tile(aes(fill = intensity))
 output$waterIntensityGrid = renderPlot({ 
    waterIntGird.gg
  })
         # geom_point(aes(color = intensity))
}
#
# save grid
save(Intensity_distribution_interpolated, file = paste0(dir_output, "Intensity_distribution_interpolated.rData"))
#
message("done.")
#########################################
## Run infiltration model
#########################################
#
# combine data for config file
configfile = data.frame(dir_input,
                        dir_output,
                        precip_time,
                        water_vol_min,
                        IntensityDistribution_file = "Intensity_distribution_interpolated.rData",
                        gcomp_file = "gravity_component_grid3d.rData",
                        gravityObservations_file,
                        data_tsf,
                        spatial_col_x = spatial_col[1],
                        spatial_col_y = spatial_col[2],
                        spatial_col_z = spatial_col[3],
                        data_col,
                        discr_x = grid3d_discr$x,
                        discr_y = grid3d_discr$y,
                        discr_z = grid3d_discr$z,
                        mb_permitted_error,
                        use_macro,
                        two_macro,
                        model_runs,
                        plot_interval,
                        plot_transect_loc,
                        stringsAsFactors=FALSE)
save(configfile, file=paste0(dir_output, "configfile.rdata"))
#
## run model in conversion or inversion mode
print("Running infiltration model..")
#

if(!inverse){
  print("Model is run in conversion mode.")
st = Sys.time()
  model_result = run_model_conversion(
              dtheta_macro = dtheta_macro_start,
              dtheta_macro2 = dtheta_macro2_start,
              mdepth = mdepth_start,
              mdepth2 = mdepth2_start,
              dtheta_other = dtheta_other_start,
              latflow_fac = latflow_fac_start,
              inf_dynamics = inf_dynamics_start,
              pdepth = pdepth_start,
              output_dir = dir_output
              )
en = Sys.time()
time_elapsed = en - st
#
}else{
## Run optimization algorithm
  print("Model is run in inversion mode.")
  print("Run optimization algorithm..")
  print("..this will take some time..")
 # 
  ## run model within optimization algorithm
  # for changing optimization function additional parameters, please see ?optim_dds
st = Sys.time()
  model_result = run_model_inversion(
              dtheta_macro_min = dtheta_macro_min,
              dtheta_macro_max = dtheta_macro_max,
              dtheta_macro2_min = dtheta_macro2_min,
              dtheta_macro2_max = dtheta_macro2_max,
              dtheta_other_min = dtheta_other_min,
              dtheta_other_max = dtheta_other_max,
              mdepth_min = mdepth_min,
              mdepth_max = mdepth_max,
              mdepth2_min = mdepth2_min,
              mdepth2_max = mdepth2_max,
              latflow_fac_min = latflow_fac_min,
              latflow_fac_max = latflow_fac_max,
              inf_dynamics_min = inf_dynamics_min,
              inf_dynamics_max = inf_dynamics_max,
              pdepth_min = pdepth_min,
              pdepth_max = pdepth_max,
              dtheta_macro_start = dtheta_macro_start,
              dtheta_macro2_start = dtheta_macro2_start,
              mdepth_start = mdepth_start,
              mdepth2_start = mdepth2_start,
              dtheta_other_start = dtheta_other_start,
              latflow_fac_start = latflow_fac_start,
              inf_dynamics_start = inf_dynamics_start,
              pdepth_start = pdepth_start,
              input_dir = dir_input,
              output_dir = dir_output,
              inner_inum = 1,
              del_prev = delete_prevRuns
  )
#
en = Sys.time()
time_elapsed = en - st

  print("Finished optimization.")
# end of inversion / conversion mode infiltration model runs
}
# print time necessary for individual setup model run
  print(paste0("Model run needed: ", time_elapsed))
#
# save model data (parameters in- and output)
save(model_result, file=paste0(dir_output, "Model_stats.rdata"))
write.table(model_result, file=paste0(dir_output, "Model_stats.csv"), sep="\t", dec=".", row.names = F, col.names = T, append = F)

## render table for shiny output
# select parameters to show
model_result_output = model_result %>%
  dplyr::select(inf_process_tested, model_kge,
              contains("dtheta_macro"),
              contains("dtheta_macro2"),
              contains("dtheta_other"),
              contains("mdepth"),
              contains("mdepth2"),
              contains("pdepth"),
              contains("latflow_fac"))

output$parameter_output = renderTable({
      model_result_output
})
##
print("done.")

message("All calculations have finished.")
message(paste0("Please have a look at the output file, located at: ", dir_output))

message("If you use this software in your publication, please cite this package.
        Information can be obtained using citation()")

####################
## output generation
####################

#########################################
## Plot: Gravity reponse (observed and modeled)
#########################################

if(plot_data){
  message("Plotting modeled and observed gravity signal..")
#
if(!inverse){
   gravity_response_run.gg = plot_gravity_responses(
              gravity_obs = gravityObservations_file,
              gravity_mod = paste0("model_output/GravityResponse_Infiltration_model_1.rData"),
              input_dir = dir_input,
              output_dir = dir_output
  )
  output$plot_gravity_responses_run = renderPlot({
    gravity_response_run.gg
  })
}else{
# plot LAST (optimized) model run scenario
   gravity_response_run.gg = plot_gravity_responses(
              gravity_obs = gravityObservations_file,
              gravity_mod = paste0("model_output/GravityResponse_Infiltration_model_", (n_param - 1), ".rData"),
              # gravity_mod = paste0("model_output/GravityResponse_Infiltration_model_9.rData"),
              input_dir = dir_input,
              output_dir = dir_output
  )
  output$plot_gravity_responses_run = renderPlot({
    gravity_response_run.gg
  })
}
  message("done.")
}else{
  message("No plotting desired.")
}

#########################################
## Plot: model output along a 2d transect
#########################################

if(plot_data){
  message("Plotting 2d transect of modeled soil moisture data..")

if(!inverse){
  SM_transect_run.gg = plot_transects_2d(
              soilmoisture_mod = paste0("model_output/Infiltration_model_output_1.rData"),
              plot_int = plot_interval,
              y_pos = SG_y,
              vert_limit = NA,
              output_dir = dir_output
  )
  output$plot_soilMoisture_transect_run = renderPlot({
      SM_transect_run.gg
  })
}else{
  SM_transect_run.gg = plot_transects_2d(
              soilmoisture_mod = paste0("model_output/Infiltration_model_output_", (n_param - 1), ".rData"),
              # soilmoisture_mod = paste0("model_output/Infiltration_model_output_9.rData"),
              plot_int = plot_interval,
              y_pos = SG_y,
              vert_limit = NA,
              # input_dir = dir_input,
              output_dir = dir_output
  )
  output$plot_soilMoisture_transect_run = renderPlot({
      SM_transect_run.gg
  })
}


#########################################
## Plot: modeled soil moisture over time and depth
#########################################

  message("Plotting modeled soil moisture over depth and time..")

  SM_run.gg = plot_waterDistributionOverTime(
              soilmoisture_mod = paste0("model_output/Infiltration_model_output_", (n_param - 1), ".rData"),
              output_dir = dir_output
            )
  output$plot_soilMoisture_run = renderPlot({
    SM_run.gg
  })

  message("done.")
}else{
  message("No plotting desired.")
}

# remove iteration parameter for inner optimization function calls
# rm(n_param)

## END: tracking R console output and pass it to shiny text output
      },
        message = function(m) {
          shinyjs::html(id = "text", html = m$message, add = TRUE)
      })

#########################################
## end of "Running model button"
      }
    )
#########################################

#########################################
## showing RESULTS
#########################################

####################
## plotting ABOVE
 observeEvent(
      eventExpr = input[["plot_above"]],
      handlerExpr = {
    # input values
    dir_input = input$plot_inputDir
    dir_output = input$plot_outputDir
    gravityObservations_file = input$plot_inputFile_gObs
    plot_interval = input$plot_plot_interval
    y_loc = input$plot_plot_transect_loc
    n_param_1 = input$plot_nparam_1
    #
# start plotting, when button is pressed
## process indicator, showing time
withProgress(message = 'Still plotting..', value = 0, {

  gravity_response_1.gg = plot_gravity_responses(
              gravity_obs = gravityObservations_file,
              gravity_mod = paste0("model_output/GravityResponse_Infiltration_model_", n_param_1, ".rData"),
              input_dir = dir_input,
              output_dir = dir_output
  )
  output$plot_gravity_responses_1 = renderPlot({
    gravity_response_1.gg
  })

  SM_1.gg = plot_waterDistributionOverTime(
              soilmoisture_mod = paste0("model_output/Infiltration_model_output_", n_param_1, ".rData"),
              output_dir = dir_output
            )
  output$plot_soilMoisture_1 = renderPlot({
    SM_1.gg
  })

  SM_transect_1.gg = plot_transects_2d(
              soilmoisture_mod = paste0("model_output/Infiltration_model_output_", n_param_1, ".rData"),
              plot_int = plot_interval,
              y_pos = y_loc,
              vert_limit = NA,
              # input_dir = dir_input,
              output_dir = dir_output
  )
  output$plot_soilMoisture_transect_1 = renderPlot({
      SM_transect_1.gg
  })

## render table for shiny output
# select parameters to show
# read in model_result data
load(file=paste0(dir_output, "Model_stats.rdata"))

model_result_output = model_result %>%
  dplyr::select(inf_process_tested, model_kge,
              contains("dtheta_macro"),
              contains("dtheta_macro2"),
              contains("dtheta_other"),
              contains("mdepth"),
              contains("mdepth2"),
              contains("pdepth"),
              contains("latflow_fac"))

output$plot_parameter_output = renderTable({
      model_result_output
})
## end progress indicator
})

# end of "Plot above button"
}
)

####################
## End of shiny server
####################
})
