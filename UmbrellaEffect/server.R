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
# library(grid)
# library(gridExtra)
# library(scales)
# library(hydroGOF)
# library(data.table)
# library(ppso)
# # for kriging
# library(spacetime)
# library(sp)
# library(automap)
# library(akima)

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
      eventExpr = input[["run_reduction"]],
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
Building_x = c(input$modelDomain_xmin, input$modelDomain_xmax) # min, max
Building_y = c(input$modelDomain_ymin, input$modelDomain_ymax) # min, max
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
## Parameters for foundation of building
# these include baseplate, walls, SG pillar(s)
# please use same units as in DEM and model domain
Building_walls_x = input$building_walls_x # extension
Building_walls_y = input$building_walls_y # extension 
Building_walls_z = input$building_walls_z # extension
# local grid
Building_baseplate_z = c(input$building_baseplate_z_min, input$building_baseplate_z_max) # min, max
Building_SGpillar_x = c(input$pillar_x_min, input$pillar_x_max) # min, max
Building_SGpillar_y = c(input$pillar_y_min, input$pillar_y_max) # min, max
Building_SGpillar_z = c(input$pillar_z_min, input$pillar_z_max) # min, max
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
## Soil moisture data time series (observed or modelled)
soilMoisture_input_file = input$inputFile_SM
#
## Observed gravity data time series
# this is optional and can be left empty if no automatized reduction is desired
gravityObservations_file = input$inputFile_gObs
# gravityObservations_input_file = "iGrav006_obs_60sec.rData"
#
## SG position
# options are: Center, Corner, Side, Trans, Wettzell
SG_position = input$SGposition
#
## Hydrological site condition
# options are: 
# "Soil type sandy loam" , "Soil type silty loam", "Soil type clay loam" ; for soil type scenarios
# "Anisotropy sandy loam", "", "" ; for anisotropy of hydraulic conductivity scenarios
# "Climate AI=0.35", "Climate AI=0.575", "Climate AI=1.0", "Climate AI=2.0" ; for climate scenarios
Hydro_condition = input$Hydro_condition
#
## Vertical extent of reduction
# indicates the vertical extent until which a reduction factor should be calculated and applied
# this is valid for both Hydroglical site condition & SG position reduction as also building dimension reduction factors
# units in [m]
# maximum value: 5
VerticalExt_reduction = input$Vertical_extent
#
#
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
            grid_domain_x = Building_x,
            grid_domain_y = Building_y,
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
            range_coords_x = Building_x,
            range_coords_y = Building_y
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
## depending on selection,
# either rectangular or circle like pillar
if(pillar_correction_rect){
  message("Correcting with rectangular pillar")
  gravity_component_grid3d = correct_SGbuilding_foundation(
              gravity_gcomp = gravity_component_grid3d,
              Bdwall_ext_x = Building_walls_x,
              Bdwall_ext_y = Building_walls_y,
              Bdwall_ext_z = Building_walls_z,
              Bdbase_x = Building_x,
              Bdbase_y = Building_y,
              Bdbase_z = Building_baseplate_z,
              Pillar_x = Building_SGpillar_x,
              Pillar_y = Building_SGpillar_y,
              Pillar_z = Building_SGpillar_z,
              grid_discretization = grid3d_discr
  )
}else{
  message("Correcting with circle pillar")
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
}
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
## Extrapolate soil moisture time series data (observed or modelled) to gravity grid domain
#########################################
message("Extrapolate soil moisture time series data (observed or modelled) to gravity grid domain..")

SMgrid3d_outside = SoilMoisture_grid3d(
            grid_domain = gravity_component_grid3d,
            soilMoisture_input = soilMoisture_input_file,
            grid_discretization = grid3d_discr,
            grid_depth = grid3d_depth,
            input_dir = dir_input
            # , sep = "a", etc..
)

message("done.")
#########################################
## Calculate gravity response (from outside of building)
#########################################
message("Calculate gravity response (from outside of building)..")

gravity_response_outside_building = calculate_gravity_response(
            gcomp_grid = gravity_component_grid3d,
            mass_input = SMgrid3d_outside
)

message("done.")
#########################################
## Calculate mean soil moisture within model domain for each timestep
#########################################
message("Calculate mean soil moisture within model domain for each timestep..")

SoilMoisture_mean_ts = dplyr::group_by(SMgrid3d_outside, datetime) %>%
                           dplyr::summarize(value = mean(value, na.rm=T))

message("done.")
#########################################
## Convert gravity response (outside) to gravity response below SG building
#########################################
message("Convert gravity response (outside) to gravity response below SG building..")

# reduction factor corresponding to chosen dominant hydrological scenario and SG position
reduction_factor_hydSen_SGloc = reduction_hydScen_SGloc(
            Scenario == Hydro_condition,
            SGlocation == SG_position,
            VertLimit = VerticalExt_reduction,
            MeanSoilMoisture = SoilMoisture_mean_ts
)

# calculate SG building size
buildingSize = BdSize(
           Bd_x = Building_x,
           Bd_y = Building_y
)

# reduction factor corresponding to size (area) of the SG building
reduction_factor_BdSize = reduction_BdSize(
            SG_BdSize = buildingSize,
            VertLimit = VerticalExt_reduction
)

## convert gravity response from next to building
gravity_response_below_building = convert_gravity_response(
            gravity_input = gravity_response_outside_building,
            factor_hydScen_SGloc = reduction_factor_hydSen_SGloc,
            factor_BdSize = reduction_factor_BdSize
)

message("done.")
#########################################
## Save gravity response of mass variations, occuring below SG building
#########################################
message("Save gravity response of mass variations, occuring below SG building..")

if(output_type == "csv"){
    save(gravity_response_below_building, file=paste0(dir_output, "gravity_response_below_building.rData"))
    write.table(...)
}else{
    save(gravity_response_below_building, file=paste0(dir_output, "gravity_response_below_building.rData"))
}

message("done.")
#########################################
## Correct gravity observation data (if supplied)
#########################################

if(gravityObservations_input_file == ""){
  message("No reduction of observed gravity data desired.")
}else{
  message("Reducing gravity observation data..")
  
  gravity_data_reduced = reduce_gravity(
            gravity_obs = gravityObservations_input_file,
            gravity_below = gravity_response_below_building,
            input_dir = dir_input,
            dat_tsf = data_tsf
  )
  
  if(output_type == "csv"){
      save(gravity_data_reduced, file=paste0(dir_output, "gravity_data_reduced.rData"))
      write.table(gravity_data_reduced, file=paste0(dir_output, "gravity_data_reduced.csv"), row.names = F)
  }else{
      save(gravity_data_reduced, file=paste0(dir_output, "gravity_data_reduced.rData"))
  }
  
  message("done.")
}

####################
## output generation
####################

#########################################
## Plot all time series
#########################################

if(plot_data){
  message("Plotting time series and saving plot to output directory..")
  if(gravityObservations_input_file == ""){
    plot_ts_data.gg = plot_ts_data(
            gravity_outside = gravity_response_outside_building,
            gravity_below = gravity_response_below_building,
            input_dir = dir_input,
            output_dir = dir_output
  )
  }else{
    plot_ts_data.gg = plot_ts_data(
            gravity_obs = gravityObservations_input_file,
            gravity_outside = gravity_response_outside_building,
            gravity_below = gravity_response_below_building,
            gravity_reduced = gravity_data_reduced,
            input_dir = dir_input,
            output_dir = dir_output,
            dat_tsf = data_tsf
  )
  }

  output$plot_all_gravity_data = renderPlot({
    plot_ts_data.gg
  })

  message("done.")
}else{
  message("No plotting desired.")
}

## end CALCULATIONS
#########################################

message("ALL calculations have finished.")
message("Please have a look at the output file, located at: ")
message(dir_output)

message("If gravity observation data was supplied, the data has been recuded automatically by the UmbrellaEffect results,
and stored as well in the output directory.")


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
