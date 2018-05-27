####################################
# 		  Floodmapper v1.0
# ----------------------------------
# Cross-scale flooding predicition 
# under heavy precipitation
# ----------------------------------
# Created by Xander Wang
# Date: April 11, 2018
# Where: Lafayette, LA, USA
# Email: xiuquan.wang@gmail.com
# Copyright @ Xander Wang
####################################

################################
## [0] Required libraries
################################
library(sp)
library(raster)
library(rgdal)
library(magick)

################################
## [1] Initialization
################################
#: create a new environment to host all parameters
FM.env <- new.env()

## Constants
FM.env$s_versionname <- "FloodMapper v1.0" 	#: name of this version
FM.env$n_g <- 9.8						#: gravity acceleration (g, unit: m/s^2)
FM.env$s_land_param_file <- "land.csv"	#: a csv file stores the parameter settings for land cover
FM.env$s_soil_param_file <- "soil.csv"	#: a csv file stores the parameter settings for soil texture
FM.env$s_readme_file <- "readme.txt"	#: a txt file provides a brief description about the folder structure
FM.env$n_integer_factor <- 1e+10		#: a factor to be applied for integer conversion, for comparison purpose
FM.env$n_max_water_velocity <- 3		#: 3 m/s, to make sure the integral time step is small enough to perform a stable run
									
## [1.0] Run parameters
FM.env$s_run_name <- ""					#: a folder with this name under your current R work directory will be created
										#: under the folder: input + output
FM.env$s_run_folder_path <- ""			#: full path to the current run folder
FM.env$s_worddir <- ""					#: full path to user's work directory
FM.env$o_run_start_datetime <- NA		#: start datetime of your model run, format: "2016-06-09 12:00:00"
FM.env$o_run_end_datetime <- NA			#: end datetime of your model run, format: "2016-06-09 13:00:00"
FM.env$n_run_output_interval <- 1800	#: output time step for the simulation results, unit: seconds, default: 30 min = 30 * 60 = 1800 sec
FM.env$n_run_internal_timestep <- 30	#: internal integral time step, unit: seconds, default: 30 dec
FM.env$n_run_total_internalruns <- 0	#: indicate how many internal model runs are required
FM.env$b_debug <- TRUE					#: TRUE = a log file is generated under the run folder, FALSE = no log file is generated
FM.env$s_logfile <- ""					#: full file path to the log file

## [1.1] DEM raster 
FM.env$s_dem_raster_file <- ""			#: full file path to DEM raster
FM.env$o_dem_nonNA_index <- NA			#: a vector to store the indices in the "values()" vector of the original DEM raster
FM.env$n_dem_nrow <- 0					#: total rows in DEM raster
FM.env$n_dem_ncol <- 0					#: total columns in DEM raster
FM.env$n_dem_cellsize <- 0				#: cell size of DEM raster (unit: meters)
FM.env$s_dem_unit <- "m"				#: unit of DEM raster, MUST be in meters!
FM.env$n_dem_cellsize_conversion <- 1	#: RESERVED, 1 = unit is meters; otherwise, specify the conversion factor from your unit to meters
FM.env$o_dem_value_matrix <- NA			#: a matrix to store the DEM data
FM.env$o_waterdepth_rasfile <- ""		#: full file path to surface water depth raster
FM.env$o_waterdepth_matrix <- NA		#: a matrix to store the status variable: surface water depth, unit: m ==> to be updated very time step

## [1.2] Land cover raster
FM.env$s_land_raster_file <- ""			#: full file path to land cover raster
FM.env$o_land_Cd_matrix <- NA			#: a matrix to store the drag coefficient, Cd
FM.env$o_land_alpha_matrix <- NA		#: a matrix to store the infiltration adjustment coefficient, alpha

## [1.3] Soil texture raster, depth raster (OPTIONAL), and moisture content raster (OPTIONAL)
FM.env$s_soil_raster_file <- ""			#: full file path to soil texture raster
FM.env$o_soil_porosity_matrix <- NA		#: a matrix to store the soil porosity, n
FM.env$o_soil_beta_matrix <- NA			#: a matrix to store the damping coefficient, beta
FM.env$s_soil_depth_file <- ""			#: full file path to soil depth raster
FM.env$o_soil_depth_matrix <- NA		#: a matrix to store the depth of surface soil, unit: m
FM.env$s_soil_moistcontent_file <- ""	#: full file path to soil moisture content raster
FM.env$o_soil_moistcontent_matrix <- NA	#: a matrix to store the current soil moisture content, unit: % ==> to be updated very time step
FM.env$o_soil_maxwaterVol_matrix <- NA	#: a matrix to store the maximum water volume for the saturated surface soil, unit: m3

## [1.4] Street drain inlet
FM.env$n_draininlet_flag <- 0			#: 1 = street drain inlets will be considered in the model, default: 0 (do NOT consider)
FM.env$s_draininlet_raster_file <- ""	#: full file path to the drain inlet raster file, unit: mm/hr, a positive value indicates a drain inlet, 0 = not an drain inlet
										#: contains the maximum rainfall intensity that can be handled by each drain inlet
										#: the maximum rainfall intensity can be determined according IDF curves once return period and duration of an event are known
FM.env$o_draininlet_maxI_matrix <- NA	#: a matrix to store the max precipitation intensity accepted by the inlet per second, unit: m/sec
										#: a positive value indicates a drain inlet, 0 = not a drain inlet
										#: should be calculated based on Duration & Max Intensity, as follows:
										#: total water depth (mm) during the given duration = Duration * MaxIntensity
										#: duration is expressed as in seconds = Duration * 3600
										#: maxI in m/sec = (Duration * MaxIntensity) / (Duration * 3600) = MaxIntensity / 3600
										#: although Duration is cancelled, we should NOTE that different Durations lead to different MaxIntensitys
										#: it is the users' responsibility to determine the duration of an expected precipitation event
										#: in other words, you need to determine how long this precipitation event will last
										#: then, find the corresponding max intensity from the IDF curves.
										
## [1.5] Precipitation data
FM.env$n_precip_datatype <- 0			#: 0 = you will provide a text file containing all weather stations, default value
										#: the format of this text file is given as follows:
										#: ------------------------------------------
										#: ID	Lon			Lat			Datafile
										#: 1	-92.0296	30.1800		p1.txt		(stores obs data for station 1)
										#: 2	-92.0260	30.1779		p2.txt		(stores obs data for station 2)
										#: ------------------------------------------
										#: 1 = you will provide a text file containing the filenames of precip rasters for all time steps
										#: the format of this text file is given as follows:
										#: ------------------------------------------
										#: ID	Rasterfile
										#: 1	p_2016-06-09_12-00-00.tif	(although a user-defined filename is allowed)
										#: 2	p_2016-06-09_13-00-00.tif	(it is suggested to name it with date and time)
										#: ------------------------------------------
FM.env$s_precip_txt_file <- ""			#: full file path to one of the above text files
FM.env$o_precip_start_datetime <- NA	#: start datetime of your precipitation data records, format: "2016-06-09 12:00:00"
FM.env$o_precip_end_datetime <- NA		#: end datetime of your precipitation data records, format: "2016-06-09 13:00:00"
FM.env$n_precip_record_interval <- 3600	#: interval of your precipitation data records, unit: seconds, e.g.: 1 hr = 60 min * 60 sec = 3600 sec
FM.env$n_precip_IDW_power <- 2			#: the power of IDW method, default: 2, suggested options: {1, 2, 3}

## [1.6] Inflow from the outside of the domain (unit: m3/s)
FM.env$n_inflow_datatype <- 0			#: 0 = you will provide a text file containing all monitoring stations, default value
										#: the format of this text file is given as follows:
										#: ------------------------------------------
										#: ID	Lon			Lat			Datafile
										#: 1	-92.0296	30.1800		in1.txt		(stores obs data for station 1)
										#: 2	-92.0260	30.1779		in2.txt		(stores obs data for station 2)
										#: ------------------------------------------
FM.env$s_inflow_txt_file <- ""			#: full file path to one of the above text files
FM.env$o_inflow_start_datetime <- NA	#: start datetime of your inflow data records, format: "2016-06-09 12:00:00"
FM.env$o_inflow_end_datetime <- NA		#: end datetime of your inflow data records, format: "2016-06-09 13:00:00"
FM.env$n_inflow_record_interval <- 3600	#: interval of your inflow data records, unit: seconds, e.g.: 1 hr = 60 min * 60 sec = 3600 sec

## [1.7] Output settings
FM.env$o_png_aerial_mar <- c(1, 1, 2, 6)
FM.env$o_png_dem_mar <- c(1, 1, 2, 2.5)
FM.env$o_png_mgp <- c(5, 5, 5)
FM.env$n_png_pointsize <- 24
FM.env$n_png_resolution <- 50
FM.env$n_png_maxwidth <- 20
FM.env$n_png_maxheight <- 20
FM.env$n_png_width <- 20				#: to be adjusted according to the size of DEM raster or background aerial raster
FM.env$n_png_height <- 20				#: to be adjusted according to the size of DEM raster or background aerial raster
FM.env$s_png_png_bg_color <- "transparent"
FM.env$s_png_border_color <- "#232323"
FM.env$s_png_ras_bg_color <- "#232323"
FM.env$o_png_dem_palette <- colorRampPalette(c("#E1EBB7", "#207B2D", "#DEB53D", "#792302", "#8C443B", "#DFDDE0"), space = "rgb")
FM.env$n_png_dem_legend_breaks <- 40
FM.env$o_png_wd_breakpoints <- seq(0, 1000, 50) #: break points for surface water depth (unit: mm)
FM.env$o_png_wd_palette <- colorRampPalette(c("#FFFFFF00", "#BEFFFF9B", "#5AFFFF9B", "#32B4FF9B", "#0050FF9B", "#01008C9B", "#00004F9B"), alpha = TRUE)

################################
## [2] Functions
################################

## [2.1] model initialization
## -------------------------------------------------
## Name:	FM.init()
## Inputs:	
##			@runname:			name of the new model run
##			@startdatetime:		start datetime of the run, format: "2016-06-09 12:00:00"
##			@enddatetime:		end datetime of the run, format: "2016-06-09 13:00:00"
##			@outputinterval:	time step for the model outputs, unit: seconds
##			@internaltimestep:	OPTIONAL, internal integral time step, unit: seconds, default: 30 dec
##			@debug:				OPTIONAL, if it is in debug model, a log file will be generated under the run folder
##			@rerun:				OPTIONAL, FALSE = new run, TRUE = rerun (no need to create folders)
##			@wdbreakpoints:		OPTIONAL, a vector defines how to break the surface water depth (unit: mm)
## -------------------------------------------------
FM.init <- function(runname, startdatetime, enddatetime, outputinterval, internaltimestep = 30, debug = TRUE, rerun = FALSE, wdbreakpoints = NA, workdir = "")
{
	cat("*********************************************\n")
	cat("********  Model Initiazation: START  ********\n")
	cat("---------------------------------------------\n")
	cat("Initializing Model...")
	flush.console()

	FM.env$s_run_name <- runname
	#: determine the work directory
	if (workdir != "") {
		FM.env$s_worddir <- workdir
		FM.env$s_run_folder_path <- file.path(workdir, runname)
	} else {
		FM.env$s_worddir <- tempdir()
		FM.env$s_run_folder_path <- file.path(tempdir(), runname)
	}
	FM.env$o_run_start_datetime <- as.POSIXct(startdatetime)
	FM.env$o_run_end_datetime <- as.POSIXct(enddatetime)
	FM.env$n_run_output_interval <- as.integer(outputinterval)
	FM.env$n_run_internal_timestep <- as.integer(internaltimestep)
	if (FM.env$o_run_end_datetime <= FM.env$o_run_start_datetime) {
		return("ERROR: the end datetime should be after the start datetime!\n")
	} else if (FM.env$n_run_output_interval %% FM.env$n_run_internal_timestep != 0) {
		return(paste("ERROR: the output interval MUST be an integer multiple of integeral time step: ", internaltimestep, " sec!\n", sep = ""))
	}
	
	#: calculate total internal runs
	n_tmp_totalseconds <- as.numeric(difftime(FM.env$o_run_end_datetime, FM.env$o_run_start_datetime, units = "secs"))
	if (n_tmp_totalseconds %% FM.env$n_run_internal_timestep != 0) {
		return(paste("ERROR: the length of your simulation period MUST be an integer multiple of integeral time step: ", internaltimestep, " sec!\n", sep = ""))
	}
	FM.env$n_run_total_internalruns <- n_tmp_totalseconds / FM.env$n_run_internal_timestep
	
	#: if this is a new run
	if (!rerun) {
		#: create folders
		if (dir.exists(FM.env$s_run_folder_path)) {
			return(paste("ERROR: the folder of '", FM.env$s_run_name, "' already exists!\n", sep = "")) #: ERROR: folder already exists
		} else {
			dir.create(FM.env$s_run_folder_path)
			#: create sub folders: input, output
			dir.create(file.path(FM.env$s_run_folder_path, "input"))
			dir.create(file.path(FM.env$s_run_folder_path, "output"))
		}
		
		#: copy the internal csv parameter files for land cover and soil to the input folder
		#: change the parameters' values in these two csv files for sensitivity analysis
		#: FM.readLandcover() and FM.readSoil() must be re-executed after the changes are made
		################################## ==> to be updated for the package!!!
		file.copy(system.file("extdata", FM.env$s_land_param_file, package = "FloodMapper", mustWork = TRUE), paste(FM.env$s_run_folder_path, "/input/", FM.env$s_land_param_file, sep = ""), overwrite = TRUE)
		file.copy(system.file("extdata", FM.env$s_soil_param_file, package = "FloodMapper", mustWork = TRUE), paste(FM.env$s_run_folder_path, "/input/", FM.env$s_soil_param_file, sep = ""), overwrite = TRUE)
		
		#: copy the readme txt file to the root folder of the new model run
		file.copy(system.file("extdata", FM.env$s_readme_file, package = "FloodMapper", mustWork = TRUE), paste(FM.env$s_run_folder_path, "/", FM.env$s_readme_file, sep = ""), overwrite = TRUE)
	}

	#: create a log file if debug model is enabled
	FM.env$b_debug <- debug
	if (FM.env$b_debug) {
		FM.env$s_logfile <- paste(FM.env$s_run_folder_path, "/log.txt", sep = "")
		cat("*********************************************\n", file = FM.env$s_logfile, sep = "") #: clear all previous logs
		cat("********  Model Initiazation: START  ********\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		cat("---------------------------------------------\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		cat("Initializing Model...", file = FM.env$s_logfile, sep = "", append = TRUE)
		flush.console()
	}
	
	#: other configurations
	FM.env$o_draininlet_maxI_matrix <- matrix(0, FM.env$n_dem_nrow, FM.env$n_dem_ncol, byrow = TRUE)
	if (!is.na(wdbreakpoints[1])) FM.env$o_png_wd_breakpoints <- wdbreakpoints
	
	cat("\t>>> OK\n")
	cat("---------------------------------------------\n")
	cat("*****  Model Initiazation: SUCCESSFUL!  *****\n")
	cat("*********************************************\n")
	cat(paste("Folder to check your results:\n", FM.env$s_run_folder_path, "\n", sep = ""))
	if (FM.env$b_debug) {
		cat("\t>>> OK\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		cat("---------------------------------------------\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		cat("*****  Model Initiazation: SUCCESSFUL!  *****\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		cat("*********************************************\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		cat(paste("Folder to check your results:\n", FM.env$s_run_folder_path, "\n", sep = ""), file = FM.env$s_logfile, sep = "", append = TRUE)
	}
	flush.console()
	
	return("") #: SUCCESS
}

## [2.2] read DEM raster
## -------------------------------------------------
## Name:	FM.readDEM()
## Inputs:	
##			@demfile:			full file path of the DEM raster
##			@waterdepth:		(OPTIONAL) initial surface water depth, applied to all grid cells
##			@waterdepthfile:	(OPTIONAL) full file path of the initial surface water depth raster, default = no surface water
## -------------------------------------------------
FM.readDEM <- function(demfile, waterdepth = 0, waterdepthfile = "", rerun = FALSE)
{
	cat("*********************************************\n")
	cat("********          DEM: START         ********\n")
	cat("---------------------------------------------\n")
	cat("Processing DEM...")
	if (FM.env$b_debug) {
		cat("*********************************************\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		cat("********          DEM: START         ********\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		cat("---------------------------------------------\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		cat("Processing DEM...", file = FM.env$s_logfile, sep = "", append = TRUE)
	}
	flush.console()
	
	#: if it is a rerun, the existing raster files under the input folder will be used
	#: no need to perform copy operation
	if (rerun) {
		demfile <- FM.env$s_dem_raster_file
		waterdepthfile <- FM.env$o_waterdepth_rasfile
	}

	#: check if it is a valid DEM raster
	if (!file.exists(demfile)) {
		if (FM.env$b_debug) cat("ERROR: DEM file does NOT exist!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		return("ERROR: DEM file does NOT exist!\n") #: ERROR: DEM file does NOT exist!
	} else {
		#: check unit of DEM raster, MUST be in meters
		o_tmp_ras <- raster(demfile)
		s_tmp_str <- toString(o_tmp_ras@crs)
		n_tmp_pos <- regexpr("+units=", s_tmp_str)[1]
		s_tmp_unit <- substr(s_tmp_str, n_tmp_pos + 6, n_tmp_pos + 6)
		if (s_tmp_unit != FM.env$s_dem_unit) {
			if (FM.env$b_debug) cat("ERROR: the unit of DEM raster MUST be meters!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
			return("ERROR: the unit of DEM raster MUST be meters!\n") #: ERROR: the unit of DEM raster is NOT meters!
		} else {
			#: save DEM info for later use
			FM.env$o_dem_nonNA_index <- which(!is.na(o_tmp_ras[]))
			FM.env$n_dem_nrow <- dim(o_tmp_ras)[1]
			FM.env$n_dem_ncol <- dim(o_tmp_ras)[2]
			FM.env$n_dem_cellsize <- res(o_tmp_ras)[1] * FM.env$n_dem_cellsize_conversion #: convert into meters
			
			#: deal with negative DEM by lifting all grid cells above msl by 1.0 m
			o_tmp_dem_vector <- o_tmp_ras[]
			n_tmp_dem_min <- min(o_tmp_dem_vector, na.rm = TRUE)
			if (ceiling(n_tmp_dem_min * FM.env$n_integer_factor) < 0) {
				o_tmp_dem_vector <- o_tmp_dem_vector + abs(n_tmp_dem_min) + 1.0
			}
		
			FM.env$o_dem_value_matrix <- matrix(o_tmp_dem_vector, FM.env$n_dem_nrow, FM.env$n_dem_ncol, byrow = TRUE)
			
			#: check if the dt is small enough to meet: <= dx / FM.env$n_max_water_velocity
			if (ceiling(FM.env$n_run_internal_timestep * FM.env$n_max_water_velocity * 1.0) > floor(FM.env$n_dem_cellsize)) {
				if (FM.env$b_debug) cat("ERROR: the internal time step is too big and should be adjusted!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
				return("ERROR: the internal time step is too big and should be adjusted!\n")
			}
		}
		
		#: release tmp variables
		rm(o_tmp_ras)
		rm(s_tmp_str)
		rm(n_tmp_pos)
		rm(s_tmp_unit)
		
		#: copy this DEM to the input folder and name it as dem.tif (only if this is a new run)
		if (!rerun) {
			FM.env$s_dem_raster_file <- paste(FM.env$s_run_folder_path, "/input/dem.tif", sep = "")
			file.copy(demfile, FM.env$s_dem_raster_file, overwrite = TRUE)
		}

		#: update the png output settings
		if (FM.env$n_dem_nrow == FM.env$n_dem_ncol) {
			FM.env$n_png_width <- FM.env$n_png_maxwidth
			FM.env$n_png_height <- FM.env$n_png_maxheight
		} else if (FM.env$n_dem_nrow < FM.env$n_dem_ncol) {
			FM.env$n_png_width <- FM.env$n_png_maxwidth
			FM.env$n_png_height <- ceiling(FM.env$n_png_maxheight * FM.env$n_dem_nrow / FM.env$n_dem_ncol)
		} else if (FM.env$n_dem_nrow > FM.env$n_dem_ncol) {
			FM.env$n_png_width <- ceiling(FM.env$n_png_maxwidth * FM.env$n_dem_ncol / FM.env$n_dem_nrow)
			FM.env$n_png_height <- FM.env$n_png_maxheight
		}
	}
	cat("\t>>> OK\n")
	if (FM.env$b_debug) cat("\t>>> OK\n", file = FM.env$s_logfile, sep = "", append = TRUE)
	flush.console()
	
	cat("Initializing Water Depth...")
	if (FM.env$b_debug) cat("Initializing Water Depth...", file = FM.env$s_logfile, sep = "", append = TRUE)
	flush.console()
	
	#: initial surface water depth
	if (waterdepthfile != "") {
		#: a raster file is provided
		if (!file.exists(waterdepthfile)) {
			if (FM.env$b_debug) cat("ERROR: the raster file for surface water depth does NOT exist!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
			return("ERROR: the raster file for surface water depth does NOT exist!\n") #: ERROR: raster file does NOT exist!
		} else {
			#: check unit of the raster, MUST be in meters
			o_tmp_wrtepth_ras <- raster(waterdepthfile)
			s_tmp_wrtepth_str <- toString(o_tmp_wrtepth_ras@crs)
			n_tmp_wrtdepth_pos <- regexpr("+units=", s_tmp_wrtepth_str)[1]
			s_tmp_wrtdepth_unit <- substr(s_tmp_wrtepth_str, n_tmp_wrtdepth_pos + 6, n_tmp_wrtdepth_pos + 6)
			n_tmp_wrtdepth_nrow <- dim(o_tmp_wrtepth_ras)[1]
			n_tmp_wrtdepth_ncol <- dim(o_tmp_wrtepth_ras)[2]
			n_tmp_wrtdepth_cellsize <- res(o_tmp_wrtepth_ras)[1] * FM.env$n_dem_cellsize_conversion #: convert into meters
			if (s_tmp_wrtdepth_unit != FM.env$s_dem_unit) {
				if (FM.env$b_debug) cat("ERROR: the unit of surface water depth raster MUST be meters!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
				return("ERROR: the unit of surface water depth raster MUST be meters!\n") #: ERROR: the raster unit is NOT meters!
			} else if (n_tmp_wrtdepth_nrow != FM.env$n_dem_nrow || n_tmp_wrtdepth_ncol != FM.env$n_dem_ncol || n_tmp_wrtdepth_cellsize != FM.env$n_dem_cellsize) {
				if (FM.env$b_debug) cat("ERROR: the layout of grid cells in surface water depth raster MUST match that of DEM raster!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
				return("ERROR: the layout of grid cells in surface water depth raster MUST match that of DEM raster!\n") #: ERROR: DOES NOT match!
			} else {
				#: initialize surface water depth with the input raster file
				FM.env$o_waterdepth_matrix <- matrix(o_tmp_wrtepth_ras[], FM.env$n_dem_nrow, FM.env$n_dem_ncol, byrow = TRUE)
			}
			#: release tmp variables
			rm(o_tmp_wrtepth_ras)
			rm(s_tmp_wrtepth_str)
			rm(n_tmp_wrtdepth_pos)
			rm(s_tmp_wrtdepth_unit)
			rm(n_tmp_wrtdepth_nrow)
			rm(n_tmp_wrtdepth_ncol)
			rm(n_tmp_wrtdepth_cellsize)
		}
		
		#: copy this raster to the input folder and name it as waterdepth.tif (only if this is a new run)
		if (!rerun) {
			FM.env$o_waterdepth_rasfile <- paste(FM.env$s_run_folder_path, "/input/waterdepth.tif", sep = "")
			file.copy(waterdepthfile, FM.env$o_waterdepth_rasfile, overwrite = TRUE)
		}
		
	} else {
		#: initialize surface water depth as zero
		FM.env$o_waterdepth_matrix <- matrix(waterdepth, FM.env$n_dem_nrow, FM.env$n_dem_ncol, byrow = TRUE)
		FM.env$o_waterdepth_rasfile <- ""
	}
	
	cat("\t>>> OK\n")
	cat("---------------------------------------------\n")
	cat("*****          DEM: SUCCESSFUL!         *****\n")
	cat("*********************************************\n")
	if (FM.env$b_debug) {
		cat("\t>>> OK\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		cat("---------------------------------------------\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		cat("*****          DEM: SUCCESSFUL!         *****\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		cat("*********************************************\n", file = FM.env$s_logfile, sep = "", append = TRUE)
	}
	flush.console()
	
	return("") #: SUCCESS
}

## [2.3] read landcover raster
## -------------------------------------------------
## Name:	FM.readLandcover()
## Inputs:	
##			@lcfile:			full file path of the land cover raster
## -------------------------------------------------
FM.readLandcover <- function(lcfile, rerun = FALSE)
{
	cat("*********************************************\n")
	cat("********      Land Cover: START      ********\n")
	cat("---------------------------------------------\n")
	cat("Processing Land Cover...")
	if (FM.env$b_debug) {
		cat("*********************************************\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		cat("********      Land Cover: START      ********\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		cat("---------------------------------------------\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		cat("Processing Land Cover...", file = FM.env$s_logfile, sep = "", append = TRUE)
	}
	flush.console()

	#: if this is a rerun, the existing land raster and land txt files will be used
	if (rerun) {
		lcfile <- FM.env$s_land_raster_file
	}
	
	#: check if it is a valid raster
	if (!file.exists(lcfile)) {
		if (FM.env$b_debug) cat("ERROR: the raster file for land cover does NOT exist!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		return("ERROR: the raster file for land cover does NOT exist!\n") #: ERROR: raster file does NOT exist!
	} else {
		#: check unit of the raster, MUST be in meters
		o_tmp_ras <- raster(lcfile)
		s_tmp_str <- toString(o_tmp_ras@crs)
		n_tmp_pos <- regexpr("+units=", s_tmp_str)[1]
		s_tmp_unit <- substr(s_tmp_str, n_tmp_pos + 6, n_tmp_pos + 6)
		n_tmp_nrow <- dim(o_tmp_ras)[1]
		n_tmp_ncol <- dim(o_tmp_ras)[2]
		n_tmp_cellsize <- res(o_tmp_ras)[1] * FM.env$n_dem_cellsize_conversion #: convert into meters
		if (s_tmp_unit != FM.env$s_dem_unit) {
			if (FM.env$b_debug) cat("ERROR: the unit of land cover raster MUST be meters!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
			return("ERROR: the unit of land cover raster MUST be meters!\n") #: ERROR: the unit of land cover raster is NOT meters!
		} else {
			#: check if the layout of grid cells in land cover raster matches that of DEM raster
			if (n_tmp_nrow != FM.env$n_dem_nrow || n_tmp_ncol != FM.env$n_dem_ncol || n_tmp_cellsize != FM.env$n_dem_cellsize) {
				if (FM.env$b_debug) cat("ERROR: the layout of grid cells in land cover raster MUST match that of DEM raster!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
				return("ERROR: the layout of grid cells in land cover raster MUST match that of DEM raster!\n") #: ERROR: DOES NOT match!
			}
			
			#: generate two raster files for parameters of Cd & alpha
			#: read the internal csv file for land cover
			o_tmp_landparams <- read.csv(paste(FM.env$s_run_folder_path, "/input/", FM.env$s_land_param_file, sep = ""), header = TRUE, sep = ",")
			o_tmp_Cd_values <- o_tmp_ras[]
			o_tmp_alpha_values <- o_tmp_ras[]
			for (ir in 1:nrow(o_tmp_landparams))
			{
				o_tmp_Cd_values <- replace(o_tmp_Cd_values, o_tmp_Cd_values == as.numeric(o_tmp_landparams[ir, 1]), as.numeric(o_tmp_landparams[ir, 3]))
				o_tmp_alpha_values <- replace(o_tmp_alpha_values, o_tmp_alpha_values == as.numeric(o_tmp_landparams[ir, 1]), as.numeric(o_tmp_landparams[ir, 4]))
			}
			FM.env$o_land_Cd_matrix <- matrix(o_tmp_Cd_values, FM.env$n_dem_nrow, FM.env$n_dem_ncol, byrow = TRUE)
			FM.env$o_land_alpha_matrix <- matrix(o_tmp_alpha_values, FM.env$n_dem_nrow, FM.env$n_dem_ncol, byrow = TRUE)
			
			#: save them into two raster files under the input folder
			o_tmp_ras[] <- o_tmp_Cd_values
			writeRaster(o_tmp_ras, filename = paste(FM.env$s_run_folder_path, "/input/Cd.tif", sep = ""), format = "GTiff", overwrite = TRUE)
			o_tmp_ras[] <- o_tmp_alpha_values
			writeRaster(o_tmp_ras, filename = paste(FM.env$s_run_folder_path, "/input/alpha.tif", sep = ""), format = "GTiff", overwrite = TRUE)
			
			#: release tmp variables
			rm(o_tmp_landparams)
			rm(o_tmp_Cd_values)
			rm(o_tmp_alpha_values)
		}

		#: release tmp variables
		rm(o_tmp_ras)
		rm(s_tmp_str)
		rm(n_tmp_pos)
		rm(s_tmp_unit)
		rm(n_tmp_nrow)
		rm(n_tmp_ncol)
		rm(n_tmp_cellsize)
		
		#: copy this land cover raster to the input folder and name it as land.tif
		if (!rerun) {
			FM.env$s_land_raster_file <- paste(FM.env$s_run_folder_path, "/input/land.tif", sep = "")
			file.copy(lcfile, FM.env$s_land_raster_file, overwrite = TRUE)
		}
		
		cat("\t>>> OK\n")
		cat("---------------------------------------------\n")
		cat("*****      Land Cover: SUCCESSFUL!      *****\n")
		cat("*********************************************\n")
		if (FM.env$b_debug) {
			cat("\t>>> OK\n", file = FM.env$s_logfile, sep = "", append = TRUE)
			cat("---------------------------------------------\n", file = FM.env$s_logfile, sep = "", append = TRUE)
			cat("*****      Land Cover: SUCCESSFUL!      *****\n", file = FM.env$s_logfile, sep = "", append = TRUE)
			cat("*********************************************\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		}
		flush.console()
		
		return("") #: SUCCESS
	}
}

## [2.4] read rasters for soil texture, depth, and initial moisture content
## -------------------------------------------------
## Name:	FM.readSoil()
## Inputs:	
##			@soilfile:				full file path of the soil texture raster
##			@soildepth:				depth of surface soil (unsaturated soil), unit: m
##			@soildepthfile:			a raster file to store the soil depth of all grid cells
##			@soilmoistcontent:		initial soil moisture content, range: 0.0 - 1.0
##			@soilmoistcontentfile:	a raster file to store the initial soil moisture content of all grid cells
## -------------------------------------------------
FM.readSoil <- function(soilfile, soildepth, soildepthfile = "", soilmoistcontent, soilmoistcontentfile = "", rerun = FALSE)
{
	cat("*********************************************\n")
	cat("************     Soil: START     ************\n")
	cat("---------------------------------------------\n")
	if (FM.env$b_debug) {
		cat("*********************************************\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		cat("************     Soil: START     ************\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		cat("---------------------------------------------\n", file = FM.env$s_logfile, sep = "", append = TRUE)
	}
	flush.console()
	
	#: if this is rerun, the existing soil raster, soil depth raster, and soil moisture files will be used
	if (rerun) {
		soilfile <- FM.env$s_soil_raster_file
		soildepthfile <- FM.env$s_soil_depth_file
		soilmoistcontentfile <- FM.env$s_soil_moistcontent_file
	}
	
	#: check if it is a valid raster
	if (!file.exists(soilfile)) {
		if (FM.env$b_debug) cat("ERROR: the raster file for soil texture does NOT exist!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		return("ERROR: the raster file for soil texture does NOT exist!\n") #: ERROR: raster file does NOT exist!
	} else {
		#: check unit of the raster, MUST be in meters
		o_tmp_ras <- raster(soilfile)
		s_tmp_str <- toString(o_tmp_ras@crs)
		n_tmp_pos <- regexpr("+units=", s_tmp_str)[1]
		s_tmp_unit <- substr(s_tmp_str, n_tmp_pos + 6, n_tmp_pos + 6)
		n_tmp_nrow <- dim(o_tmp_ras)[1]
		n_tmp_ncol <- dim(o_tmp_ras)[2]
		n_tmp_cellsize <- res(o_tmp_ras)[1] * FM.env$n_dem_cellsize_conversion #: convert into meters
		if (s_tmp_unit != FM.env$s_dem_unit) {
			if (FM.env$b_debug) cat("ERROR: the unit of soil texture raster MUST be meters!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
			return("ERROR: the unit of soil texture raster MUST be meters!\n") #: ERROR: the unit of soil texture raster is NOT meters!
		} else {
			#: check if the layout of grid cells in soil texture raster matches that of DEM raster
			if (n_tmp_nrow != FM.env$n_dem_nrow || n_tmp_ncol != FM.env$n_dem_ncol || n_tmp_cellsize != FM.env$n_dem_cellsize) {
				if (FM.env$b_debug) cat("ERROR: the layout of grid cells in soil texture raster MUST match that of DEM raster!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
				return("ERROR: the layout of grid cells in soil texture raster MUST match that of DEM raster!\n") #: ERROR: DOES NOT match!
			}
			
			cat("Processing Soil Texture...")
			if (FM.env$b_debug) cat("Processing Soil Texture...", file = FM.env$s_logfile, sep = "", append = TRUE)
			flush.console()
			
			#: generate two raster files for parameters of porosity & beta
			#: read the internal csv file for soil texture
			o_tmp_soilparams <- read.csv(paste(FM.env$s_run_folder_path, "/input/", FM.env$s_soil_param_file, sep = ""), header = TRUE, sep = ",")
			o_tmp_porosity_values <- o_tmp_ras[]
			o_tmp_beta_values <- o_tmp_ras[]
			for (ir in 1:nrow(o_tmp_soilparams))
			{
				o_tmp_porosity_values <- replace(o_tmp_porosity_values, o_tmp_porosity_values == as.numeric(o_tmp_soilparams[ir, 1]), as.numeric(o_tmp_soilparams[ir, 3]))
				o_tmp_beta_values <- replace(o_tmp_beta_values, o_tmp_beta_values == as.numeric(o_tmp_soilparams[ir, 1]), as.numeric(o_tmp_soilparams[ir, 4]))
			}
			FM.env$o_soil_porosity_matrix <- matrix(o_tmp_porosity_values, FM.env$n_dem_nrow, FM.env$n_dem_ncol, byrow = TRUE)
			FM.env$o_soil_beta_matrix <- matrix(o_tmp_beta_values, FM.env$n_dem_nrow, FM.env$n_dem_ncol, byrow = TRUE)
			
			#: save them into two raster files under the input folder
			o_tmp_ras[] <- o_tmp_porosity_values
			writeRaster(o_tmp_ras, filename = paste(FM.env$s_run_folder_path, "/input/porosity.tif", sep = ""), format = "GTiff", overwrite = TRUE)
			o_tmp_ras[] <- o_tmp_beta_values
			writeRaster(o_tmp_ras, filename = paste(FM.env$s_run_folder_path, "/input/beta.tif", sep = ""), format = "GTiff", overwrite = TRUE)
			
			cat("\t>>> OK\n")
			if (FM.env$b_debug) cat("\t>>> OK\n", file = FM.env$s_logfile, sep = "", append = TRUE)
			flush.console()
			
			#: release tmp variables
			rm(o_tmp_soilparams)
			rm(o_tmp_porosity_values)
			rm(o_tmp_beta_values)
			
			#: soil depth 
			cat("Processing Soil Depth...")
			if (FM.env$b_debug) cat("Processing Soil Depth...", file = FM.env$s_logfile, sep = "", append = TRUE)
			flush.console()
			
			if (soildepthfile != "") {
				#: a raster file is provided
				if (!file.exists(soildepthfile)) {
					if (FM.env$b_debug) cat("ERROR: the raster file for soil depth does NOT exist!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
					return("ERROR: the raster file for soil depth does NOT exist!\n") #: ERROR: raster file does NOT exist!
				} else {
					#: check unit of the raster, MUST be in meters
					o_tmp_sldepth_ras <- raster(soildepthfile)
					s_tmp_sldepth_str <- toString(o_tmp_sldepth_ras@crs)
					n_tmp_sldepth_pos <- regexpr("+units=", s_tmp_sldepth_str)[1]
					s_tmp_sldepth_unit <- substr(s_tmp_sldepth_str, n_tmp_sldepth_pos + 6, n_tmp_sldepth_pos + 6)
					n_tmp_sldepth_nrow <- dim(o_tmp_sldepth_ras)[1]
					n_tmp_sldepth_ncol <- dim(o_tmp_sldepth_ras)[2]
					n_tmp_sldepth_cellsize <- res(o_tmp_sldepth_ras)[1] * FM.env$n_dem_cellsize_conversion #: convert into meters
					if (s_tmp_sldepth_unit != FM.env$s_dem_unit) {
						if (FM.env$b_debug) cat("ERROR: the unit of soil depth raster MUST be meters!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
						return("ERROR: the unit of soil depth raster MUST be meters!\n") #: ERROR: the raster unit is NOT meters!
					} else if (n_tmp_sldepth_nrow != FM.env$n_dem_nrow || n_tmp_sldepth_ncol != FM.env$n_dem_ncol || n_tmp_sldepth_cellsize != FM.env$n_dem_cellsize) {
						if (FM.env$b_debug) cat("ERROR: the layout of grid cells in soil depth raster MUST match that of DEM raster!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
						return("ERROR: the layout of grid cells in soil depth raster MUST match that of DEM raster!\n") #: ERROR: DOES NOT match!
					} else {
						#: save the soil texture matrix
						FM.env$o_soil_depth_matrix <- matrix(o_tmp_sldepth_ras[], FM.env$n_dem_nrow, FM.env$n_dem_ncol, byrow = TRUE)
					}
					#: release tmp variables
					rm(o_tmp_sldepth_ras)
					rm(s_tmp_sldepth_str)
					rm(n_tmp_sldepth_pos)
					rm(s_tmp_sldepth_unit)
					rm(n_tmp_sldepth_nrow)
					rm(n_tmp_sldepth_ncol)
					rm(n_tmp_sldepth_cellsize)
				}
				
				#: copy this soil depth raster to the input folder and name it as soildepth.tif (only if this is a new run)
				if (!rerun) {
					FM.env$s_soil_depth_file <- paste(FM.env$s_run_folder_path, "/input/soildepth.tif", sep = "")
					file.copy(soildepthfile, FM.env$s_soil_depth_file, overwrite = TRUE)
				}
				
			} else {
				#: apply the value of soildepth to all grid cells
				FM.env$o_soil_depth_matrix <- matrix(soildepth, FM.env$n_dem_nrow, FM.env$n_dem_ncol, byrow = TRUE)
				o_tmp_ras[] <- soildepth
				FM.env$s_soil_depth_file <- paste(FM.env$s_run_folder_path, "/input/soildepth.tif", sep = "")
				writeRaster(o_tmp_ras, filename = FM.env$s_soil_depth_file, format = "GTiff", overwrite = TRUE)
			}
			
			cat("\t>>> OK\n")
			if (FM.env$b_debug) cat("\t>>> OK\n", file = FM.env$s_logfile, sep = "", append = TRUE)
			flush.console()
			#: calculate the total water volume while the surface soil is saturated
			FM.env$o_soil_maxwaterVol_matrix <- FM.env$n_dem_cellsize ^ 2 * FM.env$o_soil_depth_matrix * FM.env$o_soil_porosity_matrix																													   

			#: soil moisture content 
			cat("Processing Soil Moisture...")
			if (FM.env$b_debug) cat("Processing Soil Moisture...", file = FM.env$s_logfile, sep = "", append = TRUE)
			flush.console()
			
			if (soilmoistcontentfile != "") {
				#: a raster file is provided
				if (!file.exists(soilmoistcontentfile)) {
					if (FM.env$b_debug) cat("ERROR: the raster file for soil moisture content does NOT exist!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
					return("ERROR: the raster file for soil moisture content does NOT exist!\n") #: ERROR: raster file does NOT exist!
				} else {
					#: check unit of the raster, MUST be in meters
					o_tmp_slmoist_ras <- raster(soilmoistcontentfile)
					s_tmp_slmoist_str <- toString(o_tmp_slmoist_ras@crs)
					n_tmp_slmoist_pos <- regexpr("+units=", s_tmp_slmoist_str)[1]
					s_tmp_slmoist_unit <- substr(s_tmp_slmoist_str, n_tmp_slmoist_pos + 6, n_tmp_slmoist_pos + 6)
					n_tmp_slmoist_nrow <- dim(o_tmp_slmoist_ras)[1]
					n_tmp_slmoist_ncol <- dim(o_tmp_slmoist_ras)[2]
					n_tmp_slmoist_cellsize <- res(o_tmp_slmoist_ras)[1] * FM.env$n_dem_cellsize_conversion #: convert into meters
					if (s_tmp_slmoist_unit != FM.env$s_dem_unit) {
						if (FM.env$b_debug) cat("ERROR: the unit of soil moisture content raster MUST be meters!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
						return("ERROR: the unit of soil moisture content raster MUST be meters!\n") #: ERROR: the raster unit is NOT meters!
					} else if (n_tmp_slmoist_nrow != FM.env$n_dem_nrow || n_tmp_slmoist_ncol != FM.env$n_dem_ncol || n_tmp_slmoist_cellsize != FM.env$n_dem_cellsize) {
						if (FM.env$b_debug) cat("ERROR: the layout of grid cells in soil moisture content raster MUST match that of DEM raster!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
						return("ERROR: the layout of grid cells in soil moisture content raster MUST match that of DEM raster!\n") #: ERROR: DOES NOT match!
					} else {
						#: save the soil texture matrix
						FM.env$o_soil_moistcontent_matrix <- matrix(o_tmp_slmoist_ras[], FM.env$n_dem_nrow, FM.env$n_dem_ncol, byrow = TRUE)																					
						#: if impermeable surface or water surface, the max water volume is zero, no infiltration is accepted!
						#: update the soil moisture to be saturated
						FM.env$o_soil_moistcontent_matrix[which(ceiling(FM.env$o_soil_maxwaterVol_matrix * FM.env$n_integer_factor) <= 0)] <- 1
						o_tmp_ras[] <- as.vector(t(FM.env$o_soil_moistcontent_matrix))
						FM.env$s_soil_moistcontent_file <- paste(FM.env$s_run_folder_path, "/input/soilmoist.tif", sep = "")
						writeRaster(o_tmp_ras, filename = FM.env$s_soil_moistcontent_file, format = "GTiff", overwrite = TRUE)																		
					}
					#: release tmp variables
					rm(o_tmp_slmoist_ras)
					rm(s_tmp_slmoist_str)
					rm(n_tmp_slmoist_pos)
					rm(s_tmp_slmoist_unit)
					rm(n_tmp_slmoist_nrow)
					rm(n_tmp_slmoist_ncol)
					rm(n_tmp_slmoist_cellsize)
				}
				
			} else {
				#: apply the value of soilmoistcontent to all grid cells
				FM.env$o_soil_moistcontent_matrix <- matrix(soilmoistcontent, FM.env$n_dem_nrow, FM.env$n_dem_ncol, byrow = TRUE)																						  
				#: if impermeable surface or water surface, the max water volume is zero, no infiltration is accepted!
				#: update the soil moisture to be saturated
				FM.env$o_soil_moistcontent_matrix[which(ceiling(FM.env$o_soil_maxwaterVol_matrix * FM.env$n_integer_factor) <= 0)] <- 1
				o_tmp_ras[] <- as.vector(t(FM.env$o_soil_moistcontent_matrix))
				FM.env$s_soil_moistcontent_file <- paste(FM.env$s_run_folder_path, "/input/soilmoist.tif", sep = "")
				writeRaster(o_tmp_ras, filename = FM.env$s_soil_moistcontent_file, format = "GTiff", overwrite = TRUE)
			}
			cat("\t>>> OK\n")
			if (FM.env$b_debug) cat("\t>>> OK\n", file = FM.env$s_logfile, sep = "", append = TRUE)	
			flush.console()
		}
		
		#: release tmp variables
		rm(o_tmp_ras)
		rm(s_tmp_str)
		rm(n_tmp_pos)
		rm(s_tmp_unit)
		rm(n_tmp_nrow)
		rm(n_tmp_ncol)
		rm(n_tmp_cellsize)
		
		#: copy this land cover raster to the input folder and name it as land.tif (only if this is a new run)
		if (!rerun) {
			FM.env$s_soil_raster_file <- paste(FM.env$s_run_folder_path, "/input/soil.tif", sep = "")
			file.copy(soilfile, FM.env$s_soil_raster_file, overwrite = TRUE)
		}
		
		cat("---------------------------------------------\n")
		cat("**********    Soil: SUCCESSFUL!     *********\n")
		cat("*********************************************\n")
		if (FM.env$b_debug) {
			cat("---------------------------------------------\n", file = FM.env$s_logfile, sep = "", append = TRUE)
			cat("**********    Soil: SUCCESSFUL!     *********\n", file = FM.env$s_logfile, sep = "", append = TRUE)
			cat("*********************************************\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		}
		flush.console()
		
		return("") #: SUCCESS
	}
}
									
## [2.5] read drain inlets (OPTIONAL)
## -------------------------------------------------
## Name:	FM.readDraininlets()
## Inputs:	
##			@draininletfile:		full file path of the drain inlet raster
## -------------------------------------------------
FM.readDraininlets <- function(draininletfile, rerun = FALSE)
{
	cat("*********************************************\n")
	cat("********     Drain Inlets: START     ********\n")
	cat("---------------------------------------------\n")
	cat("Processing...")
	if (FM.env$b_debug) {
		cat("*********************************************\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		cat("********     Drain Inlets: START     ********\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		cat("---------------------------------------------\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		cat("Processing...", file = FM.env$s_logfile, sep = "", append = TRUE)
	}
	flush.console()
	
	#: if this is a rerun, the existing drain inlet raster will be used
	if (rerun) draininletfile <- FM.env$s_draininlet_raster_file
	
	#: check if it is a valid raster
	if (!file.exists(draininletfile)) {
		if (FM.env$b_debug) cat("ERROR: the raster file for drain inlet does NOT exist!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		return("ERROR: the raster file for drain inlet does NOT exist!\n") #: ERROR: raster file does NOT exist!
	} else {
		#: check unit of the raster, MUST be in meters
		o_tmp_ras <- raster(draininletfile)
		s_tmp_str <- toString(o_tmp_ras@crs)
		n_tmp_pos <- regexpr("+units=", s_tmp_str)[1]
		s_tmp_unit <- substr(s_tmp_str, n_tmp_pos + 6, n_tmp_pos + 6)
		n_tmp_nrow <- dim(o_tmp_ras)[1]
		n_tmp_ncol <- dim(o_tmp_ras)[2]
		n_tmp_cellsize <- res(o_tmp_ras)[1] * FM.env$n_dem_cellsize_conversion #: convert into meters
		if (s_tmp_unit != FM.env$s_dem_unit) {
			if (FM.env$b_debug) cat("ERROR: the unit of drain inlet raster MUST be meters!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
			return("ERROR: the unit of drain inlet raster MUST be meters!\n") #: ERROR: the unit is NOT meters!
		} else {
			#: check if the layout of grid cells in drain inlet raster matches that of DEM raster
			if (n_tmp_nrow != FM.env$n_dem_nrow || n_tmp_ncol != FM.env$n_dem_ncol || n_tmp_cellsize != FM.env$n_dem_cellsize) {
				if (FM.env$b_debug) cat("ERROR: the layout of grid cells in drain inlet raster MUST match that of DEM raster!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
				return("ERROR: the layout of grid cells in drain inlet raster MUST match that of DEM raster!\n") #: ERROR: DOES NOT match!
			}
			FM.env$n_draininlet_flag <- 1
			#: unit conversion: mm/hr -> m/sec
			o_tmp_ras[] <- o_tmp_ras[] * 0.001 / 3600
			#: save this into the matrix
			FM.env$o_draininlet_maxI_matrix <- matrix(o_tmp_ras[], FM.env$n_dem_nrow, FM.env$n_dem_ncol, byrow = TRUE)
			cat("\t>>> OK\n")
			if (FM.env$b_debug) cat("\t>>> OK\n", file = FM.env$s_logfile, sep = "", append = TRUE)
			flush.console()
		}

		#: release tmp variables
		rm(o_tmp_ras)
		rm(s_tmp_str)
		rm(n_tmp_pos)
		rm(s_tmp_unit)
		rm(n_tmp_nrow)
		rm(n_tmp_ncol)
		rm(n_tmp_cellsize)
		
		#: copy this drain inlet raster to the input folder and name it as draininlet.tif (only if this is a new run)
		if (!rerun) {
			FM.env$s_draininlet_raster_file <- paste(FM.env$s_run_folder_path, "/input/draininlet.tif", sep = "")
			file.copy(draininletfile, FM.env$s_draininlet_raster_file, overwrite = TRUE)
		}
		
		cat("---------------------------------------------\n")
		cat("*****     Drain Inlets: SUCCESSFUL!     *****\n")
		cat("*********************************************\n")
		if (FM.env$b_debug) {
			cat("---------------------------------------------\n", file = FM.env$s_logfile, sep = "", append = TRUE)
			cat("*****     Drain Inlets: SUCCESSFUL!     *****\n", file = FM.env$s_logfile, sep = "", append = TRUE)
			cat("*********************************************\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		}
		flush.console()

		return("") #: SUCCESS
	}
}

## [2.6] read precipitation data
## -------------------------------------------------
## Name:	FM.readPrecip()
## Inputs:	
##			@datatype:			0 = weather stations, 1 = raster data from NWP models or reanalysis datasets
##			@txtfile:			full file path of the text file required for precipitation data
##			@startdatetime:		start datetime of your precipitation data records, format: "2016-06-09 12:00:00"
##			@enddatetime:		end datetime of your precipitation data records, format: "2016-06-09 13:00:00"
##			@interval:			interval of your precipitation data records, unit: seconds, e.g.: 1 hr = 60 min * 60 sec = 3600 sec
##			@idwpower:			the power of IDW method, default: 2, suggested options: {1, 2, 3}
## -------------------------------------------------
FM.readPrecip <- function(datatype = 0, txtfile, startdatetime, enddatetime, interval, idwpower = 2)
{
	cat("*********************************************\n")
	cat("********    Precipitation: START     ********\n")
	cat("---------------------------------------------\n")
	if (FM.env$b_debug) {
		cat("*********************************************\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		cat("********    Precipitation: START     ********\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		cat("---------------------------------------------\n", file = FM.env$s_logfile, sep = "", append = TRUE)
	}
	flush.console()
	
	#: if this is a valid inteval
	if (interval < FM.env$n_run_internal_timestep) {
		if (FM.env$b_debug) cat("ERROR: invalid precipitation interval!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		return("ERROR: invalid precipitation interval!\n")
	}
	
	FM.env$n_precip_datatype <- datatype
	FM.env$s_precip_txt_file <- txtfile
	FM.env$o_precip_start_datetime <- as.POSIXct(startdatetime)
	FM.env$o_precip_end_datetime <- as.POSIXct(enddatetime)
	FM.env$n_precip_record_interval <- as.integer(interval)
	FM.env$n_precip_IDW_power <- as.integer(idwpower)
	#: check datetime
	if (FM.env$o_precip_end_datetime <= FM.env$o_precip_start_datetime) {
		if (FM.env$b_debug) cat("ERROR: the end datetime should be after the start datetime!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		return("ERROR: the end datetime should be after the start datetime!\n")
	}
	n_tmp_totalseconds <- as.numeric(difftime(FM.env$o_precip_end_datetime, FM.env$o_precip_start_datetime, units = "secs"))
	if (n_tmp_totalseconds %% FM.env$n_precip_record_interval != 0) {
		if (FM.env$b_debug) cat("ERROR: the interval of data records does NOT match the end datetime!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		return("ERROR: the interval of data records does NOT match the end datetime!\n")
	}
	n_tmp_expected_records <- n_tmp_totalseconds / FM.env$n_precip_record_interval
	
	#: process precipitation data
	if (datatype == 0){
		#: data from weather stations
		#: the format of this text file is given as follows:
		#: ------------------------------------------
		#: ID	Lon			Lat			Datafile
		#: 1	-92.0296	30.1800		p1.txt		(stores obs data for station 1)
		#: 2	-92.0260	30.1779		p2.txt		(stores obs data for station 2)
		#: ------------------------------------------
		o_tmp_P_stations <- read.table(FM.env$s_precip_txt_file, header = TRUE)
		n_tmp_totalstations <- nrow(o_tmp_P_stations)
		
		#: check its format
		if (ncol(o_tmp_P_stations) != 4) {
			if (FM.env$b_debug) cat("ERROR: the format of precipitation text file is NOT correct!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
			return("ERROR: the format of precipitation text file is NOT correct!\n") #: ERROR: format is not correct.
		}
		
		#: get the longitude and latitude of all grid cells
		o_tmp_dem_ras <- raster(FM.env$s_dem_raster_file)
		o_tmp_ras_to_points <- rasterToPoints(o_tmp_dem_ras, spatial = TRUE)
		o_tmp_ras_to_points <- spTransform(o_tmp_ras_to_points, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
		o_tmp_points_attrtable <- as.data.frame(o_tmp_ras_to_points)
		
		#: calculate distances
		o_tmp_distances <- matrix(0, nrow(o_tmp_points_attrtable), n_tmp_totalstations)
		o_tmp_all_Pdata <- matrix(0, n_tmp_expected_records, n_tmp_totalstations)
		o_tmp_lons_x <- as.numeric(o_tmp_points_attrtable$x)
		o_tmp_lats_y <- as.numeric(o_tmp_points_attrtable$y)
		for (ip in 1:n_tmp_totalstations) {
			#: save all precip records into a matrix
			
			#: for package release, updated on May 14, 2018, Xiamen, China
			s_tmp_tmp_tmp_stationfile <- toString(o_tmp_P_stations[ip, 4])
			#: check if it is a valid file
			if (!file.exists(s_tmp_tmp_tmp_stationfile)) {
				s_tmp_tmp_tmp_stationfile <- system.file("extdata", s_tmp_tmp_tmp_stationfile, package = "FloodMapper", mustWork = TRUE)
			}
			
			o_tmp_tmp_data <- read.table(s_tmp_tmp_tmp_stationfile, header = FALSE)
			if (nrow(o_tmp_tmp_data) != n_tmp_expected_records) {
				if (FM.env$b_debug) cat("ERROR: the number of precipitation records does NOT match your input!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
				return("ERROR: the number of precipitation records does NOT match your input!\n")
			}
			o_tmp_all_Pdata[ , ip] <- as.numeric(o_tmp_tmp_data[ , 1])

			#: read the lon and lat of each station, the range of lat must be [-90, +90], the range of lon must be [-180, +180]
			n_tmp_p_lon <- as.numeric(o_tmp_P_stations[ip, 2])
			n_tmp_p_lat <- as.numeric(o_tmp_P_stations[ip, 3])
			if (n_tmp_p_lon > 180.0) {
				n_tmp_p_lon <- n_tmp_p_lon - 360.0
			} else if (n_tmp_p_lon < -180.0) {
				n_tmp_p_lon <- n_tmp_p_lon + 360.0
			}

			#: calculate the distances of all grid cells to each station
			o_tmp_distances[ , ip] <- sqrt((o_tmp_lons_x - n_tmp_p_lon) ^ 2 + (o_tmp_lats_y - n_tmp_p_lat) ^ 2)
		}
		#: in case that the distance is zero, the IDW method won't work! Add a very small number to avoid zero.
		o_tmp_distances <- o_tmp_distances + 0.000000001  

		#: generate one raster for each time step of precipitation data records
		for (ir in 1:n_tmp_expected_records) {
			o_tmp_dem_ras[] <- 0
			n_tmp_denominator <- rep(0, nrow(o_tmp_points_attrtable))
			n_tmp_numerator <- rep(0, nrow(o_tmp_points_attrtable))
			for (ipp in 1:n_tmp_totalstations)
			{
				n_tmp_denominator <- n_tmp_denominator + o_tmp_all_Pdata[ir, ipp] / as.numeric(o_tmp_distances[ , ipp]) ^ FM.env$n_precip_IDW_power
				n_tmp_numerator <- n_tmp_numerator + 1 / as.numeric(o_tmp_distances[ , ipp]) ^ FM.env$n_precip_IDW_power
			}
			o_tmp_dem_ras[FM.env$o_dem_nonNA_index] <- n_tmp_denominator / n_tmp_numerator
			#: print the current datetime
			n_tmp_currentdatetime <- FM.env$o_precip_start_datetime + (ir - 1) * FM.env$n_precip_record_interval
			cat(format(n_tmp_currentdatetime, "%Y-%m-%d %H:%M:%S"))
			if (FM.env$b_debug) cat(format(n_tmp_currentdatetime, "%Y-%m-%d %H:%M:%S"), file = FM.env$s_logfile, sep = "", append = TRUE)
			flush.console()
			s_tmp_ras_filename <- format(n_tmp_currentdatetime, "%Y-%m-%d %H:%M:%S")
			s_tmp_ras_filename <- gsub(" ", "_", s_tmp_ras_filename)
			s_tmp_ras_filename <- gsub(":", "-", s_tmp_ras_filename)
			s_tmp_ras_filename <- paste(FM.env$s_run_folder_path, "/input/p_", s_tmp_ras_filename, ".tif", sep = "")
			writeRaster(o_tmp_dem_ras, filename = s_tmp_ras_filename, format = "GTiff", overwrite = TRUE)
			cat("\t>>> OK\n")
			if (FM.env$b_debug) cat("\t>>> OK\n", file = FM.env$s_logfile, sep = "", append = TRUE)
			flush.console()
		}

		#: release tmp variables
		rm(o_tmp_P_stations)
		rm(n_tmp_totalstations)
		rm(o_tmp_dem_ras)
		rm(o_tmp_ras_to_points)
		rm(o_tmp_points_attrtable)
		rm(o_tmp_distances)
		rm(o_tmp_all_Pdata)
		rm(o_tmp_lons_x)
		rm(o_tmp_lats_y)
		
	} else {
		#: data from rasters
		#: the format of this text file is given as follows:
		#: ------------------------------------------
		#: ID	Rasterfile
		#: 1	p_2016-06-09_12-00-00.tif	(although a user-defined filename is allowed)
		#: 2	p_2016-06-09_13-00-00.tif	(it is suggested to name it with date and time)
		#: ------------------------------------------
		o_tmp_P_rasters <- read.table(FM.env$s_precip_txt_file, header = TRUE)
		n_tmp_P_totalrecords <- nrow(o_tmp_P_rasters)
		
		#: check if the record length matches with the defined one
		if (n_tmp_expected_records != n_tmp_P_totalrecords) {
			if (FM.env$b_debug) cat("ERROR: the number of precipitation raster files does NOT match your input!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
			return("ERROR: the number of precipitation raster files does NOT match your input!\n")
		}
		
		#: check if each raster file exists and if its format is correct
		for (ipr in 1:n_tmp_P_totalrecords) {
			s_tmp_p_ras_file <- toString(o_tmp_P_rasters[ipr, 2])
			
			#: check if it is a valid file
			if (!file.exists(s_tmp_p_ras_file)) {
				if (FM.env$b_debug) cat("ERROR: the precipitation raster file does NOT exist!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
				return("ERROR: the precipitation raster file does NOT exist!\n")
			}
			
			#: check unit of the raster, MUST be in meters
			o_tmp_ras <- raster(s_tmp_p_ras_file)
			s_tmp_str <- toString(o_tmp_ras@crs)
			n_tmp_pos <- regexpr("+units=", s_tmp_str)[1]
			s_tmp_unit <- substr(s_tmp_str, n_tmp_pos + 6, n_tmp_pos + 6)
			n_tmp_nrow <- dim(o_tmp_ras)[1]
			n_tmp_ncol <- dim(o_tmp_ras)[2]
			n_tmp_cellsize <- res(o_tmp_ras)[1] * FM.env$n_dem_cellsize_conversion #: convert into meters
			if (s_tmp_unit != FM.env$s_dem_unit) {
				if (FM.env$b_debug) cat("ERROR: the unit of precipitation raster MUST be meters!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
				return("ERROR: the unit of precipitation raster MUST be meters!\n")
			}
			
			#: check if the layout of grid cells matches that of DEM raster
			if (n_tmp_nrow != FM.env$n_dem_nrow || n_tmp_ncol != FM.env$n_dem_ncol || n_tmp_cellsize != FM.env$n_dem_cellsize) {
				if (FM.env$b_debug) cat("ERROR: the layout of grid cells in precipitation raster MUST match that of DEM raster!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
				return("ERROR: the layout of grid cells in precipitation raster MUST match that of DEM raster!\n")
			}
			
			#: copy this raster to the input folder
			#: print the current datetime
			n_tmp_currentdatetime <- FM.env$o_precip_start_datetime + (ipr - 1) * FM.env$n_precip_record_interval
			cat(format(n_tmp_currentdatetime, "%Y-%m-%d %H:%M:%S"))
			if (FM.env$b_debug) cat(format(n_tmp_currentdatetime, "%Y-%m-%d %H:%M:%S"), file = FM.env$s_logfile, sep = "", append = TRUE)
			flush.console()
			s_tmp_ras_filename <- format(n_tmp_currentdatetime, "%Y-%m-%d %H:%M:%S")
			s_tmp_ras_filename <- gsub(" ", "_", s_tmp_ras_filename)
			s_tmp_ras_filename <- gsub(":", "-", s_tmp_ras_filename)
			s_tmp_ras_filename <- paste(FM.env$s_run_folder_path, "/input/p_", s_tmp_ras_filename, ".tif", sep = "")
			file.copy(s_tmp_p_ras_file, s_tmp_ras_filename, overwrite = TRUE)
			cat("\t>>> OK\n")
			if (FM.env$b_debug) cat("\t>>> OK\n", file = FM.env$s_logfile, sep = "", append = TRUE)
			flush.console()
		}
	}

	cat("---------------------------------------------\n")
	cat("*****    Precipitation: SUCCESSFUL!     *****\n")
	cat("*********************************************\n")
	if (FM.env$b_debug) {
		cat("---------------------------------------------\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		cat("*****    Precipitation: SUCCESSFUL!     *****\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		cat("*********************************************\n", file = FM.env$s_logfile, sep = "", append = TRUE)
	}
	flush.console()
	
	return("") #: SUCCESS
}

## [2.7] read infow data
##		 this allows users to add inflows to the domain
## -------------------------------------------------
## Name:	FM.readInflow()
## Inputs:	
##			@datatype:			0 = monitoring stations, 1 = raster data (NOT supported in this version)
##			@txtfile:			full file path of the text file required for inflow data
##			@startdatetime:		start datetime of your inflow data records, format: "2016-06-09 12:00:00"
##			@enddatetime:		end datetime of your inflow data records, format: "2016-06-09 13:00:00"
##			@interval:			interval of your inflow data records, unit: seconds, e.g.: 1 hr = 60 min * 60 sec = 3600 sec
## -------------------------------------------------
FM.readInflow <- function(datatype = 0, txtfile, startdatetime, enddatetime, interval)
{
	cat("*********************************************\n")
	cat("***********    Inflow: START     ************\n")
	cat("---------------------------------------------\n")
	if (FM.env$b_debug) {
		cat("*********************************************\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		cat("***********    Inflow: START     ************\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		cat("---------------------------------------------\n", file = FM.env$s_logfile, sep = "", append = TRUE)
	}
	flush.console()
	
	#: if this is a valid inteval
	if (interval < FM.env$n_run_internal_timestep) {
		if (FM.env$b_debug) cat("ERROR: invalid inflow interval!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		return("ERROR: invalid inflow interval!\n")
	}
	
	FM.env$n_inflow_datatype <- datatype
	FM.env$s_inflow_txt_file <- txtfile
	FM.env$o_inflow_start_datetime <- as.POSIXct(startdatetime)
	FM.env$o_inflow_end_datetime <- as.POSIXct(enddatetime)
	FM.env$n_inflow_record_interval <- as.integer(interval)
	#: check datetime
	if (FM.env$o_inflow_end_datetime <= FM.env$o_inflow_start_datetime) {
		if (FM.env$b_debug) cat("ERROR: the end datetime should be after the start datetime!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		return("ERROR: the end datetime should be after the start datetime!\n")
	}
	n_tmp_totalseconds <- as.numeric(difftime(FM.env$o_inflow_end_datetime, FM.env$o_inflow_start_datetime, units = "secs"))
	if (n_tmp_totalseconds %% FM.env$n_inflow_record_interval != 0) {
		if (FM.env$b_debug) cat("ERROR: the interval of data records does NOT match the end datetime!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		return("ERROR: the interval of data records does NOT match the end datetime!\n")
	}
	n_tmp_expected_records <- n_tmp_totalseconds / FM.env$n_inflow_record_interval
	
	#: process inflow data
	if (datatype == 0){
		#: data from monitoring stations
		#: the format of this text file is given as follows:
		#: ------------------------------------------
		#: ID	Lon			Lat			Datafile
		#: 1	-92.0296	30.1800		in1.txt		(stores obs data for station 1)
		#: 2	-92.0260	30.1779		in2.txt		(stores obs data for station 2)
		#: ------------------------------------------
		o_tmp_inflow_stations <- read.table(FM.env$s_inflow_txt_file, header = TRUE)
		n_tmp_totalstations <- nrow(o_tmp_inflow_stations)
		
		#: check its format
		if (ncol(o_tmp_inflow_stations) != 4) {
			if (FM.env$b_debug) cat("ERROR: the format of inflow text file is NOT correct!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
			return("ERROR: the format of inflow text file is NOT correct!\n") #: ERROR: format is not correct.
		}
		
		#: get the longitude and latitude of all grid cells
		o_tmp_dem_ras <- raster(FM.env$s_dem_raster_file)
		o_tmp_ras_to_points <- rasterToPoints(o_tmp_dem_ras, spatial = TRUE)
		o_tmp_ras_to_points <- spTransform(o_tmp_ras_to_points, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
		o_tmp_points_attrtable <- as.data.frame(o_tmp_ras_to_points)

		#: find the nearest grid cell for each station
		o_tmp_all_inflow_data <- matrix(0, n_tmp_expected_records, n_tmp_totalstations)
		o_tmp_lons_x <- as.numeric(o_tmp_points_attrtable$x)
		o_tmp_lats_y <- as.numeric(o_tmp_points_attrtable$y)
		o_tmp_nearest_vectorindex <- rep(0, n_tmp_totalstations)
		for (ip in 1:n_tmp_totalstations) {
			#: save all inflow records into a matrix
			
			#: for package release, updated on May 23, 2018, Lafayette, LA, USA
			s_tmp_tmp_tmp_stationfile <- toString(o_tmp_inflow_stations[ip, 4])
			#: check if it is a valid file
			if (!file.exists(s_tmp_tmp_tmp_stationfile)) {
				s_tmp_tmp_tmp_stationfile <- system.file("extdata", s_tmp_tmp_tmp_stationfile, package = "FloodMapper", mustWork = TRUE)
			}
			
			o_tmp_tmp_data <- read.table(s_tmp_tmp_tmp_stationfile, header = FALSE)
			if (nrow(o_tmp_tmp_data) != n_tmp_expected_records) {
				if (FM.env$b_debug) cat("ERROR: the number of inflow records does NOT match your input!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
				return("ERROR: the number of inflow records does NOT match your input!\n")
			}
			o_tmp_all_inflow_data[ , ip] <- as.numeric(o_tmp_tmp_data[ , 1])

			#: read the lon and lat of each station, the range of lat must be [-90, +90], the range of lon must be [-180, +180]
			n_tmp_p_lon <- as.numeric(o_tmp_inflow_stations[ip, 2])
			n_tmp_p_lat <- as.numeric(o_tmp_inflow_stations[ip, 3])
			if (n_tmp_p_lon > 180.0) {
				n_tmp_p_lon <- n_tmp_p_lon - 360.0
			} else if (n_tmp_p_lon < -180.0) {
				n_tmp_p_lon <- n_tmp_p_lon + 360.0
			}

			#: calculate the distances of all grid cells to each station
			o_tmp_distance <- ceiling(sqrt((o_tmp_lons_x - n_tmp_p_lon) ^ 2 + (o_tmp_lats_y - n_tmp_p_lat) ^ 2) * 1000000) #: convert an integer for comparison
			#: find the nearest grid cell
			o_tmp_nearest_vectorindex[ip] <- which(o_tmp_distance == min(o_tmp_distance))[1] #: in the case of a tie, choose the first one
		}
		
		#: generate one raster for each time step of inflow data records
		for (ir in 1:n_tmp_expected_records) {
			o_tmp_dem_ras[] <- 0
			o_tmp_dem_ras[o_tmp_nearest_vectorindex] <- o_tmp_all_inflow_data[ir, 1:n_tmp_totalstations]
			#: print the current datetime
			n_tmp_currentdatetime <- FM.env$o_inflow_start_datetime + (ir - 1) * FM.env$n_inflow_record_interval
			cat(format(n_tmp_currentdatetime, "%Y-%m-%d %H:%M:%S"))
			if (FM.env$b_debug) cat(format(n_tmp_currentdatetime, "%Y-%m-%d %H:%M:%S"), file = FM.env$s_logfile, sep = "", append = TRUE)
			flush.console()
			s_tmp_ras_filename <- format(n_tmp_currentdatetime, "%Y-%m-%d %H:%M:%S")
			s_tmp_ras_filename <- gsub(" ", "_", s_tmp_ras_filename)
			s_tmp_ras_filename <- gsub(":", "-", s_tmp_ras_filename)
			s_tmp_ras_filename <- paste(FM.env$s_run_folder_path, "/input/in_", s_tmp_ras_filename, ".tif", sep = "")
			writeRaster(o_tmp_dem_ras, filename = s_tmp_ras_filename, format = "GTiff", overwrite = TRUE)
			cat("\t>>> OK\n")
			if (FM.env$b_debug) cat("\t>>> OK\n", file = FM.env$s_logfile, sep = "", append = TRUE)
			flush.console()
		}

		#: release tmp variables
		rm(o_tmp_inflow_stations)
		rm(n_tmp_totalstations)
		rm(o_tmp_dem_ras)
		rm(o_tmp_ras_to_points)
		rm(o_tmp_points_attrtable)
		rm(o_tmp_all_inflow_data)
		rm(o_tmp_lons_x)
		rm(o_tmp_lats_y)
	}

	cat("---------------------------------------------\n")
	cat("********    Inflow: SUCCESSFUL!     *********\n")
	cat("*********************************************\n")
	if (FM.env$b_debug) {
		cat("---------------------------------------------\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		cat("********    Inflow: SUCCESSFUL!     *********\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		cat("*********************************************\n", file = FM.env$s_logfile, sep = "", append = TRUE)
	}
	flush.console()
	
	return("") #: SUCCESS
}

## [2.8] calculation of outflow direction & instantaneous slope [INTERNAL]
## -------------------------------------------------
## Name:	FM.calOutflowDir()
## Inputs:	
##			@o_instHeights:		a 9-element vector containing 9 instantaneous neights of a given cell and its immediate neighbors
##								format: c = (L0, L1, L2, L3, L4, L5, L6, L7, L8)
## Outputs:
##			@o_outputs:			a 4-element vector containing: outflow direction, direction flag, instantaneous slope, and flag of an edge/corner cell
##								1> outflow direction: {1, 2, 3, 4, 5, 6, 7, 8}
##								2> direction flag: 1 = to parallel cell, 2 = to diagonal cell
##								3> instantaneous slope: calculated from instantaneous heights rather than elevation, unit: %
##								4> flag of an edge/corner cell: 0 = normal cell, a positive value (1 to 8) = index of the highest neighbor within the domain
## -------------------------------------------------
FM.calOutflowDir <- function(o_instHeights)
{
	o_outputs <- c(0, 0, 0)

	#: find the grid cell with the lowest instantaneous height
	n_tmp_minindex <- which(o_instHeights == min(o_instHeights, na.rm = TRUE))
	n_tmp_maxindex <- which(o_instHeights == max(o_instHeights, na.rm = TRUE))
	if (length(n_tmp_minindex) > 1) n_tmp_minindex <- n_tmp_minindex[1] #: might be more than one cells with the same height
	if (length(n_tmp_maxindex) > 1) n_tmp_maxindex <- n_tmp_maxindex[1]
	
	if (n_tmp_maxindex == n_tmp_minindex){
		#: completely flat, no outflow
		o_outputs[1] <- 0
		o_outputs[2] <- 0
		o_outputs[3] <- 0
	} else if (n_tmp_minindex == 1) {
		#: the center cell is the lowest one
		#: check if this is an edge or corner cell by calculating mean(o_instHeights) without removing NA
		if(is.na(mean(o_instHeights))) {
			#: yes, it is an edge or corner cell, use the following rule to determine its exit flow direction
			#: >> find the highest cell and use its direction to the center cell to represent the exit flow direction:
			#		6  7  8
			#		5  0  1
			#       4  3  2
			# highest cell -> center cell, then center cell -> out-domain cell
			# 1 -> 0 then 0 -> 5
			# 2 -> 0 then 0 -> 6
			# 3 -> 0 then 0 -> 7
			# 4 -> 0 then 0 -> 8
			# 5 -> 0 then 0 -> 1
			# 6 -> 0 then 0 -> 2
			# 7 -> 0 then 0 -> 3
			# 8 -> 0 then 0 -> 4
			o_exit_directions <- c(5, 6, 7, 8, 1, 2, 3, 4)
			o_outputs[1] <- o_exit_directions[n_tmp_maxindex - 1]
			
			if (o_outputs[1] == 1 || o_outputs[1] == 3 || o_outputs[1] == 5 || o_outputs[1] == 7) {
				#: to parallel cells
				o_outputs[2] <- 1
			} else {
				#: to diagonal cells
				o_outputs[2] <- 2
			}
			o_outputs[3] <- n_tmp_maxindex - 1
		} else {
			#: this cell is a pit, no outflow
			o_outputs[1] <- 0
			o_outputs[2] <- 0
			o_outputs[3] <- 0
		}
	} else {
		#: normal case
		o_outputs[1] <- n_tmp_minindex - 1
		if (o_outputs[1] == 1 || o_outputs[1] == 3 || o_outputs[1] == 5 || o_outputs[1] == 7) {
			#: to parallel cells
			o_outputs[2] <- 1
		} else {
			#: to diagonal cells
			o_outputs[2] <- 2
		}
		
		o_outputs[3] <- 0
	}
	
	return(o_outputs)
}

## [2.9] determine the row and col indices of downstream cell 
##		 according to the outflow direction of current grid cell  [INTERNAL]
## -------------------------------------------------
## Name:	FM.convertOutflowDirtoRowCol()
## Inputs:	
##			@rowindex:			row index of current grid cell
##			@colindex:			col index of current grid cell
##			@outflowdir:		outflow direction of current grid cell
## Outputs:
##			@o_outputs:			a 2-element vector containing the row and col indices of the downstream cell
##								format: c(rowindex, colindex)
## -------------------------------------------------
FM.convertOutflowDirtoRowCol <- function(rowindex, colindex, outflowdir)
{
	o_outputs <- c(0, 0)
	
	#: row index
	if (outflowdir == 1 || outflowdir == 5) {
		o_outputs[1] <- rowindex
	} else if (outflowdir == 2 || outflowdir == 3 || outflowdir == 4) {
		o_outputs[1] <- rowindex + 1
	} else {
		o_outputs[1] <- rowindex - 1
	}
	
	#: col index
	if (outflowdir == 3 || outflowdir == 7) {
		o_outputs[2] <- colindex
	} else if (outflowdir == 1 || outflowdir == 2 || outflowdir == 8) {
		o_outputs[2] <- colindex + 1
	} else {
		o_outputs[2] <- colindex - 1
	}
	
	return(o_outputs)
}

## [2.10] start the model simulation
## -------------------------------------------------
## Name:	FM.start()
## Inputs:	
##			@animation:			FALSE = no animation png files will be generated
##								TRUE = png files for animationa will be generated
##			@bgtype:			0 = transparent background
##								1 = use DEM raster as the background
##								2 = use user-defined aerial raster as the background (must specify the file name in the next parameter)
##			@aerialraster:		full file path to an aerial image raster for your domain
##								if this is specified, this aerial image will be used as the background of the annimation png files
##			@pdfoutput:			FALSE = no pdf files will be generated
##								TRUE = yes, each time step will come with a pdf file showing surface water depth
## -------------------------------------------------
FM.start <- function(animation = FALSE, bgtype = 0, aerialraster = "", pdfoutput = FALSE)
{
	cat("*********************************************\n")
	cat("*********    FloodMapper: START     *********\n")
	cat("---------------------------------------------\n")
	if (FM.env$b_debug) {
		cat("*********************************************\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		cat("*********    FloodMapper: START     *********\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		cat("---------------------------------------------\n", file = FM.env$s_logfile, sep = "", append = TRUE)
	}
	flush.console()
	
	#: check if the aerial raster exists
	if (animation && bgtype == 2) {
		if (aerialraster == "") {
			if (FM.env$b_debug) cat("ERROR: the path to an aerial raster MUST be specified!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
			return("ERROR: the path to an aerial raster MUST be specified!\n")
		} else if (!file.exists(aerialraster)) {
			if (FM.env$b_debug) cat("ERROR: the aerial raster does NOT exist!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
			return("ERROR: the aerial raster does NOT exist!\n")
		} else {
			#: update the png output settings according this aerial image
			o_tmp_tmp_ras <- raster(aerialraster)
			n_aerial_nrow = dim(o_tmp_tmp_ras)[1]
			n_aerial_ncol = dim(o_tmp_tmp_ras)[2]
			rm(o_tmp_tmp_ras)
			if (n_aerial_nrow == n_aerial_ncol) {
				FM.env$n_png_width <- FM.env$n_png_maxwidth
				FM.env$n_png_height <- FM.env$n_png_maxheight
			} else if (n_aerial_nrow < n_aerial_ncol) {
				FM.env$n_png_width <- FM.env$n_png_maxwidth
				FM.env$n_png_height <- ceiling(FM.env$n_png_maxheight * n_aerial_nrow / n_aerial_ncol)
			} else if (n_aerial_nrow > n_aerial_ncol) {
				FM.env$n_png_width <- ceiling(FM.env$n_png_maxwidth * n_aerial_ncol / n_aerial_nrow)
				FM.env$n_png_height <- FM.env$n_png_maxheight
			}
		}
	}
	
	#: check if start and end datetimes of model run match those of precipitation data
	if (FM.env$o_run_start_datetime < FM.env$o_precip_start_datetime || FM.env$o_run_end_datetime > FM.env$o_precip_end_datetime) {
		if (FM.env$b_debug) cat("ERROR: precipitation data are not fully available for the simulation period!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		return("ERROR: precipitation data are not fully available for the simulation period!\n")
	}
	
	#: store the start time
	time_start <- proc.time()
	
	#: =======================
	for (irun in 1:FM.env$n_run_total_internalruns)
	{
		#: [+] for each internal model run, do the following
		n_tmp_currentdatetime <- FM.env$o_run_start_datetime + (irun - 1) * FM.env$n_run_internal_timestep
		
		#: if it reaches the end datetime, no need to run this
		if (n_tmp_currentdatetime >= FM.env$o_run_end_datetime) next
		
		cat(format(n_tmp_currentdatetime, "%Y-%m-%d %H:%M:%S"))
		if (FM.env$b_debug) cat(format(n_tmp_currentdatetime, "%Y-%m-%d %H:%M:%S"), file = FM.env$s_logfile, sep = "", append = TRUE)
		flush.console()

		#: check if the results of this run need to be saved
		n_tmp_output_yes <- 0
		n_tmp_current_totalseconds <- irun * FM.env$n_run_internal_timestep
		if (n_tmp_current_totalseconds %% FM.env$n_run_output_interval == 0) n_tmp_output_yes <- 1
			
		#: temporary matrixes for each internal run
		o_tmp_Outflow_matrix <- matrix(0, FM.env$n_dem_nrow, FM.env$n_dem_ncol)
		o_tmp_Inflow_matrix <- matrix(0, FM.env$n_dem_nrow, FM.env$n_dem_ncol)
		o_tmp_Watertankheight_matrix <- FM.env$o_waterdepth_matrix
		o_tmp_Instheight_matrix <- matrix(0, FM.env$n_dem_nrow, FM.env$n_dem_ncol)
		o_tmp_PrecipIntensity_matrix <- matrix(0, FM.env$n_dem_nrow, FM.env$n_dem_ncol)

		#: calculate the instantaneous height = elevation + surface water depth
		#: convert to an integer for comparison
		o_tmp_Instheight_matrix <- ceiling((FM.env$o_dem_value_matrix + FM.env$o_waterdepth_matrix) * FM.env$n_integer_factor)
		
		#: store the start time
		time_start_A <- proc.time()
		
		#: [A] calcualte the surface outflow for each grid cell
		#: NOTE: the inflow to each grid cell (outflows of its neighbors) will be calculated in the meantime
		for (irow in 1:FM.env$n_dem_nrow)
		{
			for (icol in 1:FM.env$n_dem_ncol)
			{
				#: note that there might be some background cells (NA) in the original DEM ==> no need to handle it
				if (is.na(o_tmp_Instheight_matrix[irow, icol])) {
					#: if it is a background cell, skip it.
					next
				}
				
				#: 1) determine outflow direction and instantaneous slope for current grid cell
				o_insthgts <- rep(NA, 9)
				o_insthgts[1] <- o_tmp_Instheight_matrix[irow, icol]
				if (icol + 1 <= FM.env$n_dem_ncol) o_insthgts[2] <- o_tmp_Instheight_matrix[irow, icol + 1] #: 1
				if (irow + 1 <= FM.env$n_dem_nrow && icol + 1 <= FM.env$n_dem_ncol) o_insthgts[3] <- o_tmp_Instheight_matrix[irow + 1, icol + 1] #: 2
				if (irow + 1 <= FM.env$n_dem_nrow) o_insthgts[4] <- o_tmp_Instheight_matrix[irow + 1, icol] #: 3
				if (irow + 1 <= FM.env$n_dem_nrow && icol - 1 >= 1) o_insthgts[5] <- o_tmp_Instheight_matrix[irow + 1, icol - 1] #: 4
				if (icol - 1 >= 1) o_insthgts[6] <- o_tmp_Instheight_matrix[irow, icol - 1] #: 5
				if (irow - 1 >= 1 && icol - 1 >= 1) o_insthgts[7] <- o_tmp_Instheight_matrix[irow - 1, icol - 1] #: 6
				if (irow - 1 >= 1) o_insthgts[8] <- o_tmp_Instheight_matrix[irow - 1, icol] #: 7
				if (irow - 1 >= 1 && icol + 1 <= FM.env$n_dem_ncol) o_insthgts[9] <- o_tmp_Instheight_matrix[irow - 1, icol + 1] #: 8
				
				o_tmps <- FM.calOutflowDir(o_insthgts)
				n_outflowdirection <- o_tmps[1]
				n_outflowdir_flag <- o_tmps[2]
				n_edge_corner_flag <- o_tmps[3]
				rm(o_tmps)
				
				#: if the outflow direction is 0, water will stay within this cell, no need to handle it
				if (n_outflowdirection == 0) {
					#: set the water tank height as zero
					o_tmp_Watertankheight_matrix[irow, icol] <- 0
					next
				}
				
				#: 2) determine the downstream grid cell's location according to the outflow direction
				o_downstream_indices <- FM.convertOutflowDirtoRowCol(irow, icol, n_outflowdirection)
				o_maxheightcell_indices <- FM.convertOutflowDirtoRowCol(irow, icol, n_edge_corner_flag) #: needed for edge and corner cells.

				#: 3) calculate the outflow velocity
				n_h1 <- 0 #: upstream
				n_l1 <- 0 #: upstream
				n_h2 <- 0 #: downstream
				n_l2 <- 0 #: downstream
				if (n_edge_corner_flag > 0) {
					#: this is an edge or corner cell
					#: upstream = max-height cell, downstream = current grid cell
					n_h1 <- FM.env$o_waterdepth_matrix[o_maxheightcell_indices[1], o_maxheightcell_indices[2]]
					n_l1 <- FM.env$o_dem_value_matrix[o_maxheightcell_indices[1], o_maxheightcell_indices[2]]
					n_h2 <- FM.env$o_waterdepth_matrix[irow, icol]
					n_l2 <- FM.env$o_dem_value_matrix[irow, icol]
				} else {
					#: this is a normal cell
					#: upstream = current grid cell, downstream = the lowest one
					n_h1 <- FM.env$o_waterdepth_matrix[irow, icol]
					n_l1 <- FM.env$o_dem_value_matrix[irow, icol]
					n_h2 <- FM.env$o_waterdepth_matrix[o_downstream_indices[1], o_downstream_indices[2]]
					n_l2 <- FM.env$o_dem_value_matrix[o_downstream_indices[1], o_downstream_indices[2]]
				}
				
				#: calculate dx factor ==> sqrt(2)
				n_dx_factor <- 1.0
				if (n_outflowdir_flag == 2) {
					#: diagonal flow to 2, 4, 6, 8
					n_dx_factor <- sqrt(2)
				}
				
				#**********************
				#: NOTE that the water can only flow into its immediate neighbors, cannot move into any of its non-immediate neighbors
				#: thus a threshold for the velocity should be specified according to the maximum allowable moving distance: dt * u = dx
				#**********************
				n_dx_max <- n_dx_factor * FM.env$n_dem_cellsize
				
				#: calculate the internal part in Eq. (9) or (10)
				n_inner_part <- abs((n_l1 - n_l2) / sqrt((n_l1 - n_l2) ^ 2 + n_dx_max ^ 2) - (n_h2 - n_h1) / n_dx_max)
				
				#: calculate outflow velocity
				#: convert h1 and h2 to integer
				n_h1 <- ceiling(n_h1 * FM.env$n_integer_factor)
				n_h2 <- ceiling(n_h2 * FM.env$n_integer_factor)
				n_velocity <- sqrt(FM.env$n_g * n_inner_part * (max(n_h1, n_h2) / FM.env$n_integer_factor) / FM.env$o_land_Cd_matrix[irow, icol])
				
				#: if the outflow velocity is 0, no outflow
				if (ceiling(n_velocity * FM.env$n_integer_factor) <= 0) next

				#: check if the water will move into the non-immediate neighbors
				#: UPDATE: only allow 1/2 of the water tank to be moved at a maximum speed to avoid creating a inverse water depth gradient!
				#: This is because the spirit of water flow is to find a stable gradient in water depth
				n_real_dx <- n_velocity * FM.env$n_run_internal_timestep
				if (ceiling(n_real_dx * FM.env$n_integer_factor) > ceiling(0.5 * n_dx_max * FM.env$n_integer_factor)) {
					#: limit the outflow only to the immediate neighbors by adjusting the velocity
					#: thus this requires that the integral time step should be set as small as possible
					#: but meanwhile, a small integral time step will lead to high computation requirements
					n_real_dx <- 0.5 * n_dx_max
					n_velocity <- n_real_dx / FM.env$n_run_internal_timestep
				}

				#: 4) calculate the volume of the flowing-out water tank for the current grid
				n_height_watertank <- 0
				n_vol_watertank <- 0
				if (n_edge_corner_flag > 0) {
					#: if this is an edge or corner cell, all the surface water in the cell will be included in the water tank
					n_height_watertank <- FM.env$o_waterdepth_matrix[irow, icol]
				} else {
					#: if this is a normal cell
					#: upstream = the current grid cell, downstream = the lowest one
					n_height_watertank <- (o_tmp_Instheight_matrix[irow, icol] - o_tmp_Instheight_matrix[o_downstream_indices[1], o_downstream_indices[2]]) * 1.0 / FM.env$n_integer_factor
				}
				#: the height of the water tank SHOULD NOT be higher than its water depth ==> IMPORTANT!!!
				if (ceiling(n_height_watertank * FM.env$n_integer_factor) > ceiling(FM.env$o_waterdepth_matrix[irow, icol] * FM.env$n_integer_factor)) n_height_watertank <- FM.env$o_waterdepth_matrix[irow, icol]
				#: convert it to integer
				n_height_watertank <- ceiling(n_height_watertank * FM.env$n_integer_factor)
				if (n_height_watertank < 0) n_height_watertank <- 0
				
				#: save the water tank height for vertical mass balance calculation later on
				o_tmp_Watertankheight_matrix[irow, icol] <- n_height_watertank * 1.0 / FM.env$n_integer_factor

				#: calculate the volume of the water tank
				#: convert this to integer
				n_vol_watertank <- ceiling(n_height_watertank * 1.0 * FM.env$n_dem_cellsize ^ 2)
				
				#: 5) calculate the volume of outflow to each immediate neighbor
				n_vol_to_downstream <- 0
				n_vol_to_each_side <- 0
				if (n_outflowdir_flag == 1) {
					#: parallel outflow
					n_vol_to_downstream <- ceiling(FM.env$n_dem_cellsize * n_real_dx * n_height_watertank * 1.0)
					n_vol_to_each_side <- 0
				} else if (n_outflowdir_flag == 2) {
					#: diagonal outflow
					n_vol_to_downstream <- ceiling((n_real_dx / sqrt(2)) ^ 2 * n_height_watertank * 1.0)
					n_vol_to_each_side <- ceiling((n_real_dx / sqrt(2)) * (FM.env$n_dem_cellsize - n_real_dx / sqrt(2)) * n_height_watertank * 1.0)
				} else {
					#: no outflow ==> need to calculate outflow for this grid cell
					n_vol_to_downstream <- 0
					n_vol_to_each_side <- 0
					next
				}

				#: 6) save outflow and infow information
				#: NOTE that the outflow matrix below is used to store the ***remaining*** water volume in the current grid cell
				o_tmp_Outflow_matrix[irow, icol] <- n_vol_watertank - (n_vol_to_downstream + 2 * n_vol_to_each_side)
				if (o_tmp_Outflow_matrix[irow, icol] < 0) o_tmp_Outflow_matrix[irow, icol] <- 0
				
				#: update inflow to its neighbors
				if (n_outflowdir_flag == 1) {
					#: parallel outflow ==> only flow into its downstream neighbor
					n_dscell_row_index <- o_downstream_indices[1]
					n_dscell_col_index <- o_downstream_indices[2]
					if (n_dscell_row_index < 1 || n_dscell_row_index > FM.env$n_dem_nrow) n_dscell_row_index <- 0
					if (n_dscell_col_index < 1 || n_dscell_col_index > FM.env$n_dem_ncol) n_dscell_col_index <- 0
					if (n_dscell_row_index != 0 && n_dscell_col_index != 0) {
						#: flow into the downstream cell (without flow to two sides)
						o_tmp_Inflow_matrix[n_dscell_row_index, n_dscell_col_index] <- o_tmp_Inflow_matrix[n_dscell_row_index, n_dscell_col_index] + n_vol_to_downstream
					}
				} else if (n_outflowdir_flag == 2) {
					#: diagonal outflow ==> update both the downstream and the side neighbors
					#: determine two side neighbors
					n_cell1_row_index <- 0
					n_cell1_col_index <- 0
					n_cell2_row_index <- 0
					n_cell2_col_index <- 0
					if (n_outflowdirection == 2) {
						n_cell1_row_index <- irow
						n_cell1_col_index <- icol + 1
						n_cell2_row_index <- irow + 1
						n_cell2_col_index <- icol
					} else if (n_outflowdirection == 4) {
						n_cell1_row_index <- irow + 1
						n_cell1_col_index <- icol
						n_cell2_row_index <- irow
						n_cell2_col_index <- icol - 1
					} else if (n_outflowdirection == 6) {
						n_cell1_row_index <- irow
						n_cell1_col_index <- icol - 1
						n_cell2_row_index <- irow - 1
						n_cell2_col_index <- icol
					} else {
						n_cell1_row_index <- irow - 1
						n_cell1_col_index <- icol
						n_cell2_row_index <- irow
						n_cell2_col_index <- icol + 1
					}
					if (n_cell1_row_index < 1 || n_cell1_row_index > FM.env$n_dem_nrow) n_cell1_row_index <- 0
					if (n_cell1_col_index < 1 || n_cell1_col_index > FM.env$n_dem_ncol) n_cell1_col_index <- 0
					if (n_cell2_row_index < 1 || n_cell2_row_index > FM.env$n_dem_nrow) n_cell2_row_index <- 0
					if (n_cell2_col_index < 1 || n_cell2_col_index > FM.env$n_dem_ncol) n_cell2_col_index <- 0
					
					#: recalculate the inflow volume to each cell
					o_new_vol_to_downstream <- n_vol_to_downstream

					#: side cell 1
					if (n_cell1_row_index != 0 && n_cell1_col_index != 0) {
						#: get the instantaneous height of this side cell
						#: NOTE you need to convert this to an integer for comparison
						n_cell1_insthgt_diff <- o_tmp_Instheight_matrix[n_cell1_row_index, n_cell1_col_index] - o_tmp_Instheight_matrix[irow, icol]
						if (!is.na(n_cell1_insthgt_diff)) {
							#: if this side cell is NA, it is also treated as an out-of-domain cell
							#: thus, the water will flow out of the domain
							if (n_cell1_insthgt_diff < 0) {
								#: flow into this side cell only if it is lower than the current grid cell
								o_tmp_Inflow_matrix[n_cell1_row_index, n_cell1_col_index] <- o_tmp_Inflow_matrix[n_cell1_row_index, n_cell1_col_index] + n_vol_to_each_side
							} else {
								#: flow into the downstream cell
								o_new_vol_to_downstream <- o_new_vol_to_downstream + n_vol_to_each_side
							}
						}
					}
						
					#: side cell 2
					if (n_cell2_row_index != 0 && n_cell2_col_index != 0) {
						#: get the instantaneous height of this side cell
						n_cell2_insthgt_diff <- o_tmp_Instheight_matrix[n_cell2_row_index, n_cell2_col_index] - o_tmp_Instheight_matrix[irow, icol]
						if (!is.na(n_cell2_insthgt_diff)) {
							#: if this side cell is NA, it is also treated as an out-of-domain cell
							#: thus, the water will flow out of the domain
							if (n_cell2_insthgt_diff < 0) {
								#: flow into this side cell only if it is lower than the current grid cell
								o_tmp_Inflow_matrix[n_cell2_row_index, n_cell2_col_index] <- o_tmp_Inflow_matrix[n_cell2_row_index, n_cell2_col_index] + n_vol_to_each_side
							} else {
								#: flow into the downstream cell
								o_new_vol_to_downstream <- o_new_vol_to_downstream + n_vol_to_each_side
							}
						}
					}
					
					#: update the inflow of the downstream cell
					n_dscell_row_index <- o_downstream_indices[1]
					n_dscell_col_index <- o_downstream_indices[2]
					if (n_dscell_row_index < 1 || n_dscell_row_index > FM.env$n_dem_nrow) n_dscell_row_index <- 0
					if (n_dscell_col_index < 1 || n_dscell_col_index > FM.env$n_dem_ncol) n_dscell_col_index <- 0
					if (n_dscell_row_index != 0 && n_dscell_col_index != 0) {
						#: flow into the downstream cell (without flow to two sides)
						o_tmp_Inflow_matrix[n_dscell_row_index, n_dscell_col_index] <- o_tmp_Inflow_matrix[n_dscell_row_index, n_dscell_col_index] + o_new_vol_to_downstream
					}
				}
			}
		}
		
		#: convert outflow and inflow back to float numbers
		o_tmp_Outflow_matrix <- o_tmp_Outflow_matrix * 1.0 / FM.env$n_integer_factor
		o_tmp_Inflow_matrix <- o_tmp_Inflow_matrix * 1.0 / FM.env$n_integer_factor
		
		# : calculate the total time used
		time_end_A <- (proc.time() - time_start_A)[[3]]
		#cat("\n[A] Time Used:\t", time_end_A, "\n")
		if (FM.env$b_debug) cat(paste("\n[A] Time Used:\t", time_end_A, " s\n", sep = ""), file = FM.env$s_logfile, sep = "", append = TRUE)
		flush.console()
		
		#: ------------------------------------------
		#: [B] calculate the precipitation intensity (m/s) for this time step

		#: store the start time
		time_start_B <- proc.time()
		
		#: dertermine the precip raster file
		n_tmp_current_pdata_records <- floor(n_tmp_current_totalseconds * 1.0 / FM.env$n_precip_record_interval)
		#: NOTE the model start datetime should match on of the precip datetimes
		if (n_tmp_current_pdata_records > 0) n_tmp_current_pdata_records <- n_tmp_current_pdata_records - 1
		n_tmp_current_pdata_datetime <- FM.env$o_run_start_datetime + n_tmp_current_pdata_records * FM.env$n_precip_record_interval
		s_tmp_current_pdata_filename <- format(n_tmp_current_pdata_datetime, "%Y-%m-%d %H:%M:%S")
		s_tmp_current_pdata_filename <- gsub(" ", "_", s_tmp_current_pdata_filename)
		s_tmp_current_pdata_filename <- gsub(":", "-", s_tmp_current_pdata_filename)
		s_tmp_current_pdata_filename <- paste(FM.env$s_run_folder_path, "/input/p_", s_tmp_current_pdata_filename, ".tif", sep = "")
		#: check if the raster file exists
		if (!file.exists(s_tmp_current_pdata_filename)) {
			if (FM.env$b_debug) cat(paste("\nERROR: NO precipitation raster file exists!\nFile path: ", s_tmp_current_pdata_filename, "\n", sep = ""), file = FM.env$s_logfile, sep = "", append = TRUE)
			return(paste("\nERROR: NO precipitation raster file exists!\nFile path: ", s_tmp_current_pdata_filename, "\n", sep = ""))
		}
		
		#: read precip data (unit: mm)
		o_tmp_current_pdata_raster <- raster(s_tmp_current_pdata_filename)
		o_tmp_PrecipIntensity_matrix <- matrix(o_tmp_current_pdata_raster[], FM.env$n_dem_nrow, FM.env$n_dem_ncol, byrow = TRUE)
		#: convert precip depth (mm) to intensity (m/s) by assuming a constant intensity for the entire period
		o_tmp_PrecipIntensity_matrix <- o_tmp_PrecipIntensity_matrix * 0.001 / FM.env$n_precip_record_interval
		
		rm(o_tmp_current_pdata_raster)
		rm(n_tmp_current_pdata_records)
		rm(n_tmp_current_pdata_datetime)
		rm(s_tmp_current_pdata_filename)

		#: dertermine the inflow raster file
		n_tmp_current_indata_records <- floor(n_tmp_current_totalseconds * 1.0 / FM.env$n_inflow_record_interval)
		#: NOTE the model start datetime should match on of the precip datetimes
		if (n_tmp_current_indata_records > 0) n_tmp_current_indata_records <- n_tmp_current_indata_records - 1
		n_tmp_current_indata_datetime <- FM.env$o_run_start_datetime + n_tmp_current_indata_records * FM.env$n_inflow_record_interval
		s_tmp_current_indata_filename <- format(n_tmp_current_indata_datetime, "%Y-%m-%d %H:%M:%S")
		s_tmp_current_indata_filename <- gsub(" ", "_", s_tmp_current_indata_filename)
		s_tmp_current_indata_filename <- gsub(":", "-", s_tmp_current_indata_filename)
		s_tmp_current_indata_filename <- paste(FM.env$s_run_folder_path, "/input/in_", s_tmp_current_indata_filename, ".tif", sep = "")
		#: update the precipitation intensity only if the raster file exists
		if (file.exists(s_tmp_current_indata_filename)) {
			#: read inflow data (unit: m3/s)
			o_tmp_current_indata_raster <- raster(s_tmp_current_indata_filename)
			#: convert inflow rate (m3/s) to intensity (m/s) by dividing the area of the grid cell
			o_tmp_inflow_raster <- matrix(o_tmp_current_indata_raster[], FM.env$n_dem_nrow, FM.env$n_dem_ncol, byrow = TRUE)
			#: update precipitation intensity
			o_tmp_PrecipIntensity_matrix <- o_tmp_PrecipIntensity_matrix + o_tmp_inflow_raster / FM.env$n_dem_cellsize ^ 2
			rm(o_tmp_current_indata_raster)
			rm(o_tmp_inflow_raster)
		}
		rm(n_tmp_current_indata_records)
		rm(n_tmp_current_indata_datetime)
		rm(s_tmp_current_indata_filename)
		
		# : calculate the total time used
		time_end_B <- (proc.time() - time_start_B)[[3]]
		#cat("[B] Time Used:\t", time_end_B, "\n")
		if (FM.env$b_debug) cat(paste("[B] Time Used:\t", time_end_B, " s\n", sep = ""), file = FM.env$s_logfile, sep = "", append = TRUE)
		flush.console()
		
		#: ------------------------------------------
		#: [C] calculate vertical water mass balance
		
		#: store the start time
		time_start_C <- proc.time()
	
		#: for each grid cell
		for (irow in 1:FM.env$n_dem_nrow)
		{
			for (icol in 1:FM.env$n_dem_ncol)
			{
				#: note that there might be some background cells (NA) in the original DEM ==> no need to handle it
				if (is.na(o_tmp_Instheight_matrix[irow, icol]))
				{
					#: if it is a background cell, skip it.
					next
				}
				
				#: 1) assume no infiltration, calculate the expected surface water
				n_expected_surfwater <- (FM.env$o_waterdepth_matrix[irow, icol] - o_tmp_Watertankheight_matrix[irow, icol]) * FM.env$n_dem_cellsize ^ 2 + 
						o_tmp_PrecipIntensity_matrix[irow, icol] * FM.env$n_run_internal_timestep * FM.env$n_dem_cellsize ^ 2  + o_tmp_Outflow_matrix[irow, icol] + o_tmp_Inflow_matrix[irow, icol]

				#: update
				if (floor(FM.env$o_soil_moistcontent_matrix[irow, icol] * 100) >= 100) {
					#: if it is already saturated, all expected surface water will become surface water
					#: NOTE that, during two extreme precipitation events (no precip), the soil moisture content can decrease to be less than 100%!
					#:            this situation is not considered in the current version
					FM.env$o_waterdepth_matrix[irow, icol] <- n_expected_surfwater / FM.env$n_dem_cellsize ^ 2
				} else {
					#: 2) infiltration rate: k
					n_infiltration_rate <- exp(-1 * FM.env$o_soil_beta_matrix[irow, icol] * FM.env$o_soil_moistcontent_matrix[irow, icol])

					#: 3) adjusted infiltration rate: k'
					n_adjusted_infiltration_rate <- n_infiltration_rate * FM.env$o_land_alpha_matrix[irow, icol]
					
					#: 4) expected infiltrated water
					n_expected_infiltrated_water <- n_expected_surfwater * n_adjusted_infiltration_rate

					#: 5) maximum infiltratable water
					o_max_infiltratable_water <- FM.env$o_soil_maxwaterVol_matrix[irow, icol] * (1 - FM.env$o_soil_moistcontent_matrix[irow, icol])

					#: 6) recalculate remaining surface water and update soil moisture content
					o_current_surfacewater <- 0.0
					if (ceiling(n_expected_infiltrated_water * FM.env$n_integer_factor) > ceiling(o_max_infiltratable_water * FM.env$n_integer_factor)) {
						#: only partial of the water can be infiltrated as the soil will become saturated
						o_current_surfacewater <- n_expected_surfwater - o_max_infiltratable_water
						#: update the moisture soil content to 100%
						FM.env$o_soil_moistcontent_matrix[irow, icol] <- 1
					} else {
						#: all will be infiltrated
						o_current_surfacewater <- n_expected_surfwater * (1 - n_adjusted_infiltration_rate)
						#: update the moisture soil content
						FM.env$o_soil_moistcontent_matrix[irow, icol] <- FM.env$o_soil_moistcontent_matrix[irow, icol] + n_expected_infiltrated_water / FM.env$o_soil_maxwaterVol_matrix[irow, icol]
						#: cannot be greater than 100%
						if (floor(FM.env$o_soil_moistcontent_matrix[irow, icol] * 100) > 100) FM.env$o_soil_moistcontent_matrix[irow, icol] <- 1
					}
					if (o_current_surfacewater < 0) o_current_surfacewater <- 0
					#: 7) update surface water depth
					FM.env$o_waterdepth_matrix[irow, icol] <- o_current_surfacewater / FM.env$n_dem_cellsize ^ 2
				}

				#: check if there is a drain inlet within this grid cell, to be dealt with differently
				if (FM.env$n_draininlet_flag == 1 && ceiling(FM.env$o_draininlet_maxI_matrix[irow, icol] * FM.env$n_integer_factor) > 0) {
					#: use another approach to calculate the vertical water mass balance
					#: only if the drain inlet flag is turned on and this grid cell has a drain inlet
					n_expected_draininlet_capacity <- FM.env$o_draininlet_maxI_matrix[irow, icol] * FM.env$n_run_internal_timestep * FM.env$n_dem_cellsize ^ 2
					if (ceiling(n_expected_surfwater * FM.env$n_integer_factor) <= ceiling(n_expected_draininlet_capacity * FM.env$n_integer_factor)) {
						#: the drain inlet is able to accommodate all surface water
						FM.env$o_waterdepth_matrix[irow, icol] <- 0.0
					} else {
						#: the drain inlet is only able to accommodate a part of surface water
						FM.env$o_waterdepth_matrix[irow, icol] <- (n_expected_surfwater - n_expected_draininlet_capacity) / FM.env$n_dem_cellsize ^ 2
					}
				}

			}
		}

		# : calculate the total time used
		time_end_C <- (proc.time() - time_start_C)[[3]]
		#cat("[C] Time Used:\t", time_end_C, "\n")
		if (FM.env$b_debug) cat(paste("[C] Time Used:\t", time_end_C, " s\n", sep = ""), file = FM.env$s_logfile, sep = "", append = TRUE)
		flush.console()
		
		#: ------------------------------------------
		#: [D] save the results to the output folder
		#:     only do this if the current timing of the model run matches the output interval
		
		#: store the start time
		time_start_D <- proc.time()
		
		if (n_tmp_output_yes == 1) {
			#: determine the suffix of filenames
			n_tmp_output_datetime <- FM.env$o_run_start_datetime + n_tmp_current_totalseconds
			s_tmp_output_filename <- format(n_tmp_output_datetime, "%Y-%m-%d %H:%M:%S")
			s_tmp_output_filename <- gsub(" ", "_", s_tmp_output_filename)
			s_tmp_output_filename <- gsub(":", "-", s_tmp_output_filename)
			
			#: save rasters
			#: construct template rasters to store output variables
			o_tmp_output_wd_raster <- raster(FM.env$s_dem_raster_file) #: surface water depth
			o_tmp_output_mc_raster <- raster(FM.env$s_dem_raster_file) #: soil moisture content
			o_tmp_output_wd_raster[] <- as.vector(t(FM.env$o_waterdepth_matrix)) * 1000 #: convert m to mm
			o_tmp_output_mc_raster[] <- as.vector(t(FM.env$o_soil_moistcontent_matrix))
			s_tmp_output_wd_rasterfile <- paste(FM.env$s_run_folder_path, "/output/wd_", s_tmp_output_filename, ".tif", sep = "")
			writeRaster(o_tmp_output_wd_raster, filename = s_tmp_output_wd_rasterfile, format = "GTiff", overwrite = TRUE)
			s_tmp_output_mc_rasterfile <- paste(FM.env$s_run_folder_path, "/output/mc_", s_tmp_output_filename, ".tif", sep = "")
			writeRaster(o_tmp_output_mc_raster, filename = s_tmp_output_mc_rasterfile, format = "GTiff", overwrite = TRUE)

			#: process extremely high values
			o_tmp_tmp_wd_values <- o_tmp_output_wd_raster[]
			o_tmp_tmp_wd_values[which(floor(o_tmp_tmp_wd_values) >= max(FM.env$o_png_wd_breakpoints))] <- max(FM.env$o_png_wd_breakpoints)
			o_tmp_output_wd_raster[] <- o_tmp_tmp_wd_values
			rm(o_tmp_tmp_wd_values)
			
			#: infos of the pdf or png imge
			s_tmp_pdfpng_datetime <- format(n_tmp_output_datetime, "%Y-%m-%d %H:%M:%S")
			s_tmp_pdfpng_varname <- "Surface Water Depth (Unit: mm)"
			
			#: save pdf for each step
			o_tmp_ani_background_raster <- NA
			if (pdfoutput) {
				s_tmp_output_pdffile <- paste(FM.env$s_run_folder_path, "/output/", s_tmp_output_filename, ".pdf", sep = "")
			
				#: generate gradual colors
				o_tmp_wd_png_colors <- FM.env$o_png_wd_palette(length(FM.env$o_png_wd_breakpoints) - 1)
				#o_tmp_wd_png_colors[1] <- NA 

				if (bgtype == 2) {
					#: use the defined aerial raster as the background image
					o_tmp_ani_background_raster <- stack(aerialraster)
					pdf(s_tmp_output_pdffile, width = FM.env$n_png_width, height = FM.env$n_png_height, 
						pointsize = FM.env$n_png_pointsize, bg = FM.env$s_png_png_bg_color)
					par(mar = FM.env$o_png_aerial_mar, mgp = FM.env$o_png_mgp)
					plotRGB(o_tmp_ani_background_raster, r = 1, g = 2, b = 3, stretch = "lin", axes = TRUE)
					mtext(s_tmp_pdfpng_datetime, side = 3, line = 0.2, cex = 1.05, adj = 0)
					mtext(s_tmp_pdfpng_varname, side = 3, line = 0.2, cex = 1.05)
					mtext(FM.env$s_versionname, side = 3, line = 0.2, cex = 1.05, adj = 1)
					u <- par("usr") # The coordinates of the plot area
					rect(u[1], u[3], u[2], u[4], col = FM.env$s_png_ras_bg_color, border = FM.env$s_png_border_color)
					plotRGB(o_tmp_ani_background_raster, r = 1, g = 2, b = 3, stretch = "lin", axes = TRUE, add = TRUE)
					plot(o_tmp_output_wd_raster, axes = FALSE, frame.plot = FALSE, bty = "n", box = FALSE, legend = TRUE,
							legend.width = 1, legend.shrink = 0.65,
							breaks = FM.env$o_png_wd_breakpoints, col = o_tmp_wd_png_colors, add = TRUE)
					dev.off()
				} else if (bgtype == 1) {
					#: use the dem as the background image
					o_tmp_ani_background_raster <- raster(FM.env$s_dem_raster_file)
					pdf(s_tmp_output_pdffile, width = FM.env$n_png_width, height = FM.env$n_png_height, 
						pointsize = FM.env$n_png_pointsize, bg = FM.env$s_png_png_bg_color)
					par(mar = FM.env$o_png_dem_mar, mgp = FM.env$o_png_mgp)
					plot(o_tmp_ani_background_raster, axes = FALSE, frame.plot = FALSE, bty = "n", box = FALSE, 
						legend = FALSE, col = FM.env$o_png_dem_palette(FM.env$n_png_dem_legend_breaks))
					mtext(s_tmp_pdfpng_datetime, side = 3, line = 0.2, cex = 1.05, adj = 0)
					mtext(s_tmp_pdfpng_varname, side = 3, line = 0.2, cex = 1.05)
					mtext(FM.env$s_versionname, side = 3, line = 0.2, cex = 1.05, adj = 1)
					u <- par("usr") # The coordinates of the plot area
					rect(u[1], u[3], u[2], u[4], col = FM.env$s_png_ras_bg_color, border = FM.env$s_png_border_color)
					plot(o_tmp_ani_background_raster, axes = FALSE, frame.plot = FALSE, bty = "n", box = FALSE, 
						legend = FALSE, col = FM.env$o_png_dem_palette(FM.env$n_png_dem_legend_breaks), add = TRUE)
					plot(o_tmp_output_wd_raster, axes = FALSE, frame.plot = FALSE, bty = "n", box = FALSE, legend = TRUE,
							legend.width = 1, legend.shrink = 0.65,
							breaks = FM.env$o_png_wd_breakpoints, col = o_tmp_wd_png_colors, add = TRUE)
					dev.off()
				} else {
					#: transparent background
					o_tmp_ani_background_raster <- raster(FM.env$s_dem_raster_file)
					pdf(s_tmp_output_pdffile, width = FM.env$n_png_width, height = FM.env$n_png_height, 
						pointsize = FM.env$n_png_pointsize, bg = FM.env$s_png_png_bg_color)
					par(mar = FM.env$o_png_dem_mar, mgp = FM.env$o_png_mgp)
					plot(o_tmp_ani_background_raster, axes = FALSE, frame.plot = FALSE, bty = "n", box = FALSE, 
						legend = FALSE, col = NA)
					mtext(s_tmp_pdfpng_datetime, side = 3, line = 0.2, cex = 1.05, adj = 0)
					mtext(s_tmp_pdfpng_varname, side = 3, line = 0.2, cex = 1.05)
					mtext(FM.env$s_versionname, side = 3, line = 0.2, cex = 1.05, adj = 1)
					u <- par("usr") # The coordinates of the plot area
					rect(u[1], u[3], u[2], u[4], col = NA, border = FM.env$s_png_border_color)
					plot(o_tmp_ani_background_raster, axes = FALSE, frame.plot = FALSE, bty = "n", box = FALSE, 
						legend = FALSE, col = NA, add = TRUE)
					plot(o_tmp_output_wd_raster, axes = FALSE, frame.plot = FALSE, bty = "n", box = FALSE, legend = TRUE,
							legend.width = 1, legend.shrink = 0.65,
							breaks = FM.env$o_png_wd_breakpoints, col = o_tmp_wd_png_colors, add = TRUE)
					dev.off()
				}
			}
			
			#: save pngs for animation
			o_tmp_ani_background_raster <- NA
			if (animation) {
				s_tmp_output_pngfile <- paste(FM.env$s_run_folder_path, "/output/", s_tmp_output_filename, ".png", sep = "")
				
				#: generate gradual colors
				o_tmp_wd_png_colors <- FM.env$o_png_wd_palette(length(FM.env$o_png_wd_breakpoints) - 1)
				#o_tmp_wd_png_colors[1] <- NA 
				
				if (bgtype == 2) {
					#: use the defined aerial raster as the background image
					o_tmp_ani_background_raster <- stack(aerialraster)
					png(filename = s_tmp_output_pngfile, units = "in", width = FM.env$n_png_width, height = FM.env$n_png_height, 
						pointsize = FM.env$n_png_pointsize, res = FM.env$n_png_resolution, bg = FM.env$s_png_png_bg_color)
					par(mar = FM.env$o_png_aerial_mar, mgp = FM.env$o_png_mgp)
					plotRGB(o_tmp_ani_background_raster, r = 1, g = 2, b = 3, stretch = "lin", axes = TRUE)
					mtext(s_tmp_pdfpng_datetime, side = 3, line = 0.2, cex = 1.05, adj = 0)
					mtext(s_tmp_pdfpng_varname, side = 3, line = 0.2, cex = 1.05)
					mtext(FM.env$s_versionname, side = 3, line = 0.2, cex = 1.05, adj = 1)
					u <- par("usr") # The coordinates of the plot area
					rect(u[1], u[3], u[2], u[4], col = FM.env$s_png_ras_bg_color, border = FM.env$s_png_border_color)
					plotRGB(o_tmp_ani_background_raster, r = 1, g = 2, b = 3, stretch = "lin", axes = TRUE, add = TRUE)
					plot(o_tmp_output_wd_raster, axes = FALSE, frame.plot = FALSE, bty = "n", box = FALSE, legend = TRUE,
							legend.width = 1, legend.shrink = 0.65,
							breaks = FM.env$o_png_wd_breakpoints, col = o_tmp_wd_png_colors, add = TRUE)
					dev.off()
				} else if (bgtype == 1) {
					#: use the dem as the background image
					o_tmp_ani_background_raster <- raster(FM.env$s_dem_raster_file)
					png(filename = s_tmp_output_pngfile, units = "in", width = FM.env$n_png_width, height = FM.env$n_png_height, 
						pointsize = FM.env$n_png_pointsize, res = FM.env$n_png_resolution, bg = FM.env$s_png_png_bg_color)
					par(mar = FM.env$o_png_dem_mar, mgp = FM.env$o_png_mgp)
					plot(o_tmp_ani_background_raster, axes = FALSE, frame.plot = FALSE, bty = "n", box = FALSE, 
						legend = FALSE, col = FM.env$o_png_dem_palette(FM.env$n_png_dem_legend_breaks))
					mtext(s_tmp_pdfpng_datetime, side = 3, line = 0.2, cex = 1.05, adj = 0)
					mtext(s_tmp_pdfpng_varname, side = 3, line = 0.2, cex = 1.05)
					mtext(FM.env$s_versionname, side = 3, line = 0.2, cex = 1.05, adj = 1)
					u <- par("usr") # The coordinates of the plot area
					rect(u[1], u[3], u[2], u[4], col = FM.env$s_png_ras_bg_color, border = FM.env$s_png_border_color)
					plot(o_tmp_ani_background_raster, axes = FALSE, frame.plot = FALSE, bty = "n", box = FALSE, 
						legend = FALSE, col = FM.env$o_png_dem_palette(FM.env$n_png_dem_legend_breaks), add = TRUE)
					plot(o_tmp_output_wd_raster, axes = FALSE, frame.plot = FALSE, bty = "n", box = FALSE, legend = TRUE,
							legend.width = 1, legend.shrink = 0.65,
							breaks = FM.env$o_png_wd_breakpoints, col = o_tmp_wd_png_colors, add = TRUE)
					dev.off()
				} else {
					#: transparent background
					o_tmp_ani_background_raster <- raster(FM.env$s_dem_raster_file)
					png(filename = s_tmp_output_pngfile, units = "in", width = FM.env$n_png_width, height = FM.env$n_png_height, 
						pointsize = FM.env$n_png_pointsize, res = FM.env$n_png_resolution, bg = FM.env$s_png_png_bg_color)
					par(mar = FM.env$o_png_dem_mar, mgp = FM.env$o_png_mgp)
					plot(o_tmp_ani_background_raster, axes = FALSE, frame.plot = FALSE, bty = "n", box = FALSE, 
						legend = FALSE, col = NA)
					mtext(s_tmp_pdfpng_datetime, side = 3, line = 0.2, cex = 1.05, adj = 0)
					mtext(s_tmp_pdfpng_varname, side = 3, line = 0.2, cex = 1.05)
					mtext(FM.env$s_versionname, side = 3, line = 0.2, cex = 1.05, adj = 1)
					u <- par("usr") # The coordinates of the plot area
					rect(u[1], u[3], u[2], u[4], col = NA, border = FM.env$s_png_border_color)
					plot(o_tmp_ani_background_raster, axes = FALSE, frame.plot = FALSE, bty = "n", box = FALSE, 
						legend = FALSE, col = NA, add = TRUE)
					plot(o_tmp_output_wd_raster, axes = FALSE, frame.plot = FALSE, bty = "n", box = FALSE, legend = TRUE,
							legend.width = 1, legend.shrink = 0.65,
							breaks = FM.env$o_png_wd_breakpoints, col = o_tmp_wd_png_colors, add = TRUE)
					dev.off()
				}
			}
			
			#: release tmp variables
			rm(o_tmp_output_wd_raster)
			rm(o_tmp_output_mc_raster)
			rm(o_tmp_ani_background_raster)
		}
		
		# : calculate the total time used
		time_end_D <- (proc.time() - time_start_D)[[3]]
		#cat("[D] Time Used:\t", time_end_D, "\n")
		if (FM.env$b_debug) cat(paste("[D] Time Used:\t", time_end_D, " s\n", sep = ""), file = FM.env$s_logfile, sep = "", append = TRUE)
		flush.console()
		
		cat("\t>>> OK\n")
		if (FM.env$b_debug) cat("\t>>> OK\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		flush.console()
	}

	# : calculate the total time used
	time_end <- (proc.time() - time_start)[[3]]
	Hours <- time_end %/% (60*60)
	Minutes <- (time_end %% 3600) %/% 60
	Seconds <- time_end %% 60
	time_used <- paste(Hours, " h ", Minutes, " m ", Seconds, " s.", sep = "")
	cat("Time Used:\t", time_used, "\n")
	cat("---------------------------------------------\n")
	cat("*********  FloodMapper: SUCCESSFUL!  ********\n")
	cat("*********************************************\n")
	if (FM.env$b_debug) {
		cat(paste("Time Used:\t", time_used, "\n", sep = ""), file = FM.env$s_logfile, sep = "", append = TRUE)
		cat("---------------------------------------------\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		cat("*********  FloodMapper: SUCCESSFUL!  ********\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		cat("*********************************************\n", file = FM.env$s_logfile, sep = "", append = TRUE)
	}
	flush.console()
	
	return("") #: SUCCESS
}

## [2.11] create an animiation gif file for the outputs
## -------------------------------------------------
## Name:	FM.animate()
## Inputs:	
##			@runname:		(OPTIONAL) the name of a model run, if it is empty, the current model run will be selected
##			@workdir:		(OPTIONAL) the full path to user's work directory, if it is empty, the current work dir will be selected.
## -------------------------------------------------
FM.animate <- function(runname = "", workdir = "")
{
	cat("*********************************************\n")
	cat("**********     Animation: START     *********\n")
	cat("---------------------------------------------\n")
	cat("Processing...")
	if (FM.env$b_debug) {
		cat("*********************************************\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		cat("**********     Animation: START     *********\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		cat("---------------------------------------------\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		cat("Processing...", file = FM.env$s_logfile, sep = "", append = TRUE)
	}
	flush.console()
	
	s_tmp_runname <- runname
	if (runname == "") {
		s_tmp_runname <- FM.env$s_run_name
	}
	
	#: check if the folder exists
	s_tmp_folder <- file.path(FM.env$s_worddir, s_tmp_runname)
	if (workdir != "") s_tmp_folder <- file.path(workdir, s_tmp_runname)
	if (!dir.exists(s_tmp_folder)) {
		if (FM.env$b_debug) cat(paste("ERROR: the model run of '", s_tmp_runname, "' does NOT exist!\n", sep = ""), file = FM.env$s_logfile, sep = "", append = TRUE)
		return(paste("ERROR: the model run of '", s_tmp_runname, "' does NOT exist!\n", sep = ""))
	}
	
	#: get all png files and order them by their creation datetime
	o_tmp_filedetails <- file.info(list.files(path = paste(s_tmp_folder, "/output", sep = ""), pattern = "*.png"))
	o_tmp_filedetails <- o_tmp_filedetails[with(o_tmp_filedetails, order(as.POSIXct(ctime))), ] 
	o_tmp_pngfiles <- rownames(o_tmp_filedetails)
	
	#: create a gif file
	if (length(o_tmp_pngfiles) <= 1) {
		cat("\nWARNING: there are no enough png files to create an animation!\n")
		if (FM.env$b_debug) cat("\nWARNING: there are no enough png files to create an animation!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		flush.console()
	} else {
		o_tmp_pngimgs <- c(image_read(paste(s_tmp_folder, "/output/", o_tmp_pngfiles[1], sep = "")))
		for (ip in 2:length(o_tmp_pngfiles))
		{
			o_tmp_pngimgs <- c(o_tmp_pngimgs, image_read(paste(s_tmp_folder, "/output/", o_tmp_pngfiles[ip], sep = "")))
		}
		o_animation <- image_animate(o_tmp_pngimgs, fps = 1, dispose = "previous")
		s_output_filename <- paste(s_tmp_folder, "/output/", "animation.gif", sep = "")
		image_write(o_animation, s_output_filename)
		rm(o_tmp_pngimgs)
		cat("\t>>> OK\n")
		cat("Animation file: ", s_output_filename, "\n")
		if (FM.env$b_debug) cat("\t>>> OK\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		if (FM.env$b_debug) cat(paste("Animation file: ", s_output_filename, "\n", sep = ""), file = FM.env$s_logfile, sep = "", append = TRUE)
		flush.console()
	}
	
	cat("---------------------------------------------\n")
	cat("*******     Animation: SUCCESSFUL!     ******\n")
	cat("*********************************************\n")
	if (FM.env$b_debug) {
		cat("---------------------------------------------\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		cat("*******     Animation: SUCCESSFUL!     ******\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		cat("*********************************************\n", file = FM.env$s_logfile, sep = "", append = TRUE)
	}
	flush.console()
	
	return("") #: SUCCESS
}

## [2.12] re-run the model with new settings
##		 this allows users to change the raster files and txt files under the [input] folder and then re-run the model.
##		 no need to specify the raster files, the model will read all existing files from the [input] folder.
##		 if users need to update the raster files for land, soil, precip, and drain inlet separately,
##		 users can run the corresponding FM.read***() function separately.
##		 NOTE: FM.readPrecip() and FM.readInflow() MUST be run separately in order to update the settings!
## -------------------------------------------------
## Name:	FM.rerun()
## Inputs:	
##			@@@@@@:			refer to the parameters of FM.init() and FM.start()
## -------------------------------------------------
FM.rerun <- function(runname, startdatetime, enddatetime, outputinterval, internaltimestep = 30, debug = TRUE, wdbreakpoints = NA, workdir = "", animation = FALSE, bgtype = 0, aerialraster = "", pdfoutput = FALSE)
{
	#: initialize the model with the new settings
	s_return <- FM.init(runname = runname, startdatetime = startdatetime, enddatetime = enddatetime, outputinterval = outputinterval, 
				internaltimestep = internaltimestep, debug = debug, rerun = TRUE, wdbreakpoints = wdbreakpoints, workdir = workdir)
	if (s_return != "0") {
		if (FM.env$b_debug) cat("Rerun -> Initialization... FAILED!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		return("Rerun -> Initialization... FAILED!\n")
	}
	
	#: update initial conditions and parameters' settings
	s_return <- FM.readDEM(demfile = "", waterdepthfile = "", rerun = TRUE)
	if (s_return != "0") {
		if (FM.env$b_debug) cat("Rerun -> DEM... FAILED!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		return("Rerun -> DEM... FAILED!\n")
	}
	s_return <- FM.readLandcover(lcfile = "", rerun = TRUE)
	if (s_return != "0") {
		if (FM.env$b_debug) cat("Rerun -> Land cover... FAILED!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		return("Rerun -> Land cover... FAILED!\n")
	}
	s_return <- FM.readSoil(soilfile = "", soildepth = 0, soildepthfile = "", soilmoistcontent = 0, soilmoistcontentfile = "", rerun = TRUE)
	if (s_return != "0") {
		if (FM.env$b_debug) cat("Rerun -> Soil... FAILED!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		return("Rerun -> Soil... FAILED!\n")
	}
	if (FM.env$n_draininlet_flag == 1 && FM.env$s_draininlet_raster_file != "") {
		s_return <- FM.readDraininlets(draininletfile = "", rerun = TRUE)
		if (s_return != "0") {
			if (FM.env$b_debug) cat("Rerun -> Drain inlet... FAILED!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
			return("Rerun -> Drain inlet... FAILED!\n")
		}
	}
	
	#: start model run
	s_return <- FM.start(animation = animation, bgtype = bgtype, aerialraster = aerialraster, pdfoutput = pdfoutput)
	if (s_return != "0") {
		if (FM.env$b_debug) cat("Rerun -> Start... FAILED!\n", file = FM.env$s_logfile, sep = "", append = TRUE)
		return("Rerun -> Start... FAILED!\n")
	}
	
	#: SUCCESS
	return("")
}
