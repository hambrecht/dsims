#' Method to run a simulation
#'
#' Runs the simulation and returns the simulation object with results. If
#' running in parallel and max.cores is not specified it will default to using
#' one less than the number of cores / threads on your machine. For example
#' code see \code{\link{make.simulation}}
#'
#' @param simulation \code{\link{Simulation-class}} object
#' @param run.parallel logical option to use multiple processors
#' @param max.cores integer maximum number of cores to use, if not specified then
#' one less than the number available will be used.
#' @param counter logical indicates if you would like to see the progress counter.
#' @param transect.path character gives the pathway to a folder of shapefiles or
#' the path to a single shapefile (.shp file) which give the transects which should
#' be used for the simulations. If a folder of transects a new shapefile will be
#' used for each repetition. If a path specifying a single shapefile then the same
#' transects will be used for each repetition.
#' @param progress.file character path with filename to output progress to file
#' for Distance for Windows progress counter. Not to be used when running directly
#' in R.
#' @return the \code{\link{Simulation-class}} object which now includes
#' the results
#' @export
#' @importFrom parallel detectCores makeCluster clusterEvalQ stopCluster parLapply
#' @importFrom rstudioapi versionInfo
#' @rdname run.simulation-methods
#' @seealso \code{\link{make.simulation}}
run.simulation <- function(simulation, run.parallel = FALSE, max.cores = NA, 
                         counter = TRUE, transect.path = character(0), 
                         progress.file = character(0)) {
  # Input validation
  if (!inherits(simulation, "Simulation")) {
    stop("Input must be a Simulation-class object")
  }
  
  # Initialize variables
  save.data <- load.data <- FALSE
  data.path <- character()
  
  # Process transect path
  if (length(transect.path) > 0) {
    transect.path <- process.transect.path(transect.path, simulation@reps)
  }
  
  # Reset results and warnings
  simulation@results <- create.results.arrays(simulation@reps,
                                            simulation@design@region,
                                            simulation@ds.analysis,
                                            simulation@population.description)
  simulation@warnings <- list(message = list(),
                            counter = list(),
                            index = list())
  
  # Check parallel capabilities
  run.parallel <- check.parallel.capability(run.parallel, max.cores)
  
  # Main execution
  if (run.parallel) {
    # Initialize cluster
    myCluster <- parallel::makeCluster(nCores)
    on.exit(stopCluster(myCluster))
    
    parallel::clusterEvalQ(myCluster, {
      require(dsims)
      gc()
    })
    
    # Run parallel simulation
    simulation <- run.parallel.simulation(simulation, myCluster, counter,
                                        save.data, load.data, data.path,
                                        transect.path, progress.file)
    
  } else {
    # Run serial simulation
    simulation <- run.serial.simulation(simulation, counter,
                                      save.data, load.data, data.path,
                                      transect.path, progress.file)
  }
  
  # Process and display warnings
  if (length(simulation@warnings$message) > 0) {
    display.warnings(simulation@warnings)
  }
  
  return(simulation)
}

# Helper functions
process.transect.path <- function(transect.path, reps) {
  # Remove trailing slash if present
  if (substr(transect.path, nchar(transect.path), nchar(transect.path)) == "/") {
    transect.path <- substr(transect.path, 0, nchar(transect.path)-1)
  }
  
  # Check if single file or directory
  check.path <- unlist(strsplit(transect.path, split = "[.]"))
  if (utils::tail(check.path, 1) == "shp") {
    return(rep(transect.path, reps))
  } else {
    paths <- file.path(transect.path, list.files(path = transect.path, pattern = "shp"))
    if (length(paths) == 0) {
      stop(sprintf("No shapefiles found at %s", transect.path))
    }
    if (length(paths) < reps) {
      warning(sprintf("Insufficient transect files (%d supplied, %d required)", 
                     length(paths), reps))
      mult.factor <- ceiling(reps/length(paths))
      return(rep(paths, mult.factor)[1:reps])
    }
    return(paths)
  }
}

run.parallel.simulation <- function(simulation, myCluster, counter,
                                  save.data, load.data, data.path,
                                  transect.path, progress.file) {
  # Process in batches
  batch.size <- 10
  batches <- split(1:simulation@reps, 
                  ceiling(seq_along(1:simulation@reps)/batch.size))
  
  sim.results <- list()
  sim.warnings <- list()
  
  for (batch_idx in seq_along(batches)) {
    current_batch <- batches[[batch_idx]]
    
    batch_results <- tryCatch({
      if (counter) {
        pbapply::pblapply(
          X = as.list(current_batch),
          FUN = single.sim.loop,
          simulation = simulation,
          save.data = save.data,
          load.data = load.data,
          data.path = data.path,
          transect.path = transect.path,
          save.transects = FALSE,
          progress.file = progress.file,
          cl = myCluster,
          counter = FALSE
        )
      } else {
        parallel::parLapply(
          myCluster,
          X = as.list(current_batch),
          fun = single.sim.loop,
          simulation = simulation,
          save.data = save.data,
          load.data = load.data,
          data.path = data.path,
          counter = FALSE,
          transect.path = transect.path,
          save.transects = FALSE,
          progress.file = progress.file
        )
      }
    }, error = function(e) {
      warning(sprintf("Error in batch %d: %s", batch_idx, e$message))
      return(NULL)
    })
    
    if (!is.null(batch_results)) {
      for (i in seq_along(batch_results)) {
        sim.results[[current_batch[i]]] <- batch_results[[i]]$results
        sim.warnings[[current_batch[i]]] <- batch_results[[i]]$warnings
      }
    }
    
    rm(batch_results)
    parallel::clusterEvalQ(myCluster, gc())
    gc()
  }
  
  simulation <- accumulate.PP.results(simulation = simulation, 
                                    results = sim.results)
  simulation@warnings <- accumulate.warnings(sim.warnings)
  
  rm(sim.results, sim.warnings)
  gc()
  
  return(simulation)
}

run.serial.simulation <- function(simulation, counter,
                                save.data, load.data, data.path,
                                transect.path, progress.file) {
  temp_dir <- file.path(tempdir(), "dsims_temp")
  dir.create(temp_dir, showWarnings = FALSE)
  on.exit(unlink(temp_dir, recursive = TRUE))
  
  if (counter) {
    pb <- txtProgressBar(min = 0, max = simulation@reps, style = 3)
    on.exit(close(pb), add = TRUE)
  }
  
  batch.size <- 40
  n.batches <- ceiling(simulation@reps/batch.size)
  
  accumulated_results <- create.results.arrays(simulation@reps,
                                            simulation@design@region,
                                            simulation@ds.analysis,
                                            simulation@population.description)
  accumulated_warnings <- list(message = list(),
                             counter = list(),
                             index = list())
  
  for (batch in 1:n.batches) {
    process.serial.batch(batch, batch.size, simulation, counter, pb,
                        save.data, load.data, data.path,
                        transect.path, progress.file, temp_dir,
                        accumulated_results, accumulated_warnings)
  }
  
  simulation@results <- accumulated_results
  simulation@warnings <- accumulated_warnings
  
  return(simulation)
}