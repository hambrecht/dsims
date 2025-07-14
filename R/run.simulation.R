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
run.simulation <- function(simulation, run.parallel = FALSE, max.cores = NA, counter = TRUE, transect.path = character(0), progress.file = character(0)){
  save.data <- load.data <- FALSE
  data.path <- character()
  # Check if the transect path ends in / if so remove - hard to do the same for windows!
  if(length(transect.path) > 0){
    if(substr(transect.path, nchar(transect.path), nchar(transect.path)) == "/"){
      transect.path <- substr(transect.path, 0, nchar(transect.path)-1)
    }
  }
  transect.path.master <- transect.path
  # Check if it is a single transect set or a folder
  if(length(transect.path) > 0){
    # Check if a folder or file have been specified
    check.path <- unlist(strsplit(transect.path, split = "[.]"))
    index <- length(check.path)
    if(check.path[index] == "shp"){
      # Set transect path to repeat the same value for each rep (for running in parallel)
      transect.path <- rep(transect.path, simulation@reps)
    }else{
      # Make up the vector of shapefile paths
      transect.path <- file.path(transect.path, list.files(path = transect.path, pattern = "shp"))
      if(length(transect.path) == 0){
        stop(paste("No shapefiles found at ", transect.path.master,", cannot run simulation.", sep = ""), call. = FALSE)
      }
      if(length(transect.path) < simulation@reps){
        warning(paste("Insufficient transect shapefiles supplied (", length(transect.path)," supplied, ", simulation@reps, " required). Simulation will use some sets of transects more than once, this may influence the results.", sep = ""), immediate. = TRUE, call. = FALSE)
        mult.factor <- ceiling(simulation@reps/length(transect.path))
        transect.path <- rep(transect.path, mult.factor)[1:simulation@reps]
      }
    }
  }
  #Reset results arrays
  simulation@results <- create.results.arrays(simulation@reps,
                                              simulation@design@region,
                                              simulation@ds.analysis,
                                              simulation@population.description)
  #reset the error/warning message
  simulation@warnings$message <- list()
  simulation@warnings$counter <- list()
  simulation@warnings$index <- list()
  #check the data.path ends in "/"
  if(length(data.path) > 0){
    temp.path <- strsplit(data.path, split = "")
    if(temp.path[length(temp.path)] != "/"){
      # if not add it
      data.path <- paste(data.path, "/", sep = "")
    }
    rm(temp.path)
  }
  # Two if(run.parallel) checks as if libraries not present will change to
  # false before running in parallel and then run in serial
  if(run.parallel){
    if(!requireNamespace('parallel', quietly = TRUE) | !requireNamespace('pbapply', quietly = TRUE)){
      warning("Could not run in parallel, check pbapply library is installed.", immediate. = TRUE, call. = FALSE)
      run.parallel = FALSE
    }else{
      # counts the number of cores you have
      nCores <- getOption("cl.cores", detectCores()) - 1
      if(!is.na(max.cores)){
        nCores <- min(nCores, max.cores)
      }
      if(nCores <= 1){
        warning("Could not run in parallel only one core available/requested (dsims limits running in parallel to 1 less than the number of cores on the machine).", immediate. = TRUE, call. = FALSE)
        run.parallel = FALSE
      }
    }
  }
if(run.parallel){
  # Initialize cluster with memory management
  myCluster <- parallel::makeCluster(nCores)
  parallel::clusterEvalQ(myCluster, {
    require(dsims)
    gc() # Initial garbage collection
  })
  on.exit(stopCluster(myCluster))
  
  # Process in batches to manage memory
  batch.size <- 10  # Adjust based on system memory
  batches <- split(1:simulation@reps, 
                  ceiling(seq_along(1:simulation@reps)/batch.size))
  
  # Initialize accumulated results
  sim.results <- list()
  sim.warnings <- list()
  
  # Process each batch
  for(batch_idx in seq_along(batches)) {
    current_batch <- batches[[batch_idx]]
    
    if(counter){
      batch_results <- pbapply::pblapply(
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
      batch_results <- parLapply(
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
    
    # Extract and accumulate batch results
    for(i in seq_along(batch_results)) {
      sim.results[[current_batch[i]]] <- batch_results[[i]]$results
      sim.warnings[[current_batch[i]]] <- batch_results[[i]]$warnings
    }
    
    # Clean up batch data
    rm(batch_results)
    parallel::clusterEvalQ(myCluster, gc())
    gc()
  }
  
  # Accumulate final results
  simulation <- accumulate.PP.results(simulation = simulation, 
                                    results = sim.results)
  simulation@warnings <- accumulate.warnings(sim.warnings)
  
  # Clean up cluster
  stopCluster(myCluster)
  on.exit()
  
  # Clean up remaining objects
  rm(sim.results, sim.warnings)
  gc()
  if(!run.parallel){
  # Create temp directory for intermediate results
  temp_dir <- file.path(tempdir(), "dsims_temp")
  dir.create(temp_dir, showWarnings = FALSE)
  
  # Initialize progress counter if needed
  if(counter) {
    pb <- txtProgressBar(min = 0, max = simulation@reps, style = 3)
  }
  
  # Process in smaller batches
  batch.size <- 40  # Adjust based on your system
  n.batches <- ceiling(simulation@reps/batch.size)
  
  # Initialize accumulated results
  accumulated_results <- create.results.arrays(simulation@reps,
                                            simulation@design@region,
                                            simulation@ds.analysis,
                                            simulation@population.description)
  accumulated_warnings <- list(message = list(),
                             counter = list(),
                             index = list())
  
  for(batch in 1:n.batches) {
    start.idx <- ((batch-1) * batch.size) + 1
    end.idx <- min(batch * batch.size, simulation@reps)
    
    batch_results <- list()
    for(i in start.idx:end.idx) {
      # Run single simulation
      result <- single.sim.loop(i = i,
                              simulation = simulation,
                              save.data = save.data,
                              load.data = load.data,
                              data.path = data.path,
                              counter = FALSE,  # Handle progress separately
                              transect.path = transect.path,
                              save.transects = FALSE,
                              progress.file = progress.file)
      
      # Save to temp file
      saveRDS(result, file = file.path(temp_dir, paste0("rep_", i, ".rds")))
      
      # Update progress bar if needed
      if(counter) {
        setTxtProgressBar(pb, i)
      }
      
      # Clean up
      rm(result)
      gc()
    }
    
    # Process batch results
    for(i in start.idx:end.idx) {
      result <- readRDS(file.path(temp_dir, paste0("rep_", i, ".rds")))
      
      # Accumulate results
      accumulated_results <- accumulate.PP.results(
        simulation = list(results = accumulated_results),
        results = list(result$results)
      )$results
      
      # Accumulate warnings
      if(length(result$warnings$message) > 0) {
        accumulated_warnings$message <- c(accumulated_warnings$message, result$warnings$message)
        accumulated_warnings$counter <- c(accumulated_warnings$counter, result$warnings$counter)
        accumulated_warnings$index <- c(accumulated_warnings$index, result$warnings$index)
      }
      
      # Clean up
      rm(result)
      unlink(file.path(temp_dir, paste0("rep_", i, ".rds")))
      gc()
    }
  }
  
  # Close progress bar if needed
  if(counter) {
    close(pb)
  }
  
  # Update simulation object
  simulation@results <- accumulated_results
  simulation@warnings <- accumulated_warnings
  
  # Clean up temp directory
  unlink(temp_dir, recursive = TRUE)
  
  # Process warnings
  if(length(simulation@warnings$message) > 0) {
    message("Summary of warnings and errors:")
    for(i in seq(along = simulation@warnings$message)) {
      rep.info <- ifelse(is.null(simulation@warnings$index), "",
                        paste(" in repetition(s): ",
                              paste(simulation@warnings$index[[i]], collapse = ", ")))
      message(paste(simulation@warnings$message[[i]],
                   " (occurred ", simulation@warnings$counter[[i]],
                   " time(s)", rep.info, ")", sep = ""))
    }
    message("-----")
  }
  
  return(simulation)
}

