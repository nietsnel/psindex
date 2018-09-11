#' Parameter Stability Index 
#' @description A function used to calculate fungible parameter estimates.
#' @keywords FPE, fungible parameter estimates
#' @author Jordan L. Prendez, \email{joradanprendez@@gmail.com}
#' @param model user-specified SEM that using \pkg{lavaan} syntax.
#' @param data_set data frame of measured variables.
#' @param RMSEA_pert numeric. Defines the maximum size of the perturbation in data-model fit (in the scale of RMSEA). For example, a RMSEA_pert = .01 indicates that all estimates will be stored that have a RMSEA value within .01 RMSEA of the MLE RMSEA value. 
#' @param genanneal_max_iters integer. The maximum number of iterations of the simulated annealing algorithm. 
#' @param plot_fpe a logical value indicating whether FPEs should be graphed alongside the maximum likelihood estimates. Defaults to \code{FALSE}.
#' @param frac_plot the fraction of FPEs that are graphed. Defaults to 1. 
#' @param iterations_bin numeric. Represents the maximum number of fungible estimates that can be stored. Defaults to 40,000
#' @examples 


#' @export
ps_index <- function(model              =  NULL,
                     data_set           =  NULL,
                     # control_args       =  NULL,
                     group              =  NULL,
                     RMSEA_pert         =  0,
                     meanstructure      =  NULL,
                     genanneal_max_iters    =  100,
                     plot_fpe           =  FALSE,
                     frac_plot          =  1, 
                     iterations_bin     =  40000,
                     control_genSA      =  list(threshold.stop=global.min+tol, verbose=TRUE, temperature=6, 
                          trace.mat = FALSE) 
                                          ){
  first_iteration_indicator <- as.integer(1)
  secondary_optimization_iterations <- 1
  control_genSA <- control_genSA
  
  
  
  library(data.table)
  library(tidyr)
  library(ggplot2)
  library(dplyr)
  # browser()
  # new_counter <- function() {
  #   i <- 0
  #   function() {
  #     i <<- i + 1
  #     i
  #   }
  # }
  # 
  # counter_one <- new_counter()
  
  
  
  # browser()
  
  iters.env <- new.env()

  iters.env$RMSEA_pert2 <- RMSEA_pert


  # assign(x = "iters.env", value = iters.env, envir = parent.frame())
  # assign(x = "iters.env", value = iters.env, envir = parent.frame())
  
  
  # control.args <- control_args
  study <-list(NA)
  iters <- iterations_bin
  reps <- replications
  # assign('iii', value=1, envir=as.environment(iters.env), inherits=TRUE) --- deprecated
  # assign('sd_of_perturbation', value=0, envir=globalenv())##vistigial--used pre GenAnneal -- deprecated
  # assign('rep_fx_mat', value=1, envir=as.environment(iters.env), inherits=TRUE)
  # assign('bandwidth', value=bandwidth, envir=globalenv()) --deprecated

  # fit.mle <<- lavaan:::sem(model=model, data=data_set, verbose = FALSE, debug = FALSE, estimator="ML", group=group)
  fit.mle <<- psindex::sem(model=model, data=data_set, verbose = FALSE, debug = FALSE,  estimator="ML", control=list(optim.method = "NLMINB0"), group=group)
  
  
  variables  <- length(fit.mle@ParTable$est)
  M = matrix(0, nrow = length(fit.mle@optim$x) + 1, ncol = iterations_bin)
  iters_assign <- 1
  dt_speed = as.data.table(M) ###data.table for speed. 
  assign('fpe_wide', value=dt_speed, envir = globalenv())
  # iters_assign <- 1 
  # assign('iters_assign', value=1, envir = globalenv()) ## -- depcrecated
  # browser()
  iters2.env <- new.env()
  iters2.env$iters_assign <- 1
  
  
  
  vars <- (variables+1)
  # answer.array <- matrix(NA, nrow = iters, ncol = vars+1) #bin size #used to be an array
  # fit_mat <- as.data.frame(matrix(NA, nrow=iters, ncol=2)) --deprecated
  secondary_optimization_iterations <- genanneal_max_iters  ###GENANNEAL steps.
  # for(ii in 1:reps){ test
  for(ii in 1:1){

    
    # browser()
    # >> fit1 <- psindex(model = mod1, data = dat, meanstructure = TRUE, control = 
    #                     > +                list(optim.method = "BGFS", eval.max = 40000, iter.max = 200000)) 
    
    # fit <- psindex::sem(model=model, data=data_set, verbose = TRUE, debug = TRUE, partrace = TRUE, estimator="ML", control=list(control.args), group=group)
    fit <- psindex::sem(model=model, data=data_set, verbose = FALSE, debug = FALSE,  estimator="ML", control=list(optim.method = "GENSA"), group=group)
    
    # results<- get('cool.output2', envir=globalenv()) ##skipping this -- now using dt_speed
    results <- get('fpe_wide', envir = globalenv())
    # results[1] <- NULL ##remove first column of dataframe
    results <- as.data.frame(t(results))
    
    
    # browser()
    # results <- as.data.frame(cool.output2)
    fixed_params <- which(fit@ParTable$free == 0)
    # fixed_params <- fixed_params + 1 #shifted 1 to accomodate fit function col.
    estim_params <- which(fit@ParTable$free !=0)
    # estim_params <- estim_params + 1 #shifted 1 to accomodate fit function col.
    fixed2<- length(fixed_params)
    iters_emp <- dim(results)[1]
    fx_vals <- results[,1]
    # results_main <- as.data.frame(matrix(NA, nrow = iters_emp, ncol=(vars-1)))
    
    
    
    results_main <- as.data.frame(matrix(NA, nrow = iters_emp, ncol=(vars)))
    for(i in 1:(vars-length(fixed_params)-1)){ #changed to 2.  # to one...test
      results_main[,(estim_params[i])] <- results[,i+1] #2 extra columns 1. Counter, 2 Function value
    }
    results_main[,dim(results_main)[2]] <- fx_vals
    # browser()
    ##experimental section ###attempting to save initial fit values.
    results_mle <- as.data.frame(fit.mle@ParTable$est)
    fixed_params_mle <- which(fit@ParTable$free == 0)
    # fixed_params <- fixed_params + 1 #shifted 1 to accomodate fit function col.
    estim_params_mle <- which(fit@ParTable$free !=0)
    # estim_params <- estim_params + 1 #shifted 1 to accomodate fit function col.
    fixed2_mle<- length(fixed_params)
    # iters_emp <- dim(results)[1]
    
    results_main_mle <- as.data.frame(matrix(NA, nrow = 1, ncol=(vars-1)))
    results_main_mle[1,] <- results_mle[,1]  #2 extra columns 1. Counter, 2 Function value
    results_main$counter <- as.numeric(iters_emp:1)
    # results_main <- cbind(results[,1],results_main) ##test.
    # results_
    results_main[is.na(results_main)] <- 1
    # results_main <- cbind(results_main,results[,length(results)]) ##adds function value to results
    results_main <- results_main %>%
      arrange(counter)
    
    answer.array <- results_main
    # assign('iii', value=1, envir=as.environment(iters.env), inherits=TRUE)--deprecated
  }
  # assign("fit_mat", value=fit_mat, envir=globalenv()) ##assigns to current environment  --deprecated
  
  # assign("results_array", value=answer.array, envir=globalenv()) ##assigns to current environment
  # assign("results_array_mle", value=results_main_mle, envir=globalenv()) ##assigns to current environment
  
  names1 <- as.data.frame(matrix(NA, nrow=variables, ncol=2))
  names1[,1] <- fit@ParTable$lhs
  names1[,2] <- fit@ParTable$rhs
  names2 <- names1 %>%
    unite(new, V1, V2, sep = "_", remove = TRUE)
  # assign("names", value=names2, envir=globalenv())
  
  # browser()
  ##Prepare for graphing::

  results_array_subset <- answer.array
  data_iterations <- as.data.frame(as.matrix(results_array_subset[,]))

  data_iterations <- data_iterations %>%
    drop_na()
  names2 <- rbind(names2, "discrepancy fx", "count")

  names(data_iterations) <- names2$new
  data_iterations <- subset(data_iterations, `discrepancy fx`!=0) ##removes unused rows.
  data_iterations <- sample_frac(data_iterations, size = frac_plot, replace = FALSE)


  results_main_long <- data_iterations %>%
    # arrange(`discrepancy fx`)
    # mutate(mle = ifelse(`discrepancy fx`))
    gather(key=variable, value=estimate, -count, -`discrepancy fx`, (1:(dim(data_iterations)[2]-2))) %>%
    arrange(variable)

  mle_dat <- results_main_long %>%
    filter(count == 1)
  names(results_main_mle)<- names2$new[1:(dim(results_main_mle)[2])]

  results_array_mle_long <- results_main_mle %>%
    # arrange(`discrepancy fx`)
    # mutate(mle = ifelse(`discrepancy fx`))
    gather(key=variable, value=estimate) %>%
    arrange(variable)

  # fit_meas <- round(as.data.frame(fitMeasures(fit.mle, c("cfi","rmsea","srmr"))),3)

  # assign("fpe_long", value=results_main_long, envir=globalenv()) ##assigns to current environment
  # assign("mle_est_long", value=results_array_mle_long, envir=globalenv()) ##assigns to current environment
  assign("fpe_wide", value=data_iterations, envir=globalenv()) ###Need to address this. -- should not use global variables.
  

  fpe_long <- results_main_long
  mle_est_long <- results_array_mle_long
   
  
  if (plot_fpe==TRUE) {
    
    ggplot() +
      geom_point(data=fpe_long, aes(x=variable, y=estimate), stat="identity") +
      geom_point(data = mle_est_long, aes(x=variable, y=estimate, fill="#FF3300", colour = "#FF3300")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="none") +
      coord_cartesian(ylim = c(-1, 7.5)) 
      # theme_bw()
    
  }
  
}
