#' Parameter Stability Index
#' @description A function used to calculate fungible parameter estimates.
#' @keywords FPE, fungible parameter estimates
#' @author Jordan Yee Prendez, \email{jordanyeeprendez@@gmail.com}
#' @param model user-specified SEM that using \pkg{lavaan} syntax.
#' @param data_set data frame of measured variables.
#' @param RMSEA_pert numeric. Defines the maximum size of the perturbation in data-model fit (in the scale of RMSEA). For example, a RMSEA_pert = .01 indicates that all estimates will be stored that have a RMSEA value within .01 RMSEA of the MLE RMSEA value.
#' @param genanneal_max_iters integer. The maximum number of iterations of the simulated annealing algorithm.
#' @param plot_fpe a logical value indicating whether FPEs should be graphed alongside the maximum likelihood estimates. Defaults to \code{FALSE}.
#' @param output_long a logical value that indicates whether the FPEs and the MLE should be output in long data format. Defaults to \code{FALSE}
#' @param frac_plot the fraction of FPEs that are graphed. Defaults to 1.
#' @param iterations_bin numeric. Represents the maximum number of fungible estimates that can be stored. Defaults to 40,000
#' @param control_genSA list. used to control the simulated annealing algorithm (see \link[GenSA]{GenSA} documentation)
#' @param index_method Character. Specify the index used to calculate FPEs. Either RMSEA or AIC. Defaults to RMSEA
#' @param lower Numeric. Lower bound of the parameter search space to be optimized over. Defaults to \code{-10}
#' @param upper Numeric. Upper bound of the parameter search space to be optimized over. Defaults to \code{10}
#' @param starting_val_MLE a logical value indicating whether the MLEs should be used as starting values for the SA algorithm. Defaults to \code{FALSE}
#' @examples 
#' \dontrun{
#' model <- '
#' measurement model
#' ind60 =~ x1 + x2 + x3
#' dem60 =~ y1 + y2 + y3 + y4
#' dem65 =~ y5 + y6 + y7 + y8
#' # regressions
#' dem60 ~ ind60
#' dem65 ~ ind60 + dem60
#' # residual correlations
#' y1 ~~ y5
#' y2 ~~ y4 + y6
#' y3 ~~ y7
#' y4 ~~ y8
#' y6 ~~ y8
#' '
#
#'ps_index(model = model, data_set = PoliticalDemocracy,
#'         RMSEA_pert = .01,
#'         plot_fpe   = TRUE,
#'         index_method = "rmsea",
#'        frac_plot = .3, iterations_bin = 100000,
#'         control_genSA = list(maxit = 500))
#'}




#' @export
ps_index <- function(model               =  NULL,
                     data_set            =  NULL,
                     # control_args      =  NULL,
                     group               =  NULL,
                     RMSEA_pert          =  0,
                     meanstructure       =  NULL,
                     genanneal_max_iters =  100,
                     plot_fpe            =  FALSE,
                     output_long         =  FALSE,
                     frac_plot           =  1,
                     iterations_bin      =  40000,
                     suppress_message    = FALSE,
                     control_genSA       = NULL,
                     index_method        = "rmsea",
                     lower               = -10,
                     upper               = 10,
                     starting_val_MLE    = FALSE
                     # control_genSA      =  list(threshold.stop=global.min+tol, verbose=TRUE, temperature=6,
                     #      trace.mat = FALSE)
                                          ){
    # new_counter2 <- function(){
    #     i <- 5000000
    #     print(i)
    #     function(){
    #         i <- i - 1
    #         if(i<(5000000-1)){
    #             i<-1
    #         }
    #         i <<- i
    #         i
    #
    #     }
    # }
    # Gen_SA_controller <- new_counter2()
    # fpe_sample_satisfied <- Gen_SA_controller()


  # first_iteration_indicator <<- as.integer(1)
  # secondary_optimization_iterations <- 1

  # control_genSA = control_genSA
  # control_genSA <- c(control_genSA, max.call = fpe_sample_satisfied)
  # control_genSA = list(max.call=fpe_sample_satisfied, max.time = 400, maxit=300)
  # control_genSA <<- control_genSA
  # index_method <<- index_method
  # lower <<- lower
  # upper <<- upper
  # starting_val_MLE <<- starting_val_MLE  

  # control_genSA <<- control_genSA
  # control=list(max.call=fpe_sample_satisfied))  ##temp disabled.



  
  
  assign(x = "first_iteration_indicator", value = as.integer(1), envir = cacheEnv)
  assign(x = "secondary_optimization_iterations", value = 1, envir = cacheEnv)
  assign(x = "control_genSA", value = control_genSA, envir = cacheEnv)
  assign(x = "index_method", value = index_method, envir = cacheEnv)
  assign(x = "lower", value = lower, envir = cacheEnv)
  assign(x = "upper", value = upper, envir = cacheEnv)
  assign(x = "starting_val_MLE", value = starting_val_MLE, envir = cacheEnv)
  assign(x = "suppress_message", value = suppress_message, envir = cacheEnv)
  assign(x = "RMSEA_pert", value = RMSEA_pert, envir = cacheEnv)
  
  
  
  


  # suppress_message <<-suppress_message

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
  # RMSEA_pert <<- RMSEA_pert

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
  fit.mle <- psindex::sem(model=model, data=data_set, verbose = FALSE, debug = FALSE,  estimator="ML", control=list(optim.method = "NLMINB0"), group=group)
  assign(x = "fit.mle", value = fit.mle, envir = cacheEnv)
  

  variables  <- length(fit.mle@ParTable$est)
  M = matrix(0, nrow = length(fit.mle@optim$x) + 1, ncol = iterations_bin)
  iters_assign <- 1
  dt_speed = as.data.table(M) ###data.table for speed.
  
  assign('fpe_wide', value=dt_speed, envir = cacheEnv)
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
    results <- get('fpe_wide', envir = cacheEnv)
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

  ##Prepare for graphing::

  results_array_subset <- answer.array
  data_iterations <- as.data.frame(as.matrix(results_array_subset[,]))

  data_iterations <- data_iterations %>%
    drop_na()
  names2 <- rbind(names2, "discrepancy fx", "count")

  names(data_iterations) <- names2$new
  data_iterations_temp <- subset(data_iterations, `discrepancy fx`!=0) ##removes unused rows.
  data_iterations_temp <- data_iterations_temp %>% ###removes any estimate that is equal to another
    distinct(.keep_all = TRUE)
  assign("fpe_wide", value=data_iterations_temp, envir=globalenv()) ###Need to address this. -- should not use global variables.

  data_iterations <- sample_frac(data_iterations_temp, size = frac_plot, replace = FALSE)


  results_main_long <- data_iterations %>%
    gather(key=variable, value=estimate, -count, -`discrepancy fx`, (1:(dim(data_iterations)[2]-2))) %>%
    arrange(variable)
  names(results_main_mle)<- names2$new[1:(dim(results_main_mle)[2])]

  results_array_mle_long <- results_main_mle %>%
    # arrange(`discrepancy fx`)
    # mutate(mle = ifelse(`discrepancy fx`))
    gather(key=variable, value=estimate) %>%
    arrange(variable)

  if (output_long==TRUE) {
  results_main_long1 <- data_iterations_temp %>%
    gather(key=variable, value=estimate, -count, -`discrepancy fx`, (1:(dim(data_iterations)[2]-2))) %>%
    arrange(variable)
  names(results_main_long1)<- names2$new[1:(dim(results_main_long1)[2])]


  assign("fpe_long", value=results_main_long1, envir=globalenv())
  assign("mle_long", value=results_array_mle_long, envir=globalenv())

  rm(data_iterations_temp, results_main_long1)
  }


  # fit_meas <- round(as.data.frame(fitMeasures(fit.mle, c("cfi","rmsea","srmr"))),3)

  # assign("fpe_long", value=results_main_long, envir=globalenv()) ##assigns to current environment
  # assign("mle_est_long", value=results_array_mle_long, envir=globalenv()) ##assigns to current environment

  fpe_long <- results_main_long
  mle_est_long <- results_array_mle_long


  if (plot_fpe==TRUE) {

    ggplot() +
      geom_point(data=fpe_long, aes(x=variable, y=estimate), stat="identity") +
      geom_point(data = mle_est_long, aes(x=variable, y=estimate, fill="#FF3300", colour = "#FF3300")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="none")
    # +
      # coord_cartesian(ylim = c(-1, 7.5))
      # theme_bw()

  }

}



cacheEnv <- new.env()

