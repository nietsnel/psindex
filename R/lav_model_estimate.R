# model estimation
lav_model_estimate <- function(lavmodel       = NULL,
                               lavsamplestats = NULL,
                               lavdata        = NULL,
                               lavoptions     = NULL,
                               lavcache       = list(),
                               do.fit         = TRUE) {

    estimator     <- lavoptions$estimator
    verbose       <- lavoptions$verbose
    debug         <- lavoptions$debug
    ngroups       <- lavsamplestats@ngroups

    if(lavsamplestats@missing.flag || estimator == "PML") {
        group.weight <- FALSE
    } else {
        group.weight <- TRUE
    }

    # temp test
    if(lavoptions$partrace) {
        # fx + parameter values
        PENV <- new.env()
        PENV$PARTRACE <- matrix(NA, nrow=0, ncol=lavmodel@nx.free + 1L)
    }

    # function to be minimized
    minimize.this.function <- function(x, verbose=FALSE, infToMax=FALSE) {
      
        #cat("DEBUG: x = ", x, "\n")

        # current strategy: forcePD is by default FALSE, except
        # if missing patterns are used
        #if(any(lavsamplestats@missing.flag)) {
        #    forcePD <- TRUE
        #} else {
            forcePD <- FALSE
        #}

        # transform variances back
        #x[lavmodel@x.free.var.idx] <- tan(x[lavmodel@x.free.var.idx])

        # update GLIST (change `state') and make a COPY!
        if(lavmodel@eq.constraints) {
            x <- as.numeric(lavmodel@eq.constraints.K %*% x) + 
                            lavmodel@eq.constraints.k0
        }
        GLIST <- lav_model_x2GLIST(lavmodel, x=x)

        fx <- lav_model_objective(lavmodel       = lavmodel, 
                                  GLIST          = GLIST, 
                                  lavsamplestats = lavsamplestats, 
                                  lavdata        = lavdata,
                                  lavcache       = lavcache,
                                  verbose        = verbose,
                                  forcePD        = forcePD)

        # only for PML: divide by N (to speed up convergence)
        if(estimator == "PML") {
            fx <- fx / lavsamplestats@ntotal
        }



        if(debug || verbose) { 
            cat("Objective function  = ", sprintf("%18.16f", fx), "\n", sep="") 
        }
        if(debug) {
            #cat("Current unconstrained parameter values =\n")
            #tmp.x <- lav_model_get_parameters(lavmodel, GLIST=GLIST, type="unco")
            #print(tmp.x); cat("\n")
            cat("Current free parameter values =\n"); print(x); cat("\n")
        }

        if(lavoptions$partrace) {
            PENV$PARTRACE <- rbind(PENV$PARTRACE, c(fx, x))
        }

        # for L-BFGS-B
        #if(infToMax && is.infinite(fx)) fx <- 1e20
        if(!is.finite(fx)) {
            fx <- 1e20
        }

        fx
    }

    first.derivative.param <- function(x, verbose=FALSE, infToMax=FALSE) {

        # transform variances back
        #x[lavmodel@x.free.var.idx] <- tan(x[lavmodel@x.free.var.idx])

        # update GLIST (change `state') and make a COPY!
        if(lavmodel@eq.constraints) {
            x <- as.numeric(lavmodel@eq.constraints.K %*% x) +
                            lavmodel@eq.constraints.k0
        }
        GLIST <- lav_model_x2GLIST(lavmodel, x=x)

        dx <- lav_model_gradient(lavmodel       = lavmodel, 
                                 GLIST          = GLIST, 
                                 lavsamplestats = lavsamplestats,
                                 lavdata        = lavdata,
                                 lavcache       = lavcache,
                                 type           = "free", 
                                 group.weight   = group.weight, ### check me!!
                                 verbose        = verbose,
                                 forcePD        = TRUE)

        if(debug) {
            cat("Gradient function (analytical) =\n"); print(dx); cat("\n")
        }

        #print( dx %*% lavmodel@eq.constraints.K )
        #stop("for now")

        # handle linear equality constraints
        if(lavmodel@eq.constraints) {
            dx <- as.numeric( dx %*% lavmodel@eq.constraints.K )
        }

        # only for PML: divide by N (to speed up convergence)
        if(estimator == "PML") {
            dx <- dx / lavsamplestats@ntotal
        }

        if(debug) {
            cat("Gradient function (analytical, after eq.constraints.K) =\n"); print(dx); cat("\n")
        }

        dx
    } 

    first.derivative.param.numerical <- function(x, verbose=FALSE) {

        # transform variances back
        #x[lavmodel@x.free.var.idx] <- tan(x[lavmodel@x.free.var.idx])

        # numerical approximation using the Richardson method
        npar <- length(x)
        h <- 10e-6
        dx <- numeric( npar )
 
        ## FIXME: call lav_model_objective directly!!
        for(i in 1:npar) {
            x.left <- x.left2 <- x.right <- x.right2 <- x
            x.left[i]  <- x[i] - h; x.left2[i]  <- x[i] - 2*h
            x.right[i] <- x[i] + h; x.right2[i] <- x[i] + 2*h
            fx.left   <- minimize.this.function(x.left)
            fx.left2  <- minimize.this.function(x.left2)
            fx.right  <- minimize.this.function(x.right)
            fx.right2 <- minimize.this.function(x.right2)
            dx[i] <- (fx.left2 - 8*fx.left + 8*fx.right - fx.right2)/(12*h)
        }

        #dx <- lavGradientC(func=minimize.this.function, x=x)
        # does not work if pnorm is involved... (eg PML)

        if(debug) {
            cat("Gradient function (numerical) =\n"); print(dx); cat("\n")
        }        

        dx
    }
 
    # starting values
    start.x <- lav_model_get_parameters(lavmodel)
    if(lavmodel@eq.constraints) {
        start.x <- as.numeric( (start.x - lavmodel@eq.constraints.k0) %*% 
                                lavmodel@eq.constraints.K )
    }

    if(debug) {
        #cat("start.unco = ", lav_model_get_parameters(lavmodel, type="unco"), "\n")
        cat("start.x = ", start.x, "\n")
    }

    # check if the initial values produce a positive definite Sigma
    # to begin with -- but only for estimator="ML"
    #if(estimator %in% c("ML","PML","FML","MML")) {
    if(estimator %in% c("ML","FML","MML")) {
        Sigma.hat <- computeSigmaHat(lavmodel, extra=TRUE, debug=lavoptions$debug)
        for(g in 1:ngroups) {
            if(!attr(Sigma.hat[[g]], "po")) {
                group.txt <- ifelse(ngroups > 1, 
                                    paste(" in group ",g,".",sep=""), ".")
                if(debug) print(Sigma.hat[[g]])
                stop("lavaan ERROR: initial model-implied matrix (Sigma) is not positive definite;\n  check your model and/or starting parameters", group.txt)
                # FIXME: should we stop here?? or try anyway?
                x <- start.x
                fx <- as.numeric(NA)
                attr(fx, "fx.group") <- rep(as.numeric(NA), ngroups)
                attr(x, "converged")  <- FALSE
                attr(x, "iterations") <- 0L
                attr(x, "control")    <- lavoptions@control
                attr(x, "fx")         <- fx
                return(x)
            }
        }
    }

    # scaling factors
    # FIXME: what is the best way to set the scale??
    # current strategy: if startx > 1.0, we rescale by using
    # 1/startx
    SCALE <- rep(1.0, length(start.x))
    #idx <- which(abs(start.x) > 10.0)
    idx <- which(abs(start.x) > 1.0)
    if(length(idx) > 0L) SCALE[idx] <- abs(1.0/start.x[idx])
    #idx <- which(abs(start.x) < 1.0 & start.x != 0.0)
    #if(length(idx) > 0L) SCALE[idx] <- abs(1.0/start.x[idx])
    if(debug) {
        cat("SCALE = ", SCALE, "\n")
    }

    # transforming variances using atan (or another sigmoid function?)
    # FIXME: better approach?
    #start.x[lavmodel@x.free.var.idx] <- atan(start.x[lavmodel@x.free.var.idx])


    # first some nelder mead steps? (default = FALSE)
    INIT_NELDER_MEAD <- lavoptions$optim.init_nelder_mead

    # gradient: analytic, numerical or NULL?
    if(is.character(lavoptions$optim.gradient)) {
        if(lavoptions$optim.gradient %in% c("analytic","analytical")) {
            GRADIENT <- first.derivative.param
        } else if(lavoptions$optim.gradient %in% c("numerical", "numeric")) {
            GRADIENT <- first.derivative.param.numerical
        } else if(lavoptions$optim.gradient %in% c("NULL", "null")) {
            GRADIENT <- NULL
        } else {
            warning("lavaan WARNING: gradient should be analytic, numerical or NULL")
        }
    } else if(is.logical(lavoptions$optim.gradient)) {
        if(lavoptions$optim.gradient) {
            GRADIENT <- first.derivative.param
        } else {
            GRADIENT <- NULL
        }
    } else if(is.null(lavoptions$optim.gradient)) {
        GRADIENT <- first.derivative.param
    }


    # optimizer
    if(length(lavmodel@ceq.nonlinear.idx) == 0L &&
       length(lavmodel@cin.linear.idx)    == 0L &&
       length(lavmodel@cin.nonlinear.idx) == 0L) {
        if(is.null(lavoptions$optim.method)) {
            OPTIMIZER <- "NLMINB"
            #OPTIMIZER <- "BFGS"  # slightly slower, no bounds; better scaling!
            #OPTIMIZER <- "L-BFGS-B"  # trouble with Inf values for fx!
        } else {
            OPTIMIZER <- toupper(lavoptions$optim.method)
            stopifnot(OPTIMIZER %in% c("NLMINB0", "NLMINB1", "NLMINB2", 
                      "NLMINB", "BFGS", "L-BFGS-B", "NONE"))
            if(OPTIMIZER == "NLMINB1") {
                OPTIMIZER <- "NLMINB"
            }
        }
    } else {
        OPTIMIZER <- "NLMINB.CONSTR"
    }

    if(INIT_NELDER_MEAD) {
        if(verbose) cat("Initial Nelder-Mead step:\n")
        trace <- 0L; if(verbose) trace <- 1L
        optim.out <- optim(par=start.x,
                           fn=minimize.this.function,
                           method="Nelder-Mead",
                           #control=list(maxit=10L, 
                           #             parscale=SCALE,
                           #             trace=trace),
                           hessian=FALSE,
                           verbose=verbose)
        cat("\n")
        start.x <- optim.out$par
    }

    

    if(OPTIMIZER == "NLMINB0") {
        if(verbose) cat("Quasi-Newton steps using NLMINB0 (no analytic gradient):\n")
        #if(debug) control$trace <- 1L;
        control.nlminb <- list(eval.max=20000L,
                               iter.max=10000L,
                               trace=0L,
                               #abs.tol=1e-20, ### important!! fx never negative
                               abs.tol=(.Machine$double.eps * 10),
                               rel.tol=1e-10,
                               #step.min=2.2e-14, # in =< 0.5-12
                               step.min=1.0, # 1.0 in < 0.5-21
                               step.max=1.0,
                               x.tol=1.5e-8,
                               xf.tol=2.2e-14)
        control.nlminb <- modifyList(control.nlminb, lavoptions$control)
        control <- control.nlminb[c("eval.max", "iter.max", "trace", 
                                    "step.min", "step.max",
                                    "abs.tol", "rel.tol", "x.tol", "xf.tol")]
        #cat("DEBUG: control = "); print(str(control.nlminb)); cat("\n")
        optim.out <- nlminb(start=start.x,
                            objective=minimize.this.function,
                            gradient=NULL,
                            control=control,
                            scale=SCALE,
                            verbose=verbose) 
        if(verbose) {
            cat("convergence status (0=ok): ", optim.out$convergence, "\n")
            cat("nlminb message says: ", optim.out$message, "\n")
            cat("number of iterations: ", optim.out$iterations, "\n")
            cat("number of function evaluations [objective, gradient]: ",
                optim.out$evaluations, "\n")
        }

        iterations <- optim.out$iterations
        x          <- optim.out$par
        if(optim.out$convergence == 0) {
            converged <- TRUE
        } else {
            converged <- FALSE
        }

    } else if(OPTIMIZER == "NLMINB") {
        if(verbose) cat("Quasi-Newton steps using NLMINB:\n")
        #if(debug) control$trace <- 1L;
        control.nlminb <- list(eval.max=20000L,
                               iter.max=10000L,
                               trace=0L,
                               #abs.tol=1e-20, ### important!! fx never negative
                               abs.tol=(.Machine$double.eps * 10),
                               rel.tol=1e-10,
                               #step.min=2.2e-14, # in =< 0.5-12
                               step.min=1.0, # 1.0 in < 0.5-21
                               step.max=1.0,
                               x.tol=1.5e-8,
                               xf.tol=2.2e-14)
        control.nlminb <- modifyList(control.nlminb, lavoptions$control)
        control <- control.nlminb[c("eval.max", "iter.max", "trace", 
                                    "step.min", "step.max",
                                    "abs.tol", "rel.tol", "x.tol", "xf.tol")]
        #cat("DEBUG: control = "); print(str(control.nlminb)); cat("\n")
        optim.out <- nlminb(start=start.x,
                            objective=minimize.this.function,
                            gradient=GRADIENT,
                            control=control,
                            scale=SCALE,
                            verbose=verbose) 
        if(verbose) {
            cat("convergence status (0=ok): ", optim.out$convergence, "\n")
            cat("nlminb message says: ", optim.out$message, "\n")
            cat("number of iterations: ", optim.out$iterations, "\n")
            cat("number of function evaluations [objective, gradient]: ",
                optim.out$evaluations, "\n")
        }

        iterations <- optim.out$iterations
        x          <- optim.out$par
        if(optim.out$convergence == 0) {
            converged <- TRUE
        } else {
            converged <- FALSE
        }
    } else if(OPTIMIZER == "BFGS") {

        # warning: Bollen example with estimator=GLS does NOT converge!
        # (but WLS works!)
        # - BB.ML works too

        control.bfgs <- list(trace=0L, fnscale=1, 
                             parscale=SCALE, ## or not?
                             ndeps=1e-3,
                             maxit=10000,
                             abstol=1e-20,
                             reltol=1e-10,
                             REPORT=1L)
        control.bfgs <- modifyList(control.bfgs, lavoptions$control)
        control <- control.bfgs[c("trace", "fnscale", "parscale", "ndeps",
                                  "maxit", "abstol", "reltol", "REPORT")]
        #trace <- 0L; if(verbose) trace <- 1L
        optim.out <- optim(par=start.x,
                           fn=minimize.this.function,
                           gr=GRADIENT,
                           method="BFGS",
                           control=control,
                           hessian=FALSE,
                           verbose=verbose)
        if(verbose) {
            cat("convergence status (0=ok): ", optim.out$convergence, "\n")
            cat("optim BFGS message says: ", optim.out$message, "\n")
            #cat("number of iterations: ", optim.out$iterations, "\n")
            cat("number of function evaluations [objective, gradient]: ",
                optim.out$counts, "\n")
        }

        #iterations <- optim.out$iterations
        iterations <- optim.out$counts[1]
        x          <- optim.out$par
        if(optim.out$convergence == 0L) {
            converged <- TRUE
        } else {
            converged <- FALSE
        }
    } else if(OPTIMIZER == "L-BFGS-B") {

        # warning, does not cope with Inf values!!

        control.lbfgsb <- list(trace=0L, fnscale=1,
                               parscale=SCALE, ## or not?
                               ndeps=1e-3,
                               maxit=10000,
                               REPORT=1L,
                               lmm=5L,
                               factr=1e7,
                               pgtol=0)
        control.lbfgsb <- modifyList(control.lbfgsb, lavoptions$control)
        control <- control.lbfgsb[c("trace", "fnscale", "parscale", 
                                    "ndeps", "maxit", "REPORT", "lmm", 
                                    "factr", "pgtol")]
        optim.out <- optim(par=start.x,
                           fn=minimize.this.function,
                           gr=GRADIENT,
                           method="L-BFGS-B",
                           control=control,
                           hessian=FALSE,
                           verbose=verbose,
                           infToMax=TRUE)
        if(verbose) {
            cat("convergence status (0=ok): ", optim.out$convergence, "\n")
            cat("optim L-BFGS-B message says: ", optim.out$message, "\n")
            #cat("number of iterations: ", optim.out$iterations, "\n")
            cat("number of function evaluations [objective, gradient]: ",
                optim.out$counts, "\n")
        }

        #iterations <- optim.out$iterations
        iterations <- optim.out$counts[1]
        x          <- optim.out$par
        if(optim.out$convergence == 0L) {
            converged <- TRUE
        } else {
            converged <- FALSE
        }
    } else if(OPTIMIZER == "NLMINB.CONSTR") {

        ocontrol <- list(verbose=verbose)
        if(!is.null(lavoptions$control$control.outer)) {
            ocontrol <- c(lavoptions$control$control.outer, verbose=verbose)
        }
        control.nlminb <- list(eval.max=20000L,
                               iter.max=10000L,
                               trace=0L,
                               #abs.tol=1e-20,
                               abs.tol=(.Machine$double.eps * 10),
                               rel.tol=1e-9) # 1e-10 seems 'too strict'
        control.nlminb <- modifyList(control.nlminb, lavoptions$control)
        control <- control.nlminb[c("eval.max", "iter.max", "trace",
                                    "abs.tol", "rel.tol")]
        cin <- cin.jac <- ceq <- ceq.jac <- NULL
        if(!is.null(body(lavmodel@cin.function))) cin     <- lavmodel@cin.function
        if(!is.null(body(lavmodel@cin.jacobian))) cin.jac <- lavmodel@cin.jacobian
        if(!is.null(body(lavmodel@ceq.function))) ceq     <- lavmodel@ceq.function
        if(!is.null(body(lavmodel@ceq.jacobian))) ceq.jac <- lavmodel@ceq.jacobian
        trace <- FALSE; if(verbose) trace <- TRUE
        optim.out <- nlminb.constr(start = start.x,
                                   objective=minimize.this.function,
                                   gradient=GRADIENT,
                                   control=control,
                                   scale=SCALE,
                                   verbose=verbose,
                                   cin = cin, cin.jac = cin.jac,
                                   ceq = ceq, ceq.jac = ceq.jac,
                                   control.outer = ocontrol
                                  )
        if(verbose) {
            cat("convergence status (0=ok): ", optim.out$convergence, "\n")
            cat("nlminb.constr message says: ", optim.out$message, "\n")
            cat("number of outer iterations: ", optim.out$outer.iterations, "\n")
            cat("number of inner iterations: ", optim.out$iterations, "\n")
            cat("number of function evaluations [objective, gradient]: ",
                optim.out$evaluations, "\n")
        }

        iterations <- optim.out$iterations
        x          <- optim.out$par
        if(optim.out$convergence == 0) {
            converged <- TRUE
        } else {
            converged <- FALSE
        }
    } else if(OPTIMIZER == "NONE") {
        x <- start.x
        iterations <- 0L
        converged <- TRUE
        control <- list()
        optim.out <- list()
    }

    fx <- minimize.this.function(x) # to get "fx.group" attribute
    # transform back
    if(lavmodel@eq.constraints) {
        x <- as.numeric(lavmodel@eq.constraints.K %*% x) +
                        lavmodel@eq.constraints.k0
    }

    # transform variances back
    #x[lavmodel@x.free.var.idx] <- tan(x[lavmodel@x.free.var.idx])

    attr(x, "converged")  <- converged
    attr(x, "iterations") <- iterations
    attr(x, "control")    <- control
    attr(x, "fx")         <- fx
    if(!is.null(optim.out$con.jac)) attr(x, "con.jac")    <- optim.out$con.jac
    if(!is.null(optim.out$lambda))  attr(x, "con.lambda") <- optim.out$lambda
    if(lavoptions$partrace) {
        attr(x, "partrace") <- PENV$PARTRACE
    }

    x
}

# backwards compatibility
# estimateModel <- lav_model_estimate
