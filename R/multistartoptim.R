
multiStartoptim <-  function(start0 = NULL, objectivefn,  gradient = NULL, ..., hessian = NULL, 
            lower = -Inf, upper = Inf, control = list(),
             method = c("L-BFGS-B", "Nelder-Mead","nlminb"), nbtrials = NULL, 
             typerunif = c("runifbase", "runifantithetics", "sobol", "torus", "niederreiterlodisp"), localsearch = c("exhaustive", "median"), 
              verb=FALSE, nbclusters = NULL)
{  
  if(missing(objectivefn)) stop("Objective function must be provided")
  
  nbpar <- length(lower)       
  if (nbpar != length(upper) || lower > upper) stop("parameter lower must be lower than parameter upper componentwise, and they must have the same length")
  
  if (missing(method) || !(method %in% c("L-BFGS-B", "Nelder-Mead", "nlminb")))  
  {
    warning("optimization method is missing or inadequate, default is set to nlminb")
    themethod <- "nlminb"
  }
  else  
  {
    themethod <- match.arg(method)
  }
   
  if (missing(nbclusters) || is.null(nbclusters) || nbclusters <= 1)
  {      
      if (missing(start0) && (missing(nbtrials) || nbtrials <= 0)) stop("When nbtrials is missing or equal to 0, starting parameters 'start0' should be provided")
                  
      if (missing(nbtrials) || is.null(nbtrials) || nbtrials <= 0)
      {        
          if (!missing(typerunif)) warning("unused argument typerunif")
          
          if (missing(hessian) && themethod %in% c("L-BFGS-B", "Nelder-Mead")) hessian <- FALSE
          
          
          return(call_localoptim(start0 = start0, objectivefn = objectivefn, gradient = gradient, ..., hessian = hessian, 
                          lower = lower, upper = upper, control = control, themethod = themethod, nbtrials = nbtrials , 
                          verb = verb, cl = NULL))          
      } 
      else
      {          
            if (!missing(start0)) warning("starting value is unused when nbtrials is provided")
            
            if (missing(typerunif) || !(typerunif %in% c("runifbase", "runifantithetics", "sobol", "torus", "niederreiterlodisp")))
            {
              warning("typerunif is either missing or inadequate, and was set to default 'runifbase'")
              typerunif <- "runifbase"
            }
            else 
            {
              typerunif <- match.arg(typerunif)
            }
            
            if (!is.finite(lower) || !is.finite(upper))
              stop("Both the lower and upper bounds for parameters should be provided and finite")
            
            localsearch <- match.arg(localsearch)
            
            if (length(localsearch) == 0) localsearch <- "exhaustive"
                        
            if (localsearch == "exhaustive")
            {
              start0 <- rUnif(nbtrials = nbtrials, nbpar = nbpar, lower = lower, upper = upper, 
                             method=typerunif)
            }
            
            if (localsearch == "median")
            {                                    
                  if (nbpar == 1)
                {  
                  U <-  rUnif(nbtrials = nbtrials, nbpar = nbpar, lower = lower, upper = upper, 
                              method=typerunif)
                  fwU <- objectivefn(U, ...)
                  start0 <- U[fwU < median(fwU)] 
                  nbtrials <- length(start0)
                }
                else
                { 
                  U <-  rUnif(nbtrials = nbtrials, nbpar = nbpar, lower = lower, upper = upper, 
                              method=typerunif)
                  fwU <- apply(U, 1, function(x) objectivefn(x, ...))
                  start0 <- U[fwU < median(fwU), ]
                  nbtrials <- dim(start0)[1]
                }    
            }
              
            if (missing(hessian) && themethod %in% c("L-BFGS-B", "Nelder-Mead")) hessian <- FALSE
            return(call_localoptim(start0 = start0, objectivefn = objectivefn, gradient = gradient, ..., hessian = hessian, 
                                   lower = lower, upper = upper, control = control, themethod = themethod, nbtrials = nbtrials , 
                                   verb = verb, cl = NULL))                 
        } 
  } 
  else
  {
    if (nbclusters < 0 || floor(nbclusters) != nbclusters || is.array(nbclusters) || !is.numeric(nbclusters)) stop("nbclusters must be a positive integer")
        
    if(missing(lower) || missing(upper)) stop("lower and upper are missing and must be provided")
    
    if (!is.finite(lower) || !is.finite(upper)) stop("Both the lower and upper bounds for parameters should be provided and finite")
    
    if(verb) warning("argument verb is unused in parallel computation")
        
    if (missing(nbtrials) || is.null(nbtrials) || nbtrials <= 0) 
    {
      warning("nbtrials is not provided (or negative), default is set to 50")
      nbtrials <- 50
    }
    
    if (missing(typerunif) || !(typerunif %in% c("runifbase", "runifantithetics", "sobol", "torus", "niederreiterlodisp")))
    {
      warning("typerunif is either missing or inadequate, and was set to default 'runifbase'")
      typerunif <- "runifbase"
    } else 
    {
      typerunif <- match.arg(typerunif)
    }
    
    localsearch <- match.arg(localsearch)
    
    if (length(localsearch) == 0) localsearch <- "exhaustive"
    
    if (localsearch == "exhaustive")
    {
      start0 <- rUnif(nbtrials = nbtrials, nbpar = nbpar, lower = lower, upper = upper, 
                      method=typerunif)
    }
    
    if (localsearch == "median")
    {                                    
      if (nbpar == 1)
      {  
        U <-  rUnif(nbtrials = nbtrials, nbpar = nbpar, lower = lower, upper = upper, 
                    method=typerunif)
        fwU <- objectivefn(U, ...)
        start0 <- U[fwU < median(fwU)] 
        nbtrials <- length(start0)
      }
      else
      { 
        U <-  rUnif(nbtrials = nbtrials, nbpar = nbpar, lower = lower, upper = upper, 
                    method=typerunif)
        fwU <- apply(U, 1, function(x) objectivefn(x, ...))
        start0 <- U[fwU < median(fwU), ]
        nbtrials <- dim(start0)[1]
      }    
    }     
        
    if (missing(hessian) && themethod %in% c("L-BFGS-B", "Nelder-Mead")) hessian <- FALSE
    
    packageStartupMessage("Processing...", appendLF = FALSE)  
    cat("\n")    
    return(call_localoptim(start0 = start0, objectivefn = objectivefn, gradient = gradient, ..., hessian = hessian, 
                                  lower = lower, upper = upper, control = control, themethod = themethod, nbtrials = nbtrials , 
                                  verb = verb, cl = nbclusters))
  }
  
}

rUnif <- function(nbtrials, nbpar, lower, upper, 
                  method=c("runifbase", "runifantithetics", "sobol", "torus", "niederreiterlodisp"))
{
        callrunifmat <- function(nbtrials, nbpar, lower, upper)
        { 
          if (nbpar == 1) 
          {
            return(lower + (upper - lower) * randtoolbox::SFMT(nbtrials))
          }
          else 
          {  
            Umin <- matrix(rep.int(lower, nbtrials), nrow = nbtrials, ncol = nbpar, byrow=T)
            etendue <- upper - lower
            return(Umin + matrix(rep.int(etendue, nbtrials), nrow = nbtrials, ncol=nbpar, byrow=T)*randtoolbox::SFMT(nbtrials,nbpar))
          }
        }
        
        callrunifmatantith <- function(nbtrials, nbpar, lower, upper)
        { 
          
          nbtrials_2 <- 0.5*nbtrials
          if (nbpar == 1) 
          {
            Uni <- randtoolbox::SFMT(nbtrials_2)
            return(lower + (upper - lower) * c(Uni, 1-Uni))
          }
          else 
          {  
            Umin <- matrix(rep.int(lower, nbtrials), nrow = nbtrials, ncol = nbpar, byrow=T)
            etendue <- upper - lower
            Uni <- randtoolbox::SFMT(nbtrials_2,nbpar)
            return(Umin + matrix(rep.int(etendue, nbtrials), nrow = nbtrials, ncol = nbpar, byrow = T)*rbind(Uni, 1-Uni))
          }
        }
        
        callxlodisp <- function(nbtrials, lower, upper)
        {
          nbpar <- length(lower)
          
          if (nbpar == 1) 
          { 
            ylodisp <- log(2*seq_len(nbtrials)[-1] - 3)/(log(2)) 
            xlodisp <- c(1, ylodisp - floor(ylodisp))     
            return(lower + (upper - lower) * xlodisp)
          }
          else
          {
            ylodisp <- log(2*seq_len(nbtrials)[-1] - 3)/(log(2)) 
            xlodisp <- c(1, ylodisp - floor(ylodisp))
            Umin <- matrix(rep.int(lower, nbtrials), nrow = nbtrials, ncol = nbpar, byrow=T)
            etendue <- upper - lower
            return(Umin + matrix(rep.int(etendue, nbtrials), nrow = nbtrials, ncol = nbpar, byrow=T)*
                     sapply(seq_len(nbpar), function (x) sample(xlodisp, nbtrials, replace = FALSE)))
          }
        }
        
  if (method == "runifbase")
  {
    return (callrunifmat(nbtrials = nbtrials, nbpar = nbpar, lower = lower, upper = upper))
  }

  if (method == "runifantithetics")
  {
    if(floor(nbtrials) != nbtrials) 
    {
      stop("To achieve this, nbtrials should be an even number")
    }    
      return(callrunifmatantith(nbtrials = nbtrials, nbpar = nbpar, lower = lower, upper = upper))
  }
  
  if (method == "sobol")
  {
    if (nbpar == 1) 
    {
      return(lower + (upper - lower) * randtoolbox::sobol(n = nbtrials, dim = nbpar, scrambling = 2))
    }
    else 
    {  
      Umin <- matrix(rep.int(lower, nbtrials), nrow = nbtrials, ncol=nbpar, byrow=T)
      etendue <- upper - lower
      return(Umin + matrix(rep.int(etendue, nbtrials), nrow = nbtrials, ncol=nbpar, byrow=T)*randtoolbox::sobol(n = nbtrials, dim = nbpar, scrambling = 2))
    }
  }
  
  if (method == "torus")
  {
    if (nbpar == 1) 
    {
      return(lower + (upper - lower) * randtoolbox::torus(n = nbtrials, dim = nbpar))
    }
    else 
    {  
      Umin <- matrix(rep.int(lower, nbtrials), nrow = nbtrials, ncol=nbpar, byrow=T)
      etendue <- upper - lower
      return(Umin + matrix(rep.int(etendue, nbtrials), nrow = nbtrials, ncol=nbpar, byrow=T)*randtoolbox::torus(n = nbtrials, dim = nbpar))
    }
  }
  
  if (method == "niederreiterlodisp")
  {
    return(callxlodisp(nbtrials = nbtrials, lower = lower, upper = upper))
  }
}

 
 call_localoptim <- function(start0, objectivefn,  gradient = NULL, ..., hessian = NULL, 
                 lower, upper, control = list(), themethod, nbtrials = NULL, 
                 verb=FALSE, cl = NULL)
{ 
   nbpar <- length(lower)

   if (missing(cl) || is.null(cl))
   {
         if (missing(nbtrials) || is.null(nbtrials) || nbtrials <= 0)
         {            
               if (themethod == "L-BFGS-B")
               {
                 if (missing(hessian)) hessian <- FALSE
                 return(optim(par = start0, fn = objectivefn, gr = gradient,  ..., 
                              method = themethod, lower = lower, upper = upper, 
                              control = control, hessian = hessian))            
               } 
               
               if (themethod == "Nelder-Mead")
               {
                 if (missing(hessian)) hessian <- FALSE
                 return(optim(par = start0, fn = objectivefn, gr = gradient,  ..., 
                              method = themethod, control = control, hessian = hessian))
               }
               
               if (themethod == "nlminb") 
               { 
                 
                 return(nlminb(start = start0, objective = objectivefn, gradient = gradient, 
                               hessian = hessian, ..., scale = 1, control = control, 
                               lower = lower, upper = upper))                    
               }
          } 
         else 
          {        
                    ObjectifNSScour <- 1000000
                    ObjectifNSSi <- 0
                    res_cour <- 0 
                    nbpar <- length(lower)
                    if (length(upper) != nbpar || lower > upper) stop("lower and upper are incorrect")
                                    
                    if (verb == TRUE)
                    {            
                      verbtrace <- list(par = matrix(NA, nrow=nbtrials, ncol=nbpar), 
                                        objective = rep.int(NA, nbtrials),
                                        iteration = rep.int(NA, nbtrials))
                    }
                    
                    howFar <- txtProgressBar(min=0,max=nbtrials,style=3)
                    
                    for(i in seq_len(nbtrials))
                    {          
                          if(nbpar == 1)
                          {
                            start_0 <- start0[i]
                          } else              
                          {
                            start_0 <- start0[i,]
                          }
                          
                          if (themethod == "L-BFGS-B")
                          { 
                            if (missing(hessian)) hessian <- FALSE
                            res_cour <- suppressWarnings(try(optim(par = start_0, fn = objectivefn, gr = gradient,  ..., 
                                                               method = themethod, lower = lower, upper = upper, 
                                                               control = control, hessian = hessian), silent = TRUE))
                            
                            try(ObjectifNSSi <- res_cour$value, silent = TRUE)
                            
                          } 
                          
                          if (themethod == "Nelder-Mead")
                          { 
                            if (missing(hessian)) hessian <- FALSE
                            res_cour <- suppressWarnings(try(optim(par = start_0, fn = objectivefn, gr = gradient,  ..., 
                                                               method = themethod, 
                                                               control = control, hessian = hessian), silent = TRUE))
                            
                            
                            ObjectifNSSi <- res_cour$value
                            
                          }
                          
                          if (themethod == "nlminb") 
                          { 
                            res_cour <- suppressWarnings(try(nlminb(start = start_0, objective = objectivefn, gradient = gradient, hessian = hessian,
                                                                ..., scale = 1, control = control, lower = lower, upper = upper), silent = TRUE))
                            
                            
                            ObjectifNSSi <- res_cour$objective
                            
                          }
                          
                                if(ObjectifNSSi < ObjectifNSScour) 
                              {
                                
                                res <- res_cour
                                ObjectifNSScour <- ObjectifNSSi
                                if (verb == TRUE)
                                {
                                  verbtrace$par[i, ] <- res_cour$par
                                  verbtrace$objective[i] <- ObjectifNSSi
                                  verbtrace$iteration[i] <- i
                                }
                              }
                             
                          
                          setTxtProgressBar(howFar, value=i)
                    }    
                    close(howFar)
                    
                    if (verb == TRUE)
                    { 
                      booltrace <- !is.na(verbtrace$iteration)
                      if (nbpar == 1)
                      {startingparams <- start0[booltrace]}
                      else startingparams <- start0[booltrace,]
                      foundparams <- verbtrace$par[booltrace, ]
                      objective_val <- verbtrace$objective[booltrace]
                      iteration_no <- verbtrace$iteration[booltrace]
                      return(list(res=res, 
                                  iteration_no = iteration_no,
                                  startingparams_sequence = startingparams,
                                  foundparams_sequence = foundparams, 
                                  objective_val_sequence = objective_val))  
                    } else return(res)        
          }
 } 
  else 
   {
    
    force(objectivefn)     
    if (!missing(gradient)) force(gradient)    
    
    if (length(upper) != nbpar || lower > upper) stop("lower and upper are incorrect")
    
    cl <- makeCluster(cl, type="SOCK")    
    .env <- environment()
    ldots <- list(...)  
    clusterExport(cl, "ldots", envir = .env)
     
     if (themethod == "L-BFGS-B")
     {       
       if (missing(hessian)) hessian <- FALSE
         if (nbpar <= 1)
         {
           objectivemin <- suppressWarnings(parSapply(cl, start0, function (x) try(optim(x, fn = objectivefn, gr = gradient, ..., 
                                                                                     method = themethod,  
                                                                                     lower = lower, upper = upper,
                                                                                     control = control, hessian = hessian), silent=TRUE)))
           on.exit(stopCluster(cl))           
           packageStartupMessage("             Done !")        
           
           return(objectivemin[ , which.min(objectivemin[2,])])
         } 
          else 
        {
           objectivemin <- suppressWarnings(parRapply(cl, start0, function (x) try(optim(x, fn = objectivefn, gr = gradient, ..., 
                                                                                     method = themethod,  
                                                                                     lower = lower, upper = upper,
                                                                                     control = control, hessian = hessian), silent=TRUE)))
           on.exit(stopCluster(cl))
           packageStartupMessage("             Done !")        
           
           return(objectivemin[[unique(which.min(parSapply(cl, objectivemin, '[[', "value")))]])
           
         }
     }
   
   if (themethod == "Nelder-Mead")
   {   
     
     if (missing(hessian)) hessian <- FALSE
         if (nbpar <= 1)
         {
           objectivemin <- suppressWarnings(parSapply(cl, start0, function (x) try(optim(x, fn = objectivefn, gr = gradient, ..., 
                                                                                     method = themethod,
                                                                                     control = control, hessian = hessian), silent=TRUE)))
           on.exit(stopCluster(cl))
           packageStartupMessage("             Done !")        
           
           return(objectivemin[ , which.min(objectivemin[2,])])
         } 
          else 
        {
           objectivemin <- suppressWarnings(parRapply(cl, start0, function (x) try(optim(x, fn = objectivefn, gr = gradient, ..., 
                                                                                     method = themethod,
                                                                                     control = control, hessian = hessian), silent=TRUE)))
           on.exit(stopCluster(cl))
           packageStartupMessage("             Done !")        
           
           return(objectivemin[[unique(which.min(parSapply(cl, objectivemin, '[[', "value")))]])       
         }
   }
   
   if (themethod == "nlminb")
   {      
     if (!missing(hessian)) force(hessian)
     
           if (nbpar <= 1)
           {
             objectivemin <- suppressWarnings(parSapply(cl, start0, function (x) try(nlminb(x, objective = objectivefn,
                                                                                        gradient = gradient, hessian = hessian, ...
                                                                                        , scale = 1, control = control,
                                                                                        lower = lower, upper = upper), silent=TRUE)))
             on.exit(stopCluster(cl))  
             packageStartupMessage("             Done !")        
             
             return(objectivemin[ , which.min(objectivemin[2,])])
           } 
            else
           {
             objectivemin <- suppressWarnings(parRapply(cl, start0, function (x) try(nlminb(x, objective = objectivefn,
                                                                                        gradient = gradient, hessian = hessian, ...
                                                                                        , scale = 1, control = control,
                                                                                        lower = lower, upper = upper), silent=TRUE)))
             on.exit(stopCluster(cl))
             packageStartupMessage("             Done !")        
             
             return(objectivemin[[unique(which.min(parSapply(cl, objectivemin, '[[', "objective")))]])
             
           }    
   }
    
  }
 }
