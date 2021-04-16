sglfitpath <- function(x, y, nlam, flmin, ulam, isd, intr, nf, eps, peps, dfmax, pmax, jd, 
                    pf, gindex, ngroups,  maxit, gamma, nobs, nvars, vnames) {
    #################################################################################
    # data setup
    storage.mode(y) <- "double"
    storage.mode(x) <- "double"
    
    # gamma setup
    if (gamma < 0 || gamma > 1) stop("gamma must be in [0,1]")  
    
	# gindex should contain only the variable index at the end of each group - no overlap is allowed 
    if (any(diff(gindex)>1)) 
        stop("only adjacent group memberships are allowed")
    
    gindex <- which(c(diff(gindex),1)==1) 
    #################################################################################
    # call Fortran
    if (nf == 0) {
        fit <- .Fortran("sglfitF", as.double(gamma), as.integer(ngroups), as.integer(gindex),
                        as.integer(nobs), as.integer(nvars), as.matrix(x), y, pf, dfmax, pmax, nlam, flmin, ulam, 
                        eps, as.double(peps), isd, intr, maxit, nalam = integer(1), b0 = double(nlam), 
                        beta = double(pmax * nlam), ibeta = integer(pmax), nbeta = integer(nlam), 
                        alam = double(nlam), npass = integer(1), jerr = integer(1), 
                        PACKAGE = "midasml")
    } else {
        fit <- .Fortran("panelsglfitF", as.integer(nf), as.integer(nobs/nf), as.double(gamma), as.integer(ngroups), as.integer(gindex),
                        as.integer(nobs), as.integer(nvars), as.matrix(x), y, pf, dfmax, pmax, nlam, flmin, ulam, 
                        eps, as.double(peps), isd, intr, maxit, nalam = integer(1), b0 = double(nlam), a0 = double(nf * nlam),
                        beta = double(pmax * nlam), ibeta = integer(pmax), nbeta = integer(nlam), 
                        alam = double(nlam), npass = integer(1), jerr = integer(1), 
                        PACKAGE = "midasml")
        
    }
    #################################################################################
    # output
    if (nf == 0){
        nf <- intr
    }
    fit$nf <- nf
    outlist <- getoutput(fit, maxit, pmax, nvars, vnames)
    outlist <- c(outlist, list(npasses = fit$npass, jerr = fit$jerr))
    outlist$dimx <- c(nobs, nvars)
    class(outlist) <- c("sglpath")
    outlist
} 
