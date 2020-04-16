forecast <- function(mcout,            ## object of class mcmc or mcmc.list
                     nFuture,          ## Future years
                     M,                ## Can be different than M in model
                     harvest=NULL,     ## How many bears to harvest/year
                     stochastic=FALSE, ## Should harvest be stochastic?
                     season=c("postrepro","prerepro"), ## When to harvest
                     maxRecruit=Inf,   ## Maximum per-capita recruitment
                     NbarSims=1,       ## nSims to for each MCMC iteration
                                       ## to compute expected values of N
                     report=Inf) {
    nYears <- nFuture+1
    if(length(harvest)==1L & !stochastic) {
        harvest <- rep(harvest, nYears)
    }
    if(length(harvest)>1L & stochastic) {
        stop("If stochastic=TRUE, provide 'harvest' as a single value, taken to be the expected number of individuals harvested")
    }
    if(length(harvest)==1L & stochastic)
        meanHarvest <- harvest
    if(NbarSims<1 || !is.finite(NbarSims))
        NbarSims <- 1
    season <- match.arg(season[1], c("prerepro","postrepro"))
    mcmat <- as.matrix(mcout)
    nIter <- nrow(mcmat)
    v.names <- colnames(mcmat)
    z.names <- grep("z[", v.names, fixed=TRUE, value=TRUE)
    if(length(z.names)<1)
        stop("mcout must include posterior samples for z")
    z.ind0 <- sapply(strsplit(z.names, "[", fixed=TRUE), "[", 2)
    z.ind <- sapply(strsplit(z.ind0, "]", fixed=TRUE), "[", 1)
    Mout <- max(as.integer(sapply(strsplit(z.ind, ","), "[", 1)))
    nPast <- max(as.integer(sapply(strsplit(z.ind, ","), "[", 2)))
    ## Extract posterior samples
    z.last.names <- paste("z[", 1:Mout, ",", nPast, "]", sep="")
    z.last <- t(mcmat[,z.last.names])
    a.last.names <- paste("a[", 1:Mout, ",", nPast, "]", sep="")
    a.last <- t(mcmat[,a.last.names])
    phi <- matrix(NA_real_, nYears, nIter)
    if("gamma1" %in% v.names) {
        DD <- TRUE
        gamma1 <- mcmat[,"gamma1"]
    } else DD <- FALSE
    if("igamma0" %in% v.names) {
        rmod <- "gaus"
        igamma0 <- mcmat[,"igamma0"]
    } else {
        rmod <- "loglin"
        gamma0 <- mcmat[,"gamma0"]
    }
    gamma <- matrix(NA_real_, nYears, nIter)
    if(!DD) {
        mu.lgamma <- mcmat[,"mu.lgamma"]
        sig.lgamma <- mcmat[,"sig.lgamma"]
    }
    random.phi <- "mu.lphi" %in% v.names
    if(random.phi) {
        mu.lphi <- mcmat[,"mu.lphi"]
        sig.lphi <- mcmat[,"sig.lphi"]
    }
    ## Could save a lot of memory by writing over z and a in each iter
    z <- array(0L, c(M, nYears, nIter))
    a <- array(1L, c(M, nYears, nIter))
    Nbar <- array(NA_integer_, c(nYears, NbarSims, nIter))
    reportit <- report>0
    badM <- rep(FALSE, nIter)
    for(iter in 1:nIter) {
        if(reportit) {
            if(iter %% report == 0)
                cat("  iter", iter, "of", nIter, "\n")
        }
        z[1:Mout,1,iter] <- z.last[,iter]
        a[1:Mout,1,iter] <- a.last[,iter]
        ## The q loop is used to simulate trajectories of abundance
        ## for each posterior draw. Allows for Monte Carlo integration to
        ## compute expected values of quantities such as extinction risk
        ## For all variables other than Nbar, only last value in the q loop
        ## will be stored and returned
        for(q in 1:NbarSims) {
        Nbar[1,q,iter] <- sum(z[,1,iter]) ##N[1,iter]
        if(stochastic) {
            harvest <- rpois(nYears, meanHarvest)
        }
        for(t in 2:nYears) {
            NafterHarvest <- max(Nbar[t-1,q,iter]-harvest[t],0)
            if(season=="prerepro") {
                Nrepro <- NafterHarvest ##max(Nbar[t-1,q,iter]-harvest[t],0)
            } else {
                Nrepro <- Nbar[t-1,q,iter]
            }
            if(DD) {
                if(rmod=="loglin") {
                    ## Standardize N (as is done in the model)
                    gamma[t-1,iter] <- min(exp(gamma0[iter] -
                                               gamma1[iter]*(Nrepro-150)/20),
                                           maxRecruit)
                } else if(rmod=="gaus") {
                    gamma[t-1,iter] <- igamma0[iter]*exp(-gamma1[iter]*Nrepro^2)
                } else
                    stop("huh?")
            } else {
                gamma[t-1,iter] <- exp(rnorm(1, mu.lgamma[iter],
                                             sqrt(sig.lgamma[iter])))
            }
            if(random.phi) {
                phi[t-1,iter] <- plogis(rnorm(1, mu.lphi[iter],
                                              sqrt(sig.lphi[iter])))
            } else {
                phi[t-1,iter] <- mcmat[iter, "phi"]
            }
            ## Even if harvest occurs after reproduction, highly unlikely
            ## that a sow's cubs will survive. Harvest 'season'
            ## therefore has an impact on the per-capita recruitment
            ## rate, but not on the number of females that contribute
            ## to recruitment 
            ER <- NafterHarvest*gamma[t-1,iter]
            Avail <- sum(a[,t-1,iter])
            if(Avail<1) {
                psi <- 1
            } else {
                psi <- ER/Avail
            }
            if(psi >= 1) {
                psi <- 1
                phi[t-1,iter] <- 1
                badM[iter] <- TRUE
                warning("Entrance probability >1. Increase M. For now, N fixed at max.")
            }
            ## Harvest could be spatially explicit
            not.harvested <- rep(1, M)
            if(!is.null(harvest) && (Nbar[t-1,q,iter]>0) && (season=="postrepro")) {
                harvested <- sample(M, min(harvest[t-1], Nbar[t-1,q,iter]),
                                    prob=z[,t-1,iter]/Nbar[t-1,q,iter])
                not.harvested[harvested] <- 0
            } else if(season=="prerepro") {
                warning("This hasn't been implemented yet")
            }
            Ez <- z[,t-1,iter]*phi[t-1,iter]*not.harvested +
                a[,t-1,iter]*psi
            z[,t,iter] <- rbinom(M, 1, Ez)
            a[,t,iter] <- a[,t-1,iter]*(1-z[,t,iter])
            Nbar[t,q,iter] <- sum(z[,t,iter]) ##N[t,iter]
        }
        }
    }

    return(list(z=z, ##N=N,
                Nbar=Nbar,
                gamma=gamma, phi=phi, badM=badM)) 
}
