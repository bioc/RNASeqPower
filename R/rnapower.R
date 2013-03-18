rnapower <- function(depth, n, n2=n, cv, cv2=cv, effect, alpha, power) {
    if (missing(depth)) stop("Depth must be specified")
    ismiss <- c(missing(n), missing(cv), missing(alpha), missing(power),
                missing(effect))
    if (sum(ismiss) != 1) 
        stop("Exactly one of n, cv, effect, alpha, or power should be unspecified")
    
    if (!missing(n)) {
        n1 <- n
        if (length(n) != length(n2))
            stop("n and n2 must be the same length")
    }
    if (!missing(cv)){
        cv1 <- cv
        if (length(cv) != length(cv2))
            stop("cv and cv2 must be the same length")
    }
    if (!missing(effect))  e2 <- log(effect)  # the scale for effect size
    
    if (missing(n)) {  #solve for the number of samples needed
        out <- array(0., dim=c(length(depth), length(cv1), length(effect),
                               length(alpha), length(power)))
        dimnames(out) <- list(depth, cv1, effect, alpha, power)
        temp <- outer(-qnorm(alpha/2), qnorm(power), "+")
        for (i1 in seq(along=depth)) {
            for (i2 in seq(along=cv1)) {
                for (i3 in seq(along=effect)) {
                    vtemp <-  (1/depth[i1] + cv1[i2]^2) + 
                              (1/depth[i1] + cv2[i2]^2)
                    out[i1, i2,i3,,] <- vtemp * (temp/e2[i3])^2
                }
            }
        }
    }
    
    else if (missing(cv1)) { # solve for the cv (odd choice)
        out <- array(0., dim=c(length(depth), length(n1), length(effect),
                               length(alpha), length(power)))
        dimnames(out) <- list(depth, n1, effect, alpha, power)
        temp <- outer(-qnorm(alpha/2), qnorm(power), "+")
         for (i1 in seq(along=depth)) {
            for (i2 in seq(along=n1)) {
                for (i3 in seq(along=effect)) {
                    out[i1, i2,i3,,] <- (e2[i3]/temp)^2/(1/n1[i2] + 1/n2[i2]) -
                        1/depth[i1]
                }
            }
        }
    }

    else if (missing(effect)) { # solve for the effect size
        out <- array(0., dim=c(length(depth), length(n1), length(cv1),
                               length(alpha), length(power)))
        dimnames(out) <- list(depth, n1, cv1, alpha, power)
        temp <- outer(-qnorm(alpha/2), qnorm(power), "+")
         for (i1 in seq(along=depth)) {
            for (i2 in seq(along=n1)) {
                for (i3 in seq(along=cv1)) {
                    out[i1, i2,i3,,] <- temp * (1/n1[i2] + 1/n2[i2]) * 
                        (2/depth[i1] + cv1[i3]^2 + cv2[i3]^2)
                }
            }
        }
        out <- exp(out)
    }
      
    else if (missing(alpha)) { # solve for the false pos rate
        out <- array(0., dim=c(length(depth), length(n1), length(cv1),
                               length(effect), length(power)))
        dimnames(out) <- list(depth, n1, cv1, effect, power)
        temp1 <- outer(1/depth, cv1^2, '+')
        temp2 <- outer(1/depth, cv2^2, '+')
         for (i1 in seq(along=n1)) {
            for (i2 in seq(along=effect)) {
                for (i3 in seq(along=power)) {
                    za <- abs(e2[i2])/sqrt(temp1/n1[i1] + temp2/n2[i1])
                out[,i1,, i2,i3] <- pnorm(qnorm(power[i3]) - za)*2
                }
            }
        }
    }

    else if (missing(power)) { # solve for the power
        out <- array(0., dim=c(length(depth), length(n1), length(cv1),
                               length(effect), length(alpha)))
        dimnames(out) <- list(depth, n1, cv1, effect, alpha) 
        temp1 <- outer(1/depth, cv1^2, '+')
        temp2 <- outer(1/depth, cv2^2, '+')
         for (i1 in seq(along=n1)) {
            for (i2 in seq(along=effect)) {
                for (i3 in seq(along= alpha)) {
                    za <- abs(e2[i2])/sqrt(temp1/n1[i1] + temp2/n2[i1])
                    out[, i1,, i2,i3] <- pnorm(za+ qnorm(alpha[i3]/2))
                }
            }
        }
    }
    drop(out)
}
