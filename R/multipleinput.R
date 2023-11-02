simulatehook_multiple_inputs <- function(x,data,...) {
    minp <- x$attributes$multiple.inputs
    if (length(minp)>0) {
        for (i in seq_along(minp)) {            
            outcome <- names(minp[i])
            inp <- minp[[i]]$input
            fun <- minp[[i]]$fun
            data[,outcome] <- fun(x, data, inp)
        }
    }
    return(data)
}

addhook("simulatehook_multiple_inputs","sim.hooks")

printhook_multiple_inputs <- function(x,...) {
    minp <- x$attributes$multiple.inputs
    if (length(minp)>0) {
        outcomes <- names(minp)
        for (i in seq_along(minp)) {
            cat(minp[[i]]$type, ":\n\n")
            st <- paste0(outcomes[i]," ~ ", paste0(minp[[i]]$input,collapse=" | "))
            cat("  ", st, "\n")
            cat("\n")
        }
    }
    return(NULL)
}

addhook("printhook_multiple_inputs","print.hooks")
