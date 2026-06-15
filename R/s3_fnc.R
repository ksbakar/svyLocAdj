
print.svyLocAdj <- function(x, ...){
 cat("Survey spatial location displacement model\n")
 cat(paste0("Model family: ",x$results$family,"\n"))
 waic <- data.frame(lppd = x$results$lppd,
                    waic_approximation = x$results$waic_approx)
 row.names(waic) <- c("Est.")
 print(x$results$fixed_parameters)
 cat(paste0("Watanabe-Akaike information criterion: \n"))
 print(t(waic))
}

plot.svyLocAdj <- function(x, ...){
   print(x$sp_data$plot)
}

summary.svyLocAdj <- function(object, ...){
  cat("Survey spatial location displacement model\n")
  cat(paste0("Model family: ",object$results$family,"\n"))
  cat(paste0("----------------------------\n"))
  print(object$results$fixed_parameters)
  cat(paste0("----------------------------\n"))
  cat(paste0("Model variability parameters\n"))
  print(object$results$variability_parameters)
  cat(paste0("----------------------------\n"))
  waic <- data.frame(lppd = object$results$lppd,
          waic_approximation = object$results$waic_approx)
  row.names(waic) <- c("Est.")
  cat(paste0("Watanabe–Akaike information criterion \n"))
  print(t(waic))
  cat(paste0("----------------------------\n"))
}
