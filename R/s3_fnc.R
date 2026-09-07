
print.svyLocAdj <- function(x, ...){
 cat("Survey spatial location displacement model\n")
 cat(paste0("Model family: ",x$results$family,"\n"))
 waic <- data.frame(waic_approximation = x$results$waic_approx)
 row.names(waic) <- c("Est.")
 print(x$results$fixed_parameters)
 cat(paste0("Watanabe-Akaike information criterion: \n"))
 print(t(waic))
}

plot.svyLocAdj <- function(x, type = "map",
                           title = "Forest Plot",
                           drop_intercept = TRUE,
                           ...){
  # type=="para" or "map"
  if(type=="map"){
    print(x$sp_data$plot)
  }
  if(type=="para"){
   # forest plot
    if(x$results$family%in%"binomial"){
      df <- x$results$fixed_parameters
      df <- data.frame(
        variable = rownames(df),
        est = df[, "OR_median"],
        lower = df[, "OR_lower_95"],
        upper = df[, "OR_upper_95"]
      )
      if(isTRUE(drop_intercept)){ df = df[-1,]}
      # forest plot
      p <- ggplot(df, aes(x = est, y = variable)) +
        geom_point(size = 3.5, shape = 21, fill = "#2C7FB8", color = "black", stroke = 0.4) +
        geom_errorbarh(aes(xmin = lower, xmax = upper),
                       width = 0.15,
                       linewidth = 0.8,
                       color = "gray35") +
        geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.8, color = "gray50") +
        #scale_x_log10(
        #  breaks = scales::pretty_breaks(n = 6),
        #  labels = scales::label_number(accuracy = 0.01)
        #) +
        labs(
          x = "Odds Ratio",
          y = NULL,
          title = title
        ) +
        theme_minimal(base_size = 13) +
        theme(
          plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
          axis.text.y = element_text(size = 11),
          axis.text.x = element_text(size = 10),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          plot.title.position = "plot"
        )
      print(p)
    }
    if(x$results$family%in%"gaussian"){
      df <- x$results$fixed_parameters
      df <- data.frame(
        variable = rownames(df),
        est = df[, "median"],
        lower = df[, "lower_95"],
        upper = df[, "upper_95"]
      )
      if(isTRUE(drop_intercept)){ df = df[-1,]}
      # forest plot
      p <- ggplot(df, aes(x = est, y = variable)) +
        geom_point(size = 3.5, shape = 21, fill = "#2C7FB8", color = "black", stroke = 0.4) +
        geom_errorbarh(aes(xmin = lower, xmax = upper),
                       width = 0.15,
                       linewidth = 0.8,
                       color = "gray35") +
        geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.8, color = "gray50") +
        labs(
          x = "Estimate",
          y = NULL,
          title = title
        ) +
        theme_minimal(base_size = 13) +
        theme(
          plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
          axis.text.y = element_text(size = 11),
          axis.text.x = element_text(size = 10),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          plot.title.position = "plot"
        )
      print(p)
    }
    if(x$results$family%in%"gaussian_sp"){
      df <- x$results$spatial_parameters
      df <- data.frame(
        variable = rownames(df),
        est = df[, "median"],
        lower = df[, "lower_95"],
        upper = df[, "upper_95"]
      )
      #if(isTRUE(drop_intercept)){ df = df[-1,]}
      df <- df %>%
        mutate(
          significance = case_when(
            lower > 0 & upper > 0 ~ "Significant",
            lower < 0 & upper < 0 ~ "Significant",
            TRUE ~ "Not significant"
          )
        )
      # forest plot
      df <- df %>%
        mutate(
          variable_label = ifelse(
            significance == "Significant",
            as.character(variable),
            ""
          )
        )
      p <- ggplot(df, aes(x = est, y = variable)) +
        # All variables
        geom_point(
          size = 2.5,
          shape = 21,
          fill = "gray75",
          color = "gray50",
          stroke = 0.3
        ) +
        # Significant variables
        geom_errorbar(
          data = subset(df, significance == "Significant"),
          aes(xmin = lower, xmax = upper),
          width = 0.15,
          linewidth = 0.9,
          color = "#2C7FB8",
          orientation = "y"
        ) +
        geom_point(
          data = subset(df, significance == "Significant"),
          size = 3.5,
          shape = 21,
          fill = "#2C7FB8",
          color = "black",
          stroke = 0.4
        ) +
        geom_vline(
          xintercept = 0,
          linetype = "dashed",
          linewidth = 0.8,
          color = "gray50"
        ) +
        scale_y_discrete(
          labels = setNames(df$variable_label, df$variable)
        ) +
        labs(
          x = "Estimate",
          y = NULL,
          title = title
        ) +
        theme_minimal(base_size = 13) +
        theme(
          plot.title = element_text(
            face = "bold",
            size = 14,
            hjust = 0.5
          ),
          axis.text.y = element_text(
            size = 9,
            color = "black"
          ),
          axis.text.x = element_text(size = 10),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank()
        )
      print(p)
    }
    if(x$results$family%in%"poisson"){
      df <- x$results$fixed_parameters
      df <- data.frame(
        variable = rownames(df),
        est = df[, "IRR_median"],
        lower = df[, "IRR_lower_95"],
        upper = df[, "IRR_upper_95"]
      )
      if(isTRUE(drop_intercept)){ df = df[-1,]}
      # forest plot
      p <- ggplot(df, aes(x = est, y = variable)) +
        geom_point(size = 3.5, shape = 21, fill = "#2C7FB8", color = "black", stroke = 0.4) +
        geom_errorbarh(aes(xmin = lower, xmax = upper),
                       width = 0.15,
                       linewidth = 0.8,
                       color = "gray35") +
        geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.8, color = "gray50") +
        labs(
          x = "Risk Ratio",
          y = NULL,
          title = title
        ) +
        theme_minimal(base_size = 13) +
        theme(
          plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
          axis.text.y = element_text(size = 11),
          axis.text.x = element_text(size = 10),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          plot.title.position = "plot"
        )
      print(p)
    }
  }
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
  waic <- data.frame(waic_approximation = object$results$waic_approx)
  row.names(waic) <- c("Est.")
  cat(paste0("Watanabe–Akaike information criterion \n"))
  print(t(waic))
  cat(paste0("----------------------------\n"))
}
