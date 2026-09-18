sim_ac_area_model <- function(Smax, I, E, V1, V2, L, G, A, time, plot_SI = FALSE, plot_SA = FALSE, plot_SC = FALSE, plot_ST = FALSE, plot_ALL = FALSE) {
  # Smax: richness of species in the sps pool
  # I: immigration per linage per million years
  # E: extinction per linage per million years
  # V1: anagenesis rate per linage per million years
  # V2: cladogenesis rate per linage per million years
  # L: parameter that affects E
  # G: parameter that affects V2
  # A: area of the insular system
  # time: time of the simulation
  time <- 1:time
  eq_community <- data.frame("Immigrant_sps_eq" = as.numeric(), "Anagenetic_sps_eq" = as.numeric(), "Cladogenetic_sps_eq" = as.numeric(), "Total_sps_eq" = as.numeric())
  plot_data <- data.frame("Richness" = as.numeric(), "Type" = as.character(), "Time" = as.numeric())

  SI_pred <- (I * Smax) / (I + E * A^(-L) + V1 + V2 * A^G) * (1 - exp(-(I + E * A^(-L) + V1 + V2 * A^G) * time))
  
  if (any(!is.finite(SI_pred))) {
    warning("Combination of parameters produces non-finite values in immigrant species richness, try another")
  }
  if (any(!is.finite(SI_pred)) & (plot_SI == T | plot_ST == T | plot_ALL = T  )) {
    stop("Combination of parameters produces non-finite values in imigrabt species richness. Try another. GRAPH CANT BE PRODUCED")
  }
  if (length(unique(round(tail(SI_pred, 3), 3))) == 1) {
    Immigrant_sps_eq <- 1
  } else {
    Immigrant_sps_eq <- 0
  }
  
  plot_data <- rbind(plot_data, data.frame("Richness" = SI_pred, "Type" = "Immigrants", "Time" = time))
  
  SA_pred <- (V1 * I * Smax / (I + E * A^(-L) + V1 + V2 * A^G)) * ((1 / (E * A^(-L) + V2 * A^G)) * (1 - exp(-(E * A^(-L) + V2 * A^G) * time)) + (1 / (I + V1)) * (exp(-(I + E * A^(-L) + V1 + V2 * A^G) * time) - exp(-(E * A^(-L) + V2 * A^G) * time)))

  if (any(!is.finite(SA_pred))) {
    warning("Combination of parameters produces non-finite values in anagenetic species richness, try another")
  }
  if (any(!is.finite(SA_pred)) & (plot_SA == T | plot_ST == T | plot_ALL = T )) {
    stop("Combination of parameters produces non-finite values in anagenetic species richness. Try another. GRAPH CANT BE PRODUCED")
  }
  if (length(unique(round(tail(SA_pred, 3), 3))) == 1) {
    Anagenetic_sps_eq <- 1
  } else {
    Anagenetic_sps_eq <- 0
  }
  
  plot_data <- rbind(plot_data, data.frame("Richness" = SA_pred, "Type" = "Anagenetic", "Time" = time))
  
  SC_pred <-  (V2 * A^G * I * Smax / (I + E * A^(-L) + V1 + V2 * A^G)) * ((2 / (E * A^(-L) - V2 * A^G)) * (1 + V1 / (E * A^(-L) + V2 * A^G)) * (1 - exp((V2 * A^G - E * A^(-L)) * time))
                                                                          + (2 / (2 * V2 * A^G + I + V1)) * (1 - V1 / (I + V1)) * (exp((V2 * A^G - E * A^(-L)) * time) - exp(-(I + E * A^(-L) + V1 + V2 * A^G) * time))
                                                                          + V1 /(V2 * A^G) * (1 / (E * A^(-L) + V1) + 1 / (V1 + I)) * (exp(-(V2 * A^G + E * A^(-L)) * time) - exp((V2 * A^G - E * A^(-L)) * time)))
  
    
  if (any(!is.finite(SC_pred))) {
    warning("Combination of parameters produces non-finite values in cladogenetic species richness, try another")
  }  
  if (any(!is.finite(SC_pred)) & (plot_SC == T | plot_ST == T | plot_ALL = T )) {
    stop("Combination of parameters produces non-finite values in cladogenetic species richness. Try another. GRAPH CANT BE PRODUCED")
  }
  if (length(unique(round(tail(SC_pred, 3), 3))) == 1) {
    Cladogenetic_sps_eq <- 1
  } else {
    Cladogenetic_sps_eq <- 0
  }
  
  plot_data <- rbind(plot_data, data.frame("Richness" = SC_pred, "Type" = "Cladogenetic", "Time" = time))
  
  ST_pred <- SI_pred + SA_pred + SC_pred
  
  if (length(unique(round(tail(ST_pred, 3), 3))) == 1) {
    Total_sps_eq <- 1
  } else {
    Total_sps_eq <- 0
  }
  eq_community <- data.frame("Immigrant_sps_eq" = Immigrant_sps_eq, "Anagenetic_sps_eq" = Anagenetic_sps_eq, "Cladogenetic_sps_eq" = Cladogenetic_sps_eq, "Total_sps_eq" = Total_sps_eq)
  
  
  plot_data <- rbind(plot_data, data.frame("Richness" = ST_pred, "Type" = "Total", "Time" = time))
  
  output_data <- data.frame("Immigrant_sps" = SI_pred, "Anagenetic_sps" = SA_pred, "Cladogenetic_sps" = SC_pred, "Total_sps" = ST_pred, "Time" = time)
  
  
  if (plot_SI == T | plot_ALL == T) {
    SI_plot <- ggplot(output_data, aes(x = Time, y = Immigrant_sps)) +
      theme_classic() +
      labs(
        x = "Time (My)", y = "Immigrant species richness",
        title = "Dynamics of\nimmigrant species"
      ) +
      geom_line(size = 1, colour = "#CD3278")
  } else {
    SI_plot <- NULL
  }
  
  if (plot_SA == T | plot_ALL == T) {
    SA_plot <- ggplot(output_data, aes(x = Time, y = Anagenetic_sps)) +
      geom_line(size = 1, colour = "#CD3278") +
      theme_classic() +
      labs(
        x = "Time (My)", y = "Anagenetic species richness",
        title = "Dynamics of\nanagenetic species"
      )
  } else {
    SA_plot <- NULL
  }
  
  if (plot_SC == T | plot_ALL == T) {
    SC_plot <- ggplot(output_data, aes(x = Time, y = Cladogenetic_sps)) +
      geom_line(size = 1, colour = "#CD3278") +
      theme_classic() +
      labs(
        x = "Time (My)", y = "Cladogenetic species richness",
        title = "Dynamics of\ncladogenetic species"
      )
  } else {
    SC_plot <- NULL
  }
  
  if (plot_ST == T | plot_ALL == T) {
    ST_plot <- ggplot(output_data, aes(x = Time, y = Total_sps)) +
      geom_line(size = 1, colour = "#CD3278") +
      theme_classic() +
      labs(
        x = "Time (My)", y = "Total species richness",
        title = "Dynamics of\ntotal community"
      )
  } else {
    ST_plot <- NULL
  }
  if (plot_SI == FALSE & plot_SA == FALSE & plot_SC == FALSE & plot_ST == TRUE){
    print(ST_plot)
  }else{
    all_graphs <- (SI_plot | SA_plot | SC_plot | ST_plot)
    print(all_graphs)}
  
  return(list(output_data, eq_community))
}
