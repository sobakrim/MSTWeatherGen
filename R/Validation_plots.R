utils::globalVariables(c("y", "type"))
#' Plot Observed vs. Simulated Density
#'
#' Compares the density of observed and simulated data for specific weather variables across seasons.
#' It generates density plots for each variable, distinguishing between observed and simulated data,
#' and arranges them by season. This function relies on `ggplot2` for plotting and `patchwork` for
#' combining plots.
#'
#' @param sim A 3-dimensional array of simulated weather data with dimensions [time, location, variable].
#' @param observed A 3-dimensional array of observed weather data with dimensions [time, location, variable].
#' @param dates A vector of dates corresponding to the first dimension of `sim` and `observed`.
#' @param seasons A list where each element defines the start and end days and months of a season.
#' @param location The index of the location for which to plot the data.
#' @param names A vector of strings indicating the names of the variables to be plotted.
#' @param names_seasons (Optional) A vector of season names corresponding to the `seasons` list.
#' If `NULL`, the names are derived from the `seasons` list itself.
#'
#' @return A `ggplot` object representing the combined density plots for observed and simulated data
#' across the specified seasons and variables.
#'
#' @import ggplot2
#' @import patchwork
#' @examples
#' \dontrun{
#'   # Assuming 'sim', 'observed', 'dates', 'seasons', and 'coordinates' are already defined:
#'   
#'   # Define the seasons (example)
#'   seasons <- list(
#'     s1 = list(min_day = 1, max_day = 29, min_month = 12, max_month = 2),
#'     s2 = list(min_day = 1, max_day = 31, min_month = 3, max_month = 5),
#'     s3 = list(min_day = 1, max_day = 31, min_month = 6, max_month = 8),
#'     s4 = list(min_day = 1, max_day = 30, min_month = 9, max_month = 11)
#'   )
#'   
#'   # Define variable names to be plotted
#'   names <- c("Precipitation", "Wind", "Temp_max")
#'   
#'   # Specify the location index (for example, 1 for the first location)
#'   j <- 1
#'   
#'   # Run the plot function (this example assumes that the necessary data is loaded)
#'   plot_result <- plot_observed_vs_simulated_density(sim, observed, dates, seasons, j, names)
#'   
#'   # Print or view the plot
#'   print(plot_result)
#' }
#' @export
plot_observed_vs_simulated_density = function(sim, observed, dates, seasons, location, names, names_seasons=NULL) {
  
  id = which(sim==0, arr.ind = T)
  id = id[id[,3]==5,]
  sim[id] = NaN
  
  id = which(observed==0, arr.ind = T)
  id = id[id[,3]==5,]
  observed[id] = NaN
  
  combined_plot = NULL
  names_seasons = if(is.null(names_seasons)){
    names(seasons)
  }else{names_seasons}
  for (s in 1:length(seasons)) {
    season_indice = season_indices(dates, seasons[[s]])
    
    combined_data = lapply(names, function(v) {
      sim_data = data.frame(y = sim[season_indice, location, v], type = 'Simulated', variable = v)
      obs_data = data.frame(y = observed[season_indice, location, v], type = 'Observed', variable = v)
      rbind(sim_data, obs_data)
    })
    
    df = do.call(rbind, combined_data)
    
    p = ggplot2::ggplot(df, ggplot2::aes(x = y, color = type)) +
      ggplot2::geom_density() +
      ggplot2::facet_wrap(~variable, scales = "free") + 
      ggplot2::theme_bw() +
      ggplot2::labs(x = "Value", y = "Density", color = "Type", title = names_seasons[s]) +
      ggplot2::scale_color_manual("",values = c("Simulated" = "blue", "Observed" = "red")) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    
    if(is.null(combined_plot)) {
      combined_plot = p
    } else {
      combined_plot = combined_plot + p
    }
  }
  
  combined_plot = combined_plot + patchwork::plot_layout(guides = 'collect')
  
  return(combined_plot)
  
}

utils::globalVariables(c("r", "y", "v", "lon", "lat"))
#' Plot Dry and Wet Spells Maps
#'
#' Generates maps showing the frequency of consecutive wet days for observed and simulated precipitation data. This function uses precipitation data to classify days as wet or dry based on a threshold and calculates the length of consecutive wet day spells. The results are then plotted on a map, with points colored by the number of spells and faceted by observation type and spell length category.
#'
#' @param sim A 3-dimensional array of simulated weather data with dimensions [time, location, variable], where "Precipitation" is expected to be one of the variables.
#' @param observed A 3-dimensional array of observed weather data with the same dimensions and variable expectations as `sim`.
#' @param coordinates A matrix or data frame of geographic coordinates for the locations represented in `sim` and `observed`, with longitude in the first column and latitude in the second.
#' @param dates A vector of dates corresponding to the time dimension in `sim` and `observed`. This parameter is currently not used in the function but may be intended for future use or filtering within the provided data.
#'
#' @return A `ggplot` object representing the maps of dry and wet spells for observed and simulated precipitation data across the specified locations.
#'
#' @import ggplot2
#' @import patchwork
#' @import dplyr
#' @import viridis
#' @examples
#' \dontrun{
#'   # Assuming `sim`, `observed`, and `coordinates` are already defined and correctly formatted:
#'   plot_dry_wet_spells_maps(sim, observed, coordinates, dates)
#' }
#' @export
plot_dry_wet_spells_maps = function(sim, observed, coordinates, dates){

  df = lapply(1:nrow(coordinates), function(j){
    
    prec_obs = observed[,j,"Precipitation"]
    prec_sim = sim[,j,"Precipitation"]
    qx = 0
    
    prec_obs[prec_obs>qx] = 1
    prec_sim[prec_sim>qx] = 1
    
    ro = rle(prec_obs)
    rs = rle(prec_sim)
    if(length(ro$lengths[ro$values==0])==0){
      dfo=NULL
    }else{
      dfo = data.frame(r = ro$lengths[ro$values==1], y = "Observed")
    }
    if(length(rs$lengths[rs$values==0])==0){
      dfs=NULL
    }else{
      dfs = data.frame(r = rs$lengths[rs$values==1], y = "Simulated")
    }
    df = rbind(dfo,dfs)
    df_sum <- df %>%
      dplyr::group_by(r, y) %>%
      dplyr::tally() %>%
      dplyr::mutate(v = dplyr::case_when(
        r > 1 & r < 4 ~ 1,
        r >= 4 & r < 7 ~ 2,
        r >= 7 & r < 10 ~ 3,
        r >= 10 ~ 4,
        TRUE ~ 0
      )) %>%
      dplyr::group_by(y, v) %>%
      dplyr::summarise(n = sum(n), .groups = 'drop') %>%
      dplyr::mutate(v_label = dplyr::case_when(
        v == 1 ~ "1 < NC < 4",
        v == 2 ~ "4 <= NC < 7",
        v == 3 ~ "7 <= NC < 10",
        v == 4 ~ "NC >= 10",
        TRUE ~ "Other"
      )
      )
    df_sum$lon = coordinates[j,1]
    df_sum$lat = coordinates[j,2]
    return(df_sum)
  })
  df = do.call(rbind, df)
  p = ggplot2::ggplot(df[df$v>0,], ggplot2::aes(lon, lat)) +
    ggplot2::borders("world", colour="black", fill= "grey", xlim=range(df$lon), ylim = range(df$lat)) +
    ggplot2::coord_cartesian(xlim=range(df$lon), ylim = range(df$lat)) +
    ggplot2::geom_point(ggplot2::aes(color = n), size=7, shape=15) + ggplot2::xlab("Longitude (degree)") + 
    ggplot2::ylab("Latitude (degree)") + 
    ggplot2::scale_color_gradientn("Number of consecutive wet days (NC)", colours = rev(viridis::viridis(10, option="magma"))) + 
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "top") + 
    ggplot2::facet_grid(y ~ v_label) 
  
  return(p)
}
utils::globalVariables(c("mu", "x"))
#' Plot Monthly Mean Comparisons
#'
#' Compares the monthly mean values of observed and simulated weather data for specified variables and locations. The function generates boxplots to visually compare the means for each month, facilitating an easy comparison between observed and simulated data across different variables and locations.
#'
#' @param sim A 3-dimensional array of simulated weather data with dimensions [time, location, variable].
#' @param observed A 3-dimensional array of observed weather data with the same dimensions as `sim`.
#' @param places A vector of indices representing the locations for which the data will be plotted. If `places` has only one element, the function generates plots for that specific location; otherwise, it generates plots for each location specified.
#' @param names_places A vector of names corresponding to the `places` indices, used for labeling in the plots.
#' @param names A vector of variable names to be included in the comparison. These should correspond to the third dimension of the `sim` and `observed` arrays.
#' @param dates A vector of dates corresponding to the first dimension of the `sim` and `observed` arrays. This is used to determine the months and years for grouping the data.
#'
#' @return A `ggplot` object containing the boxplots comparing the monthly mean values of the observed and simulated data for the specified variables and locations.
#'
#' @examples
#' \dontrun{
#'   # Assuming `sim`, `observed`, `dates`, `places`, and `names` are predefined:
#'   plot <- plot_mean_by_month(sim, observed, places, names_places, names, dates)
#'   print(plot)
#' }
#'
#' @import ggplot2
#' @importFrom lubridate month year
#' @export
plot_mean_by_month = function(sim, observed,places,names_places, names, dates){
  if(length(places) == 1){
    df = lapply(names, function(variable){
      months = unique(month(dates))
      years = unique(year(dates))
      df = lapply(months, function(m){
        df = lapply(years, function(y){
          month_indices = which(month(dates) %in% m & year(dates) %in% y)
          xt = observed[month_indices,places,variable]
          xs = sim[month_indices,places,variable]
          df = data.frame(month = m, mu = mean(xt), x = "Observed", variable = variable, 
                          som = "mean", location = names_places)
          df = rbind(df,data.frame(month = m, mu = mean(xs), x = "Simulated", variable= variable,
                                   som = "mean", location = names_places))
          return(df)
        })
        return(do.call(rbind, df))
      })
      df = do.call(rbind, df)
    })
    df = do.call(rbind, df)
    p = ggplot2::ggplot(df, ggplot2::aes(x=factor(month), y=mu, fill=x)) + ggplot2::scale_fill_discrete("")+
      ggplot2::geom_boxplot(position="dodge", width=0.4)+ggplot2::theme_bw()+ ggplot2::xlab("Month")+ ggplot2::ylab("Mean")+
      ggplot2::facet_wrap(~variable,  scales = "free")
  }else{
    df = lapply(1:length(places), function(j){
      df = lapply(names, function(variable){
        months = unique(month(dates))
        years = unique(year(dates))
        df = lapply(months, function(m){
          df = lapply(years, function(y){
            month_indices = which(month(dates) %in% m & year(dates) %in% y)
            xt = observed[month_indices,places[j],variable]
            xs = sim[month_indices,places[j],variable]
            df = data.frame(month = m, mu = mean(xt), x = "Observed", variable = variable, 
                            som = "mean", location = names_places[j])
            df = rbind(df,data.frame(month = m, mu = mean(xs), x = "Simulated", variable= variable,
                                     som = "mean", location = names_places[j]))
            return(df)
          })
          return(do.call(rbind, df))
        })
        df = do.call(rbind, df)
      })
      df = do.call(rbind, df)
    })
    df = do.call(rbind, df)
    p = ggplot2::ggplot(df, ggplot2::aes(x=factor(month), y=mu, fill=x)) + ggplot2::scale_fill_discrete("")+ 
      ggplot2::geom_boxplot(position="dodge", width=0.4)+ggplot2::theme_bw()+ ggplot2::xlab("Month")+ ggplot2::ylab("Mean") +
      ggplot2::facet_wrap(~location+variable,  scales = "free")
  }
  return(p)
}

utils::globalVariables(c("v1", "v2", "lagtime", "vgmsim"))
#' Plot multivariate space-time covariances
#'
#' Function plot multivariate space-time covariances.
#'
#' @import ggplot2
#' @import patchwork
#' @import ggpubr
#' @keywords internal
plot_multivariate_st_covariances = function(vgm, names){
  
  # Create an empty plot placeholder for upper triangle
  empty_plot <- ggplot2::ggplot() + ggplot2::theme_void()
  
  # Initialize an empty plot matrix
  plot_matrix <- matrix(list(empty_plot), nrow = length(names), ncol = length(names))
  
  # Fill the matrix with plots, lower triangle and diagonal
  for (i in 1:length(names)) {
    for (j in 1:i) {
      variable1 <- names[i]
      variable2 <- names[j]
      filtered_vgm <- subset(vgm, v1 == variable1 & v2 == variable2 | v1 == variable2 & v2 == variable1)
      
      p <- ggplot2::ggplot(filtered_vgm, ggplot2::aes(x = dist, color = factor(lagtime))) +
        ggplot2::geom_line(ggplot2::aes(y = vgmsim)) +
        ggplot2::geom_point(ggplot2::aes(y = cov)) +
        ggplot2::labs(title = paste(variable1, "vs", variable2), x = "", y = "") +
        ggplot2::theme_minimal() +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 8),
              legend.position = "none")
      last_row = length(names)
      mid_col = round(length(names)/2)
      # Add xlab to the bottom middle plot
      if (i == last_row && j == mid_col) {
        p <- p + ggplot2::labs(x = "Distance (degree)") + 
          ggplot2::theme(axis.text.x = ggplot2::element_text(), axis.ticks.x = ggplot2::element_line())
      }
      
      # Add ylab to the middle left plot
      if (i == mid_col && j == 1) {
        p <- p + ggplot2::labs(y = "Covariance") + 
          ggplot2::theme(axis.text.y = ggplot2::element_text(), axis.ticks.y = ggplot2::element_line())
      }
      plot_matrix[j, i] <- list(p)
    }
  }
  # Print the arranged plot layout
  
  legend_plot <- ggplot2::ggplot(filtered_vgm, ggplot2::aes(x = dist, color = factor(lagtime))) +
    ggplot2::geom_line(ggplot2::aes(y = vgmsim)) +
    ggplot2::geom_point(ggplot2::aes(y = cov)) +
    ggplot2::labs(title = paste(variable1, "vs", variable2), x = "", y = "") +
    ggplot2::theme_minimal() + ggplot2::scale_color_discrete("Time lag (day)") + 
    ggplot2::theme(plot.title = ggplot2::element_text(size = 8))
  leg = ggpubr::get_legend(legend_plot)
  num_vars <- length(names)
  
  plot_matrix[length(names),1] <- list(leg)
  
  # Convert matrix to patchwork layout
  plot_layout <- patchwork::wrap_plots(plot_matrix)
  
  return(plot_layout)
}

utils::globalVariables(c("lon", "lat", "Frequency"))
#' Plot Wet Day Frequency for Observed and Simulated Data
#'
#' Generates plots comparing the frequency of wet days between observed and simulated precipitation data across specified seasons and geographic locations. Each season is visualized in a separate plot, with observed and simulated data presented side-by-side for direct comparison. The function uses a threshold to classify days as wet or dry and calculates the percentage of wet days for each location.
#'
#' @param sim A 3-dimensional array of simulated precipitation data with dimensions [time, location, variable].
#' @param observed A 3-dimensional array of observed precipitation data with the same dimensions as `sim`.
#' @param dates A vector of dates corresponding to the first dimension of the `sim` and `observed` arrays.
#' @param seasons A list where each element defines the start and end days and months of a season.
#' @param coordinates A data frame or matrix containing the geographic coordinates (longitude and latitude) for each location in the `sim` and `observed` arrays.
#' @param names_seasons (Optional) A vector of names corresponding to the seasons defined in the `seasons` list. If `NULL`, the names are derived from the `seasons` list itself.
#'
#' @return A combined `ggplot` object that includes the plots for all seasons, arranged side-by-side, each comparing the wet day frequency between observed and simulated data for the specified locations.
#'
#' @examples
#' \dontrun{
#'   # Assuming `sim`, `observed`, `dates`, `seasons`, and `coordinates` are predefined:
#'   plot_wet_frequency(sim, observed, dates, seasons, coordinates, c("Spring", "Summer", "Fall", "Winter"))
#' }
#'
#' @import ggplot2
#' @import viridis
#' @import ggpubr
#' @import dplyr
#' @export

plot_wet_frequency = function(sim, observed, dates, seasons, coordinates, names_seasons){
  
  names_seasons = if(is.null(names_seasons)){
    names(seasons)
  }else{names_seasons}
  plot_list = list()
  
  for (s in 1:length(seasons)) {
    season_indices = season_indices(dates, seasons[[s]])
    
    # Preprocess observed data
    prec_obs = observed[season_indices, , "Precipitation"]
    prec_obs = ifelse(prec_obs <= 0, 0, 1)
    
    # Preprocess simulated data
    prec_sim = sim[season_indices, , "Precipitation"]
    prec_sim = ifelse(prec_sim <= 0, 0, 1)
    
    # Combine into one data frame
    df_obs = data.frame(Type = "Observed", lon = coordinates$longitude, lat = coordinates$latitude, Frequency = 100 * colSums(prec_obs) / nrow(prec_obs))
    df_sim = data.frame(Type = "Simulated", lon = coordinates$longitude, lat = coordinates$latitude, Frequency = 100 * colSums(prec_sim) / nrow(prec_sim))
    df = rbind(df_obs, df_sim)
    
    # Create the plot
    p = ggplot2::ggplot(df, ggplot2::aes(lon, lat )) + 
      ggplot2::borders("world", colour="black",fill= "grey",xlim=range(df$lon), ylim = range(df$lat)) +
      ggplot2::coord_cartesian(xlim=range(df$lon), ylim = range(df$lat))+
      ggplot2::scale_color_gradientn(name = "Frequency of wet days (%)", colours = rev(viridis::viridis(10, option = "magma"))) +
      ggplot2::geom_point(ggplot2::aes(color = Frequency),size=7, shape=15)+
      ggplot2::facet_wrap(~Type, scales = "free", ncol = 1) +
      ggplot2::theme_light() +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 15, hjust = 0.5), 
            panel.spacing = ggplot2::unit(1, "lines"), aspect.ratio = 0.9, legend.position = "top") +
      ggplot2::xlab("Longitude (degree)") + ggplot2::ylab("Latitude (degree)") +
      ggplot2::ggtitle(names_seasons[s])
    
    plot_list[[s]] = p
  }
  
  # Combine all plots into one with a common legend
  combined_plot = ggpubr::ggarrange(plotlist = plot_list, ncol = length(seasons), common.legend = TRUE, legend = "top")
  return(combined_plot)
}

