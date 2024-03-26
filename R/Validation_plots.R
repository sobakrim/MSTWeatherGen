
########### plots #############
observed_vs_simulated_density = function(sim, observed, dates, seasons, j, names, names_seasons=NULL) {
  
  require(ggplot2)
  require(patchwork)
  
  id = which(sim==0, arr.ind = T)
  id = id[id[,3]==5,]
  sim[id] = NaN
  
  id = which(observed==0, arr.ind = T)
  id = id[id[,3]==5,]
  observed[id] = NaN
  
  combined_plot = NULL
  names_seasons = if(is.null(names_seasons)){
    names(seasons)
  }
  for (s in 1:length(seasons)) {
    season_indices = season_indices(dates, seasons[[s]])
    
    combined_data = lapply(names, function(v) {
      sim_data = data.frame(y = sim[season_indices, j, v], type = 'Simulated', variable = v)
      obs_data = data.frame(y = observed[season_indices, j, v], type = 'Observed', variable = v)
      rbind(sim_data, obs_data)
    })
    
    df = do.call(rbind, combined_data)
    
    p = ggplot(df, aes(x = y, color = type)) +
      geom_density() +
      facet_wrap(~variable, scales = "free") + 
      theme_bw() +
      labs(x = "Value", y = "Density", color = "Type", title = names_seasons[s]) +
      scale_color_manual("",values = c("Simulated" = "blue", "Observed" = "red")) +
      theme(plot.title = element_text(hjust = 0.5))
    
    if(is.null(combined_plot)) {
      combined_plot = p
    } else {
      combined_plot = combined_plot + p
    }
  }
  
  combined_plot = combined_plot + plot_layout(guides = 'collect')
  
  return(combined_plot)
  
}
Plot_dry_wet_spells_maps_season = function(sim, observed, coordinates, dates, seasons, names_seasons){
  require(ggplot2)
  require(patchwork)
  
  
  seasons_assigned <- assign_seasons(dates, seasons)
  change_season_indices <- c(which(diff(seasons_assigned) != 0) + 1, length(seasons_assigned) + 1)
  names_seasons <- if(is.null(names_seasons)){
    names(seasons)
  }else{
    names_seasons
  }
  
  df = lapply(1:(length(change_season_indices)-2), function(s){
    ido = change_season_indices[s]:(change_season_indices[s+1]-1)
    df = lapply(1:nrow(coordinates), function(j){
      
      prec_obs = observed[ido,j,"Precipitation"]
      prec_sim = sim[ido,j,"Precipitation"]
      qx = 0
      
      prec_obs[prec_obs>qx] = 1
      prec_sim[prec_sim>qx] = 1
      
      ro = rle(prec_obs)
      rs = rle(prec_sim)
      if(length(ro$lengths[ro$values==0])==0){
        dfo=NULL
      }else{
        dfo = data.frame(r = ro$lengths[ro$values==0], y = "Observed", season = seasons_assigned[change_season_indices[s]])
      }
      if(length(rs$lengths[rs$values==0])==0){
        dfs=NULL
      }else{
        dfs = data.frame(r = rs$lengths[rs$values==0], y = "Simulated", season = seasons_assigned[change_season_indices[s]])
      }
      df = rbind(dfo,dfs)
      library(dplyr)
      df_sum <- df %>%
        group_by(r, y, season) %>%
        tally()
      df_sum$v = with(df_sum, ifelse(r<5, 1, ifelse(r>=5 & r<10, 2, ifelse(r>=10 & r<20, 3, ifelse(r>=20,4,0)))))
      df_sum <- df_sum %>%
        group_by(y, v) %>% summarise(n=sum(n))
      df_sum$lon = coordinates[j,1]
      df_sum$lat = coordinates[j,2]
      df_sum$season =  names_seasons[seasons_assigned[change_season_indices[s]]]
      return(df_sum)
    })
    df = do.call(rbind, df)
    df$season =  names_seasons[seasons_assigned[change_season_indices[s]]]
    return(df)
  })
  df = do.call(rbind, df)
  df$season = factor(df$season, levels = names_seasons)
  p = ggplot(df[df$v==1,], aes(lon, lat)) +
    borders("world", colour="black", fill= "grey", xlim=range(df$lon), ylim = range(df$lat)) +
    coord_cartesian(xlim=range(df$lon), ylim = range(df$lat)) +
    geom_point(aes(color = n), size=7, shape=15) + xlab("Longitude (degree)") + 
    ylab("Latitude (degree)") + 
    scale_color_gradientn("Number of up to 5 consicutive dry days", colours = rev(viridis(10, option="magma"))) + 
    theme_bw() +
    theme(legend.position = "top") + 
    facet_grid(y ~ season) 
  return(p)
}

Plot_dry_wet_spells_maps = function(sim, observed, coordinates, dates){
  require(ggplot2)
  require(patchwork)
  
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
    library(dplyr)
    df_sum <- df %>%
      group_by(r, y) %>%
      tally() %>%
      mutate(v = case_when(
        r > 1 & r < 4 ~ 1,
        r >= 4 & r < 7 ~ 2,
        r >= 7 & r < 10 ~ 3,
        r >= 10 ~ 4,
        TRUE ~ 0
      )) %>%
      group_by(y, v) %>%
      summarise(n = sum(n), .groups = 'drop') %>%
      mutate(v_label = case_when(
        v == 1 ~ "1<NC<4",
        v == 2 ~ "4≤NC<7",
        v == 3 ~ "7≤NC<10",
        v == 4 ~ "NC≥10",
        TRUE ~ "Other"
      ))
    df_sum$lon = coordinates[j,1]
    df_sum$lat = coordinates[j,2]
    return(df_sum)
  })
  df = do.call(rbind, df)
  p = ggplot(df[df$v>0,], aes(lon, lat)) +
    borders("world", colour="black", fill= "grey", xlim=range(df$lon), ylim = range(df$lat)) +
    coord_cartesian(xlim=range(df$lon), ylim = range(df$lat)) +
    geom_point(aes(color = n), size=7, shape=15) + xlab("Longitude (degree)") + 
    ylab("Latitude (degree)") + 
    scale_color_gradientn("Number of consecutive wet days", colours = rev(viridis(10, option="magma"))) + 
    theme_bw() +
    theme(legend.position = "top") + 
    facet_grid(y ~ v_label) 
  
  return(p)
}


plot_mean_sd_by_month = function(sim, observed,places,names_places, names, dates){
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
  p = ggplot(df, aes(x=factor(month), y=mu, fill=x)) + 
    geom_boxplot(position="dodge", width=0.4)+theme_bw()+ xlab("Month")+ 
    facet_wrap(~location+variable,  scales = "free")
  return(p)
}


plot_multivariate_st_covariances = function(vgm, names){
  
  library(ggplot2)
  library(patchwork)
  
  # Create an empty plot placeholder for upper triangle
  empty_plot <- ggplot() + theme_void()
  
  # Initialize an empty plot matrix
  plot_matrix <- matrix(list(empty_plot), nrow = length(names), ncol = length(names))
  
  # Fill the matrix with plots, lower triangle and diagonal
  for (i in 1:length(names)) {
    for (j in 1:i) {
      variable1 <- names[i]
      variable2 <- names[j]
      filtered_vgm <- subset(vgm, v1 == variable1 & v2 == variable2 | v1 == variable2 & v2 == variable1)
      
      p <- ggplot(filtered_vgm, aes(x = dist, color = factor(lagtime))) +
        geom_line(aes(y = vgmsim)) +
        geom_point(aes(y = cov)) +
        labs(title = paste(variable1, "vs", variable2), x = "", y = "") +
        theme_minimal() +
        theme(plot.title = element_text(size = 8),
              legend.position = "none")
      last_row = length(names)
      mid_col = round(length(names)/2)
      # Add xlab to the bottom middle plot
      if (i == last_row && j == mid_col) {
        p <- p + labs(x = "Distance (degree)") + 
          theme(axis.text.x = element_text(), axis.ticks.x = element_line())
      }
      
      # Add ylab to the middle left plot
      if (i == mid_col && j == 1) {
        p <- p + labs(y = "Covariance") + 
          theme(axis.text.y = element_text(), axis.ticks.y = element_line())
      }
      plot_matrix[j, i] <- list(p)
    }
  }
  # Print the arranged plot layout
  
  legend_plot <- ggplot(filtered_vgm, aes(x = dist, color = factor(lagtime))) +
    geom_line(aes(y = vgmsim)) +
    geom_point(aes(y = cov)) +
    labs(title = paste(variable1, "vs", variable2), x = "", y = "") +
    theme_minimal() + scale_color_discrete("Time lag (day)") + 
    theme(plot.title = element_text(size = 8))
  leg = get_legend(legend_plot)
  num_vars <- length(names)
  
  plot_matrix[num_vars-1,2] <- list(leg)
  
  # Convert matrix to patchwork layout
  plot_layout <- wrap_plots(plot_matrix)
  
  return(plot_layout)
}

wet_frequency_months = function(sim, observed, dates, seasons, coordinates, names_seasons){
  require(ggplot2)
  require(viridis)
  require(ggpubr)
  require(dplyr)
  
  names_seasons = if(is.null(names_seasons)){
    names(seasons)
  }
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
    p = ggplot(df, aes(lon, lat )) +  borders("world", colour="black",fill= "grey",xlim=range(df$lon), ylim = range(df$lat)) +
      coord_cartesian(xlim=range(df$lon), ylim = range(df$lat))+
      scale_color_gradientn(name = "Frequency of wet days (%)", colours = rev(viridis(10, option = "magma"))) +
      geom_point(aes(color = Frequency),size=7, shape=15)+
      facet_wrap(~Type, scales = "free", ncol = 1) +
      theme_light() +
      theme(plot.title = element_text(size = 15, hjust = 0.5), 
            panel.spacing = unit(1, "lines"), aspect.ratio = 0.9, legend.position = "top") +
      xlab("Longitude (degree)") + ylab("Latitude (degree)") +
      ggtitle(names_seasons[s])
    
    plot_list[[s]] = p
  }
  
  # Combine all plots into one with a common legend
  combined_plot = ggarrange(plotlist = plot_list, ncol = length(seasons), common.legend = TRUE, legend = "top")
  return(combined_plot)
}

