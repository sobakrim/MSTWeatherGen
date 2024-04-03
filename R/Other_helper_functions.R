#' @importFrom lubridate year leap_year
season_indices = function(dates, season, Year){
  
  years = unique(lubridate::year(dates))
  years = c(years, max(years) + 1)
  md = lapply(years, function(Year) {
    if(Year == years[1]){
      is_leap = lubridate::leap_year(Year)
      
      # Adjust min_day for non-leap years if necessary
      adjusted_min_day = ifelse(!is_leap && season$min_month == 2 && season$min_day == 29, 28, season$min_day)
      adjusted_max_day = ifelse(!is_leap && season$max_month == 2 && season$max_day == 29, 28, season$max_day)
      
      # Check if the season spans the end of the year
      if (season$min_month > season$max_month || (season$min_month == season$max_month && adjusted_min_day > adjusted_max_day)) {
        d1 = as.Date(paste(Year, min(month(dates[lubridate::year(dates) %in% years[1]])), adjusted_min_day, sep = "-"))
        d2 = as.Date(paste(Year, season$max_month, adjusted_max_day, sep = "-"))
      }else{
        d1 = as.Date(paste(Year, season$min_month, adjusted_min_day, sep = "-"))
        d2 = as.Date(paste(Year, season$max_month, adjusted_max_day, sep = "-"))
      }
      
      return(seq(d1, d2, by = "day"))
    }else{
      if(Year == max(years)){
        Year = Year - 1
        is_leap = lubridate::leap_year(Year)
        
        # Adjust min_day for non-leap years if necessary
        adjusted_min_day = ifelse(!is_leap && season$min_month == 2 && season$min_day == 29, 28, season$min_day)
        adjusted_max_day = ifelse(!is_leap && season$max_month == 2 && season$max_day == 29, 28, season$max_day)
        
        # Check if the season spans the end of the year
        if (season$min_month > season$max_month || (season$min_month == season$max_month && adjusted_min_day > adjusted_max_day)) {
          d1 = as.Date(paste(Year, season$min_month, adjusted_min_day, sep = "-"))
          d2 = max(dates)
          if(d1>d2){
            return(NULL)
          }else{
            return(seq(d1, d2, by = "day"))
          }
        }else{
          d1 = as.Date(paste(Year, season$min_month, adjusted_min_day, sep = "-"))
          d2 = as.Date(paste(Year, season$max_month, adjusted_max_day, sep = "-"))
          
        }
        
      }else{
        
        # Check for leap year
        is_leap = lubridate::leap_year(Year)
        
        # Adjust min_day for non-leap years if necessary
        adjusted_min_day = ifelse(!is_leap && season$min_month == 2 && season$min_day == 29, 28, season$min_day)
        d1 = as.Date(paste(Year, season$min_month, adjusted_min_day, sep = "-"))
        
        # Adjust max_day for non-leap years if necessary
        if (season$min_month > season$max_month){
          is_leap = lubridate::leap_year(Year)
          adjusted_max_day = ifelse(!is_leap && season$max_month == 2 && season$max_day == 29, 28, season$max_day)
          is_leap = lubridate::leap_year(Year-1)
          adjusted_min_day = ifelse(!is_leap && season$min_month == 2 && season$min_day == 29, 28, season$min_day)
          d2 = as.Date(paste(Year, season$max_month, adjusted_max_day, sep = "-"))
          d1 = as.Date(paste(Year-1, season$min_month, adjusted_min_day, sep = "-"))
        }else{
          is_leap = lubridate::leap_year(Year)
          adjusted_max_day = ifelse(!is_leap && season$max_month == 2 && season$max_day == 29, 28, season$max_day)
          d2 = as.Date(paste(Year, season$max_month, adjusted_max_day, sep = "-"))
        }
        
        return(seq(d1, d2, by = "day"))
      }
    }
  })
  md = do.call(c, md)
  md = format(md, "%Y-%m-%d CET")
  dates = format(dates, "%Y-%m-%d CET")
  if(missing(Year)){
    return(which(dates %in% md))
  }else{
    return(which(dates %in% md & lubridate::year(dates)==Year))
  }
}
filter_season_data <- function(data, dates, season, names) {
  # Ensure 'dates' is converted to Date class
  dates <- as.Date(dates)
  
  # Obtain indices for the specified season
  # Assuming 'season_indices' is a function that returns indices of dates within the specified season
  season_indices <- season_indices(dates, season)
  
  # Filter data and dates based on season indices
  data_filtered <- data[season_indices,,]
  dates_filtered <- dates[season_indices]
  
  # Further filter data to include only the specified variables (names)
  data_filtered <- data_filtered[,,names]
  
  # Order the filtered data by dates
  data_filtered <- data_filtered[order(dates_filtered),,, drop = FALSE]
  
  # Return the filtered and ordered data along with the corresponding dates
  return(list(data_filtered = data_filtered, dates_filtered = dates_filtered))
}

haversine <- function(point1, point2) {
  # Define the Haversine function
  # Convert degrees to radians
  lon1 <- c(point1[,1] * pi / 180)
  lat1 <- c(point1[,2] * pi / 180)
  lon2 <- c(point2[,1] * pi / 180)
  lat2 <- point2[,2] * pi / 180
  
  # Differences in coordinates
  dlon <- lon2 - lon1
  dlat <- lat2 - lat1
  
  # Haversine formula
  a <- as.numeric(sin(dlat / 2)^2 + cos(lat1) * cos(lat2) * sin(dlon / 2)^2)
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  
  # Radius of Earth in kilometers
  R <- 6371
  
  # Distance in kilometers
  d <- R * c
  
  return(d)
}

#' @importFrom geosphere distHaversine
ds = function(i,j,coordinates) {
  return(geosphere::distHaversine(coordinates[i,], coordinates[j,])/1000)
}
colMins = function(x) {
  return(apply(x, 2, min))
}
colMaxs = function(x) {
  return(apply(x, 2, max))
}

colSds = function(x) {
  return(apply(x, 2, sd))
}
calculate_ICI <- function(radiation, time) {
  ## ICI : Insolation Clearness Index
  day_of_year <- as.integer(format(time, "%j"))
  TOA_irradiance <- 1361 + 0.033 * cos((360*day_of_year)/365)
  # Calculate the Insolation Clearness Index (ICI)
  ICI <- radiation / TOA_irradiance
  
  return(ICI)
}

