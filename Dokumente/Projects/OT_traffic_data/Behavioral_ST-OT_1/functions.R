############################################## DATA PREPROC

#Process data for the weekdays
process_data_wd = function(raw_data){
  # Clean, time-stemp, use morning rush hour window (7 AM to 11 AM)
  study_window <- raw_data %>%
    mutate(
      start_time = as.POSIXct(started_at),
      end_time = as.POSIXct(ended_at),
      # Create 1-hour time bins
      T_Start = hour(start_time),
      T_End = hour(end_time),
      # other vars
      wd = strftime(as.Date(start_time), format = "%w"), #weekday: Sun = 0
      doy = as.numeric(strftime(as.Date(start_time), format = "%j"))#day of year
    ) %>%
    filter(T_Start >= 7, T_End <= 12) %>%
    # Filter out trips that start and end at the exact same station 
    filter(start_station_id != end_station_id)
  
  # Join with clustered station lookup table to the trip data
  trips_mapped <- study_window %>%
    left_join(stations_final %>% select(station_id, start_neigh = neighborhood_id),
              by = c("start_station_id" = "station_id")) %>%
    left_join(stations_final %>% select(station_id, end_neigh = neighborhood_id),
              by = c("end_station_id" = "station_id")) %>%
    drop_na(start_neigh, end_neigh)
  
  # True O-D
  true_od_matrix <- trips_mapped %>%
    count(start_neigh, end_neigh, T_Start, T_End, wd, name = "actual_trips")
  
  # Departures per Neighborhood per Hour
  supply_st <- trips_mapped %>%
    count(start_neigh, T_Start, wd, name = "S") %>%
    mutate(Node_ID = paste0(start_neigh, "_T", T_Start, "_", wd))
  
  # Arrivals per Neighborhood per Hour
  demand_st <- trips_mapped %>%
    count(end_neigh, T_End, wd, name = "D") %>%
    mutate(Node_ID = paste0(end_neigh, "_T", T_End, "_", wd))
  
  
  # Create the 10 neighborhoods x 5 hours grid
  time_bins <- 7:11 # hours to consider
  st_grid <- expand_grid(neighborhood_id = 1:10, Time = time_bins, wd = 0:6) %>%
    mutate(Node_ID = paste0(neighborhood_id, "_T", Time, "_", wd)) %>%
    left_join(supply_st %>% select(Node_ID, S), by = "Node_ID") %>%
    left_join(demand_st %>% select(Node_ID, D), by = "Node_ID") %>%
    mutate(S = replace_na(as.numeric(S), 0.1), D = replace_na(as.numeric(D), 0.1))
  
  out = list(st_grid = st_grid, trips_mapped = trips_mapped)
  return(out)
}

process_data_doy = function(raw_data25){
  # Clean, time-stemp, use morning rush hour window (7 AM to 11 AM)
  study_window <- raw_data25 %>%
    mutate(
      start_time = as.POSIXct(started_at),
      end_time = as.POSIXct(ended_at),
      # Create 1-hour time bins
      T_Start = hour(start_time),
      T_End = hour(end_time),
      # other vars
      wd = strftime(as.Date(start_time), format = "%w"), #weekday: Sun = 0
      doy = as.numeric(strftime(as.Date(start_time), format = "%j"))#day of year
    ) %>%
    filter(T_Start >= 7, T_End <= 12) %>%
    # Filter out trips that start and end at the exact same station
    filter(start_station_id != end_station_id)
  # Join with the clustered station lookup table to the trip data
  trips_mapped25 <- study_window %>%
    left_join(stations_final %>% select(station_id, start_neigh = neighborhood_id),
              by = c("start_station_id" = "station_id")) %>%
    left_join(stations_final %>% select(station_id, end_neigh = neighborhood_id),
              by = c("end_station_id" = "station_id")) %>%
    drop_na(start_neigh, end_neigh)
  
  # True O-D
  true_od_matrix25 <- trips_mapped25 %>%
    count(start_neigh, end_neigh, T_Start, T_End, wd, name = "actual_trips")
  
  # Departures per Neighborhood per Hour
  supply_st <- trips_mapped25 %>%
    count(start_neigh, T_Start, doy, wd, name = "S") %>%
    mutate(Node_ID = paste0(start_neigh, "_T", T_Start, "_", wd))
  
  # Arrivals per Neighborhood per Hour
  demand_st <- trips_mapped25 %>%
    count(end_neigh, T_End, doy, wd, name = "D") %>%
    mutate(Node_ID = paste0(end_neigh, "_T", T_End, "_", wd))
  
  # Create the 10 neighborhoods x 5 hours grid
  time_bins <- 7:11 # hours to consider
  doy <- sort(unique(trips_mapped25$doy)); doy = doy[doy<=365] # daily basis
  doy_wd <- trips_mapped25[,c("doy","wd")]
  wd <- sapply(doy, function(x){doy_wd$wd[which(doy_wd$doy == x)][1]})
  
  st_grid25 <- expand_grid(neighborhood_id = 1:10, Time = time_bins, doy = doy) %>%
    mutate(Node_ID = paste0(neighborhood_id, "_T", Time, "_", wd)) %>%
    left_join(supply_st %>% select(Node_ID, S, doy), by = c("Node_ID", "doy")) %>%
    left_join(demand_st %>% select(Node_ID, D, doy), by = c("Node_ID", "doy")) %>%
    mutate(S = replace_na(as.numeric(S), 0.1), D = replace_na(as.numeric(D), 0.1))
  
  out = list(st_grid = st_grid25, trips_mapped = trips_mapped25, doy = doy, wd = wd)
  return(out)
}



############################################## COST
compute_cost <- function(C_spatial_informed, st_grid, gamma = 1){
  n_nodes <- nrow(st_grid)
  C_st <- matrix(Inf, n_nodes, n_nodes)
  # Fill C_st using our 10x10 OSRM neighborhood travel times
  for(r in 1:n_nodes) {
    for(c in 1:n_nodes) {
      # Causal constraint: can't travel backward in time
      if(st_grid$Time[c] >= st_grid$Time[r]) {
        orig <- st_grid$neighborhood_id[r]
        dest <- st_grid$neighborhood_id[c]
        
        # Use OSRM durations (assuming C_neighborhoods is in hours)
        travel_time <- C_spatial_informed[orig, dest]
        
        # Cost = Physical Travel Time + Temporal Penalty
        C_st[r,c] <- (travel_time) + 
          (gamma * abs((st_grid$Time[c] - st_grid$Time[r]) - travel_time))
      }
    }
  }
  return(C_st)
}
############################################### Sinkhorn 
unbalanced_sinkhorn <- function(a, b, C, eps = 1, tau = 10, n_iter = 200) {
  a_mat <- matrix(as.numeric(a), ncol = 1)
  b_mat <- matrix(as.numeric(b), ncol = 1)
  C_mat <- as.matrix(C)
  
  if (nrow(C_mat) != nrow(a_mat) || ncol(C_mat) != nrow(b_mat)) stop("Dimension mismatch!")
  
  K <- exp(-C_mat / eps)
  K[is.nan(K) | is.na(K)] <- 0 
  u <- matrix(1, nrow = nrow(a_mat), ncol = 1)
  fi <- tau / (tau + eps)
  
  for (i in 1:n_iter) {
    Kv <- t(K) %*% u + 1e-10
    v <- (b_mat / Kv)^fi
    Ku <- K %*% v + 1e-10
    u <- (a_mat / Ku)^fi
  }
  
  P <- diag(as.vector(u)) %*% K %*% diag(as.vector(v))
  rownames(P) <- colnames(P) <- rownames(C_mat)
  return(P)
}

unbalanced_sinkhorn_with_prior <- function(a, b, C, Q, eps = 5, tau = 100) {
  # Incorporate the Prior
  K <- Q * exp(-as.matrix(C) / eps)
  
  # Standard Sinkhorn
  u <- rep(1, length(a))
  fi <- tau / (tau + eps)
  
  for (i in 1:100) {
    v <- (b / (as.vector(t(K) %*% u) + 1e-10))^fi
    u <- (a / (as.vector(K %*% v) + 1e-10))^fi
  }
  
  return(sweep(sweep(K, 1, u, "*"), 2, v, "*"))
}

compute_plan <- function(st_grid, C_st, prior = FALSE, Q_prior = NULL, eps = 1, tau = 1000){
  if(prior){
    P_est <- unbalanced_sinkhorn_with_prior(st_grid$S, st_grid$D, C_st, Q = Q_prior, eps = eps, tau = tau)
  } else {
    P_est <- unbalanced_sinkhorn(st_grid$S, st_grid$D, C_st, eps = eps, tau = tau)
  }
  return(P_est)
}
########################################################## true P
compute_true_P <- function(st_grid_d, trips_mapped_d){ # everything on a daily basis
  wd = trips_mapped_d$wd[1] # all should be the same day, hence the same wd
  n_nodes <- nrow(st_grid_d)
  node_lookup <- setNames(1:n_nodes, st_grid_d$Node_ID)
  true_P <- matrix(0, nrow = n_nodes, ncol = n_nodes)
  st_counts <- trips_mapped_d %>%
    mutate(
      from_node = paste0(start_neigh, "_T", T_Start, "_",wd),
      to_node = paste0(end_neigh, "_T", T_End, "_",wd)
    ) %>%
    count(from_node, to_node)
  #  Fill the matrix
  for(i in 1:nrow(st_counts)) {
    row_idx <- node_lookup[st_counts$from_node[i]]
    col_idx <- node_lookup[st_counts$to_node[i]]
    
    # Only fill if the nodes exist in our study grid
    if(!is.na(row_idx) & !is.na(col_idx)) {
      true_P[row_idx, col_idx] <- st_counts$n[i]
    }
  }
  # Assign names from your grid to the matrix
  rownames(true_P) <- colnames(true_P) <- st_grid_d$Node_ID
  return(true_P)
}
############################################################# Calibration

calibrate_st_ot_wd <- function(params, st_grid_d, C_spatial_informed, trips_mapped_d, prior = FALSE, Q_prior_wd = NULL) {
  # Set parameters
  eps_val <- params[1]
  tau_val <- params[2]
  gamma   <- params[3]
  wd = trips_mapped_d$wd[1] # all should be the same day, hence the same wd
  
  # Rebuild the cost matrix with given gamma
  C_calib = compute_cost(C_spatial_informed, st_grid_d, gamma = gamma)
  
  # Run model
  Q_prior = Q_prior_wd[[as.numeric(wd)+1]] 
  P_hat = compute_plan(st_grid_d, C_calib, prior = prior, Q_prior = Q_prior, eps = eps_val, tau = tau_val)
  
  # True P
  true_P = compute_true_P(st_grid_d, trips_mapped_d)
  
  # Aggregate and compare to true_P
  # Return RMSE and MAE and Common part of commuters (CPC)
  error <- c(sqrt(mean(((true_P - P_hat)^2))),  #RMSE
             mean(abs(true_P - P_hat)), #MAE
             2 * sum(pmin(true_P, P_hat))/ 
               (sum(true_P) + sum(P_hat)) #CPC
  )
  return(error)
}
choose_best_eps <- function(eps_min = 0, eps_max = 10000, st_grid_d, trips_mapped_d, C_spatial_informed, prior, Q_prior_wd){
  bys = c(1000,100,10,1,0.1)
  best_params = rep(NA,3)
  metrics = rep(NA,3)
  if(sum(st_grid_d$S!=0.1)>0){
    for(b in 1:length(bys)){
      bseq= bys[b]
      eps_grid = seq(eps_min, eps_max, bseq)
      search_grid <- expand_grid(eps   = eps_grid,  tau   = 1000,  gamma = 0)
      res = apply(search_grid, 1, calibrate_st_ot_wd, st_grid_d = st_grid_d, 
                  trips_mapped_d = trips_mapped_d, C_spatial_informed = C_spatial_informed,
                  prior = TRUE, Q_prior_wd = Q_prior_wd)
      best_params <- as.numeric(search_grid[which.max(res[3,]),])
      eps_min = best_params[1] - bseq
      eps_max = best_params[1] + bseq
    }
    metrics = res[,which.max(res[3,])]
  }
  out = list(best_params = best_params, metrics = metrics)
  return(out)
}

