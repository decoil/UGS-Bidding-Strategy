library(lubridate)
library(ggplot2)
library(readr)
library(tidyverse)
library(ompr)
library(ompr.roi)
library(ROI.plugin.glpk)
library(magrittr)
library(cowplot)

WGV_max <- 1000000             # Working gas volume
I_max <- 20000                 # Maximum injection rate
W_max <- 30000                 # Maximum withdrawal rate
Inv_thresh <- 0.5 * WGV_max    # Storage inventory threshold/half of WGV
I_r1 <- 1.0                    # Injection rate when storage <50%
I_r2 <- 0.7                    # Injection rate when storage >50%
W_int <- 0.4                   # Intercept of linear withdrawal function
W_slo <- 0.6                   # Slope of linear withdrawal function
c_var <- 0.012                 # Injection fee
Inv_0 <- 0                     # Starting inventory
Inv_N <- 0                     # Ending inventory
M <- WGV_max                   # Big M value
epsilon <- 0.001               # Small epsilon tolerance

fw_curve <- read.csv("fwcurve.csv") %>% 
  mutate(
    date = dmy(date), # Parse date to Date type
    price = as.numeric(price)
  ) %>%
  filter(date >= ymd("2026-04-01") & date <= ymd("2027-03-31")) %>% # Only include contract relevant dates
  arrange(date)

if(nrow(fw_curve)!= 365){
  stop("Incomplete forward curve sequence")
}

if(any(is.na(fw_curve$price))){
  stop("Missing price data")
}

N_days <- nrow(fw_curve)
P <- fw_curve$price

# Define Mixed Integer Linear Programming model
model <- MILPModel() %>%
  # Define variables
  add_variable(I[t], t = 1:N_days, type = "continuous", lb = 0) %>%
  add_variable(W[t], t = 1:N_days, type = "continuous", lb = 0) %>%
  add_variable(Inv[t], t = 0:N_days, type = "continuous", lb = 0, ub = WGV_max) %>%
  add_variable(b[t], t = 1:N_days, type = "binary") %>%

  # Define objective function
  set_objective(sum_over(W[t] * P[t], t=1:N_days) - sum_over(I[t] * P[t] * (1 + c_var), t = 1:N_days), "max") %>%
  
  # Define constraints
  add_constraint(Inv[0] == Inv_0) %>% 
  add_constraint(Inv[t] == Inv[t-1] + I[t] - W[t], t = 1:N_days) %>% 
  add_constraint(Inv[N_days] == Inv_N) %>%
  
  #Prevent simultaneous injection and withdrawal
  add_constraint(I[t] + W[t] <= max(I_max, W_max), t = 1:N_days) %>%
  
  # Injection rate parameters
  add_constraint(Inv[t-1] - Inv_thresh <= M * b[t] - epsilon, t = 1:N_days) %>%
  add_constraint(Inv[t-1] - Inv_thresh >= (-M) * (1 - b[t]), t = 1:N_days) %>%
  add_constraint(I[t] <= (I_max * I_r1 * (1 - b[t])) + (I_max * I_r2 * b[t]), t = 1:N_days) %>%
  
  # Withdrawal rate parameter
  add_constraint(W[t] <= (W_max * W_int) + (((W_max * W_slo) / WGV_max) * Inv[t-1]), t = 1:N_days)

# Solve model  
result <- tryCatch({
  solve_model(model, with_ROI(solver = "glpk", verbose = TRUE))
}, error = function(e) {
  message("Error solving model: ", e$message)
  NULL
})

# Process results and plots
if (!is.null(result) && result$status == "success"){
  optimal_profit <- objective_value(result)
  intrinsic_value <- optimal_profit / WGV_max # Intrinsic value per MWh
  
  message(sprintf("Solver Status: %s", result$status))
  message(sprintf("Total Intrinsic Value: \u20AC%.2f", optimal_profit))
  message(sprintf("Normalized Intrinsic Value: \u20AC%.2f per MWh WGV", intrinsic_value))
  
  # Optimal schedule values
  sol_I <- get_solution(result, I[t])
  sol_W <- get_solution(result, W[t])
  sol_Inv <- get_solution(result, Inv[t]) %>% slice(-1)
  
  # Instantiate table for optimal schedule
  optimal_schedule <- tibble(
    day_index = 1:N_days,
    date = fw_curve$date,
    price = fw_curve$price
    )
  
  # Add Injections, Withdrawals, and Inventory
  optimal_schedule %<>% mutate(Injection = sol_I$value) %>%
    mutate(Withdrawal = sol_W$value) %>%
    mutate(Inventory_End = sol_Inv$value)
  
  print(optimal_schedule, n=365) # Daily injection, withdrawal, and inventory schedule
  
  # Inventory Plot
  inventory_plot <- ggplot(optimal_schedule) +
    geom_line(aes(x = date, y = Inventory_End / WGV_max * 100), color = "red") +
    scale_y_continuous(name = "Inventory (% Full)") +
    labs(
      title = "Inventory Level (%)",
      x = "Date"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5)) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b %Y")

  # Price Plot
  price_plot <- ggplot(optimal_schedule) +
    geom_line(aes(x = date, y = price),color = "blue" ) +
    scale_y_continuous(name = "Price (\u20AC/MWh)") +
    labs(
      title = "Forward Price",
      x = "Date"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5)) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b %Y")

  plot_grid(inventory_plot, price_plot, ncol = 1, align = "v")
  
} else {
  message("Optimal solution not found or failed.")
  if (!is.null(result)) {
    message(sprintf("Solver Status: %s", result$status))
  }
}

# Monthly summary of optimal schedule
if (!is.null(result) && result$status == "success") {
  monthly_summary <- optimal_schedule %>%
    mutate(Month = floor_date(date, "month")) %>%
    group_by(Month) %>%
    summarise(
      Total_Injection = sum(Injection),
      Total_Withdrawal = sum(Withdrawal),
      Avg_Inventory = mean(Inventory_End),
      End_Inventory = last(Inventory_End),
     .groups = 'drop'
    )
  print("Monthly Summary of Optimal Schedule:")
  print(as_tibble(monthly_summary)) 
}

    
  
