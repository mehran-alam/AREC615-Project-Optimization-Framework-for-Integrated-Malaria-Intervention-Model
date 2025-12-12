README: Malaria Intervention Optimization Model (5-Year)

## Overview

This repository contains an R-based simulation and optimization model that evaluates the effectiveness of malaria interventions across multiple regions over a 5‑year horizon. The model integrates:

- **Epidemiological dynamics** (SIR-like malaria model)
- **Behavioral response dynamics**
- **Intervention effects** (LLIN, IRS, IPT, ACT, vaccine, campaign)
- **Multi-year disease trajectories**
- **Integer Linear Programming (ILP)** optimization
- **Greedy cost-effectiveness optimization**

The goal is to identify cost-effective multi-year intervention strategies that minimize malaria burden measured as **avoided person-days of infection**, subject to a budget constraint.

---

## Key Features

### 1. Epidemiological Model
Implements an ODE system using `deSolve` with climate-dependent mosquito abundance and malaria transmission parameters.

### 2. Behavioral Component
Models community behavioral effort (e.g., prevention effort), influenced by infection prevalence and campaigns.

### 3. Intervention Effects
Applies intervention-specific reductions to transmission parameters, recovery rates, and behavior.

Supports:
- LLIN
- IRS
- IPT
- ACT
- Vaccine
- Community campaign
- Combinations (LLIN+ACT, LLIN+campaign, ACT+campaign)

### 4. Multi-Year Simulation
Simulates 5 years (365×5 days) for:
- Dry regions
- Moderate regions
- Wet regions

### 5. Optimization
Includes two solvers:

#### **ILP Solver (`lpSolve`)**
Maximizes avoided person-days under a total budget constraint.

#### **Greedy Cost-Effectiveness Solver**
Prioritizes interventions with the highest avoided-person-days-per-dollar ratio, selecting at most one per district.

---

## Repository Structure

```
├── main_script.R           # Full simulation + optimization code
├── README.md               # Documentation (this file)
└── /plots                  # (optional) Output plots
```

---

## How to Run

### **1. Install Required R Packages**
```r
install.packages(c("deSolve", "ggplot2", "dplyr", "tidyr", "lpSolve"))
```

### **2. Set Budget and Execute Optimization**
```r
budget <- 1e7
run_complete_optimization(budget)
```

### **3. View Results**
The script prints:
- Total avoided person-days
- Number of selected actions
- Total cost
- Remaining budget
- Cost-effectiveness (person-days per $1,000)

---

## Customization

We can modify:

- **Intervention cost & effects** in `params_interventions`
- **Coverage levels**
- **Behavioral parameters**
- **Simulation length** (`SIM_YEARS`)
- **Regional populations**

---

## Notes

- The model includes progress bars for transition computations.
- All ODE failures are handled gracefully.
- Baseline (no-intervention) trajectories are precomputed for each climate.

## Block-by-block Explanations
# ============================================================================
# 1. PARAMETERS
# ============================================================================

# Epidemiological parameters
params_epidemiology <- list(
  a_u = 0.25,      # Biting rate (per mosquito per day)
  b_u = 0.022,     # Probability of infection per bite
  c = 0.36,        # Probability mosquito gets infected from human
  delta = 4.7895e-5, # Disease-induced mortality rate
  gamma_u = 1/180, # Natural recovery rate (1/180 days ≈ 6 months)
  m_u_dry = 5,     # Mosquito density in dry climate (mosquitoes/human)
  m_u_mod = 20,    # Mosquito density in moderate climate
  m_u_wet = 35,    # Mosquito density in wet climate
  mu = 0.095,      # Mosquito mortality rate
  omega = 274,     # Mosquito extrinsic incubation period
  tau = 10,        # Mosquito biting cycle length
  rho_u = 1/274,   # Loss of immunity rate
  N = 10000        # Population per district
)
# Behavioral parameters
params_behavior <- list(
  alpha_S = 0.02,  # Social learning rate from prevalence
  alpha_O = 0.01,  # Direct observation rate
  alpha_C = 0.03,  # Campaign influence rate
  delta_B = 0.005, # Belief decay rate
  theta_B = 0.6,   # Belief weight in effort calculation
  theta_I = 0.4,   # Prevalence weight in effort calculation
  theta_C = 0.3,   # Cost of effort coefficient
  c_effort = 0.5   # Baseline effort cost
)
# Intervention parameters
params_interventions <- list(
  cost = c(none = 0, LLIN = 1.33, IRS = 2.22, IPT = 1.13, 
           ACT = 4.82, vaccine = 20.66, campaign = 1.50),
  effects = list(
    LLIN = list(a_reduction = 0.80, m_reduction = 0.64, theta_i = 0.2),
    IRS = list(m_reduction = 0.95, theta_i = 0.15),
    IPT = list(b_reduction = 0.786, theta_i = 0.1),
    ACT = list(gamma_multiplier = 18, theta_i = 0.15),
    vaccine = list(b_reduction = 0.773, gamma_multiplier = 3.27, theta_i = 0.1),
    campaign = list(alpha_C_boost = 0.4, theta_i_boost = 0.3, compliance_boost = 0.25)
  ),
  coverage_levels = c(0.2, 0.4, 0.6, 0.8)
)
Explanation: Cost and effectiveness parameters for each intervention. Coverage levels define implementation scales. Campaigns uniquely boost behavioral parameters.
# Region definitions
regions <- expand.grid(
  climate = c("dry", "moderate", "wet"),
  dist_cost = c("low", "medium", "high"),
  stringsAsFactors = FALSE
)

# Initial conditions
initial_conditions <- list(
  dry = c(S = 0.60, I = 0.15, R = 0.25, B = 0.3, E = 0.2),
  moderate = c(S = 0.15, I = 0.15, R = 0.70, B = 0.5, E = 0.3),
  wet = c(S = 0.10, I = 0.15, R = 0.75, B = 0.4, E = 0.25)
)

SIM_YEARS <- 5
SIM_DAYS <- 365
Explanation: Defines 9 regions (3 climates × 3 cost levels) with disease state initial conditions. SIM_YEARS sets the 5-year horizon.

# Core Dynamics
calculate_force_of_infection <- function(m, a, b, I_total, E, params) {
  numerator <- m * a^2 * b * params$c * exp(-params$mu * params$tau) * I_total
  denominator <- params$mu + a * params$c * I_total
  lambda <- numerator / denominator
  return(lambda * (1 - min(E, 1)))
}
Explanation: Calculates the force of infection λ using Ross-Macdonald formula with behavioral effort (E) reducing transmission.
combined_odes <- function(t, state, params) {
  S <- state["S"]; I <- state["I"]; R <- state["R"]; B <- state["B"]
  
  # Extract parameters
  m <- params$m; a <- params$a; b <- params$b; c <- params$c
  gamma <- params$gamma; delta <- params$delta; rho <- params$rho
  mu <- params$mu; tau <- params$tau
  alpha_S <- params$alpha_S; alpha_O <- params$alpha_O
  alpha_C <- params$alpha_C; delta_B <- params$delta_B
  theta_B <- params$theta_B; theta_I <- params$theta_I
  theta_C <- params$theta_C; c_effort <- params$c_effort
  theta_i_boost <- params$theta_i_boost %||% 0
  C_intensity <- params$C_intensity %||% 0
  alpha_C_boost <- params$alpha_C_boost %||% 0
  
  # Calculate effort (E)
  theta_total <- theta_B * B + theta_I * I + theta_i_boost
  E <- min(1, max(0, theta_total - theta_C * c_effort))
  
  # Force of infection
  lambda <- calculate_force_of_infection(m, a, b, I, E, params)
  
  # Epidemiological equations
  dS <- -lambda * S + rho * R
  dI <- lambda * S - gamma * I - delta * I
  dR <- gamma * I - rho * R
  
  # Behavioral equation
  dB <- alpha_S * (1 - I) * I * (1 - B) +
    alpha_O * I * (1 - B) +
    alpha_C * (1 + alpha_C_boost) * C_intensity * (1 - B) -
    delta_B * B
  
  return(list(c(dS = dS, dI = dI, dR = dR, dB = dB), E = E))
}
Explanation: Main ODE system combining epidemiological (SIR) and behavioral (B) dynamics. Effort E mediates between beliefs and actual prevention behavior.
get_intervention_effects <- function(intervention_type, coverage, m_u, params_epi) {
  effects <- list(
    m = m_u, a = params_epi$a_u, b = params_epi$b_u,
    gamma = params_epi$gamma_u, rho = params_epi$rho_u,
    delta = params_epi$delta, mu = params_epi$mu,
    tau = params_epi$tau, c = params_epi$c,
    theta_i_boost = 0, alpha_C_boost = 0, C_intensity = 0
  )
  
  if (is.null(intervention_type) || intervention_type == "none" || intervention_type == "") {
    return(effects)
  }
  
  interventions <- strsplit(intervention_type, "_")[[1]]
  
  for (int in interventions) {
    if (int %in% c("0.2", "0.4", "0.6", "0.8")) next
    
    if (int == "LLIN" && int %in% names(params_interventions$effects)) {
      int_effects <- params_interventions$effects[[int]]
      effects$a <- effects$a * (1 - int_effects$a_reduction * coverage)
      effects$m <- effects$m * (1 - int_effects$m_reduction * coverage)
      effects$theta_i_boost <- effects$theta_i_boost + int_effects$theta_i * coverage
    }
    
    # Similar logic for other interventions...
  }
  
  return(effects)
}
Explanation: Applies intervention effects multiplicatively based on coverage. Allows combination interventions (e.g., "LLIN_ACT").
## Simulation Function
simulate_year <- function(initial_state, intervention_type, coverage,
                          region_type, days,
                          params_epi = params_epidemiology,
                          params_behav = params_behavior) {
  
  m_u <- switch(region_type,
                "dry" = params_epi$m_u_dry,
                "moderate" = params_epi$m_u_mod,
                "wet" = params_epi$m_u_wet,
                params_epi$m_u_mod)
  
  int_effects <- get_intervention_effects(intervention_type, coverage, m_u, params_epi)
  combined_params <- c(int_effects, params_behav)
  
  y0 <- c(S = as.numeric(initial_state["S"]),
          I = as.numeric(initial_state["I"]),
          R = as.numeric(initial_state["R"]),
          B = as.numeric(initial_state["B"]))
  
  times <- seq(0, days, by = 1)
  
  tryCatch({
    solution <- ode(y = y0, times = times, func = combined_odes,
                    parms = combined_params, method = "lsoda")
    
    final_idx <- nrow(solution)
    final_state <- solution[final_idx, c("S", "I", "R", "B")]
    
    # Calculate person-days of infection
    I_trajectory <- solution[, "I"]
    person_days <- sum(I_trajectory[-1] + I_trajectory[-length(I_trajectory)]) / 2 * 
      params_epi$N
    
    return(list(
      final_state = final_state,
      person_days = person_days,
      trajectories = as.data.frame(solution),
      success = TRUE,
      intervention = intervention_type,
      coverage = coverage,
      region = region_type
    ))
    
  }, error = function(e) {
    warning(paste("ODE solver failed for", intervention_type, ":", e$message))
    return(list(success = FALSE))
  })
}
Explanation: Runs 5-year simulation for specific intervention and region. Uses trapezoidal rule to integrate infection person-days.
## Optimization
generate_actions <- function() {
  interventions <- c("none", "LLIN", "IRS", "IPT", "ACT", "vaccine", "campaign")
  actions <- list()
  
  actions[["none_0"]] <- list(type = "none", coverage = 0, cost = 0)
  
  for (int in setdiff(interventions, "none")) {
    for (cov in params_interventions$coverage_levels) {
      actions[[paste0(int, "_", cov)]] <- list(
        type = int, coverage = cov,
        cost = params_interventions$cost[int] * cov
      )
    }
  }
  
  # Define key combinations
  key_combinations <- list(
    c("LLIN", "ACT"),
    c("LLIN", "campaign"),
    c("ACT", "campaign")
  )
  
  for (combo in key_combinations) {
    combo_name <- paste(combo, collapse = "_")
    for (cov in params_interventions$coverage_levels) {
      cost <- sum(sapply(combo, function(x) params_interventions$cost[x] * cov))
      actions[[paste(combo_name, cov, sep = "_")]] <- list(
        type = combo_name, coverage = cov, cost = cost
      )
    }
  }
  
  return(actions)
}
Explanation: Generates action space including single interventions and key combinations at different coverage levels.
precompute_transitions <- function(regions, initial_conditions, actions,
                                   params_epi = params_epidemiology) {
  transitions <- list()
  costs <- list()
  
  # Compute baseline (no intervention) for each climate
  baseline_person_days <- list()
  for (climate in unique(regions$climate)) {
    res_none <- simulate_year(
      initial_state = initial_conditions[[climate]],
      intervention_type = "none", coverage = 0,
      region_type = climate, days = SIM_DAYS * SIM_YEARS, params_epi = params_epi
    )
    baseline_person_days[[climate]] <- res_none$person_days
  }
  
  # Precompute all state transitions
  for (i in 1:nrow(regions)) {
    region <- regions[i, ]
    climate <- region$climate
    dist_cost <- region$dist_cost
    
    cost_multiplier <- switch(dist_cost, "low" = 0.8, "medium" = 1.0, "high" = 1.2)
    
    for (action_name in names(actions)) {
      action <- actions[[action_name]]
      
      result <- simulate_year(
        initial_state = initial_conditions[[climate]],
        intervention_type = action$type,
        coverage = action$coverage,
        region_type = climate, days = SIM_DAYS * SIM_YEARS,
        params_epi = params_epi
      )
      
      trans_key <- paste(climate, dist_cost, action_name, sep = "|")
      baseline_pd <- baseline_person_days[[climate]]
      avoided_pd <- max(0, baseline_pd - result$person_days)
      
      transitions[[trans_key]] <- list(
        initial = initial_conditions[[climate]],
        final = result$final_state,
        person_days = result$person_days,
        avoided_person_days = avoided_pd,
        success = result$success
      )
      
      costs[[trans_key]] <- (action$cost %||% 0) * params_epi$N * cost_multiplier * SIM_YEARS
    }
  }
  
  return(list(transitions = transitions, costs = costs))
}
Explanation: Precomputes all possible state transitions for optimization. Avoided person-days calculated relative to baseline.
solve_ilp_optimization <- function(transitions_data, budget, verbose = TRUE) {
  all_trans <- transitions_data$transitions
  all_costs <- transitions_data$costs
  keys <- names(all_trans)
  
  objective_fn <- sapply(keys, function(k) all_trans[[k]]$avoided_person_days %||% 0)
  cost_vector <- sapply(keys, function(k) all_costs[[k]] %||% 0)
  
  # Filter valid interventions (positive avoided person-days)
  valid_indices <- which(objective_fn > 0)
  objective_fn <- objective_fn[valid_indices]
  cost_vector <- cost_vector[valid_indices]
  keys_valid <- keys[valid_indices]
  
  # Setup ILP
  constraint_matrix <- matrix(cost_vector, nrow = 1)
  constraint_direction <- "<="
  constraint_rhs <- budget
  
  lp_model <- lpSolve::lp(
    direction = "max",
    objective.in = objective_fn,
    const.mat = constraint_matrix,
    const.dir = constraint_direction,
    const.rhs = constraint_rhs,
    all.bin = TRUE
  )
  
  solution_vector <- lp_model$solution
  selected_indices <- which(solution_vector == 1)
  selected_keys <- keys_valid[selected_indices]
  
  selected_cost <- sum(cost_vector[selected_indices])
  selected_apd <- sum(objective_fn[selected_indices])
  
  return(list(
    objective = selected_apd,
    total_cost = selected_cost,
    num_selected = length(selected_keys),
    solution = setNames(as.list(solution_vector[selected_indices]), selected_keys)
  ))
}
solve_greedy_ce_optimization <- function(transitions_data, budget, verbose = TRUE) {
  all_trans <- transitions_data$transitions
  all_costs <- transitions_data$costs
  keys <- names(all_trans)
  
  # Calculate cost-effectiveness ratios
  ce_ratios <- sapply(keys, function(key) {
    cost <- as.numeric(all_costs[[key]] %||% NA)
    avoided <- as.numeric(all_trans[[key]]$avoided_person_days %||% NA)
    
    if (!is.na(cost) && !is.na(avoided)) {
      if (cost > 0 && avoided > 0) {
        return(avoided / cost)
      } else if (cost == 0 && avoided > 0) {
        return(Inf)
      }
    }
    return(-Inf)
  }, USE.NAMES = TRUE)
  
  # Sort by cost-effectiveness
  sorted_keys <- names(sort(ce_ratios[is.finite(ce_ratios)], decreasing = TRUE))
  
  remaining_budget <- budget
  solution <- list()
  total_avoided <- 0
  total_cost <- 0
  
  # Greedy selection (one intervention per region)
  selected_regions <- character(0)
  
  for (key in sorted_keys) {
    cost <- as.numeric(all_costs[[key]])
    avoided <- as.numeric(all_trans[[key]]$avoided_person_days)
    region <- extract_region(key)
    
    if (!is.na(cost) && !is.na(avoided) &&
        cost <= remaining_budget && avoided > 0 &&
        !(region %in% selected_regions)) {
      
      solution[[key]] <- 1
      selected_regions <- c(selected_regions, region)
      remaining_budget <- remaining_budget - cost
      total_avoided <- total_avoided + avoided
      total_cost <- total_cost + cost
    }
  }
  
  return(list(
    solution = solution,
    objective = total_avoided,
    total_cost = total_cost,
    remaining_budget = remaining_budget,
    num_selected = length(solution),
    selected_regions = selected_regions
  ))
}
Explanation: Greedy algorithm selecting most cost-effective interventions, limited to one per region.
## Main Execution
run_complete_optimization <- function(budget) {
  cat("\nCOMPLETE MALARIA INTERVENTION OPTIMIZATION (5 YEARS)\n")
  cat("═══════════════════════════════════════════════════════════════════\n\n")
  
  cat("Step 1: Generating action space...\n")
  actions <- generate_actions()
  
  cat("\nStep 2: Computing state transitions...\n")
  transitions_data <- precompute_transitions(regions, initial_conditions, actions)
  
  cat("\nStep 3: Running optimization algorithms...\n")
  
  cat("\n--- METHOD 1: Integer Linear Program (True Optimum) ---\n")
  result_ilp <- solve_ilp_optimization(transitions_data, budget, verbose = TRUE)
  
  cat("\n--- METHOD 2: Cost-Effectiveness Focused (Greedy Approximation) ---\n")
  result_ce <- solve_greedy_ce_optimization(transitions_data, budget, verbose = TRUE)
  
  cat("\nStep 4: Comparing optimization methods...\n")
  
  # Print comparison table
  ilp_pd_per_k <- (result_ilp$objective / max(1, result_ilp$total_cost)) * 1000
  ce_pd_per_k <- (result_ce$objective / max(1, result_ce$total_cost)) * 1000
  
  cat(paste(rep("─", 90), collapse = ""), "\n")
  cat(sprintf("%-25s %-10s %-15s %-20s %-10s %-10s\n",
              "Method", "Actions", "Total Cost", "Avoided PD (5Y)", "Used %", "PD/$1k"))
  cat(paste(rep("─", 90), collapse = ""), "\n")
  
  cat(sprintf("%-25s %-10d $%-14s %-20s %-10.1f %-10.0f\n",
              "ILP (True Optimum)",
              result_ilp$num_selected,
              format(round(result_ilp$total_cost), big.mark = ","),
              format(round(result_ilp$objective), big.mark = ","),
              (result_ilp$total_cost / budget) * 100,
              ilp_pd_per_k))
  
  cat(sprintf("%-25s %-10d $%-14s %-20s %-10.1f %-10.0f\n",
              "CE-Focused (Greedy)",
              result_ce$num_selected,
              format(round(result_ce$total_cost), big.mark = ","),
              format(round(result_ce$objective), big.mark = ","),
              (result_ce$total_cost / budget) * 100,
              ce_pd_per_k))
  
  return(list(
    actions = actions,
    transitions_data = transitions_data,
    ilp_optimum = result_ilp,
    ce_focused = result_ce,
    timestamp = Sys.time()
  ))
}
