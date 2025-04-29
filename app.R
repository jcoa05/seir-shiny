# Load and install packages if needed
#required_pck <- c("shiny", "deSolve", "ggplot2", "tidyr", "dplyr", "scales", "shinyBS", "bslib", "httpuv")
#load_or_install <- function(packages) {
#  for (pkg in packages) {
#    if (!require(pkg, character.only = TRUE)) {
#      install.packages(pkg)
#      library(pkg, character.only = TRUE)
#    }
#  }
#}
#load_or_install(required_pck)

#Packages
library(shiny)
library(deSolve)
library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)
library(shinyBS)
library(bslib)
library(httpuv)  

# Equations for the deterministic SEI2R model with Intervention
# S = Susceptible = Healthy individuals at risk of infection.
# E = Exposed = Infected but not yet infectious (latent period).
# I1 = Infectious, Untreated = Infectious and not yet treated.
# I2 = Infectious, Treated = Treated individuals with reduced transmission.
# R = Removed = Recovered or deceased (no longer infectious).

seir_intervention_model <- function (t, x, params) {
  S  <- x[1]
  E  <- x[2]
  I1 <- x[3]
  I2 <- x[4]
  R  <- x[5]
  
  N_total <- S + E + I1 + I2 + R
  if (N_total == 0) N_total <- 1e-9 # Avoid division by zero
  
  beta <- params["beta"]    # Transmission rate
  omega <- params["omega"] # Rate of progression from E to I (1/latency_period)
  gamma <- params["gamma"] # Rate of recovery (1/infectious_period)
  theta <- params["theta"] # Transmission reduction factor (0 to 1)
  phi <- params["phi"]      # Rate of treatment initiation (1/treatment_delay)
  
  # ODEs
  dS <- -beta * S * (I1 + (1 - theta) * I2) / N_total  
  dE <- beta * S * (I1 + (1 - theta) * I2) / N_total - omega * E  

  # Handle special cases, with theta = 0 checked first
  if (theta == 0) {
    # Case where theta = 0 (no transition to I2)
    dI1 <- omega * E - gamma * I1 
    dI2 <- 0  # No change in I2 (assuming it's not part of the process)
    dR <- gamma * I1 
  } else if (is.infinite(phi)) { 
    # Special case for zero delay (phi = Inf)
    dI1 <- -gamma * I1  # Only recovery leaves I1 (should remain 0 if I1_0=0)
    dI2 <- omega * E - gamma * I2  # Incoming from E, recovery leaves I2
    dR <- gamma * I1 + gamma * I2  # Recovery from both I1 and I2
  } else {
    # Normal case: newly infectious go to I1, then transition to I2 at rate phi
    dI1 <- omega * E - phi * I1 - gamma * I1
    dI2 <- phi * I1 - (1/(1/gamma - 1/phi)) * I2  # Adjusted recovery rate for treated
    dR <- gamma * I1 + (1/(1/gamma - 1/phi)) * I2  # Recovery from both compartments
  }

  dx <- c(dS, dE, dI1, dI2, dR)
  list(dx)
}

# Define the User Interface (UI)
ui <- fluidPage(
  theme = bs_theme(version = 5),
  
  titlePanel(tags$b("Treatment effectiveness and delay in an influenza epidemic simulation")),
  p("An antiviral treatment reduces the contagiousness of influenza cases. This model simulates its impact on the dynamics of an influenza epidemic."),
  
  sidebarLayout(
    sidebarPanel(
      width = 4,
      
      wellPanel(
        h5("Scenarios: Delay at individual level"),
        fluidRow(
          column(3, actionButton("scenario_mild", "No treatment", icon = icon("exclamation"))),
          column(3, actionButton("scenario_moderate", "Optimistic", icon = icon("check"))),
          column(3, actionButton("scenario_severe", "Pessimistic", icon = icon("skull")))
        )
      ),
      
      wellPanel(
        h5(strong("Treatment Parameters")),
        h6(strong("Model A")),
        sliderInput("theta_A", "Treatment effectiveness:",
                    min = 0, max = 100, value = 50, step = 5, post = "%"),
        helpText("Reduction in transmission for treated individuals.", style = "font-size: small;"),
        sliderInput("treatment_delay_A", "Average treatment delay (days):",
                    min = 0, max = 8, value = 2, step = 1),
        helpText("Average time between becoming infectious and starting treatment.", style = "font-size: small;"),
        
        hr(),
        
        h6(strong("Model B")),
        sliderInput("theta_B", "Treatment effectiveness:",
                    min = 0, max = 100, value = 50, step = 5, post = "%"),
        helpText("Reduction in transmission for treated individuals.", style = "font-size: small;"),
        sliderInput("treatment_delay_B", "Average treatment delay (days):",
                    min = 0, max = 8, value = 2, step = 1),
        helpText("Average time between becoming infectious and starting treatment.", style = "font-size: small;")
      ),
      
      accordion(
        accordion_panel(
          "Epidemic Parameters",
          h6(strong("Model A")),
          numericInput("R0_input_A", "Basic reproduction number (R0):",
                       min = 0, max = 10, value = 3.5, step = 0.1),
          helpText("Average number of secondary infections from one infectious individual in a fully susceptible population.", style = "font-size: small;"),
          
          numericInput("N_input_A", "Population size (Millions):",
                       min = 0.1, max = 1000, value = 16, step = 1),
          helpText("Total population size for Model A.", style = "font-size: small;"),
          
          hr(),
          
          h6(strong("Model B")),
          numericInput("R0_input_B", "Basic reproduction number (R0):",
                       min = 0, max = 10, value = 3.2, step = 0.1),
          helpText("Average number of secondary infections from one infectious individual in a fully susceptible population.", style = "font-size: small;"),
          
          numericInput("N_input_B", "Population size (Millions):",
                       min = 0.1, max = 1000, value = 68, step = 1),
          helpText("Total population size for Model B.", style = "font-size: small;")
        ),
        
        accordion_panel(
          "Disease Progression",
          h6(strong("Model A")),
          sliderInput("latency_period_A", "Average latency period (days):",
                      min = 1, max = 10, value = 2, step = 1),
          helpText("Average time from infection to becoming infectious (E -> I1/I2).", style = "font-size: small;"),
          
          sliderInput("infectious_period_A", "Average infectious period (days):",
                      min = 1, max = 14, value = 8, step = 1),
          helpText("Average duration an individual remains infectious (untreated or treated).", style = "font-size: small;"),
          
          hr(),
          
          h6(strong("Model B")),
          sliderInput("latency_period_B", "Average latency period (days):",
                      min = 1, max = 10, value = 2, step = 1),
          helpText("Average time from infection to becoming infectious (E -> I1/I2).", style = "font-size: small;"),
          
          sliderInput("infectious_period_B", "Average infectious period (days):",
                      min = 1, max = 14, value = 6, step = 1),
          helpText("Average duration an individual remains infectious (untreated or treated).", style = "font-size: small;")
        ),
        
        accordion_panel(
          "Simulation Settings",
          h6(strong("Common Settings")),
          sliderInput("simtime_input_common", "Duration of simulation (days):",
                      min = 1, max = 500, value = 150, step = 1),
          helpText("Total time period for the simulation.", style = "font-size: small;"),
          
          sliderInput("E0_input_common", "Initial exposed individuals (E0):",
                      min = 1, max = 100000, value = 500, step = 1),
          helpText("Number of individuals initially exposed (infected but not yet infectious) at day 0.", style = "font-size: small;")
        ),
        
        accordion_panel(
          "Visualization Options",
          checkboxInput("use_log_scale", "Use logarithmic scale for infections plot", FALSE),
          checkboxInput("normalize_by_population", "Normalize infections by population size", FALSE)
        ),
        
        accordion_panel(
          "Key Assumptions",
          h6("1. Population Structure"),
          p("• Homogeneous mixing: Uniform random interactions among all individuals.", style = "font-size: small;"),
          p("• Closed population: No births, deaths, or migration.", style = "font-size: small;"),
          
          h6("2. Disease Progression"),
          p("• Latent phase: No infectiousness during incubation.", style = "font-size: small;"),
          p("• Constant infectiousness: Stable transmission rate during infection unless treated.", style = "font-size: small;"),
          p("• Treatment: Full coverage, instantaneous efficacy (no delays or failures).", style = "font-size: small;"),
          
          h6("3. Recovery & Immunity"),
          p("• Exponential recovery: Constant per-capita recovery rate (biologically simplistic; may misestimate timing).", style = "font-size: small;"),
          p("• Lifelong immunity: No reinfection or waning immunity.", style = "font-size: small;"),
          
          h6("4. Model Framework"),
          p("• ODE structure: Uses ordinary differential equations (no fixed delays or DDEs).", style = "font-size: small;"),
          p("• Deterministic dynamics: Ignores stochasticity (large population assumed).", style = "font-size: small;"),
          
          h6("5. Transmission Dynamics"),
          p("• Duration-independent spread: Transmission rate unaffected by time since infection.", style = "font-size: small;"),
          p("• Static behavior: No behavioral adaptations, seasonality, or co-infections.", style = "font-size: small;"),
          
          h6("6. Data & External Factors"),
          p("• Unbiased data: Epidemiological data are unaffected by external shocks or artifacts.", style = "font-size: small;"),
          p("• Pathogen isolation: No interactions with other diseases or environmental factors.", style = "font-size: small;")
        )
      )
    ),
    
    mainPanel(
      width = 8,
      tabsetPanel(
        tabPanel(strong("Epidemic curve: new infections"),
                 plotOutput("newInfectionsPlot"),
                 downloadButton("downloadNewInfectionsPlot", "Download Plot")
        ),
        tabPanel(strong("Compartment sizes over time"),
                 plotOutput("seirCurvesPlot"),
                 downloadButton("downloadSeirCurvesPlot", "Download Plot") 
        ),
        tabPanel(strong("Summary statistics"),
                 verbatimTextOutput("summaryStats")
                 #downloadButton("downloadSummaryStats", "Download Summary")
        )
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  # --- Helper function to run a single simulation ---
  run_single_simulation <- function(N_model, E0_common, simtime_common,
                                    R0, latency_period, infectious_period, theta, treatment_delay,
                                    model_name) {
    N <- N_model
    E0 <- E0_common
    simtime <- simtime_common
    I1_0 <- 0 # Fixed initial untreated infected
    I2_0 <- 0 # Fixed initial treated infected
    R_init <- 0 # Fixed initial removed
    
    # Basic Input validation with renamed parameters
    if(is.null(simtime) || is.null(latency_period) || is.null(infectious_period) ||
       is.null(N) || is.null(E0) || is.null(R0) || is.null(theta) || is.null(treatment_delay) ||
       simtime <= 0 || latency_period <= 0 || infectious_period <= 0 || N <= 0 ||
       E0 < 0 || R0 < 0 || theta < 0 || theta > 1 || treatment_delay < 0) {
      warning(paste("Invalid input parameters detected for", model_name))
      return(NULL)
    }
    
    S0 <- max(0, N - E0 - I1_0 - I2_0 - R_init) # Ensure S0 is not negative
    times <- seq(0, simtime, by = 0.25) # Simulation time grid
    xstart <- c(S = S0, E = E0, I1 = I1_0, I2 = I2_0, R = R_init)
    
    # Calculate model parameters from inputs with consistent naming
    omega <- 1/latency_period # Rate of progression from E to I
    gamma <- 1/infectious_period # Rate of recovery
    beta <- R0 * gamma # Transmission rate
    
    # Calculate phi (treatment rate) from treatment_delay
    # Handle special case: if delay is 0, set phi to Inf (immediate treatment)
    phi <- if (treatment_delay == 0) Inf else 1/treatment_delay
    
    # Bundle parameters for the ODE solver with consistent naming
    # Theta is already scaled 0-1 here
    params <- c(beta = beta, omega = omega, gamma = gamma, theta = theta, phi = phi)
    
    # Run simulation with error handling
    out_list <- tryCatch({
      as.data.frame(lsoda(xstart, times, seir_intervention_model, params))
    }, error = function(e) {
      warning(paste("ODE solver failed for", model_name, ":", e$message))
      return(NULL)
    })
    
    # Add Model identifier and return
    if (!is.null(out_list)) {
      out_list$Model <- model_name
    }
    return(out_list)
  }
  
  observeEvent(input$scenario_mild, {
    updateSliderInput(session, "R0_input_A", value = 2.5)
    updateSliderInput(session, "R0_input_B", value = 2.2)
    updateSliderInput(session, "theta_A", value = 0)
    updateSliderInput(session, "theta_B", value = 0)
    updateSliderInput(session, "treatment_delay_A", value = 1)
    updateSliderInput(session, "treatment_delay_B", value = 1)
    updateSliderInput(session, "simtime_input_common", value = 365)
  })
  
  observeEvent(input$scenario_moderate, {
    updateSliderInput(session, "R0_input_A", value = 2.5)
    updateSliderInput(session, "R0_input_B", value = 2.2)
    updateSliderInput(session, "theta_A", value = 50)
    updateSliderInput(session, "theta_B", value = 50)
    updateSliderInput(session, "treatment_delay_A", value = 3)
    updateSliderInput(session, "treatment_delay_B", value = 3)
    updateSliderInput(session, "simtime_input_common", value = 365)
  })
  
  observeEvent(input$scenario_severe, {
    updateSliderInput(session, "R0_input_A", value = 2.5)
    updateSliderInput(session, "R0_input_B", value = 2.2)
    updateSliderInput(session, "theta_A", value = 30)
    updateSliderInput(session, "theta_B", value = 30)
    updateSliderInput(session, "treatment_delay_A", value = 5)
    updateSliderInput(session, "treatment_delay_B", value = 5)
    updateSliderInput(session, "simtime_input_common", value = 365)
  })
  
  observeEvent(input$scenario_mild2, {
    updateSliderInput(session, "R0_input_A", value = 2.5)
    updateSliderInput(session, "R0_input_B", value = 2.2)
    updateSliderInput(session, "theta_A", value = 0)
    updateSliderInput(session, "theta_B", value = 0)
    updateSliderInput(session, "treatment_delay_A", value = 1)
    updateSliderInput(session, "treatment_delay_B", value = 1)
    updateSliderInput(session, "E0_common", value = 30000)
    updateSliderInput(session, "simtime_input_common", value = 365)
  })
  
  observeEvent(input$scenario_moderate2, {
    updateSliderInput(session, "R0_input_A", value = 2.5)
    updateSliderInput(session, "R0_input_B", value = 2.2)
    updateSliderInput(session, "theta_A", value = 50)
    updateSliderInput(session, "theta_B", value = 50)
    updateSliderInput(session, "treatment_delay_A", value = 3)
    updateSliderInput(session, "treatment_delay_B", value = 3)
    updateSliderInput(session, "E0_common", value = 10000)
    updateSliderInput(session, "simtime_input_common", value = 365)
  })
  
  observeEvent(input$scenario_severe2, {
    updateSliderInput(session, "R0_input_A", value = 2.5)
    updateSliderInput(session, "R0_input_B", value = 2.2)
    updateSliderInput(session, "theta_A", value = 30)
    updateSliderInput(session, "theta_B", value = 30)
    updateSliderInput(session, "treatment_delay_A", value = 5)
    updateSliderInput(session, "treatment_delay_B", value = 5)
    updateSliderInput(session, "E0_common", value = 50000)
    updateSliderInput(session, "simtime_input_common", value = 365)
  })
  
  # --- Reactive expression to run simulation for Model A with updated parameter names ---
  simulation_results_A <- reactive({
    req(input$N_input_A, input$theta_A) # Ensure inputs are available
    run_single_simulation(
      N_model = input$N_input_A * 1e6,
      E0_common = input$E0_input_common,
      simtime_common = input$simtime_input_common,
      R0 = input$R0_input_A,
      latency_period = input$latency_period_A,
      infectious_period = input$infectious_period_A,
      theta = input$theta_A / 100,
      treatment_delay = input$treatment_delay_A,
      model_name = "Model A"
    )
  })
  
  # --- Reactive expression to run simulation for Model B with updated parameter names ---
  simulation_results_B <- reactive({
    req(input$N_input_B, input$theta_B) # Ensure inputs are available
    run_single_simulation(
      N_model = input$N_input_B * 1e6,
      E0_common = input$E0_input_common,
      simtime_common = input$simtime_input_common,
      R0 = input$R0_input_B,
      latency_period = input$latency_period_B,
      infectious_period = input$infectious_period_B,
      theta = input$theta_B / 100, # Convert percentage back to 0-1 scale
      treatment_delay = input$treatment_delay_B,
      model_name = "Model B"
    )
  })
  
  # --- Reactive expression to process and combine results for daily data ---
  daily_data_combined <- reactive({
    out_A <- simulation_results_A()
    out_B <- simulation_results_B()
    
    # Helper function to process raw simulation output to daily data
    process_to_daily <- function(out, model_name) {
      if (is.null(out) || nrow(out) < 2) return(NULL) # Need at least two time points
      
      # Filter for integer times (daily snapshots) and time 0
      # Keep time 0 to calculate first day's infections correctly
      out.daily_all <- out[out$time %% 1 == 0 | out$time == 0, ]
      out.daily_all <- distinct(out.daily_all, time, .keep_all = TRUE) # Ensure unique times
      
      # Handle cases where simulation is too short
      if (nrow(out.daily_all) < 2) return(NULL)
      
      # Calculation of daily new infections (S[t-1] - S[t])
      # Use lead function from dplyr for easier calculation
      out.daily_all <- out.daily_all %>%
        arrange(time) %>% # Ensure data is sorted by time
        mutate(S_previous = lag(S, default = first(S))) %>%
        mutate(New_Infections = pmax(0, round(S_previous - S)))
      
      # Filter out time 0 *after* calculating infections for day 1
      out.processed <- out.daily_all %>% filter(time > 0)
      
      # If no daily data remains after filtering day 0
      if (nrow(out.processed) == 0) return(NULL)
      
      # Round compartment values for clarity in plots
      # Create combined Infected count
      out.processed$Infected_Untreated <- round(out.processed$I1)
      out.processed$Infected_Treated <- round(out.processed$I2)
      out.processed$Infected_Total <- round(out.processed$I1 + out.processed$I2) # Combined I
      out.processed$Exposed <- round(out.processed$E)
      out.processed$Susceptible <- round(out.processed$S)
      out.processed$Recovered <- round(out.processed$R)
      
      # Calculate normalized values if needed
      N_total <- max(out.processed$Susceptible) + max(out.processed$Infected_Total) + max(out.processed$Recovered)
      if (N_total > 0) {
        out.processed$New_Infections_Normalized <- out.processed$New_Infections / N_total * 100
      } else {
        out.processed$New_Infections_Normalized <- 0
      }
      
      out.processed <- dplyr::select(out.processed, time, Model, Susceptible, Exposed,
                                     Infected_Untreated, Infected_Treated, Infected_Total,
                                     Recovered, New_Infections, New_Infections_Normalized)
      return(out.processed)
    }
    
    # Process each simulation result
    daily_A <- process_to_daily(out_A, "Model A")
    daily_B <- process_to_daily(out_B, "Model B")
    
    # Combine the processed data frames
    combined_data <- rbind(daily_A, daily_B)
    
    # Return NULL if no valid data could be combined
    if (is.null(combined_data) || nrow(combined_data) == 0) {
      return(NULL)
    } else {
      return(combined_data)
    }
  })
  
  
  # --- Modified Daily Infections Plot with Peak Highlighting and log/normalize options ---
  generateNewInfectionsPlot <- function() {
    plot_data <- daily_data_combined()
    
    # Check if data is available
    if (is.null(plot_data) || nrow(plot_data) == 0 || all(is.na(plot_data$New_Infections)) || sum(plot_data$New_Infections, na.rm = TRUE) == 0) {
      return(ggplot() + labs(title=strong("Daily new infections: No significant outbreak data or simulation failed"), x="Time (days)", y="Daily new infections") + theme_minimal() + theme(plot.title = element_text(face = "bold")))
    }
    
    # Determine which y-value to use based on normalization choice
    y_column <- if(input$normalize_by_population) "New_Infections_Normalized" else "New_Infections"
    
    # Calculate peaks for each model for annotation
    peaks <- plot_data %>%
      filter(!is.na(!!sym(y_column))) %>% # Filter out NAs before grouping
      group_by(Model) %>%
      filter(!!sym(y_column) == max(!!sym(y_column), na.rm = TRUE)) %>%
      # In case peak value occurs on multiple days, take the first day
      summarise(time = min(time[!!sym(y_column) == max(!!sym(y_column), na.rm=TRUE)]),
                peak_value = max(!!sym(y_column), na.rm = TRUE),
                .groups = 'drop')
    
    # Create the plot using ggplot, mapping Model to color
    p <- ggplot(plot_data, aes(x = time, y = !!sym(y_column), color = Model)) +
      geom_line(linewidth = 1, na.rm = TRUE) # Use linewidth instead of size
    
    # Add peak annotations for each model if peaks were found
    if(nrow(peaks)>0){
      for (i in 1:nrow(peaks)) {
        model_name <- peaks$Model[i]
        peak_time <- peaks$time[i]
        peak_val <- peaks$peak_value[i]
        model_color <- if(model_name == "Model A") "firebrick" else "darkblue" 
        p <- p + geom_vline(xintercept = peak_time,
                            linetype = "dashed",
                            color = model_color,
                            alpha = 0.7,
                            linewidth = 0.8)
        
        # Add point at peak
        p <- p + geom_point(data = peaks[i,],
                            aes(x = time, y = peak_value),
                            color = model_color,
                            size = 3)
        
        # Add annotation text - Format peak value based on whether we're normalizing
        label_text <- if(input$normalize_by_population) {
          paste("Peak:", round(peak_val, 2), "%\nDay:", round(peak_time))
        } else {
          paste("Peak:", scales::comma(round(peak_val)), "\nDay:", round(peak_time))
        }
        
        p <- p + annotate("text",
                          x = peak_time,
                          y = peak_val,
                          label = label_text,
                          hjust = -0.1,
                          vjust = -0.5,
                          size = 4.5,
                          color = model_color)
      }
    }
    
    # Set y-axis scale based on log option
    if(input$use_log_scale) {
      p <- p + scale_y_log10(labels = scales::comma)
    } else {
      p <- p + scale_y_continuous(labels = scales::comma, limits = c(0, NA), expand = expansion(mult = c(0, 0.2)))
    }
    
    # Set appropriate y-axis label based on normalization
    y_label <- if(input$normalize_by_population) {
      "Percentage of population newly infected"
    } else {
      "Number of new infections"
    }
    
    # Complete the plot styling
    p + labs(title = "New infections over time", 
             x = "Time (days)",
             y = y_label,
             color = "Model") +
      scale_color_manual(values = c("Model A" = "firebrick", "Model B" = "darkblue")) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "bottom",
            plot.title = element_text(face = "bold")) +
      scale_x_continuous(breaks = scales::breaks_width(10))
  }
  
  # Render the plot
  output$newInfectionsPlot <- renderPlot({
    generateNewInfectionsPlot()
  })
  
  # --- Modified SEIR Curves Plot ---
  generateSeirCurvesPlot <- function() {
    plot_data_combined <- daily_data_combined() # Use the combined daily data
    
    if (is.null(plot_data_combined) || nrow(plot_data_combined) == 0) {
      return(ggplot() + labs(title=strong("Compartment curves: Simulation failed or produced no data"), x="Time (days)", y="Number of individuals") + theme_minimal() + theme(plot.title = element_text(face = "bold")))
    }
    
    # Pivot longer for plotting multiple compartments per model
    # Include I1 (Untreated) and I2 (Treated)
    plot_data_long <- tidyr::pivot_longer(plot_data_combined,
                                          cols = c("Susceptible", "Exposed", "Infected_Untreated", "Infected_Treated", "Recovered"),
                                          names_to = "Compartment",
                                          values_to = "Count")
    
    if (nrow(plot_data_long) == 0) {
      return(ggplot() + labs(title=strong("Compartment curves: No compartment data to plot"), x="Time (days)", y="Number of individuals") + theme_minimal() + theme(plot.title = element_text(face = "bold")))
    }
    
    # Create the plot: Color by Compartment, Linetype by Model
    ggplot(plot_data_long, aes(x = time, y = Count, color = Compartment, linetype = Model)) +
      geom_line(linewidth = 1, na.rm = TRUE) + # Use linewidth instead of size
      labs(title = "Compartment sizes over time", 
           x = "Time (days)",
           y = "Number of individuals",
           color = "Compartment",
           linetype = "Model") +
      scale_color_manual(values = c("Susceptible" = "grey50", 
                                    "Exposed" = "orange",
                                    "Infected_Untreated" = "magenta", 
                                    "Infected_Treated" = "purple", 
                                    "Recovered" = "darkgreen")) +
      scale_linetype_manual(values = c("Model A" = "solid", "Model B" = "dashed")) + 
      theme_minimal(base_size = 14) +
      theme(legend.position = "bottom",
            legend.box = "vertical",
            plot.title = element_text(face = "bold")) + 
      scale_y_continuous(labels = scales::comma, limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
      scale_x_continuous(breaks = scales::breaks_width(10))
  }
  
  # Render the SEIR curves plot
  output$seirCurvesPlot <- renderPlot({
    generateSeirCurvesPlot()
  })
  
  # --- Output for Summary Statistics (for both models) ---
  output$summaryStats <- renderText({
    # Fetch combined data once
    combined_data <- daily_data_combined()
    
    # Filter data for each model
    sim_data_A <- if (!is.null(combined_data)) combined_data %>% filter(Model == "Model A") else NULL
    sim_data_B <- if (!is.null(combined_data)) combined_data %>% filter(Model == "Model B") else NULL
    
    # Helper function to generate summary for one model's data
    generate_summary <- function(sim_data, model_name, N_model, simtime_common, E0_common, R0, alp, adi, theta_perc, treatment_delay) {
      # Use inputs passed to function
      N <- N_model # This is the input population size (scalar)
      simtime <- simtime_common
      E0 <- E0_common
      theta_val <- theta_perc / 100 # Convert percentage back for display/use if needed
      
      
      # --- Parameter Summary String ---
      param_summary <- paste(
        model_name, " summary:\n",
        "--------------------\n",
        "Parameters:\n",
        sprintf("  Population size (N): %s\n", format(round(N), big.mark=",")),
        sprintf("  R0: %.2f\n", R0),
        sprintf("  Avg. latency period: %d days\n", alp),
        sprintf("  Avg. infectivity period: %d days\n", adi),
        sprintf("  Simulation duration: %d days\n", simtime),
        sprintf("  Initial exposed (E0): %s\n", format(round(E0), big.mark=",")),
        "Intervention Parameters:\n",
        sprintf("  Treatment Effectiveness (theta): %.2f (%d%% reduction)\n", theta_val, theta_perc),
        sprintf("  Average treatment initiation delay: %d days\n\n", treatment_delay),
        sep=""
      )
      
      # --- Results Summary String ---
      if (is.null(sim_data) || nrow(sim_data) == 0) {
        result_summary <- paste(
          "Results:\n",
          "--------------------\n",
          "  Simulation failed or produced no results.\n\n",
          sep=""
        )
      } else {
        # Calculate peak daily infections
        peak_infections_val <- 0
        peak_infections_time <- NA
        valid_infections <- !is.na(sim_data$New_Infections) & sim_data$New_Infections > 0
        if (any(valid_infections)) {
          peak_row <- sim_data %>% filter(New_Infections == max(New_Infections[valid_infections], na.rm = TRUE)) %>% slice(1)
          peak_infections_val <- peak_row$New_Infections
          peak_infections_time <- peak_row$time
        }
        
        # Calculate peak prevalence (concurrently infected - using Infected_Total)
        peak_prevalence_val <- 0
        peak_prevalence_time <- NA
        valid_prevalence <- !is.na(sim_data$Infected_Total) & sim_data$Infected_Total > 0
        if (any(valid_prevalence)) {
          peak_row <- sim_data %>% filter(Infected_Total == max(Infected_Total[valid_prevalence], na.rm = TRUE)) %>% slice(1)
          peak_prevalence_val <- peak_row$Infected_Total
          peak_prevalence_time <- peak_row$time
        }
        
        # Calculate total infected (cumulative) = N - S(end)
        final_S <- tail(sim_data$Susceptible[!is.na(sim_data$Susceptible)], 1)
        total_infected <- if(length(final_S) > 0) round(N - final_S) else NA
        
        # --- Final State Distribution ---
        final_state <- tail(sim_data, 1)
        if (nrow(final_state) > 0 && N > 0) {
          final_S_count <- final_state$Susceptible
          final_E_count <- final_state$Exposed
          final_I1_count <- final_state$Infected_Untreated
          final_I2_count <- final_state$Infected_Treated
          final_R_count <- final_state$Recovered
          # Recalculate N_final based on sum to handle potential small discrepancies
          N_final_check = final_S_count + final_E_count + final_I1_count + final_I2_count + final_R_count
          
          final_dist_table <- paste(
            "Final State Distribution (Day ", simtime, "):\n",
            "  Compartment         Count      Percentage\n",
            "  ------------------------------------------\n",
            sprintf("  Susceptible      %10s    %6.2f%%\n", format(final_S_count, big.mark=",", trim=TRUE), (final_S_count / N) * 100),
            sprintf("  Exposed          %10s    %6.2f%%\n", format(final_E_count, big.mark=",", trim=TRUE), (final_E_count / N) * 100),
            sprintf("  Infected (Untr.) %10s    %6.2f%%\n", format(final_I1_count, big.mark=",", trim=TRUE), (final_I1_count / N) * 100),
            sprintf("  Infected (Tr.)   %10s    %6.2f%%\n", format(final_I2_count, big.mark=",", trim=TRUE), (final_I2_count / N) * 100),
            sprintf("  Recovered        %10s    %6.2f%%\n", format(final_R_count, big.mark=",", trim=TRUE), (final_R_count / N) * 100),
            "  ------------------------------------------\n",
            sprintf("  Check Total      %10s    %6.2f%%\n\n", format(round(N_final_check), big.mark=",", trim=TRUE), (N_final_check / N) * 100), # Should be close to N and 100%
            sep=""
          )
        } else {
          final_dist_table <- "  Final state distribution could not be calculated.\n\n"
        }
        
        
        # --- Combine Results Text ---
        result_summary <- paste(
          "Results:\n",
          "--------------------\n",
          "  Peak daily new infections: ", if(is.na(peak_infections_time) || peak_infections_val <= 0) "N/A" else paste0(format(round(peak_infections_val), big.mark=","), " on day ", round(peak_infections_time)), "\n",
          "  Peak number infected (prevalence I1+I2): ", if(is.na(peak_prevalence_time) || peak_prevalence_val <= 0) "N/A" else paste0(format(round(peak_prevalence_val), big.mark=","), " on day ", round(peak_prevalence_time)), "\n",
          "  Estimated total infected (cumulative): ", if(is.na(total_infected) || total_infected < 0) "N/A" else format(total_infected, big.mark=","), "\n\n",
          final_dist_table, # Add the final state table here
          sep=""
        )
      }
      
      return(paste(param_summary, result_summary, sep=""))
    } # end generate_summary function
    
    # Generate summaries for both models using their respective inputs
    # Pass model-specific N (*1e6), theta (as %), treatment_delay and common inputs
    req(input$N_input_A, input$N_input_B, input$theta_A, input$theta_B) # Ensure inputs are available
    
    summary_A <- generate_summary(sim_data_A, "Model A", input$N_input_A * 1e6, input$simtime_input_common, input$E0_input_common,
                                  input$R0_input_A, input$latency_period_A, input$infectious_period_A,
                                  input$theta_A, input$treatment_delay_A) # Pass theta_A as percentage
    summary_B <- generate_summary(sim_data_B, "Model B", input$N_input_B * 1e6, input$simtime_input_common, input$E0_input_common,
                                  input$R0_input_B, input$latency_period_B, input$infectious_period_B,
                                  input$theta_B, input$treatment_delay_B) # Pass theta_B as percentage
    
    # Combine the summaries
    paste(summary_A, summary_B, sep="")
    
  }) # end renderText
  
  # Download handler for the new infections plot
  output$downloadNewInfectionsPlot <- downloadHandler(
    filename = function() {
      paste("new-infections-plot-", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = generateNewInfectionsPlot(), device = "png", 
             width = 10, height = 7, dpi = 300)
    }
  )
  
  # Download handler for the SEIR curves plot
  output$downloadSeirCurvesPlot <- downloadHandler(
    filename = function() {
      paste("seir-curves-plot-", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = generateSeirCurvesPlot(), device = "png", 
             width = 10, height = 7, dpi = 300)
    }
  )
  
  # Download handler for the summary statistics
  output$downloadSummaryStats <- downloadHandler(
    filename = function() {
      paste("epidemic-summary-", Sys.Date(), ".txt", sep = "")
    },
    content = function(file) {
      # Create a summary function that returns the text rather than rendering it
      combined_data <- daily_data_combined()
      
      # Filter data for each model
      sim_data_A <- if (!is.null(combined_data)) combined_data %>% filter(Model == "Model A") else NULL
      sim_data_B <- if (!is.null(combined_data)) combined_data %>% filter(Model == "Model B") else NULL
      
      # Reuse your generate_summary function (already defined in renderText)
      summary_A <- generate_summary(sim_data_A, "Model A", input$N_input_A * 1e6, input$simtime_input_common, input$E0_input_common,
                                    input$R0_input_A, input$latency_period_A, input$infectious_period_A,
                                    input$theta_A, input$treatment_delay_A)
      summary_B <- generate_summary(sim_data_B, "Model B", input$N_input_B * 1e6, input$simtime_input_common, input$E0_input_common,
                                    input$R0_input_B, input$latency_period_B, input$infectious_period_B,
                                    input$theta_B, input$treatment_delay_B)
      
      # Combine and write the summaries
      writeLines(paste(summary_A, summary_B, sep=""), file)
    }
  )
} # end server

# Run the application
shinyApp(ui = ui, server = server)