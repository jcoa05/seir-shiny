# Epidemic Intervention Simulator

[![Shiny](https://img.shields.io/badge/Powered%20by-Shiny-blue?style=flat&logo=r&logoColor=white)](https://shiny.rstudio.com/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

An interactive epidemic simulation tool that demonstrates the impact of treatment timing and effectiveness on infectious disease outbreaks.

## Overview

This Shiny web application simulates the dynamics of infectious disease epidemics using a modified SE2IR (Susceptible-Exposed-Infectious[untreated]-Infectious[treated]-Recovered) compartmental model. The simulator allows users to explore how different treatment strategies affect epidemic outcomes across different populations.

### Key Features:

- **Dual Model Comparison**: Simultaneously compare epidemic trajectories between two populations with different parameters
- **Real-time Visualization**: See how changes in treatment timing and effectiveness alter epidemic curves
- **Treatment Intervention Analysis**: Model the impact of antivirals or other treatments that reduce transmission
- **Comprehensive Statistics**: Track peak infections, total cases, and final epidemic state
- **Multiple Visualization Options**: View new infections or all compartment populations over time
- **Pre-configured Scenarios**: Explore optimistic, pessimistic, and no-treatment scenarios

## The Model

The application implements a deterministic compartmental model with five states:

- **S (Susceptible)**: Healthy individuals at risk of infection
- **E (Exposed)**: Infected but not yet infectious (latent period)
- **I₁ (Infectious, Untreated)**: Infectious individuals who have not yet received treatment
- **I₂ (Infectious, Treated)**: Infectious individuals who are receiving treatment with reduced transmission
- **R (Removed)**: Recovered or deceased individuals (no longer infectious)

The model incorporates treatment timing (average delay to treatment) and effectiveness (reduction in transmission) to simulate realistic intervention scenarios.

## Interactive Features

Users can modify parameters including:

- **Treatment Effectiveness**: Reduction in transmission for treated individuals (0-100%)
- **Treatment Delay**: Average time between becoming infectious and starting treatment (0-8 days)
- **R₀**: Basic reproduction number (average secondary cases per infectious individual)
- **Population Size**: Total population for each model scenario
- **Disease Progression**: Latency and infectious period durations
- **Visualization Options**: Log scale, population normalization, and more

## Getting Started

### Prerequisites

- R 4.0.0 or higher
- Required packages: shiny, deSolve, ggplot2, tidyr, dplyr, scales, shinyBS, bslib, httpuv

### Installation

```r
# Clone this repository
git clone https://github.com/username/epidemic-intervention-simulator.git

# Navigate to the project directory
cd epidemic-intervention-simulator

# Install required packages
install.packages(c("shiny", "deSolve", "ggplot2", "tidyr", "dplyr", "scales", "shinyBS", "bslib", "httpuv"))

# Run the application
shiny::runApp()
