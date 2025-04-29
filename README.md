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
```

### Deployment

The app can be deployed to [shinyapps.io](https://www.shinyapps.io/) or [Posit Connect](https://posit.co/products/enterprise/connect/) for broader accessibility.

## Use Cases

This simulator is valuable for:

- **Education**: Teaching epidemiological concepts and intervention dynamics
- **Public Health Planning**: Evaluating the potential impact of treatment strategies
- **Research**: Exploring theoretical questions about epidemic control
- **Decision Support**: Communicating the importance of rapid diagnosis and treatment
- **Policy Development**: Demonstrating the value of healthcare system capacity and response time

## Mathematical Foundation

The application uses ordinary differential equations (ODEs) to model disease transmission:

```
dS/dt = -β * S * (I₁ + (1-θ) * I₂) / N
dE/dt = β * S * (I₁ + (1-θ) * I₂) / N - ω * E
dI₁/dt = ω * E - φ * I₁ - γ * I₁
dI₂/dt = φ * I₁ - γ * I₂
dR/dt = γ * I₁ + γ * I₂
```

Where:
- β = transmission rate
- ω = rate of progression from exposed to infectious (1/latency period)
- γ = recovery rate (1/infectious period)
- θ = treatment effectiveness (reduction in transmission)
- φ = rate of treatment initiation (1/treatment delay)
- N = total population

## Future Development

Planned enhancements include:

- Stochastic modeling options
- Age-structured populations
- Vaccination interventions
- Hospitalization and healthcare capacity modeling
- Exportable reports and data analysis
- Regional comparison tools

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add some amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- Based on epidemiological models developed for public health response planning
- Inspired by the need for accessible tools to understand intervention timing
- Special thanks to the R and Shiny communities for their excellent documentation and support

---

*This simulator is intended for educational and planning purposes. Real-world epidemics involve complex social, biological, and environmental factors not captured in this model.*
