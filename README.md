# Simulating Epigenetic and Genetic Divergence Under Migration and Environmental Influence

This repository provides a conceptual and computational framework to explore how migration, environmental variation, and epigenetic heritability shape population-level patterns of epigenetic and genetic divergence. It includes an individual-based simulation model and theoretical recursive models implemented in R.

## Overview

The code simulates two populations with distinct environmental baselines and tracks the evolution of a quantitative epigenetic trait (bounded between 0 and 1) and multilocus genotypes over multiple generations. Simulations account for:

* Epigenetic heritability (`h_epi`)
* Bidirectional migration (`m`)
* Environmental stochasticity (`env_sd`)
* Selection on trait-environment matching

Comparisons are made between simulated results and theoretical predictions to understand divergence dynamics over time.

## Directory Structure

```
├── EpiFst.R                            # Main script to run simulations and summarize results
├── output/                             # Directory for simulation output files (.rds, .txt)
├── source/                             # Contains simulation and helper functions
│   ├── Simulate.R                      # Core simulation function
│   ├── CalcFstEpi.R                    # Function to calculate FST from trait data
│   ├── CalcFstGen.R                    # Function to calculate FST from genetic genotpe data
│   ├── Simulation_plots.R              # ABM simulation output visualization
│   ├── Recursive_theory.R              # Theoretical model visualization
│   ├── Theory_simple.R                 # Theoretical model, baseline
│   ├── Theory_migration.R              # Theoretical model, migration
│   ├── Theory_full.R                   # Theoretical model, sampling error
│   ├── Theory_accumulative_migration.R # Theoretical model, accumulating ancestry
│   ├── Theory_WP.R                     # Theoretical model, approximation of Webster and Philips 2024
│   ├── TheorySim_comp.R                # Theoretical and simulation comparison 
│   └── FunctionSourcer.R               # Sourcing script that loads dependencies and functions

```

## Getting Started

Use `EpiFst.R` in R to set parameters and run the simulation.

Generated output files will be saved in the `output/` directory. Use `Plots.R` and `Theory.R` to visualize and compare results.

## Key Parameters

Set in `EpiFst.R`:

* `h_vals`: Vector of epigenetic heritability values
* `m_vals`: Vector of migration rates
* `n_reps`: Number of replicate simulations per parameter combination
* `env_A`, `env_B`: Environmental baselines for populations A and B

## Output

Each simulation run produces:

* `simulation_output_<timestamp>.rds`: Contains raw and aggregated results
* `simulation_progress_<timestamp>.txt`: Log file recording simulation progress
* `EpiFst.R`: Copy of runfile that includes parameter value run settings

