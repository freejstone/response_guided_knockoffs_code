# Response-guided knockoffs

This repository contains the manuscript and reproducibility code for **“Response-guided knockoffs for directional FDR control in linear models”** by Jack Freestone, Garth Tarr, Samuel Muller, and Uri Keich.

The project studies feature selection when discoveries are sought in a pre-specified direction. A motivating example is identifying HIV mutations that increase resistance to a particular drug. The proposed response-guided knockoff filter uses a noise-perturbed response to guide the knockoff construction while retaining finite-sample directional false discovery rate (FDR) control.

## Methods implemented

The active R implementation includes:

- `RG(signed, ME)`: response-guided knockoffs with OLS signed screening and a maximum-entropy construction;
- `RG(reg-signed, ME)`: response-guided knockoffs with elastic-net signed screening and a maximum-entropy construction;
- `S(equi)` and `S(SDP)`: signed knockoff filters using equi-correlated and SDP constructions;
- `S(ME)` and `S(MVR)`: signed knockoff filters using maximum-entropy and minimum variance-based reconstructability constructions.

The ME and MVR optimizers follow the algorithms of Spector and Janson. The R implementation supports randomized fixed-X knockoff construction through a randomized orthogonal completion. Design columns are rescaled to unit Euclidean norm without mean-centring.

## Repository structure

| Path | Purpose |
| --- | --- |
| `functions/functions.R` | Main knockoff constructions and response-guided and signed filters. |
| `functions/helpers.R` | Supporting numerical and data-handling functions. |
| `functions/simulation_functions.R` | Data-generating functions used by the simulations. |
| `naive_approaches_sim/` | Motivating comparison of naive directional procedures, the signed filter, and the response-guided filter. |
| `simulation_1_paper/` | Simulation 1, covering six covariance models with `n >= 2p`. |
| `simulation_2_paper/` | Simulation 2, covering the same six models with `p + 2 < n < 2p`. |
| `HIV_resistance/` | HIV protease-inhibitor resistance analysis and its source data. |
| `plot_ratio.R` | Combined and experiment-specific relative-power figures. |

Each paper simulation contains six numbered subdirectories corresponding to the following design families from Spector and Janson:

1. AR1;
2. correlated AR1 coefficients;
3. block equicorrelation;
4. Erdős–Rényi covariance;
5. Erdős–Rényi precision;
6. equicorrelation.

Within each numbered directory, `run_sim.R` performs one simulation run, `run_one.sh` provides its shell wrapper, and `launch.sh` runs a grid of sample sizes, sparsity levels, and Monte Carlo repetitions in parallel.

## Requirements

The code requires R. The active R dependencies are:

```r
install.packages(c(
  "knockoff", "glmnet", "pracma", "foreach", "doParallel", "doRNG",
  "EnvStats", "hdf5r", "tidyverse", "jsonlite", "withr", "remotes"
))
```

## Reproducing the results

Run the following commands from the repository root. The full simulations are computationally intensive.

### Motivating naive-procedure simulation

```bash
cd naive_approaches_sim
bash run.sh 15
```

The numeric argument is the maximum number of concurrent R workers. This command generates the shared data, runs all methods, aggregates the results, and creates the FDR and power plots.

### Simulation 1

Run all six design families:

```bash
cd simulation_1_paper
for i in 1 2 3 4 5 6; do
  (cd "simulation_${i}" && bash launch.sh -j 20)
done
Rscript plot.R
```

### Simulation 2

```bash
cd simulation_2_paper
for i in 1 2 3 4 5 6; do
  (cd "simulation_${i}" && bash launch.sh -j 20)
done
Rscript plot.R
```

### HIV drug-resistance analysis

```bash
cd HIV_resistance/run
Rscript run_PI.R
cd ..
Rscript plot.R
```

The analysis uses resistance measurements for the protease inhibitors atazanavir (`ATV`), lopinavir (`LPV`), and darunavir (`DRV`). Predictor columns indicate the presence of particular protease mutations in each viral isolate. The source data and approximate ground-truth mutation list are stored under `HIV_resistance/run/RB_data/`.

The combined relative-power visualization is produced with:

```bash
Rscript plot_ratio.R
```

Before running `plot_ratio.R`, set the `base` variable at the top of that script to the repository root for your local checkout.
