# Performance and Robustness of Numerical Optimisation Algorithms under Aggregation-Induced Likelihood Compression

This repository contains the complete production-ready R replication pipeline, dataset acquisition routines, Monte Carlo simulation engines, and figure generation code for the manuscript:

> **"Performance and Robustness of Numerical Optimisation Algorithms under Aggregation-Induced Likelihood Compression"**
> *Ahmad Hakiim Jamaluddin, Farid Zamani Che Rose, Muhammad Aslam Mohd Safari, Muhammad Jaffri Mohd Nasir, Andrea Tri Rian Dani, and Syaiful Anam.*

---

## Overview

Repeated maximum likelihood estimation (MLE) over overlapping time series windows is computationally intensive. This project introduces a histogram-symbolic MLE framework for temperature–humidity index (THI) time series, compressing raw sample evaluations from $O(n)$ to $O(B)$ bin boundary calculations per iteration.

Key methodological features of this repository include:

* **Factorial Monte Carlo Simulation Grid:** 108 core cells (3 models $\times$ 2 primary optimisers $\times$ 3 sample sizes $\times$ 3 bin resolutions $\times$ 2 bin types) over 1,000 replicates ($21,600$ total primary fits).
* **Parametric Model Families:** Normal, Logistic, and Shifted Gamma (with Shifted Lognormal sensitivity checks).
* **Numerical Optimiser Benchmarks:** Quasi-Newton routines (`BFGS`, `L-BFGS-B`), derivative-free simplex (`Nelder-Mead`), and conjugate gradients (`CG`) with analytical callbacks.
* **Empirical NASA POWER Benchmark:** Automated API ingestion of 41 years (1984–2025) of daily meteorological data across 5 contrasting climate zones (Subarctic, Tropical Rainforest, Desert, Oceanic, Humid Continental).

---

## Repository Structure

```text
├── main_thi_pipeline.R        # Master production execution script
├── README.md                  # Repository documentation and replication guide
└── results_symbolic_thi/      # Output directory (auto-generated upon execution)
    ├── cache/                 # Local CSV cache for downloaded NASA POWER raw data
    ├── run_config.json        # Execution parameter JSON record
    ├── session_info.txt       # R session environment metadata
    ├── *.csv                  # Summary and detailed Monte Carlo simulation tables
    ├── *.tex                  # Publication-ready LaTeX tabular rows
    ├── fig_speedup.png        # Simulation speed-up comparative plot
    ├── fig_relrmse.png        # Simulation relative RMSE accuracy plot
    └── fig_empirical_thi_bar_overlay.png  # Decadal overlay publication figure

```

---

## System Requirements & Dependencies

All scripts were developed and tested in the **R statistical computing environment** (version 4.5+ recommended). Parallel execution is natively supported via `future` and `future.apply`.

### Required R Packages

To install all necessary dependencies automatically, execute the following in your R console:

```R
required_pkgs <- c(
  "dplyr", "tibble", "purrr", "readr", "tidyr", "ggplot2",
  "stringr", "jsonlite", "httr2", "future", "future.apply", 
  "parallel", "R.utils", "patchwork"
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) install.packages(missing_pkgs, repos = "https://cloud.r-project.org")

```

---

## Replication Instructions

### 1. Clone the Repository

```bash
git clone https://github.com/HakiimJ/symbolic-thi-likelihood.git
cd symbolic-thi-likelihood

```

### 2. Run the Execution Pipeline

Execute the master pipeline from your terminal or R environment:

```bash
Rscript main_thi_pipeline.R

```

*Note: The pipeline automatically manages parallel core detection, caches NASA POWER daily weather streams locally in `results_symbolic_thi/cache/`, and outputs formatted LaTeX tables alongside high-resolution figures ($300\text{ DPI}$).*

---

## Citation & License

If you use or adapt this codebase for your research, please cite our manuscript:

```bibtex
@article{jamaluddin2026symbolic,
  title={Performance and Robustness of Numerical Optimisation Algorithms under Aggregation-Induced Likelihood Compression},
  author={Jamaluddin, Ahmad Hakiim and Nasir, Muhammad Jaffri Mohd and Rose, Farid Zamani Che and Safari, Muhammad Aslam Mohd and Dani, Andrea Tri Rian and Anam, Syaiful},
  journal={arXiv},
  note={Under Review},
  year={2026}
}

```

This project is licensed under the [MIT License](https://opensource.org/licenses/MIT).
