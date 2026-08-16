# Pathogen Framework — Evolutionary Pathogen/Detector Simulation

This repository contains the implementation used in the project, based on the supplied
`Pathogen Framework` Jupyter notebook.

## What is implemented

The simulation models binary-string pathogens and detectors. The main algorithmic components are:

- contiguous-match recognition (`r_contiguous_match`)
- negative selection detector generation
- pathogen mutation and one-point crossover
- genetic-algorithm-style pathogen evolution with 50% elitist survival
- Shannon diversity and mean Hamming-distance metrics
- detector coverage, pathogen fitness, and protein-score measurements
- parameter sweeps over string length, detector population size, and matching threshold

The original notebook is preserved in `notebooks/Pathogen_Framework.ipynb`.

## Repository contents

```text
.
├── README.md
├── pathogen_framework.py
├── notebooks/
│   └── Pathogen_Framework.ipynb
├── data/
│   └── README.md
├── results/
│   ├── simulation_results.csv
    └── figures/
        └── sample_run.png
```

No external biological dataset is required. Pathogens and detectors are generated
synthetically by the program.

## Installation

Python 3.11+ is recommended.

Create and activate a virtual environment:

```bash
python -m venv .venv
```

Windows:

```bash
.venv\Scripts\activate
```

macOS/Linux:

```bash
source .venv/bin/activate
```

Install the pinned dependencies:

```bash
python -m pip install --upgrade pip
pip install -r requirements.txt
```

## Run a small sample

The sample run is intentionally small so that it can be reproduced quickly:

```bash
python src/pathogen_framework.py --sample --seed 42
```

This generates:

- `results/sample_simulation_results.csv`
- `figures/sample_run.png`

A captured sample console run is included in `results/sample_run.txt`.

## Run the full experiment

The default experiment uses:

- string lengths: `l = 5, 10, 20`
- detector populations: `N = 500, 1000, 2000`
- matching thresholds: `r = 2, 3, 8` (configurations with `r > l` are skipped)
- 100 generations
- 3 repetitions per valid configuration
- 1,000 pathogens per simulation

Run:

```bash
python src/pathogen_framework.py --seed 42
```

To change the number of generations or repetitions:

```bash
python src/pathogen_framework.py --generations 100 --repetitions 3 --seed 42
```

## Jupyter notebook

The original notebook can be opened with:

```bash
jupyter notebook notebooks/Pathogen_Framework.ipynb
```

The notebook contains the simulation and plotting workflow used in the project.

## Reproducibility

The command-line implementation accepts `--seed` so that a run can be made deterministic.
The original notebook did not specify a package lockfile or fixed random seed, so the exact
original execution environment cannot be reconstructed from the notebook alone. The
`requirements.txt` versions in this repository record the versions used to prepare this
repository package.

## Sample behavior

A sample plot and sample CSV are included so that a reviewer can verify that the system
produces simulation output without first running the full experiment.

## License

Add the license required by your course/project before publishing if one is specified.
