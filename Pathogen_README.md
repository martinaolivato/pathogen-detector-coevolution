# Immune System Simulation

A genetic algorithm-based simulation implementing negative selection and pathogen evolution dynamics, inspired by biological immune system principles.

## Project Overview

This project simulates the coevolution between:
- **Detectors**: Artificial immune cells that recognize pathogens using a contiguous matching rule
- **Pathogens**: Self-replicating entities that evolve to evade detection while maintaining structural integrity

The simulation explores how immune pressure and structural constraints shape pathogen evolution.

## Installation

### Prerequisites
- Python 3.8 or higher
- pip (Python package manager)

### Setup

1. Clone the repository:
   ```bash
   git clone https://github.com/martinaolivato/immune-system-simulation.git
   cd immune-system-simulation
   ```

2. Create a virtual environment:
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

## Usage

### Running the Main Simulation

Execute the full parameter sweep:
```bash
python simulation.py
```

This will:
- Run 3 repetitions for each parameter configuration
- Test all combinations of:
  - String lengths: 5, 10, 20 bits
  - Detector counts: 500, 2000
  - Matching thresholds (r): 2, 3, 8
  - Structural weights (alpha): 0, 0.2, 0.5, 0.8, 1.0
- Generate `comprehensive_sweep_results.csv` with results

### Understanding Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `NUM_PATHOGENS` | Number of pathogens per generation | 1000 |
| `GENERATIONS` | Evolution timesteps | 100 |
| `MUTATION_DETECTOR` | Detector bit-flip probability | 0.05 |
| `MUTATION_PATHOGEN` | Pathogen bit-flip probability | 0.5 |
| `THETA` | Immune pressure scaling factor | 0.5 |
| `alpha` | Weight for structure vs. evasion (0=pure evasion, 1=pure structure) | 0.2 |

### Output

- **comprehensive_sweep_results.csv**: Results matrix containing:
  - Mean Pathogen Fitness and standard deviation
  - Mean Protein Score and standard deviation
  - Mean Detector Survival Rate

### Example Analysis

```python
import pandas as pd

# Load and analyze results
df = pd.read_csv("comprehensive_sweep_results.csv")

# Filter for high immune pressure (more detectors)
high_pressure = df[df["Initial Detectors (N)"] == 2000]
print(high_pressure[["String Length (l)", "Structural Weight (alpha)", "Mean Pathogen Fitness"]])
```

## Algorithm Description

### Negative Selection (Detectors)

1. Generate random bit strings as detector candidates
2. Remove any detector that matches the "self" pattern (all zeros) with r contiguous bits
3. Survive detectors that catch pathogens in random samples

### Pathogen Evolution (Genetic Algorithm)

1. **Fitness Evaluation**: 
   - Evasion component: Based on detection ratio and immune pressure
   - Structural component: Proportion of 1s in bitstring
   - Combined: `fitness = (1-α) * evasion + α * structure`

2. **Selection**: Keep top 50% of pathogens (elitism)

3. **Reproduction**: Randomly crossover survivors and mutate offspring

## File Structure

```
immune-system-simulation/
├── README.md
├── requirements.txt
├── simulation.py
├── sample_output.csv
└── docs/
    └── algorithm_explanation.md
```

## Results & Sample Output

See `sample_output.csv` for example results from a test run.

## References

- de Castro, L. N., & Von Zuben, F. J. (2002). Learning and optimization in the immune system. IEEE transactions on neural networks, 13(5), 1196-1209.
- Forrest, S., & Hofmeyr, S. A. (2000). Immunology as information processing. Design Principles for the Immune System and Other Distributed Autonomous Systems.

## License

MIT License - see LICENSE file for details

## Author

martinaolivato