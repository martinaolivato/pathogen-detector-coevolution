import argparse
import itertools
import random
from pathlib import Path
from collections import Counter
from itertools import combinations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import entropy


NUM_PATHOGENS = 1000
GENERATIONS = 100
MUTATION_PATHOGEN = 0.5
MUTATION_DETECTOR = 0.2
REPETITIONS = 3


def random_bitstring(n):
    return ''.join(random.choice('01') for _ in range(n))


def r_contiguous_match(s1, s2, r):
    """Return True if s1 and s2 share at least r contiguous matching symbols."""
    count = 0
    for i in range(len(s1)):
        if s1[i] == s2[i]:
            count += 1
            if count >= r:
                return True
        else:
            count = 0
    return False


def mutate_bitstring(s, mutation_rate):
    bits = list(s)
    for i in range(len(bits)):
        if random.random() < mutation_rate:
            bits[i] = '1' if bits[i] == '0' else '0'
    return ''.join(bits)


def shannon(bitstrings):
    counts = Counter(bitstrings)
    return entropy(list(counts.values()))


def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))


def mean_hamming(bitstrings):
    if len(bitstrings) < 2:
        return 0
    distances = [hamming(a, b) for a, b in combinations(bitstrings, 2)]
    return np.mean(distances)


class Detector:
    def __init__(self, bitstring, r, mutation_detector=MUTATION_DETECTOR):
        self.bitstring = bitstring
        self.r = r
        self.mutation_detector = mutation_detector
        self.fitness = 1.0

    def detect(self, pathogen):
        return r_contiguous_match(self.bitstring, pathogen.bitstring, self.r)


class Pathogen:
    def __init__(self, bitstring, mutation_pathogen=MUTATION_PATHOGEN, l=None):
        self.bitstring = bitstring
        self.mutation_pathogen = mutation_pathogen
        self.l = len(bitstring) if l is None else l
        self.fitness = 0.0

    def protein_score(self):
        return self.bitstring.count('1') / self.l

    def detection_score(self, detectors):
        if not detectors:
            return 0.0
        detections = sum(1 for d in detectors if d.detect(self))
        return detections > 0

    def evaluate_fitness(self, detectors):
        if self.detection_score(detectors):
            self.fitness = 0.0
        else:
            self.fitness = self.protein_score()

    def mutate(self):
        self.bitstring = mutate_bitstring(self.bitstring, self.mutation_pathogen)

    def crossover(self, other):
        point = random.randint(1, self.l - 1)
        child_bitstring = self.bitstring[:point] + other.bitstring[point:]
        mutation_rate = (self.mutation_pathogen + other.mutation_pathogen) / 2
        return Pathogen(child_bitstring, mutation_rate, self.l)


def generate_detectors(num_detectors, l, r):
    detectors = []
    self_string = "0" * l

    while len(detectors) < num_detectors:
        candidate = random_bitstring(l)
        if not r_contiguous_match(candidate, self_string, r):
            detectors.append(Detector(candidate, r))
    return detectors


def generate_pathogens(num_pathogens, l):
    return [Pathogen(random_bitstring(l), MUTATION_PATHOGEN, l)
            for _ in range(num_pathogens)]


def evolve_pathogens(pathogens, detectors):
    for p in pathogens:
        p.evaluate_fitness(detectors)

    pathogens.sort(key=lambda x: x.fitness, reverse=True)
    survivors = pathogens[:len(pathogens) // 2]
    children = []

    while len(children) < len(pathogens) // 2:
        parent1 = random.choice(survivors)
        parent2 = random.choice(survivors)
        child = parent1.crossover(parent2)
        child.mutate()
        children.append(child)

    return survivors + children


def run_simulation(l, n_detectors, r, generations=GENERATIONS, rep=1,
                   num_pathogens=NUM_PATHOGENS):
    detectors = generate_detectors(n_detectors, l, r)
    pathogens = generate_pathogens(num_pathogens, l)

    detector_strings = [d.bitstring for d in detectors]
    detector_shannon = shannon(detector_strings)
    detector_hamming = mean_hamming(detector_strings)
    unique_detectors = len(set(detector_strings))

    rows = []
    for generation in range(generations):
        pathogens = evolve_pathogens(pathogens, detectors)
        pathogen_strings = [p.bitstring for p in pathogens]

        mean_fitness = float(np.mean([p.fitness for p in pathogens]))
        mean_protein = float(np.mean([p.protein_score() for p in pathogens]))
        unique_pathogens = len(set(pathogen_strings))
        pathogen_shannon = shannon(pathogen_strings)
        pathogen_hamming = mean_hamming(pathogen_strings)

        active_detectors = sum(
            any(d.detect(p) for p in pathogens) for d in detectors
        )
        detector_coverages = [
            sum(d.detect(p) for p in pathogens) / len(pathogens)
            for d in detectors
        ]
        average_detector_coverage = float(np.mean(detector_coverages))

        rows.append({
            "String Length (l)": l,
            "Initial Detectors (N)": n_detectors,
            "Matching Threshold (r)": r,
            "Repetition": rep,
            "Generation": generation,
            "Mean Pathogen Fitness": mean_fitness,
            "Mean Protein Score": mean_protein,
            "Unique Detectors": unique_detectors,
            "Detector Shannon": detector_shannon,
            "Detector Hamming": detector_hamming,
            "Unique Pathogens": unique_pathogens,
            "Pathogen Shannon": pathogen_shannon,
            "Pathogen Hamming": pathogen_hamming,
            "Active Detectors": active_detectors,
            "Average Detector Coverage": average_detector_coverage,
        })
    return rows


def run_experiment(l_values, detector_grid, r_values, generations,
                   repetitions, num_pathogens, seed=None):
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)

    all_rows = []
    for l, n_det, r in itertools.product(l_values, detector_grid, r_values):
        if r > l:
            continue
        print(f"Processing Config -> l: {l}, N: {n_det}, r: {r}")
        for rep in range(1, repetitions + 1):
            all_rows.extend(
                run_simulation(l, n_det, r, generations, rep, num_pathogens)
            )
    return pd.DataFrame(all_rows)


def make_sample_plot(df, output_path):
    plot_df = (
        df.groupby(
            ["String Length (l)", "Initial Detectors (N)",
             "Matching Threshold (r)", "Generation"],
            as_index=False
        ).mean(numeric_only=True)
    )
    plt.figure(figsize=(8, 6))
    sns.lineplot(
        data=plot_df,
        x="Generation",
        y="Mean Pathogen Fitness",
        hue="Matching Threshold (r)",
        marker="o",
    )
    plt.title("Sample: Evolution of Mean Pathogen Fitness")
    plt.tight_layout()
    plt.savefig(output_path, dpi=160)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Pathogen/immune-system evolutionary simulation.")
    parser.add_argument("--generations", type=int, default=100)
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--pathogens", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--output", default="results/simulation_results.csv")
    parser.add_argument("--sample", action="store_true",
                        help="Run a small demonstration instead of the full parameter grid.")
    args = parser.parse_args()

    if args.sample:
        l_values, detector_grid, r_values = [5], [50], [2]
        generations, repetitions, num_pathogens = 10, 1, 100
    else:
        l_values, detector_grid, r_values = [5, 10, 20], [500, 1000, 2000], [2, 3, 8]
        generations, repetitions, num_pathogens = args.generations, args.repetitions, args.pathogens

    df = run_experiment(
        l_values, detector_grid, r_values,
        generations, repetitions, num_pathogens, args.seed
    )
    output = Path(args.output)
    if args.sample and args.output == 'results/simulation_results.csv':
        output = Path('results/sample_simulation_results.csv')
    output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output, index=False)
    print(f"\nSimulation process complete. Dataset exported to '{output}'.")
    print(df.head())

    if args.sample:
        make_sample_plot(df, "figures/sample_run.png")


if __name__ == "__main__":
    main()
