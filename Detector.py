import random
import numpy as np
import pandas as pd
import itertools
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
from scipy.stats import entropy

random.seed(42)
np.random.seed(42)


REPETITIONS = 3
GENERATIONS = 100
MUTATION_DETECTOR = 0.02
NUM_DETECTORS = 1000
ALPHA = 0.8
BETA = 0.2
l_values = [5, 10, 20]
r_values = [2, 3, 8]
PATHOGEN_COUNTS = [500, 1000, 2000]   


def random_bitstring(n):
    return ''.join(random.choice('01') for _ in range(n))

def r_contiguous_match(s1, s2, r):
    count = 0
    n = min(len(s1), len(s2))
    for i in range(n):
        if s1[i] == s2[i]:
            count += 1
            if count >= r:
                return True
        else:
            count = 0
    return False

def longest_contiguous_match(s1, s2):
    max_len = 0
    cur = 0
    n = min(len(s1), len(s2))
    for i in range(n):
        if s1[i] == s2[i]:
            cur += 1
            if cur > max_len:
                max_len = cur
        else:
            cur = 0
    return max_len

def mutate_bitstring(s, mutation_rate):
    bits = list(s)
    for i in range(len(bits)):
        if random.random() < mutation_rate:
            bits[i] = '1' if bits[i] == '0' else '0'
    return ''.join(bits)

def shannon(bitstrings):
    counts = Counter(bitstrings)
    return entropy(list(counts.values()))

def mean_hamming(bitstrings):
    if len(bitstrings) < 2:
        return 0
    from itertools import combinations
    distances = [sum(a != b for a, b in zip(x, y)) for x, y in combinations(bitstrings, 2)]
    return np.mean(distances)

def hamming_distance(s1, s2):
    return sum(a != b for a, b in zip(s1, s2))

def is_alternating(bitstring):
    return all(bitstring[i] != bitstring[i-1] for i in range(1, len(bitstring)))

def generate_harder_pathogen(l):
    min_zeros = min(3, l-1)
    max_zeros = min(5, l-1)
    if min_zeros > max_zeros:
        raise ValueError(f"Cannot generate harder pathogen for l={l}.")
    while True:
        num_zeros = random.randint(min_zeros, max_zeros)
        num_ones = l - num_zeros
        bits = ['0'] * num_zeros + ['1'] * num_ones
        random.shuffle(bits)
        pathogen = ''.join(bits)
        if is_alternating(pathogen) or pathogen == '0' * l:
            continue
        return pathogen


class Pathogen:
    def __init__(self, bitstring):
        self.bitstring = bitstring

class Detector:
    def __init__(self, bitstring, r):
        self.bitstring = bitstring
        self.r = r
        self.fitness = 0.0

    def detect(self, pathogen):
        return r_contiguous_match(self.bitstring, pathogen.bitstring, self.r)

    def mutate(self, rate=MUTATION_DETECTOR):
        self.bitstring = mutate_bitstring(self.bitstring, rate)

    def crossover(self, other):
        point = random.randint(1, len(self.bitstring)-1)
        child_bitstring = self.bitstring[:point] + other.bitstring[point:]
        return Detector(child_bitstring, self.r)


def generate_detectors(n, l, r, self_string='0'):
    detectors = []
    while len(detectors) < n:
        candidate = random_bitstring(l)
        if not r_contiguous_match(candidate, self_string, r):
            detectors.append(Detector(candidate, r))
    return detectors


def detection_score(detector, pathogens, l):
    if not pathogens:
        return 0.0
    total = sum(longest_contiguous_match(detector.bitstring, p.bitstring) for p in pathogens)
    return total / (l * len(pathogens))

def diversity_score(detector, population, l):
    if len(population) <= 1:
        return 0.0
    distances = []
    for other in population:
        if other is detector:
            continue
        hd = hamming_distance(detector.bitstring, other.bitstring)
        distances.append(hd / l)
    return np.mean(distances) if distances else 0.0

def evaluate_detectors(detectors, pathogens, alpha=ALPHA, beta=BETA):
    if not detectors:
        return
    l = len(detectors[0].bitstring)
    for d in detectors:
        det = detection_score(d, pathogens, l)
        div = diversity_score(d, detectors, l)
        d.fitness = alpha * det + beta * div


def evolve_detectors(detectors, pathogens, self_string):
    evaluate_detectors(detectors, pathogens)
    detectors.sort(key=lambda d: d.fitness, reverse=True)
    survivors = detectors[:len(detectors)//2]
    children = []
    while len(children) < len(detectors)//2:
        p1, p2 = random.sample(survivors, 2)
        child = p1.crossover(p2)
        child.mutate()
        if not r_contiguous_match(child.bitstring, self_string, child.r):
            children.append(child)
        else:
            candidate = random_bitstring(len(self_string))
            while r_contiguous_match(candidate, self_string, child.r):
                candidate = random_bitstring(len(self_string))
            child.bitstring = candidate
            children.append(child)

    new_population = survivors + children
    evaluate_detectors(new_population, pathogens)
    return new_population


def run_experiment_single(l, r, rep, pathogen_type, num_pathogens, generations=GENERATIONS):
    # Generate a list of pathogen bitstrings (all static)
    if pathogen_type == 'easy':
        pathogen_bitstrings = ['1' * l for _ in range(num_pathogens)]
    else:  
        pathogen_bitstrings = [generate_harder_pathogen(l) for _ in range(num_pathogens)]

    pathogens = [Pathogen(bs) for bs in pathogen_bitstrings]
    self_string = '0' * l

    detectors = generate_detectors(NUM_DETECTORS, l, r, self_string)

    results = []

    for gen in range(generations):
        detectors = evolve_detectors(detectors, pathogens, self_string)

        bitstrings = [d.bitstring for d in detectors]
        mean_detection = np.mean([
            longest_contiguous_match(d.bitstring, pathogen_bitstrings[0]) / l
            for d in detectors
        ])

        detected_count = sum(1 for p in pathogens if any(d.detect(p) for d in detectors))
        coverage = detected_count / len(pathogens)

        record = {
            'l': l,
            'r': r,
            'N': num_pathogens,
            'PathogenType': pathogen_type,
            'Repetition': rep,
            'Generation': gen,
            'Mean Fitness': np.mean([d.fitness for d in detectors]),
            'Mean Detection': mean_detection,
            'Diversity Shannon': shannon(bitstrings),
            'Diversity Hamming': mean_hamming(bitstrings),
            'Coverage': coverage,
        }
        results.append(record)

    return results



if __name__ == "__main__":
    all_results = []

    for l, r in itertools.product(l_values, r_values):
        if r > l:
            continue
        for num_pathogens in PATHOGEN_COUNTS:
            for rep in range(REPETITIONS):
                for ptype in ['easy', 'hard']:
                    print(f"Running l={l}, r={r}, N={num_pathogens}, rep={rep}, pathogen={ptype}")
                    res = run_experiment_single(l, r, rep, ptype, num_pathogens)
                    all_results.extend(res)

    df = pd.DataFrame(all_results)
    df.to_csv("detector_evolution_multi_pathogens.csv", index=False)

    sns.set_theme(style="whitegrid")

    plot_df = df.groupby(['l', 'r', 'N', 'PathogenType', 'Generation']).mean(numeric_only=True).reset_index()

    metrics = {
        'Mean Fitness': 'Fitness (detection rate)',
        'Mean Detection': 'Detection (avg longest match)',
        'Diversity Shannon': 'Shannon Diversity',
        'Diversity Hamming': 'Hamming Distance',

    }

    for ptype in ['easy', 'hard']:
        data_ptype = plot_df[plot_df['PathogenType'] == ptype]
        if data_ptype.empty:
            continue

        for N_val in PATHOGEN_COUNTS:
            data_N = data_ptype[data_ptype['N'] == N_val]
            if data_N.empty:
                continue

            fig, axes = plt.subplots(2, 2, figsize=(14, 10))
            fig.suptitle(f'{ptype.upper()} PATHOGEN-N={N_val}', fontsize=16)

            ax_list = axes.flatten()

            for idx, (metric, ylabel) in enumerate(metrics.items()):
                ax = ax_list[idx]
                for (l_val, r_val), group in data_N.groupby(['l', 'r']):
                    ax.plot(group['Generation'], group[metric],
                            label=f'l={l_val}, r={r_val}', linewidth=2)
                ax.set_xlabel('Generation')
                ax.set_ylabel(ylabel)
                ax.legend(loc='best')
                ax.grid(True)

            for j in range(len(metrics), 4):
                fig.delaxes(ax_list[j])

            plt.tight_layout()
            plt.show()