import random
import numpy as np
import pandas as pd
import itertools
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt


NUM_PATHOGENS = 1000                # Number Pathogens
GENERATIONS = 100
MUTATION_DETECTOR = 0.05            # Probability self bit flips
MUTATION_PATHOGEN = 0.5             # Probability self bit flips
THETA = 0.5                         # Immune pressure parameter

def random_bitstring(n):
    return ''.join(random.choice('01') for _ in range(n))


# Self vs non-self similarity rule
def r_contiguous_match(s1, s2, r):
    """
    Return True if s1 and s2 share at least r contiguous matching symbols
    """
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
class Detector:
    def __init__(self, bitstring, r, MUTATION_DETECTOR):
        self.bitstring = bitstring
        self.r = r
        self.MUTATION_DETECTOR = MUTATION_DETECTOR
        self.fitness = 1.0

    def detect(self, pathogen):
        # Made deterministic: directly return if a match is found, no longer using DETECTION_SUCCESS
        return r_contiguous_match(self.bitstring, pathogen.bitstring, self.r)
class Pathogen:

    def __init__(self, bitstring, THETA, MUTATION_PATHOGEN, ALPHA, l): # Added ALPHA parameter
        self.bitstring = bitstring
        self.THETA = THETA
        self.MUTATION_PATHOGEN = MUTATION_PATHOGEN
        self.ALPHA = ALPHA  # Weight for protein score (e.g., 0.2 means 20% weight on structure, 80% on evasion)
        self.l = l
        self.fitness = 0.0

    def protein_score(self):
        """
        Preserves pathogen structure: more 1s = better. 
        protein score = number of "1" / Number of bits 
        """
        return self.bitstring.count('1') / self.l

    def detection_score(self, detectors):
        """
        Less recognisable = higher score
        detection score = number of detections
        """
        if not detectors:
            return 0.0
        
        detections = sum(1 for d in detectors if d.detect(self))
        return detections / len(detectors)

    def evaluate_fitness(self, detectors):
        D = len(detectors)
        detection_ratio = self.detection_score(detectors)
        probability_evasion = np.exp(-self.THETA * (D / 1000.0)) # Scaled to prevent immediate absolute zero
        if probability_evasion > detection_ratio:
            evasion_component = probability_evasion
            protein_component = self.protein_score()
            # Blended fitness formula
            self.fitness = ((1 - self.ALPHA) * evasion_component) + (self.ALPHA * protein_component)
        else:
            self.fitness = 0.0 # Caught!
            
    # Mutation
    def mutate(self):
        self.bitstring = mutate_bitstring(self.bitstring, self.MUTATION_PATHOGEN)
    
    # Crossover
    def crossover(self, other):
        point = random.randint(1, self.l - 1)
        child_bitstring = (self.bitstring[:point] + other.bitstring[point:])

        # Inherit parameters from either parent
        theta = (self.THETA + other.THETA) / 2
        mutation_rate = (self.MUTATION_PATHOGEN + other.MUTATION_PATHOGEN) / 2
        
        # Average or inherit the alpha weight parameter
        alpha = (self.ALPHA + other.ALPHA) / 2

        return Pathogen(child_bitstring, theta, mutation_rate, alpha, self.l)
      
    def __repr__(self):
        return (f"Pathogen(bitstring='{self.bitstring}', fitness={self.fitness:.3f})")
# Negative Selection Algorithm

def generate_detectors(NUM_DETECTORS, l, R):

    detectors = []
    SELF = "0" * l   
       
    while len(detectors) < NUM_DETECTORS:
        # Random detectors are generated with a random R and Mutation Rate
        candidate = random_bitstring(l)
        # Negative Selection: detectors recognizing "SELF" are removed
        if not r_contiguous_match(candidate, SELF, R):
            detectors.append(Detector(candidate, R, MUTATION_DETECTOR))
    return detectors


def generate_pathogens(NUM_PATHOGENS, l, ALPHA):
    pathogens = []
    for _ in range(NUM_PATHOGENS):

        # start pathogen near all-ones
        s = ''.join(
            random.choice(['1', '1', '1', '0'])
            for _ in range(l)
        )
        pathogens.append(Pathogen(s, THETA, MUTATION_PATHOGEN, ALPHA, l))

    return pathogens

# Self Evolution (Genetic Algorithm) Pathogen

def evolve_pathogens(pathogens, detectors):

    # evaluate fitness
    for p in pathogens:
        p.evaluate_fitness(detectors)

    # sort pathogens based on fitness 
    pathogens.sort(key=lambda x: x.fitness, reverse=True) 

    # Elitism: it keeps only 50% of the best. (Natural Selection)
    survivors = pathogens[:NUM_PATHOGENS // 2]
    children = []

    while len(children) < NUM_PATHOGENS // 2:
        
        # Select random parent2
        parent1 = random.choice(survivors)
        parent2 = random.choice(survivors)

        # Crossover between parents 
        child = parent1.crossover(parent2)
        
        # Mutation
        child.mutate()

        children.append(child)

    return survivors + children


#  Self Evolution (Genetic Algorithm) Detectors

def evolve_detectors(detectors, pathogens):
    if not detectors: return []
    
    alive_detectors = []
    sample_size = min(10, len(pathogens))
    test_pathogens = random.sample(pathogens, sample_size)
    
    for d in detectors:
        # A detector survives if it catches at least one pathogen in the sample space
        if any(d.detect(p) for p in test_pathogens):
            alive_detectors.append(d)
            
    return alive_detectors
    
# Parameters for testing
l_values = [5, 10, 20]            
r_values = [2, 3, 8]
alpha_values = [0, 0.2, 0.5, 0.8, 1]
num_detectors_grid = [500, 2000] # Representative small vs. large immune systems
experimental_results = []
REPETITIONS = 3  # Run each configuration 3 times for statistical validation


def run_simulation(l, n_detectors, r, alpha, GENERATIONS):
    """
    Wraps your simulation logic. Returns the final performance metrics.
    Note: Ensure your Pathogen/Detector classes read these dynamic variables 
    instead of relying on hardcoded global variables!
    """
    
    detectors = generate_detectors(n_detectors, l, r)
    pathogens = generate_pathogens(NUM_PATHOGENS, l, alpha)
    
    for generation in range(GENERATIONS):
        pathogens = evolve_pathogens(pathogens, detectors)
        detectors = evolve_detectors(detectors, pathogens)
        if len(detectors) == 0:
            break  # Break early if immune system experiences absolute collapse
            
    # Calculate metric extractions at the final step
    if len(pathogens) == 0:
        final_avg_fitness = 0.0
        final_protein_score = 0.0
    else:
        final_avg_fitness = np.mean([p.fitness for p in pathogens])
        final_protein_score = np.mean([p.protein_score() for p in pathogens])
        
    detector_survival_rate = len(detectors) / n_detectors
    
    return final_avg_fitness, final_protein_score, detector_survival_rate

if __name__ == "__main__":
    experimental_results = []
    print("Beginning automated evaluation sweep...")

    # Iterate through every mathematical permutation of parameters
    for l, n_det, r, alpha in itertools.product(l_values, num_detectors_grid, r_values, alpha_values):
        
        # Enforce mathematical bound: threshold rule cannot exceed total string bits
        if r > l:
            continue
            
        print(f"Processing Config -> l: {l}, N: {n_det}, r: {r}, alpha: {alpha}")
        
        # Lists to gather statistical variance across repetitions
        rep_fitness = []
        rep_proteins = []
        rep_survival = []
        
        for rep in range(REPETITIONS):
            fit, prot, survival = run_simulation(l, n_det, r, alpha, GENERATIONS)
            rep_fitness.append(fit)
            rep_proteins.append(prot)
            rep_survival.append(survival)
            
        # Store clean means and statistical variations for report presentation
        experimental_results.append({
            "String Length (l)": l,
            "Initial Detectors (N)": n_det,
            "Matching Threshold (r)": r,
            "Structural Weight (alpha)": alpha,
            "Mean Pathogen Fitness": np.mean(rep_fitness),
            "Fitness StdDev": np.std(rep_fitness),
            "Mean Protein Score": np.mean(rep_proteins),
            "Protein StdDev": np.std(rep_proteins),
            "Mean Detector Survival Rate": np.mean(rep_survival)
        })

    # Save out structured spreadsheet data directly matching report criteria
    df = pd.DataFrame(experimental_results)
    df.to_csv("comprehensive_sweep_results.csv", index=False)
    print("\nSimulation process complete. Dataset exported to 'comprehensive_sweep_results.csv'.")
    
    # Showcase data head sample format
    print(df.head())