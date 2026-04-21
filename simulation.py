"""
Core simulation class for the Stability-Distillation Hypothesis.
"""

import numpy as np
import random
from collections import defaultdict, Counter
from .stability import get_stability
from .utils import is_trna_like, calculate_entropy


class StabilityDistillationSimulation:
    """
    Simulates the stability-distillation process for RNA sequence evolution.
    
    Key mechanisms:
    - Wet phase: Unstable sequences degrade based on thermal/acidic stress
    - Dry phase: Stable sequences serve as building blocks for ligation
    - Nucleotide supplementation: Periodic monomer input from geothermal activity
    - Harshness: Environmental stress amplifies stability differences
    """
    
    def __init__(self, 
                 initial_seq_length=6,
                 num_initial_sequences=200,
                 max_seq_length=150,
                 dry_growth_factor=1.3,
                 wet_degradation_strength=1.0,
                 ligation_base_prob=0.05,
                 supplementation_fraction=0.02,
                 initial_nucleotide_pool=1000,
                 harshness=1.0):
        """
        Initialize the simulation.
        
        Args:
            initial_seq_length: Length of randomly generated starting sequences
            num_initial_sequences: Number of initial sequences to generate
            max_seq_length: Maximum allowed sequence length
            dry_growth_factor: Multiplication factor for copy numbers during dry phase
            wet_degradation_strength: Exponent controlling degradation strength
            ligation_base_prob: Base probability for sequence ligation events
            supplementation_fraction: Fraction of initial pool to add as monomers each cycle
            initial_nucleotide_pool: Total number of ribonucleotide monomers at start
            harshness: Environmental harshness factor (amplifies stability differences)
        """
        self.max_seq_length = max_seq_length
        self.dry_growth_factor = dry_growth_factor
        self.wet_degradation_strength = wet_degradation_strength
        self.ligation_base_prob = ligation_base_prob
        self.supplementation_fraction = supplementation_fraction
        self.initial_nucleotide_pool = initial_nucleotide_pool
        self.harshness = harshness
        
        self.nucleotides = ['A', 'U', 'G', 'C']
        
        # Initialize sequence pool
        self.sequences = []
        for _ in range(num_initial_sequences):
            seq = ''.join(random.choice(self.nucleotides) 
                         for _ in range(initial_seq_length))
            self.sequences.append(seq)
        
        self.frequencies = defaultdict(int)
        for seq in self.sequences:
            self.frequencies[seq] += 1
        
        # Add initial monomer pool
        for _ in range(initial_nucleotide_pool):
            monomer = random.choice(self.nucleotides)
            self.frequencies[monomer] += 1
        
        self.current_longest_sequence = ""
        self.history = {
            'entropy': [],
            'num_unique_sequences': [],
            'avg_length': [],
            'total_copies': [],
            'max_length': [],
            'longest_sequence': [],
            'top_sequences': [],
            'trna_candidates': []
        }
    
    def get_longest_sequence(self):
        """Find the current longest RNA sequence in the pool."""
        if not self.frequencies:
            return "", 0
        
        sequences_only = [seq for seq in self.frequencies.keys() if len(seq) > 1]
        if not sequences_only:
            return "", 0
        
        longest_seq = max(sequences_only, key=len)
        return longest_seq, len(longest_seq)
    
    def add_nucleotide_supplementation(self):
        """Add fresh ribonucleotide monomers to the pool."""
        if self.supplementation_fraction <= 0:
            return 0
        
        num_monomers = max(1, int(self.initial_nucleotide_pool * self.supplementation_fraction))
        
        for _ in range(num_monomers):
            monomer = random.choice(self.nucleotides)
            self.frequencies[monomer] += 1
        
        return num_monomers
    
    def wet_phase(self):
        """Wet phase: Unstable sequences degrade."""
        new_frequencies = defaultdict(int)
        
        for seq, count in self.frequencies.items():
            stability = get_stability(seq)
            survival_prob = stability ** (self.wet_degradation_strength * self.harshness)
            survived = np.random.binomial(count, survival_prob)
            
            if survived > 0:
                new_frequencies[seq] = survived
        
        self.frequencies = new_frequencies
        
        if len(self.frequencies) == 0:
            for _ in range(50):
                seq = ''.join(random.choice(self.nucleotides) for _ in range(6))
                self.frequencies[seq] += 1
    
    def dry_phase(self):
        """Dry phase: Sequence ligation and growth."""
        total = sum(self.frequencies.values())
        if total == 0:
            return
        
        freq_ratio = {seq: count/total for seq, count in self.frequencies.items()}
        new_sequences = []
        seq_list = list(self.frequencies.keys())
        
        num_ligations = int(total * 0.3)
        
        for _ in range(num_ligations):
            if len(seq_list) < 2:
                break
            
            weights = [freq_ratio.get(seq, 0) for seq in seq_list]
            if sum(weights) == 0:
                break
            
            a = np.random.choice(seq_list, p=weights)
            b = np.random.choice(seq_list, p=weights)
            
            p_ligate = self.ligation_base_prob * (freq_ratio[a] * freq_ratio[b] * 100)
            p_ligate = min(0.8, p_ligate)
            
            if random.random() < p_ligate:
                new_seq = a + b
                if len(new_seq) <= self.max_seq_length:
                    new_sequences.append(new_seq)
        
        for seq in new_sequences:
            self.frequencies[seq] += 1
        
        for seq in list(self.frequencies.keys()):
            growth = np.random.poisson(self.dry_growth_factor - 1)
            if growth > 0:
                self.frequencies[seq] += growth
    
    def step(self, cycle_number):
        """Execute one complete wet-dry cycle."""
        self.wet_phase()
        self.add_nucleotide_supplementation()
        self.dry_phase()
        
        longest_seq, max_len = self.get_longest_sequence()
        self.current_longest_sequence = longest_seq
        
        if longest_seq:
            display_seq = longest_seq[:50] + '...' if len(longest_seq) > 50 else longest_seq
            print(f"  Cycle {cycle_number:3d}: Longest RNA = {display_seq} (length={max_len})")
    
    def run(self, num_cycles=50, record_interval=5):
        """Run the simulation for a specified number of cycles."""
        print(f"\nStarting simulation: {num_cycles} wet-dry cycles...")
        print("-" * 70)
        
        for cycle in range(num_cycles):
            self.step(cycle)
            
            if cycle % record_interval == 0 or cycle == num_cycles - 1:
                entropy = calculate_entropy(self.frequencies)
                num_unique = len(self.frequencies)
                total = sum(self.frequencies.values())
                
                seq_lengths = [len(seq) for seq in self.frequencies.keys() if len(seq) > 1]
                avg_len = np.mean(seq_lengths) if seq_lengths else 0
                _, max_len = self.get_longest_sequence()
                top_seqs = Counter(self.frequencies).most_common(5)
                
                trna_candidates = []
                for seq in list(self.frequencies.keys())[:50]:
                    if is_trna_like(seq) and self.frequencies[seq] > 1:
                        trna_candidates.append(seq)
                
                self.history['entropy'].append(entropy)
                self.history['num_unique_sequences'].append(num_unique)
                self.history['avg_length'].append(avg_len)
                self.history['total_copies'].append(total)
                self.history['max_length'].append(max_len)
                self.history['longest_sequence'].append(self.current_longest_sequence)
                self.history['top_sequences'].append(top_seqs)
                self.history['trna_candidates'].append(trna_candidates)
                
                print(f"  [Record] Cycle {cycle:3d}: Entropy={entropy:.3f} | "
                      f"Unique={num_unique:4d} | Total copies={total:6d} | "
                      f"Avg len={avg_len:.1f} | Max len={max_len} | "
                      f"tRNA candidates={len(trna_candidates)}")
                print()
        
        print("-" * 70)
        print("Simulation complete!")
