#!/usr/bin/env python3
"""
Main entry point for running the Stability-Distillation Hypothesis simulation.
"""

import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sdh_simulation import StabilityDistillationSimulation
from sdh_simulation.stability import USE_VIENNA_GLOBAL, HARSHNESS_GLOBAL
from sdh_simulation.utils import plot_results, print_summary


def run_demo(stability_method='heuristic', 
             initial_nucleotide_pool=1500,
             num_cycles=80,
             supplementation_fraction=0.02,
             ligation_prob=0.08,
             harshness=1.0):
    """
    Run the stability-distillation simulation demo.
    """
    global USE_VIENNA_GLOBAL, HARSHNESS_GLOBAL
    
    print("=" * 80)
    print("     Stability-Distillation Hypothesis Simulation")
    print("     Simulating wet-dry cycles driving RNA sequence convergence")
    print("=" * 80)
    print(f"\nPaper: https://doi.org/10.48550/arXiv.2403.17072")
    
    print(f"\nConfiguration:")
    print(f"  Stability calculation method: {stability_method.upper()}")
    print(f"  Initial nucleotide pool: {initial_nucleotide_pool} monomers")
    print(f"  Number of cycles: {num_cycles}")
    print(f"  Supplementation fraction: {supplementation_fraction * 100:.1f}% of initial pool per cycle")
    print(f"  Ligation probability (base): {ligation_prob}")
    print(f"  Environmental harshness: {harshness}")
    
    # Import ViennaRNA check
    try:
        import RNA
        VIENNA_AVAILABLE = True
    except ImportError:
        VIENNA_AVAILABLE = False
    
    if stability_method.lower() == 'vienna':
        if VIENNA_AVAILABLE:
            USE_VIENNA_GLOBAL = True
            print(f"\n  ✓ Using ViennaRNA for stability calculations")
        else:
            USE_VIENNA_GLOBAL = False
            print(f"\n  ✗ ViennaRNA requested but not available!")
            print(f"    Falling back to heuristic method.")
    else:
        USE_VIENNA_GLOBAL = False
        print(f"\n  Using heuristic stability calculation")
    
    HARSHNESS_GLOBAL = harshness
    print()
    
    sim = StabilityDistillationSimulation(
        initial_seq_length=6,
        num_initial_sequences=200,
        max_seq_length=150,
        dry_growth_factor=1.3,
        wet_degradation_strength=1.0,
        ligation_base_prob=ligation_prob,
        supplementation_fraction=supplementation_fraction,
        initial_nucleotide_pool=initial_nucleotide_pool,
        harshness=harshness
    )
    
    sim.run(num_cycles=num_cycles, record_interval=5)
    plot_results(sim)
    print_summary(sim)


if __name__ == "__main__":
    # User configurable parameters
    STABILITY_METHOD = 'heuristic'      # 'heuristic' or 'vienna'
    INITIAL_NUCLEOTIDE_POOL = 1500      # Total monomers at start
    NUM_CYCLES = 80                     # Number of wet-dry cycles
    SUPPLEMENTATION_FRACTION = 0.02     # 2% of initial pool per cycle
    LIGATION_PROB = 0.08                # Base ligation probability
    HARSHNESS = 1.0                     # Environmental harshness
    
    run_demo(
        stability_method=STABILITY_METHOD,
        initial_nucleotide_pool=INITIAL_NUCLEOTIDE_POOL,
        num_cycles=NUM_CYCLES,
        supplementation_fraction=SUPPLEMENTATION_FRACTION,
        ligation_prob=LIGATION_PROB,
        harshness=HARSHNESS
    )
