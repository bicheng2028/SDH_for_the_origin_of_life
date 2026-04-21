# SDH_for_the_origin_of_life
simulation for the origin of life
# SDH for the Origin of Life

**Stability-Distillation Hypothesis Simulation**

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![arXiv](https://img.shields.io/badge/arXiv-2403.17072-b31b1b.svg)](https://doi.org/10.48550/arXiv.2403.17072)

## Overview

This repository contains a computational simulation of the **Stability-Distillation Hypothesis (SDH)** for the origin of life, as proposed in the paper by Cheng Bi.

**Paper:** [Stability-distillation hypothesis for the origin of life](https://doi.org/10.48550/arXiv.2403.17072) (arXiv:2403.17072)

### Abstract from the Paper

> Through in-depth thinking and reasoning about the conditions required for cells to maintain unchanged material distribution, it is concluded that life metabolic reactions require high information content. However, the self-replication of a few macromolecules such as RNA and DNA can only retain information, but cannot generate information. This paper proposes that life originated through a "stability-distillation" process, which continuously changes the molecular frequency distribution in the system through periodic changes in the environment and generates information from a consequentialist perspective.

### Key Hypothesis

The hypothesis posits that life originated not from a few self-replicating macromolecules (as in the RNA World hypothesis), but from a mixture of RNA, DNA, proteins, and lipids in **terrestrial hot springs with periodic wet-dry cycles**. Through repeated cycles of synthesis and degradation:

1. **Stable RNA sequences** are selectively retained
2. Their **frequencies increase** in the population
3. They serve as **"building blocks"** for ligation
4. **Longer, functional structures** (e.g., tRNA-like molecules) emerge

## Key Mechanisms

| Mechanism | Description | Simulation Implementation |
|-----------|-------------|---------------------------|
| **Wet Phase** | Unstable RNA sequences degrade. Survival depends on structural stability. | `survival_prob = stability ** (harshness)` |
| **Dry Phase** | Frequent stable sequences ligate to form longer RNAs | `ligation_prob ∝ freq(a) × freq(b)` |
| **Nucleotide Supplementation** | Periodic monomer input from geothermal activity | Fixed fraction of initial pool added each cycle |
| **Environmental Harshness** | Amplifies stability differences (hot, acidic conditions) | `stability = 0.5 + (raw - 0.5) × harshness` |

## Simulation Outputs

| Output | Meaning |
|--------|---------|
| **Entropy decrease** | Information accumulation (negative entropy feeding) |
| **Maximum length growth** | Distillation effect - longer RNAs emerge from stable short sequences |
| **tRNA-like structure detection** | Emergence of functional molecules with ≥3 stem-loop regions |
| **Sequence frequency distribution** | Convergence toward stable, functional sequences |

## Requirements

- Python 3.8+
- numpy
- matplotlib
- ViennaRNA (optional, for accurate MFE-based stability calculation)

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/SDH_for_the_origin_of_life.git
cd SDH_for_the_origin_of_life

# Install dependencies
pip install -r requirements.txt

# (Optional) Install ViennaRNA for accurate stability calculation
pip install viennarna

## Quick Start
python scripts/run_simulation.py

## Python Script
from sdh_simulation import run_demo

# Run with default parameters
run_demo(
    stability_method='heuristic',  # 'heuristic' or 'vienna'
    initial_nucleotide_pool=1500,  # Starting monomers
    num_cycles=80,                  # Wet-dry cycles
    supplementation_fraction=0.02, # 2% of initial pool per cycle
    ligation_prob=0.08,            # Base ligation probability
    harshness=1.0                  # Environmental stress factor
)
## Example Output
================================================================================
     Stability-Distillation Hypothesis Simulation
     Simulating wet-dry cycles driving RNA sequence convergence
================================================================================

Configuration:
  Stability calculation method: HEURISTIC
  Initial nucleotide pool: 1500 monomers
  Number of cycles: 80
  Supplementation fraction: 2.0% of initial pool per cycle
  Ligation probability (base): 0.08
  Environmental harshness: 1.0

Simulation initialized
  Initial sequence count: 200
  Unique sequences: 200
  Initial nucleotide pool: 1500
  Environmental harshness: 1.0

Starting simulation: 80 wet-dry cycles...
----------------------------------------------------------------------
  Cycle   0: Longest RNA = AGGCCU (length=6)
  [Record] Cycle   0: Entropy=6.432 | Unique= 245 | Total copies=  2100 | Max len=6 | tRNA candidates=0

  Cycle   5: Longest RNA = AGGCCUAGGC... (length=12)
  [Record] Cycle   5: Entropy=5.123 | Unique= 189 | Total copies=  3450 | Max len=12 | tRNA candidates=0

  ...

  Cycle  75: Longest RNA = CGAUUCGAAUUCGAAUUCGAAUUCG... (length=64)
  [Record] Cycle  75: Entropy=3.234 | Unique=  78 | Total copies=  8900 | Max len=64 | tRNA candidates=2

----------------------------------------------------------------------
Simulation complete!

================================================================================
Top 15 Most Frequent Sequences
================================================================================
 1. [  456] CGAUUCGAAUUCGAAUUCGAAUUCG... (stability=0.723) ★ tRNA-like
 2. [  321] AGGCCUAGGCCUAGGCCU... (stability=0.687)
 3. [  234] UUCGAAUUCGAAUUCGAA... (stability=0.654)

================================================================================
Discovered tRNA-like Structures (2 sequences)
================================================================================
  Sequence: CGAUUCGAAUUCGAAUUCGAAUUCGAAUUCG...
  Length: 64 nt, Stability: 0.723
  Predicted structure: ((((....))))...((((....))))
  Minimum free energy: -23.4 kcal/mol


