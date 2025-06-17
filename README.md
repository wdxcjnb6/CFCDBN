# CFCDBN

CFCDBN: A Framework for Modeling Cross-Frequency Coupling Directed Brain Networks with Open-source Implementation

This repository provides the core code and resources used in the study:

**"Characterizing Cross-Frequency Directed Brain Networks in Adolescent Anxiety via Graph-Based Modeling"**

---

## Overview

The CFCDBN framework integrates individualized phase–amplitude coupling (PAC) features with directed brain network (DBN) modeling to analyze aberrant oscillatory communication patterns in adolescents with anxiety disorders (AD).
This repository provides the official implementation of the CFCDBN framework, designed to model dynamic, cross-frequency, and directionally-organized brain networks based on EEG data. CFCDBN integrates personalized phase–amplitude coupling (PAC) analysis with causal information flow estimation, enabling the construction of interpretable and subject-specific directed brain networks.

It supports:
- **Dynamic estimation of cross-frequency coupling (CFC)** using individualized oscillatory bands
- **Construction of directed brain networks** based on transfer entropy
- **Graph-level classification** using direction-aware graph neural networks 
- **Visualization and interpretation** of discriminative brain connectivity patterns
---

## Dependencies

### MATLAB
- **MATLAB R2021a** or later  
- [**FieldTrip Toolbox**](https://www.fieldtriptoolbox.org/) — for source modeling and connectivity estimation  
- [**EEGLAB Toolbox**](https://sccn.ucsd.edu/eeglab/index.php) — for EEG preprocessing  
- **Signal Processing Toolbox** — for time–frequency and filtering operations  
- [**MIPAC Toolbox**](https://github.com/TNTLFreiburg/MIPAC) — for robust PAC estimation  
- [**Transfer Entropy (TE) Toolbox**](https://github.com/Lobachevskyy/TE-CausalityToolbox) — for directional interaction inference  

### Python (optional, for visualization or supplementary preprocessing)
- Python 3.8+  
- [`fooof`](https://fooof-tools.github.io/fooof/) — Fitting Oscillations and One-Over-F (for separating aperiodic vs. oscillatory activity)  
- `numpy`, `scipy`, `matplotlib`, `networkx`  
- `torch` (PyTorch) — for training the GNN model

> **Note**: All core signal processing, PAC computation, and network modeling can be performed in MATLAB. Python scripts are provided for optional advanced visualization and deep learning model training.
