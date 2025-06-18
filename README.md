# CFCDBN

CFCDBN: A Framework for Modeling Cross-Frequency Coupling Directed Brain Networks with Open-source Implementation

This repository provides the core code and resources used in the study:

**"CFCDBN: Personalized Directional Brain Network Modeling of Cross-Frequency Coupling Alterations in Adolescent Anxiety Disorders"**

---

## Overview

The **CFCDBN** framework integrates individualized phase–amplitude coupling (PAC) features with directed brain network (DBN) modeling to investigate aberrant oscillatory communication patterns in adolescents with anxiety disorders (AD).

This repository provides the official implementation of the CFCDBN framework for constructing dynamic, cross-frequency, and directionally-organized brain networks from EEG data. The framework combines personalized PAC estimation with causal interaction modeling to generate interpretable, subject-specific directed brain networks, supporting both mechanistic analysis and automated classification of AD-related neural dynamics.

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
- [**MIPAC Toolbox**](https://github.com/sccn/PACTools) — for robust PAC estimation  
- [**Transfer Entropy (TE) Toolbox**](https://github.com/trentool/TRENTOOL3) — for directional interaction inference  

### Python (optional, for visualization or supplementary preprocessing)
- Python 3.8+  
- [`fooof`](https://fooof-tools.github.io/fooof/) — Fitting Oscillations and One-Over-F (for separating aperiodic vs. oscillatory activity)  
- `numpy`, `scipy`, `matplotlib`, `networkx`  
- `torch` (PyTorch) — for training the GNN model

> **Note**: All core signal processing, PAC computation, and network modeling can be performed in MATLAB. Python scripts are provided for optional advanced visualization and deep learning model training.
