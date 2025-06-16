# CFCDBN

MATLAB implementation of **CFCDBN: A framework for modeling cross-frequency coupling directed brain networks**.

This repository provides the core code and resources used in the study:

**"Characterizing Cross-Frequency Directed Brain Networks in Adolescent Anxiety via Graph-Based Modeling"**

---

## Overview

The CFCDBN framework integrates individualized phase–amplitude coupling (PAC) features with directed brain network (DBN) modeling to analyze aberrant oscillatory communication patterns in adolescents with anxiety disorders (AD).

It supports:
- **Dynamic estimation of cross-frequency coupling (CFC)** using individualized oscillatory bands
- **Construction of directed brain networks** based on transfer entropy
- **Graph-level classification** using direction-aware graph neural networks (D-GNNs)
- **Visualization and interpretation** of discriminative brain connectivity patterns

---

## Requirements

- MATLAB R2021a or later
- [FieldTrip Toolbox](https://www.fieldtriptoolbox.org/)
- [EEGLAB Toolbox](https://sccn.ucsd.edu/eeglab/index.php)
- Signal Processing Toolbox
- Deep Learning Toolbox (for GNN training)
- Python (optional, for visualization or supplementary preprocessing)
