# Ti-V BCC Dilute Alloy SQS Generator

This repository demonstrates the usage of the ATAT (Alloy Theoretic Automated Toolkit) software to create Special Quasirandom Structures (SQS) for a Ti-V BCC dilute alloy. SQS is a useful approach for generating alloy structures that approximate random distributions of atoms.

## Prerequisites

Before using the provided scripts, you need to obtain the Monte Carlo SQS (MCSQS) software from the ATAT website: [https://www.brown.edu/Departments/Engineering/Labs/avdw/atat/](https://www.brown.edu/Departments/Engineering/Labs/avdw/atat/)

If you use MCSQS, please cite the relevant paper: [https://www.sciencedirect.com/science/article/pii/S0364591613000540?via%3Dihub](https://www.sciencedirect.com/science/article/pii/S0364591613000540?via%3Dihub)

## Usage
```
bash runsqs.sh
```
## Script Overview

The `runsqs.sh` script automates the generation of Special Quasirandom Structures (SQS) for a Ti-V BCC dilute alloy using the ATAT (Alloy Theoretic Automated Toolkit) software. The script performs the following steps:

1. **Create Supercell:**
   - Generates a 5x5x5 supercell from a 2-atom BCC unit cell, resulting in a structure containing 250 atoms.

2. **Add V atoms:**
   - Randomly adds 2-20% of V atoms to the supercell in 2% increments.

3. **Run MCSQS:**
   - Executes the Monte Carlo SQS (MCSQS) algorithm on the supercell. The tolerance can be adjusted for more precision.

> **Note:**
> - The script utilizes the [`sqs2poscar.cpp`](https://github.com/c-niu/sqs2poscar) code from [c-niu/sqs2poscar](https://github.com/c-niu/sqs2poscar) to convert SQS output to POSCAR.
> - Additionally, a Python wrapper is used to convert POSCAR to a LAMMPS data file.

For any inquiries, please contact: mash@lanl.gov
