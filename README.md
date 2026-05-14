# Shapley-Shubik Attribution from Minimal Subsets

## Authors

*   **Pablo Martínez-Naredo**¹
*   **Raúl Mencía**¹
*   **Joao Marques-Silva**²
*   **Carlos Mencía**¹

## Affiliations

¹ *Universidad de Oviedo, Spain*  
² *ICREA & Universitat de Lleida, Spain*

---


## Abstract

This repository contains the source code, instances, and experimental results presented in our conference paper.

## Repository Structure

The project is organized into three main directories:

*   **`source/`**: Contains the Python tool execution scripts, internal classes, and examples.
    *   `ex/`: Directory with example formulas (`ex.cnf`, `exmus.cnf`, `exmcs.cnf`).
    *   `requirements.txt`: Python dependencies list.
    *   `attribution.py`: Shapley-Shubik attribution program (supports SAT-based and interval-based approaches).
    *   `coverage.py`: Script to compute the coverage of MUSes and MCSes formulas.
    *   `enumeration.py`: Tool to enumerate MUSes and MCSes from an unsatisfiable formula.
    *   `sh_common.py`: Shared core classes used by the other scripts.
*   **`instances/`**: Contains the benchmarks utilized in the paper.
    *   `list.txt`: Contains the instances utilized. These are unsatisfiable unweighted instances from the exact track of the MaxSAT Evaluation 2024, which can be downloaded from: https://maxsat-evaluations.github.io/2024/benchmarks.html
    *   `baseline/`: CNF formulas of MUSes and MCSes for the "Coverage and Importance Interval Convergence" section.
    *   `challenging/`: CNF formulas of MUSes and MCSes for the "Attribution in Harder Formulas" section.
*   **`results/`**: Comprehensive experimental output data.
    *   `baseline/`: Results for the "Coverage and Importance Interval Convergence" section.
    *   `challenging/`: Results for the "Attribution in Harder Formulas" section.
    
    *Both execution directories contain the following structure:*
    *   `coverage.txt`: Contains the coverage of the formula created with the MUSes and MCSes of each considered instance.
    *   `sat_based/`: Results from executing the SAT-based approach. Provides details on completed iterations, execution runtime, and the computed Shapley-Shubik values.
    *   `interval_based/`: Results from executing the interval-based approach. Provides details on completed iterations, execution runtime, and the values obtained for each metric (conservative, uniform, exponential, low, and up).

---

## Setup and Installation

Follow these steps to configure a local Python virtual environment using your web-downloaded repository files.

### 1. Create a Virtual Environment
Navigate to your project directory via terminal and run:
```bash
python -m venv venv
```

### 2. Activate the Environment
*   **Linux/macOS:**
    ```bash
    source venv/bin/activate
    ```
*   **Windows (Command Prompt):**
    ```cmd
    venv\Scripts\activate
    ```

### 3. Install Dependencies
Install the required packages listed in `requirements.txt`:
```bash
pip install -r source/requirements.txt
```

---

## Usage Instructions

To view execution parameters, detailed usage instructions, and practical examples, run the scripts directly from the `source` directory:

*   **Shapley-Shubik Attribution:**
    ```bash
    python source/attribution.py
    ```
*   **Coverage Calculation:**
    ```bash
    python source/coverage.py
    ```
*   **MUS/MCS Enumeration:**
    ```bash
    python source/enumeration.py
    ```

---

## Licensing

This repository features a multi-licensing structure to separate software logic permissions from research data and output usage:

*   **Software Code (`source/`)**: All Python scripts, core utilities, and execution frameworks are distributed under the terms of the **MIT License**. See the root [LICENSE](LICENSE) file for details.
*   **Benchmarks and Inputs (`instances/`)**: The benchmark formulas and instance lists are made available under the **Creative Commons Attribution 4.0 International License (CC BY 4.0)**. See [instances/LICENSE](instances/LICENSE) for details.
*   **Experimental Outputs (`results/`)**: All generated execution logs, metrics, and computed Shapley-Shubik values are also distributed under the **Creative Commons Attribution 4.0 International License (CC BY 4.0)**. See [results/LICENSE](results/LICENSE) for details.

