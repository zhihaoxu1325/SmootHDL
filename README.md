# SmootHDL

**SmootHDL** is an automated testing tool for FPGA synthesis compilers. It takes a Simulink model, mutates it to create functionally equivalent variants, generates HDL, and compares RTL vs post-synthesis behavior to surface compiler bugs.

## Requirements

- Two machines with different OSes (Windows + Linux).
- MATLAB 2023a or higher.
- Python 3.11 or higher.
- Vivado 2023.1 or higher.

## Project Structure

- `CPS/`: Simulink model generation and HDL generation scripts.
- `Mutant/`: EMI mutation pipeline (TPG/BDS), coverage experiments, and preprocessing.
- `Linux_test/`: RTL/netlist simulation and equivalence checks.

## Workflow Overview

1. **Model Generation (CPS)**  
   Generate seed Simulink models using `CPS/sgtest.m` and tune the generator using `CPS/cfg.m`.
2. **Mutation (Mutant/EMI)**  
   Collect coverage and create EMI mutants with TPG/BDS support.
3. **HDL Generation (CPS)**  
   Convert mutated models to HDL using `CPS/HDL_generation.m`.
4. **Equivalence Testing (Linux_test)**  
   Run RTL and post-synthesis simulations and compare results to detect bugs.

## Components (TPG / BDS / ETC)

- **TPG (Test-Program Generation)**  
  Stateflow structural mutations (state duplication, path duplication, transition expansion).  
  Config: `Mutant/+emi/cfg.m` (`TPG_ENABLE`, `TPG_MUTATIONS_PER_CHART`, `TPG_SELECTOR_NAME`).

- **BDS (Bayesian Diversity Selection)**  
  Scores candidate mutants using Stateflow graph features and HDL “pressure” metrics, then selects promising candidates (greedy or SMBO).  
  Config: `Mutant/+emi/cfg.m` (`BDS_ENABLE`, `BDS_NUM_CANDIDATES`, `BDS_SELECT_COUNT`, `BDS_USE_SMBO`, `BDS_RF_TREES`, `BDS_RF_MIN_LEAF`, `BDS_HDL_SEARCH_DIRS`).

- **ETC (Equivalent Test Check)**  
  Compares RTL and post-synthesis outputs and logs mismatches.  
  Entry point: `Linux_test/test_main.py` (uses `Linux_test/etc_checker.py`).

## End-to-End Runbook

> Paths in scripts are environment-specific. Please update them to match your machine layout.

### 1) Generate Models (Windows / MATLAB)

```matlab
sgtest
```

This generates seed models. Configure generator settings in `CPS/cfg.m`.

### 2) Mutation (Windows / MATLAB)

Copy generated models from `CPS/repoetsneo` to `Mutant/reproduce/samplecorpus`, then:

```bash
covexp.covcollect
emi.go
```

This preprocesses models and generates EMI mutants (TPG/BDS are driven by `Mutant/+emi/cfg.m`).

### 3) HDL Generation (Windows / MATLAB)

Run `CPS/HDL_generation.m` and update any hard-coded paths to the mutant corpus at
`Mutant/reproduce/samplecorpus`.

If you enable automated upload, configure FTP settings in `CPS/remote.m`.

### 4) Equivalence Testing (Linux)

On the Linux machine, run:

```bash
python test_main.py
```

The script performs RTL simulation, synthesis, netlist simulation, and ETC comparison.
Issues are logged to an `error_log` file next to the generated HDL.

## Notes

- The flow spans multiple tools and machines. If you want full automation, wrap the
  above steps into a single pipeline script with your local paths.
- For debugging, start with a small number of models and mutants and verify each stage
  (mutation → HDL → RTL simulation → netlist simulation → ETC).
