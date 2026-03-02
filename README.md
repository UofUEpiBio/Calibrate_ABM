# 📊 Calibrate_ABM


This repository contains R scripts for **epidemic simulation, parameter calibration, prediction, and evaluation** using agent-based epidemic models.

The workflow integrates:

- 🧮 Epidemic simulation via `epiworldRcalibrate`  
- 🎯 Approximate Bayesian Computation (ABC)  
- 🤖 BiLSTM-based parameter calibration  
- 📈 Epidemic curve reconstruction and statistical evaluation  

All simulations and calibration procedures are implemented in **R**.

---

# 🧠 BiLSTM Calibration Architecture

The BiLSTM model used for parameter estimation is illustrated below:
![Uploading bilstm_model.drawio (3).png…]()



## Architecture Overview

**Inputs**
- Incidence time series (T = 60)
- Additional epidemiological inputs:
  - Population size (n)
  - Recovery rate

**Model Structure**
- 3 stacked Bidirectional LSTM layers  
- 160 hidden units  
- Dropout = 0.5  
- Forward and backward hidden states concatenated  
- Epidemiological covariates concatenated with LSTM output  

**Fully Connected Layers**
- Linear: 322 → 64 (ReLU activation)  
- Linear: 64 → 3 outputs  

**Outputs**
- Transmission rate (sigmoid activation)  
- Contact rate (softplus activation)  
- \(R_0\) (softplus activation)  

An epidemiological constraint is incorporated through:

\[
R_0 = \frac{\text{contact rate} \times \text{transmission rate}}{\text{recovery rate}}
\]

This relationship is used as a penalty term during training to ensure epidemiological consistency.

---

# 📂 Repository Structure

## Core Scripts

### `00-params.R`
Generates parameter sets used for epidemic simulations.

### `01a-bilstm.R`
Implements the BiLSTM-based calibration workflow using epidemic trajectories generated via `epiworldRcalibrate`.

### `01b-abc.R`
Performs Approximate Bayesian Computation (ABC) for parameter estimation from simulated epidemic curves.

### `02-abc-bilstm-prediction.R`
Combines ABC and BiLSTM calibration outputs and produces final predicted parameter values.

### `03-epicurves-stats.R`
Generates epidemic curves from calibrated parameters and compares reconstruction accuracy across methods.

### `04-parameter-stats.R`
Evaluates calibration performance:
- Bias computation  
- Estimated vs. true parameter comparison  
- Statistical summaries  
- Visualization of parameter distributions  

---

# 🔬 Methodological Workflow

1. Generate simulation parameters  
2. Simulate epidemic trajectories using `epiworldRcalibrate`  
3. Calibrate parameters using:
   - Approximate Bayesian Computation (ABC)
   - BiLSTM regression model  
4. Compare predicted vs. true parameters  
5. Evaluate epidemic curve reconstruction accuracy  

---

# 🛠 Requirements

- R (≥ 4.x recommended)  
- `epiworldRcalibrate`  
- Required R packages for simulation, machine learning, and visualization  

Example installation:

```r
install.packages("epiworldRcalibrate")
```

---

# 🚀 Running the Pipeline

Run scripts sequentially:

```bash
Rscript 00-params.R
Rscript 01a-bilstm.R
Rscript 01b-abc.R
Rscript 02-abc-bilstm-prediction.R
Rscript 03-epicurves-stats.R
Rscript 04-parameter-stats.R
```

---

# 📈 Outputs

The pipeline produces:

- Calibrated parameter estimates  
- Predicted epidemic trajectories  
- Bias and performance summaries  
- Comparative epidemic curve visualizations  

---

# 📌 Notes

- Scripts should be executed in order for reproducibility.
- Ensure consistent random seeds when comparing calibration methods.
- Results may vary depending on simulation size and parameter ranges.
