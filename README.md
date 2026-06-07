# Identifying Monetary policy shocks using a Vector Autoregression (VAR) model

## 📌 Project Overview
This project builds a reduced-form and structural **Vector Autoregression (VAR)** system to analyse the dynamic transmission channels of macroeconomic shocks. Using 40 years of quarterly United States economic data (1971 Q1 – 2010 Q4), the model isolates the impact of an unexpected contractionary monetary policy shock, defined as a 1 percentage point increase in the Federal Funds Rate (FFR).

The system tracks 5 key endogenous variables to map out market interactions:
* **Inflation (CPI)**
* **Unemployment Rate**
* **Real GDP Growth Rate**
* **Federal Funds Rate (FFR)**
* **M2 Money Supply Growth**

---

## 🛠️ Core Skills & Technologies Used
* **Languages:** R / Python (Statsmodels)
* **Econometric Models:** Vector Autoregression (VAR), Cholesky Decomposition (Recursive Identification)
* **Statistical Testing:** Augmented Dickey-Fuller (ADF) stationarity test, Akaike Information Criterion (AIC), Bayesian Information Criterion (BIC), Portmanteau Serial Correlation test, Jarque-Bera normality test
* **Simulation:** Bootstrapping (500 runs) for 95% Confidence Interval generation

---

## 📈 Key Findings & Structural Insights

### 1. Model Specifications & Diagnostics
* **Lag Selection:** Evaluated system lag dynamics over 8 quarters. A parsimonious **VAR(2)** model was adopted using BIC criteria to avoid overfitting.
* **Stability:** Validated system equilibrium with all eigenvalues falling strictly inside the unit circle (Maximum Modulus = 0.947).
* **Policy Function:** The FFR equation achieved an adjusted R² = 0.932, proving strong empirical validation for a Taylor-rule type reaction function.

### 2. Identification Schemes & Robustness
The project tested the robustness of the economic transmission mechanism by contrasting two distinct structural ordering constraints:
* **Scheme 1 (Baseline Cholesky):** `Inflation ➔ Unemployment ➔ GDP Growth ➔ FFR ➔ M2 Growth`. This assumes real economy sticky-price rigidities where demand does not react instantly to an overnight interest rate change.
* **Scheme 2 (Alternative Monetarist):** `Inflation ➔ Unemployment ➔ GDP Growth ➔ M2 Growth ➔ FFR`. This assumes money supply is determined before interest rates within the quarter.

### 3. Impulse Response Functions (IRFs)
Both structural identification setups yielded highly consistent, matching results:
* **Output Trough:** A contractionary FFR shock drops real GDP growth to a trough of **-0.2 percentage points** within 4 to 5 quarters, showing policy operates quickly through consumer and investment demand channels.
* **Price Puzzle:** A minor initial spike in inflation occurred before shifting into a steady, long-run disinflationary path.
* **Monetary Tightening:** M2 growth drops instantly on impact, validating immediate credit condition tightening.

---

