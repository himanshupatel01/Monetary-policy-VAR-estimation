# EC3367 Project — VAR Analysis: Monetary Policy Shocks
# Candidate Number: 2606019

library(vars)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tseries)
library(gridExtra)


# 1. LOAD AND PREPARE DATA

df <- read.csv("~/Downloads/Projectdata.csv")

# Parse quarterly dates
df$Date <- as.Date(df$Date, format = "%d/%m/%Y")

cat("Dataset dimensions:", nrow(df), "rows x", ncol(df), "cols\n")
cat("Date range:", format(df$Date[1]), "to", format(df$Date[nrow(df)]), "\n")
cat("Variables:", paste(names(df), collapse = ", "), "\n\n")

# Variables selected for the VAR:
#   Inf.   - CPI inflation rate (%)
#   Unemp  - Unemployment rate (%)
#   GDPgr  - Real GDP growth rate (%)
#   FFR    - Federal Funds Rate (%) — monetary policy instrument
#   M2gr   - M2 money supply growth (%) — monetary transmission
#
# Cholesky ordering: Inf -> Unemp -> GDPgr -> FFR -> M2gr
# Rationale: slow-moving macro variables first, policy rate reacts
# to current conditions, money supply responds endogenously.

varnames <- c("Inf.", "Unemp", "GDPgr", "FFR", "M2gr")

# Create time-series object (quarterly, start 1971 Q1)
vardata <- ts(df[, varnames],
              start     = c(1971, 1),
              frequency = 4)


# 2. DESCRIPTIVE STATISTICS & TIME-SERIES PLOTS

cat("=== Descriptive Statistics ===\n")
print(summary(df[, varnames]))

# Plot all five variables
df_long <- df %>%
  select(Date, all_of(varnames)) %>%
  pivot_longer(-Date, names_to = "Variable", values_to = "Value")

p_ts <- ggplot(df_long, aes(x = Date, y = Value)) +
  geom_line(colour = "steelblue") +
  facet_wrap(~ Variable, scales = "free_y", ncol = 2) +
  theme_minimal(base_size = 11) +
  labs(title = "Figure 1: US Macroeconomic Variables, 1971 Q1 – 2010 Q4",
       x = NULL, y = NULL)

print(p_ts)


# 3. STATIONARITY TESTS

cat("\n=== ADF Unit Root Tests (H0: unit root) ===\n")
for (v in varnames) {
  test <- adf.test(na.omit(vardata[, v]))
  cat(sprintf("%-8s: statistic = %6.3f, p-value = %.4f %s\n",
              v,
              test$statistic,
              test$p.value,
              ifelse(test$p.value < 0.05, "=> Stationary", "=> Possible unit root")))
}

# All variables are growth rates or level rates, so stationarity
# is expected. Following Sims, Stock & Watson (1990), we estimate
# the VAR in levels regardless.


# 4. LAG LENGTH SELECTION

cat("\n=== Lag Length Selection (max = 8) ===\n")
lag_select <- VARselect(vardata, lag.max = 8, type = "const")
print(lag_select$selection)

# p = 2 chosen: balances AIC/BIC, consistent with quarterly macro VARs
p <- 2
cat("\nSelected lag length: p =", p, "\n")


# 5. ESTIMATE VAR

var_model <- VAR(vardata, p = p, type = "const")

cat("\n=== VAR Model Summary (FFR equation) ===\n")
summary(var_model$varresult$FFR)


# 6. DIAGNOSTIC TESTS

cat("\n=== Diagnostic Tests ===\n")

# Serial correlation
serial_test <- serial.test(var_model, lags.pt = 12, type = "PT.asymptotic")
cat("Serial correlation (Portmanteau) p-value:",
    round(serial_test$serial$p.value, 4), "\n")

# Normality
norm_test <- normality.test(var_model, multivariate.only = TRUE)
cat("Normality (JB) p-value:",
    round(norm_test$jb.mul$JB$p.value, 4), "\n")

# Stability: all eigenvalues inside the unit circle
roots_val <- roots(var_model)
cat("All roots inside unit circle:", all(roots_val < 1), "\n")
cat("Max root modulus:", round(max(roots_val), 4), "\n")


# 7. IDENTIFICATION 1: RECURSIVE (CHOLESKY) ORDERING

# Order: Inf. -> Unemp -> GDPgr -> FFR -> M2gr
# Zero restrictions: variables ordered before FFR do not respond
# to the FFR shock contemporaneously. FFR reacts to current
# inflation, unemployment, and output (Taylor Rule logic).

cat("\n=== Computing Cholesky IRFs (500 bootstrap runs) ===\n")

irf_chol <- irf(var_model,
                impulse  = "FFR",
                response = varnames,
                n.ahead  = 20,
                ortho    = TRUE,   # Cholesky orthogonalisation
                ci       = 0.95,
                boot     = TRUE,
                runs     = 500,
                seed     = 12345)

# Plot
plot(irf_chol,
     main = "Figure 2: Cholesky IRFs — Shock to FFR")


# 8. IDENTIFICATION 2: ALTERNATIVE RECURSIVE ORDERING

# Robustness check: swap M2gr and FFR in the ordering.
# New order: Inf -> Unemp -> GDPgr -> M2gr -> FFR
# This assumes money supply is set independently of the
# interest rate within the quarter (money-supply rule).

varnames_alt <- c("Inf.", "Unemp", "GDPgr", "M2gr", "FFR")

vardata_alt <- ts(df[, varnames_alt],
                  start     = c(1971, 1),
                  frequency = 4)

var_model_alt <- VAR(vardata_alt, p = p, type = "const")

cat("\n=== Computing Alternative Ordering IRFs (500 bootstrap runs) ===\n")

irf_alt <- irf(var_model_alt,
               impulse  = "FFR",
               response = varnames_alt,
               n.ahead  = 20,
               ortho    = TRUE,
               ci       = 0.95,
               boot     = TRUE,
               runs     = 500,
               seed     = 12345)

plot(irf_alt,
     main = "Figure 3: Alternative Ordering IRFs — Shock to FFR")


# 9. COMPARISON PLOT: BOTH IDENTIFICATIONS

# Helper function to extract IRF data frame
extract_irf <- function(irf_obj, var_name, id_label) {
  # irf() stores results as: irf_obj$irf$IMPULSE_VAR — a matrix with one
  # column per response variable. We grab the column matching var_name.
  impulse_name <- names(irf_obj$irf)[1]   # only one impulse variable
  
  irf_mat   <- irf_obj$irf[[impulse_name]]
  lower_mat <- irf_obj$Lower[[impulse_name]]
  upper_mat <- irf_obj$Upper[[impulse_name]]
  
  # Column names on the matrix are the response variable names
  data.frame(
    horizon  = 0:20,
    irf      = as.numeric(irf_mat[, var_name]),
    lower    = as.numeric(lower_mat[, var_name]),
    upper    = as.numeric(upper_mat[, var_name]),
    id       = id_label,
    variable = var_name
  )
}

gdp_chol <- extract_irf(irf_chol, "GDPgr", "Cholesky (M2 last)")
gdp_alt  <- extract_irf(irf_alt,  "GDPgr", "Alt. Ordering (FFR last)")
inf_chol <- extract_irf(irf_chol, "Inf.",   "Cholesky (M2 last)")
inf_alt  <- extract_irf(irf_alt,  "Inf.",   "Alt. Ordering (FFR last)")
une_chol <- extract_irf(irf_chol, "Unemp", "Cholesky (M2 last)")
une_alt  <- extract_irf(irf_alt,  "Unemp", "Alt. Ordering (FFR last)")

compare_df <- rbind(gdp_chol, gdp_alt, inf_chol, inf_alt, une_chol, une_alt)

p_compare <- ggplot(compare_df,
                    aes(x = horizon, y = irf, colour = id, fill = id)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ variable, scales = "free_y", ncol = 3) +
  scale_colour_manual(values = c("steelblue", "tomato")) +
  scale_fill_manual(values   = c("steelblue", "tomato")) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom") +
  labs(title   = "Figure 4: IRF Comparison — GDP Growth, Inflation, Unemployment",
       subtitle = "Response to a contractionary FFR shock, 95% bootstrap CI",
       x        = "Quarters after shock",
       y        = "Response (percentage points)",
       colour   = "Identification",
       fill     = "Identification")

print(p_compare)
