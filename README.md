# Monte Carlo Simulations

This repository contains R code and analysis for simulation-based statistical methods and power analyses. The work includes investigations into trimmed means, power curves, and robustness of confidence intervals.

## Contents

- MSE estimation of trimmed means from the Cauchy distribution  
- Empirical power curve for a two-sided t-test  
- Power curves for different sample sizes in t-tests  
- Coverage probability of the t-test with non-normal data (Chi-square distribution)

---

## MSE of Trimmed Means

**Objective:** Estimate the Mean Squared Error (MSE) of level-𝑘 trimmed means using samples from a standard Cauchy distribution (n = 20, m = 1000 simulations).  
**Key Insight:** The Cauchy distribution lacks a defined mean, so the center (median) is used as the target parameter. Trimming is essential to reduce the influence of outliers.

**Findings:**
- MSE was minimized when trimming 6 values from each tail.
- Low trimming (k = 1 or 2) resulted in higher MSE due to outlier influence.

---

## Empirical Power Curve for a t-Test

**Objective:** Estimate the power of a two-sided t-test under varying alternative hypotheses (µ ≠ 500), using normally distributed samples with σ = 100 and n = 20.

**Method:**
- Simulate 1000 samples for each µ ∈ [450, 650].
- Perform t-tests and calculate proportion of p-values ≤ 0.05.

**Findings:**
- Power increased as µ moved away from 500 in either direction.
- The curve confirmed expected behavior: higher power when the true mean deviates more from the null.

---

## Power Curves for Varying Sample Sizes

**Objective:** Compare the power of the two-sided t-test across different sample sizes (n = 10, 20, 30, 40, 50).

**Method:**
- Simulate tests but for each value of n.
- Plot multiple power curves on a single graph without standard error bars.

**Findings:**
- Larger sample sizes resulted in faster power gains.
- Power was significantly higher at larger deviations from the null for higher n.

---

## Coverage Probability of the t-Interval

**Objective:** Test the robustness of the t-interval on non-normal data (χ² with df = 2) to see if it retains the nominal 95% coverage rate.

**Method:**
- Simulate 1000 samples of size n = 20 from χ²(2).
- Construct t-tests and determine whether the true mean (µ = 2) was captured.

**Findings:**
- The empirical coverage probability was approximately 91.7%, lower than the expected 95%.
- Indicates reduced accuracy of t-intervals for non-normal data, but still better than variance intervals.
- Confirms robustness of t-intervals for mean estimation under mild non-normality.

---

Make sure the following packages are installed:
```R
install.packages("ggplot2")
