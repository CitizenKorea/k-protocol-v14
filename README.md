# Geometric Analysis of VLBI NGS Data

## 1. Overview
This repository provides a geometric alignment analysis of VLBI (Very Long Baseline Interferometry) NGS observation data. Unlike standard post-processing methods that rely on extensive parameter tuning and manual residual shaving, this project focuses on verifying the fundamental geometric consistency of space-time using minimal variables.

## 2. Methodology
The analysis adopts a minimalist approach to preserve the integrity of raw observation data. Instead of utilizing hundreds of atmospheric and instrumental correction parameters, we apply two core geometric constants:

* **ck (Absolute Speed of Light):** 297,880,197.6 m/s
* **S_earth (Scale Factor):** 1.006419562

The engine performs a direct geometric strike without any post-hoc manual adjustments. Our goal is to observe the raw convergence of data when aligned to a deterministic geometric frame.

## 3. Key Principles
* **No Manual Shaving:** We do not artificially "shave" residuals to force-fit the data into an idealized curve.
* **Geometric Causality:** We prioritize the causal relationship between geometry and time-delay over statistical precision achieved through overfitting.
* **System U Alignment:** The observation stations and sources are mapped into an absolute geometric coordinate system to evaluate the "Raw Truth" of the residuals.

## 4. Analysis Results
* **Dataset:** 20JAN02XE_N005.ngs
* **Median Residual:** ~108,432 ns
* **Observation:** By applying only ck and S_earth, the chaotic raw data (initially exhibiting ~20ms error) collapses into a deterministic band of approximately 0.1ms. This residual is maintained as a "honest scar"—a reflection of uncorrected local gravitational and instrumental variables—rather than being suppressed through artificial post-processing.

## 5. Conclusion
This project demonstrates that a significant portion of what is often considered "stochastic noise" can be resolved through fundamental unit redefinition and geometric alignment. While standard models achieve sub-nanosecond precision through thousands of ad-hoc corrections, this analysis suggests that the underlying geometric structure of the universe is far more consistent than complex post-diction models imply.

---
**"We do not predict the past; we analyze the causality of the future."**

### Contact
**Lead Researcher:** estake@naver.com
