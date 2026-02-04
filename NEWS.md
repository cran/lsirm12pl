# lsirm12pl 2.0.0 (2026-02-03)

## Major Changes

* **NEW**: Added Graded Response Model (GRM) for ordinal/Likert-scale data
  - `lsirmgrm()`: 1PL GRM with latent space
  - `lsirmgrm2pl()`: 2PL GRM with item discrimination parameters
  - Supports complete data, MAR, and MCAR missing data mechanisms
  - Implements threshold ordering constraint (β₁ > β₂ > ... > βₖ₋₁) for identifiability
  - Reference: De Carolis, Kang, & Jeon (2025), DOI: 10.1080/00273171.2025.2605678

* **NEW**: Adaptive MCMC algorithm
  - Automatic tuning of proposal standard deviations during burn-in
  - Parameter-specific target acceptance rates
  - Improved mixing and convergence for complex models
  - Available for all model types (1PL, 2PL, GRM)