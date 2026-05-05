# Statistics and Data Science Portfolio – Komi Roger Ayi

**Biostatistician · Data Scientist – Causal modeling & predictive analytics**

Welcome to my professional portfolio.

I am a statistician with a strong academic background in mathematics and a specialization in applied statistics. I hold a Master's degree from the *Université du Québec à Montréal* (UQAM), focused on **causal mediation modeling** and its applications in health research.

This portfolio presents a curated selection of my academic, research, and personal projects in data analysis, statistical modeling, predictive analytics, and causal inference. It also documents my technical skills, certifications, and academic distinctions.

My work reflects a commitment to **rigorous methodology**, **reproducibility**, and the practical use of data to support decision-making — particularly in public health and applied research contexts.

[![LinkedIn](https://img.shields.io/badge/linkedin-%230077B5.svg?style=for-the-badge&logo=linkedin&logoColor=white)](https://www.linkedin.com/in/komi-ayi)
[![Gmail](https://img.shields.io/badge/Gmail-D14836?style=for-the-badge&logo=gmail&logoColor=white)](mailto:ayi1rogera@gmail.com)
[![GitHub](https://img.shields.io/badge/GitHub-181717?style=for-the-badge&logo=github&logoColor=white)](https://github.com/komiayi)

---

## Master's Thesis

**[*Analyse de médiation causale pour des médiateurs non causalement liés*](https://archipel.uqam.ca/19950)**
*Causal mediation analysis for non-causally-linked mediators.*

Master's thesis in Statistics, **Université du Québec à Montréal**, defended in **2025** under the supervision of **Professor Karim Oualkacha** (Department of Mathematics, UQAM). The thesis develops two original parametric methods – CC (Constant Correlation) and CNC (Non-Constant Correlation) – for estimating natural direct and indirect effects when mediators are correlated through an unmeasured common cause, rather than through a direct causal pathway.

The methods were validated on a real-world application linking **childhood trauma**, **DNA methylation**, and **cortisol stress reactivity**, and have been operationalized through a dedicated R Shiny application (see the `rbcm` project below).

[![Read on Archipel UQAM](https://img.shields.io/badge/Archipel%20UQAM-Read%20the%20thesis-red?style=for-the-badge&logo=adobe-acrobat-reader&logoColor=white)](https://archipel.uqam.ca/19950)
[![Companion R Shiny tool](https://img.shields.io/badge/Companion%20tool-rbcm-276DC3?style=for-the-badge&logo=r&logoColor=white)](https://github.com/komiayi/rbcm)

### Associated presentations
- ESPUM Symposium – Quantitative Methods in Health Research (2025)
- Mediation Research Days, UQAM (2024)
- SSC Annual Meeting, Carleton University, Ottawa (2023)

---

## Achievements

- **Excellence Scholarship**, *Institut des Sciences Mathématiques du Québec* (ISM), for research in causal mediation analysis (2022–2023).
- **Recruitment Scholarship**, *STATQAM – Faculty Research Center in Statistics and Data Science*, UQAM (2021).
- **Excellence Scholarship**, *African Center of Excellence in Mathematical Sciences and Applications* (funded by the World Bank), for research on the controllability of integro-differential equations at the *Institute of Mathematics and Physical Sciences* (2018–2019).
- **ERMIT Program Scholarship** (Entrepreneurship, Resources, Management, Innovation, and Technologies) for the Master's degree in Statistics–Probability (2015–2017).

---

## Projects

### rbcm — R Shiny application for causal mediation

*R Shiny · Statistical software · Biostatistics*

An R Shiny application implementing the **CC** and **CNC** methods developed in my Master's thesis, designed to make these techniques accessible to applied researchers without requiring advanced R programming skills. The interface integrates data import, assumption checking, parametric estimation of natural direct and indirect effects, and interactive visualization within a single workflow.

The package follows CRAN-style structure with a full `DESCRIPTION`, `R/`, `tests/`, MIT license, and a `CITATION.cff` file enabling academic citation. Active development is currently focused on optimizing the estimators for large-dimensional matrices.

**Key competencies:** R Shiny development · S3 object-oriented design · Software packaging · Numerical optimization.

[![Source code](https://img.shields.io/badge/GitHub-Source%20code-181717?style=for-the-badge&logo=github)](https://github.com/komiayi/rbcm)

---

### Ligue 1 strategic forecasting engine – Season 2025–2026

*Predictive analytics · Machine learning · Live web application*

A high-performance predictive tool that quantifies the probabilities of match outcomes for the French Ligue 1. By moving beyond intuition-based analysis, this project leverages a data-driven approach to forecast results across the current 18-club season.

The core of the application is a **multinomial logistic regression model** developed with Scikit-Learn. It processes team-level performance metrics — such as offensive efficiency and defensive resilience – to output calibrated win/draw/loss probabilities. The interface is a custom-built **Streamlit dashboard** with light/dark theme support.

**Key competencies:** Multinomial logistic regression (Scikit-Learn) · Data normalization and pipeline management (Pandas) · Streamlit deployment.

[![Launch the Live Application](https://img.shields.io/badge/Streamlit-Launch%20the%20Live%20Application-FF4B4B?style=for-the-badge&logo=streamlit&logoColor=white)](https://ligue-1-predictor.streamlit.app/)
[![Source code](https://img.shields.io/badge/GitHub-Source%20code-181717?style=for-the-badge&logo=github)](https://github.com/komiayi/Ligue-1)

---

### Mediation analysis with unmeasured confounding

*Biostatistics · Causal modeling · Epigenetics*

A research project addressing a major challenge in multiple mediation analysis: estimating effects accurately when mediators are correlated through an unmeasured common cause. Two advanced parametric methods – CC and CNC — are applied and evaluated to analyze the influence of childhood trauma on cortisol stress reactivity via DNA methylation. The implementation successfully corrects for confounding bias and confirms the robustness of a significant direct causal effect.

This project provides the empirical application of the methods developed in my Master's thesis. The interactive implementation is available in the `rbcm` project above.

**Key competencies:** Causal inference · Advanced statistical modeling · R programming.

[![View detailed project](https://img.shields.io/badge/GitHub-View%20Detailed%20Project-181717?style=for-the-badge&logo=github)](https://github.com/komiayi/dna_mediation)

---

### Waze user churn prediction

*Machine learning · Imbalanced classification · Threshold optimization*

A supervised learning pipeline that predicts which Waze users are at risk of churning, based on behavioral features. **Three tree-based models** are compared — Decision Tree, Random Forest, and XGBoost — with explicit handling of the moderate class imbalance (18% churners) through `class_weight="balanced"` and `scale_pos_weight`.

The project includes an explicit **ethical analysis of error costs** (false negatives vs false positives) that guides the choice of recall as the primary evaluation metric, and a **decision-threshold optimization** step on the validation set to refine the precision–recall trade-off. The final XGBoost model achieves a **recall of 0.64** at an optimized threshold of 0.49.

**Key competencies:** Tree-based modeling · GridSearchCV with cross-validation · Class imbalance handling · Threshold optimization · Feature interpretability.

[![View Notebook](https://img.shields.io/badge/Jupyter-Project%20Notebook-F37626?style=for-the-badge&logo=jupyter&logoColor=white)](https://github.com/komiayi/churn-waze-prediction/blob/main/notebooks/Waze_ML.ipynb)
[![Source code](https://img.shields.io/badge/GitHub-Source%20code-181717?style=for-the-badge&logo=github)](https://github.com/komiayi/churn-waze-prediction)

---

### Regression trees and ensemble methods — MSE and stability

*Statistical learning · Monte Carlo simulation · R*

An academic project comparing regression trees, Bagging, and Random Forests on **predictive accuracy (MSE)** and **stability** under data perturbation. The study combines a controlled Monte Carlo simulation on the Friedman benchmark (1 000 replications) with a real-world application to apartment construction costs in Tehran (UCI Machine Learning Repository).

Key findings: Random Forest reduces MSE by ~53% relative to a single regression tree, and is approximately **20× more stable** under perturbation. Results consistent on both simulated and real data.

**Key competencies:** Predictive modeling · Model evaluation (MSE, stability) · Monte Carlo simulation · R (`rpart`, `randomForest`, `tidyverse`, `ipred`).

[![View Detailed Project](https://img.shields.io/badge/GitHub-View%20Detailed%20Project-181717?style=for-the-badge&logo=github)](https://github.com/komiayi/regression-trees-ensemble-methods)
[![View PDF Report](https://img.shields.io/badge/PDF-Detailed%20Report-red?style=for-the-badge&logo=adobe-acrobat-reader&logoColor=white)](https://github.com/komiayi/regression-trees-ensemble-methods/blob/main/docs/ProjetMAT8886.pdf)

---

### MixLaw — R package for mixture distributions

*R package development · S3 object-oriented programming · Documentation*

A custom R package implementing an **S3 class hierarchy** for probability distribution mixtures, with methods for random sampling, density visualization, and basic descriptive statistics. Includes a specialized subclass for mixtures of normal distributions, full roxygen2 documentation, and unit tests with `testthat`.

Developed for the graduate course **MAT8186 – Techniques avancées en programmation statistique : R** at UQAM.

**Key competencies:** R package authoring · S3 dispatch · roxygen2 · testthat.

[![MixLaw Package](https://img.shields.io/badge/R-MixLaw-276DC3?style=for-the-badge&logo=r&logoColor=white)](https://github.com/komiayi/MixLaw)

---

## Core competencies

- **Methodologies:** Causal mediation analysis · Predictive modeling · Statistical inference · Monte Carlo simulation · Imbalanced classification.
- **Languages:** Python (Pandas, NumPy, Scikit-Learn, SciPy, Matplotlib, Streamlit) · R (dplyr, tidyr, caret, ggplot2, rpart, randomForest, Shiny) · SQL.
- **Tools:** Git/GitHub · SQLiteStudio · Tableau · Jupyter · RStudio · LaTeX.

---

## Certifications

- **[Google Advanced Data Analytics — Professional Certificate](https://www.coursera.org/account/accomplishments/professional-cert/93YW7B6ODR0U)**
  - [Google Advanced Data Analytics Capstone](https://www.coursera.org/account/accomplishments/verify/HEFB7OM479ZC)
  - [The Nuts and Bolts of Machine Learning](https://www.coursera.org/account/accomplishments/verify/NHM4UFCKSGK3)
  - [Regression Analysis: Simplify Complex Data Relationships](https://coursera.org/share/ba480a8f8584086658c4aaf443656cd9)
  - [The Power of Statistics](https://www.coursera.org/account/accomplishments/verify/1VZPPBSZW8VZ)
  - [Go Beyond the Numbers: Translate Data into Insights](https://www.coursera.org/account/accomplishments/verify/KTGLU02ZLC2K)
  - [Get Started with Python](https://www.coursera.org/account/accomplishments/verify/VLJML8RZZDJA)
  - [Foundations of Data Science](https://www.coursera.org/account/accomplishments/verify/2REIFKX25FFF)
- [Introducing DAX](https://github.com/komiayi/komiayi.github.io/blob/main/Files/dax.pdf)
- [Create a Dashboard with Tableau](https://github.com/komiayi/komiayi.github.io/blob/main)

---

## Contact

I am currently open to data science, biostatistics, and applied research opportunities.

[![LinkedIn](https://img.shields.io/badge/linkedin-%230077B5.svg?style=for-the-badge&logo=linkedin&logoColor=white)](https://www.linkedin.com/in/komi-ayi)
[![Gmail](https://img.shields.io/badge/Gmail-D14836?style=for-the-badge&logo=gmail&logoColor=white)](mailto:ayi1rogera@gmail.com)
[![GitHub](https://img.shields.io/badge/GitHub-181717?style=for-the-badge&logo=github&logoColor=white)](https://github.com/komiayi)
