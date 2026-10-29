# PREVENT (Predicting Risk of Cardiovascular Disease Events)

PREVENT is the American Heart Association's comprehensive cardiovascular risk prediction model for primary prevention, replacing the previous ASCVD Pooled Cohort Equations. It incorporates additional risk factors and provides more precise risk estimates.
This implementation uses the PREVENT equations from 2026-10-14.

## Overview

PREVENT (Predicting Risk of cardiovascular disease EVENTs) was developed to improve upon the ASCVD risk calculator by incorporating additional risk factors, novel biomarkers, and contemporary data. The model provides both 10-year and 30-year risk estimates.

## Risk Factors

### Mandatory Risk Factors
- **Age**: Years
- **Sex**: Male/Female
- **Total Cholesterol**: mmol/L
- **HDL Cholesterol**: mmol/L
- **Systolic Blood Pressure**: mmHg
- **Diabetes**: Has currently diabetes
- **Smoking Status**: Current smoker
- **BMI**: Body Mass Index in kg/m²
- **eGFR**: Estimated glomerular filtration in ml/min/1.73m²
- **Antihypertensive**: Current antihypertensive medication use
- **Statin**: Current statin use (True/False/0/1)

### Optional Risk Factors
- **USA ZIP code**: USA zip code. Used to estimate the Social Deprivation Index (SDI). Ignored if SDI is specified.
- **SDI**: Social Deprivation Index tertile: 1=low, 2=medium, 3=high. Overrides USA zip code if valid (1-3).
- **uACR**: Urine albumin-to-creatinine ratio in mg/g, from 0.1-25000.
- **HbA1c**: Glycated hemoglobin in %, from 4.5-15.

## Mathematical Formulation

PREVENT uses a sophisticated modeling approach combining:

1. **Logistic regression** for 10-year risk
2. **Survival analysis** for 30-year risk
3. **Calibration factors** for different populations

## Risk Categories

### 10-Year Risk
- **Low Risk**: <5%
- **Borderline Risk**: 5-7.5%
- **Intermediate Risk**: 7.5-20%
- **High Risk**: ≥20%

### 30-Year Risk
- **Low Risk**: <20%
- **Intermediate Risk**: 20-30%
- **High Risk**: ≥30%

## Reference

Khan SS, Coresh J, Pencina MJ, et al. Novel prediction equations for absolute risk assessment of total cardiovascular disease incorporating cardiovascular-kidney-metabolic health: the PREVENT equations. *Circulation*. 2024.

DOI: [10.1161/CIRCULATIONAHA.123.067626](https://doi.org/10.1161/CIRCULATIONAHA.123.067626)

## Implementation Notes

- Comprehensive primary prevention model
- Incorporates CKM (Cardiovascular-Kidney-Metabolic) health
- Includes novel biomarkers and risk enhancers
- Provides both 10-year and 30-year risk estimates
- Represents current state-of-the-art in CVD risk prediction
- Based on R implementation (preventr) from Dr. Martin Mayer
