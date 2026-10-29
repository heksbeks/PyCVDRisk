"""
PREVENT cardiovascular risk model.

The PREVENT equations are the American Heart Association (AHA) 2023
cardiovascular risk prediction model for Total CVD, ASCVD, and Heart Failure.
This implementation uses the PREVENT equations from 2026-10-14.

This is a Python implementation based on preventr:
https://github.com/martingmayer/preventr
https://cran.r-project.org/web/packages/preventr/index.html
https://github.com/AHA-DS-Analytics/PREVENT

See preventr for detailed information.

Reference:
    Khan SS, Matsushita K, Yang C, et al. Development and Validation of the
    American Heart Association's PREVENT Equations. Circulation. 2024;149(5):388-404.
    DOI: 10.1161/CIRCULATIONAHA.123.067626
"""

import logging
import typing
from typing import Literal, Union
from pydantic import validate_call

import numpy as np
import pandas as pd

from cvd_risk.core.base import RiskModel
from cvd_risk.core.validation import PatientData, RiskResult

from . import prevent_data

logger = logging.getLogger(__name__)


class Prevent(RiskModel):
    """
    PREVENT cardiovascular risk prediction model (AHA 2023).

    Predicts 10-year and 30-year risk for:
    - Total CVD (composite of ASCVD and heart failure)
    - Heart Failure
    - ASCVD (atherosclerotic cardiovascular disease), includes coronary heart disease (CHD) and stroke as defined below
    - CHD (coronary heart disease, heart attack), includes nonfatal myocardial infarction (MI) and fatal CHD
    - Stroke

    10-year risk valid for ages 30-79.
    30-year risk valid for ages 30-59.
    """

    model_name = "PREVENT"
    model_version = "2023"
    supported_outcomes = Literal[
        "total_cvd_10", "ascvd_10", "hf_10", "chd_10", "stroke_10",
        "total_cvd_30", "ascvd_30", "hf_30", "chd_30", "stroke_30"
        ]

    @validate_call
    def __init__(self, outcome: supported_outcomes = "total_cvd_10") -> None:
        f"""
        Initialize PREVENT model.

        Parameters
        ----------
        outcome : str
            Risk outcome to calculate. One of: {typing.get_args(self.supported_outcomes)}
        """
        super().__init__()
        self.outcome = outcome
        self._debug = False

    def validate_input(self, patient: PatientData) -> None:
        super().validate_input(patient)
        if self.outcome.endswith("_10") and (patient.age < 30 or patient.age > 79):
            logger.warning(
                f"Age {patient.age} outside optimal range [30, 79] for {self.outcome}. "
                "PREVENT accuracy may be reduced."
            )
        if self.outcome.endswith("_30") and (patient.age < 30 or patient.age > 59):
            logger.warning(
                f"Age {patient.age} outside optimal range [30, 59] for {self.outcome}. "
                "PREVENT accuracy may be reduced."
            )
        if patient.egfr is None:
            raise ValueError("eGFR is required for PREVENT calculations")
        if patient.antihypertensive is None:
            raise ValueError("Antihypertensive medication status is required for PREVENT calculations")
        if patient.statin_use is None:
            raise ValueError("Statin use status is required for PREVENT calculations")
        if patient.bmi is None:
            raise ValueError("BMI is required for PREVENT calculations")

    def _transform_variables(self,
                             age: Union[float, np.ndarray],
                             total_chol_mmol: Union[float, np.ndarray],
                             hdl_chol_mmol: Union[float, np.ndarray],
                             sbp: Union[float, np.ndarray],
                             egfr: Union[float, np.ndarray],
                             bmi: Union[float, np.ndarray]) -> dict:
        """
        Transform input variables according to PREVENT methodology.

        Parameters
        ----------
        age : float or array
            Age in years
        total_chol_mmol : float or array
            Total cholesterol in mmol/L (converted to mg/dL internally)
        hdl_chol_mmol : float or array
            HDL cholesterol in mmol/L (converted to mg/dL internally)
        sbp : float or array
            Systolic blood pressure in mmHg
        egfr : float or array
            eGFR in mL/min/1.73m2
        bmi : float or array
            BMI in kg/m2

        Returns
        -------
        dict
            Transformed variables
        """

        # PREVENT expects mmol/L internally
        non_hdl_mmol = total_chol_mmol - hdl_chol_mmol

        # Transformations
        t_age = (age - 55.0) / 10.0
        t_age2 = t_age ** 2
        t_cnhdl = non_hdl_mmol - 3.5
        t_chdl = (hdl_chol_mmol - 1.3) / 0.3

        # SBP splines
        t_csbp = (np.minimum(sbp, 110.0) - 110.0) / 20.0
        t_csbp2 = (np.maximum(sbp, 110.0) - 130.0) / 20.0

        # eGFR splines
        t_cegfr = (np.minimum(egfr, 60.0) - 60.0) / -15.0
        t_cegfr2 = (np.maximum(egfr, 60.0) - 90.0) / -15.0

        t_cbmi = (np.minimum(bmi, 30.0) - 25.0) / 5.0
        t_cbmi2 = (np.maximum(bmi, 30.0) - 30.0) / 5.0

        return {
            't_age': t_age,
            't_age2': t_age2,
            't_cnhdl': t_cnhdl,
            't_chdl': t_chdl,
            't_csbp': t_csbp,
            't_csbp2': t_csbp2,
            't_cegfr': t_cegfr,
            't_cegfr2': t_cegfr2,
            't_cbmi': t_cbmi,
            't_cbmi2': t_cbmi2,
        }

    def _get_risk_thresholds(self, age: Union[int, float, np.ndarray]) -> tuple:
        """
        Get clinical risk thresholds based on age (AHA/ACC Guidelines).

        Returns
        -------
        tuple (moderate_cutoff, high_cutoff)
            >= moderate_cutoff: Moderate Risk
            >= high_cutoff: High Risk
        """
        # Vectorized logic for NumPy arrays
        if isinstance(age, (np.ndarray, pd.Series)):
            moderate_cut = np.select(
                [age < 40, (age >= 40) & (age < 60), age >= 60],
                [5.0, 7.5, 10.0],
                default=7.5
            )
            high_cut = np.select(
                [age < 40, (age >= 40) & (age < 60), age >= 60],
                [10.0, 15.0, 20.0],
                default=15.0
            )
            return moderate_cut, high_cut

        # Scalar logic for single patient
        if age < 40:
            return 5.0, 10.0
        elif age < 60:
            return 7.5, 15.0
        else:
            return 10.0, 20.0

    def calculate(self, patient: PatientData) -> RiskResult:
        """
        Calculate PREVENT risk for a single patient.

        Parameters
        ----------
        patient : PatientData
            Patient information

        Returns
        -------
        RiskResult
            Risk calculation results
        """
        self.validate_input(patient)

        # Transform variables
        tf = self._transform_variables(
            patient.age, patient.total_cholesterol, patient.hdl_cholesterol,
            patient.systolic_bp, patient.egfr, patient.bmi
        )

        t_sid_group = 0
        if patient.sdi and 1 <= patient.sdi <= 3:
            t_sid_group = patient.sdi
        elif patient.us_zip_code:
            i = patient.us_zip_code
            if i > 0 and i < len(prevent_data.SDI_GROUP):
                t_sid_group = prevent_data.SDI_GROUP[i]

        # Get coefficients

        sex = patient.sex.lower()

        c = 0
        if patient.hba1c_percent:
            c += 1
        if patient.acr:
            c += 1
        if t_sid_group:
            c += 1

        model = None
        if c > 1:
            model = 'full'
        elif c == 0:
            model = 'base'
        elif patient.hba1c_percent:
            model = 'hba1c'
        elif patient.acr:
            model = 'uacr'
        elif t_sid_group:
            model = 'sdi'
        coeffs = prevent_data.COEFFS[f'{model}_{sex}_{self.outcome}']

        # Indicator variables
        i_db = float(patient.diabetes or 0.0)
        i_smk = float(patient.smoking or 0.0)
        i_htn = float(patient.antihypertensive or 0.0)
        i_statin = float(patient.statin_use or 0.0)

        if t_sid_group:
            if t_sid_group == 1:
                sdi_term = 0.0
            elif t_sid_group == 2:
                sdi_term = coeffs[24]
            elif t_sid_group == 3:
                sdi_term = coeffs[25]
        else:
            sdi_term = coeffs[26]

        if patient.acr:
            uacr_term = np.log(patient.acr) * coeffs[27]
        else:
            uacr_term = coeffs[28]

        if patient.hba1c_percent:
            if patient.diabetes:
                hba1c_term = (patient.hba1c_percent - 5.3) * coeffs[30]
            else:
                hba1c_term = (patient.hba1c_percent - 5.3) * coeffs[29]
        else:
            hba1c_term = coeffs[31]

        # Calculate log-odds
        log_odds = \
            coeffs[0] + \
            coeffs[1] * tf['t_age'] + \
            coeffs[2] * tf['t_age'] ** 2 + \
            coeffs[3] * tf['t_cnhdl'] + \
            coeffs[4] * tf['t_chdl'] + \
            coeffs[5] * tf['t_csbp'] + \
            coeffs[6] * tf['t_csbp2'] + \
            coeffs[7] * i_db + \
            coeffs[8] * i_smk + \
            coeffs[9] * tf['t_cbmi'] + \
            coeffs[10] * tf['t_cbmi2'] + \
            coeffs[11] * tf['t_cegfr'] + \
            coeffs[12] * tf['t_cegfr2'] + \
            coeffs[13] * i_htn + \
            coeffs[14] * i_statin + \
            coeffs[15] * i_htn * tf['t_csbp2'] + \
            coeffs[16] * i_statin * tf['t_cnhdl'] + \
            coeffs[17] * tf['t_age'] * tf['t_cnhdl'] + \
            coeffs[18] * tf['t_age'] * tf['t_chdl'] + \
            coeffs[19] * tf['t_age'] * tf['t_csbp2'] + \
            coeffs[20] * tf['t_age'] * i_db + \
            coeffs[21] * tf['t_age'] * i_smk + \
            coeffs[22] * tf['t_age'] * tf['t_cbmi2'] + \
            coeffs[23] * tf['t_age'] * tf['t_cegfr'] + \
            sdi_term + \
            uacr_term + \
            hba1c_term

        # Convert to risk percentage
        e = np.exp(log_odds)
        risk_percent = e/(1+e) * 100.0
        risk_percent = round(np.clip(risk_percent, 0.0, 100.0), 16 if self._debug else 2)

        # Categorize risk
        moderate_cut, high_cut = self._get_risk_thresholds(patient.age)

        if risk_percent < moderate_cut:
            category = "low"
        elif risk_percent < high_cut:
            category = "moderate"
        else:
            category = "high"

        outcome_name = {
            "total_cvd_10": "Total CVD 10-year",
            "ascvd_10": "ASCVD 10-year",
            "hf_10": "Heart Failure 10-year",
            "chd_10": "Heart Attack 10-year",
            "stroke_10": "Stroke 10-year",
            "total_cvd_30": "Total CVD 30-year",
            "ascvd_30": "ASCVD 30-year",
            "hf_30": "Heart Failure 30-year",
            "chd_30": "Heart Attack 30-year",
            "stroke_30": "Stroke 30-year",
        }[self.outcome]

        calculation_metadata={
            "outcome": self.outcome,
            "outcome_name": outcome_name,
        }
        if self._debug:
            calculation_metadata["dbgA"] = sdi_term
            calculation_metadata["dbgB"] = uacr_term
            calculation_metadata["dbgC"] = hba1c_term

        return RiskResult(
            risk_score=risk_percent,
            risk_category=category,
            model_name=self.model_name,
            model_version=self.model_version,
            calculation_metadata=calculation_metadata
        )

    def calculate_batch(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Vectorized calculation for high-throughput analysis.

        Parameters
        ----------
        df : pd.DataFrame
            DataFrame with required columns

        Returns
        -------
        pd.DataFrame
            DataFrame with added risk columns
        """
        required = ["age", "sex", "total_cholesterol", "hdl_cholesterol",
                   "systolic_bp", "egfr", "diabetes", "smoking",
                   "antihypertensive", "statin_use"]

        if not set(required).issubset(df.columns):
            raise ValueError(f"Missing columns: {set(required) - set(df.columns)}")

        results = df.copy()
        results["risk_score"] = np.nan
        results["risk_category"] = "unknown"
        results["model_name"] = self.model_name
        results["outcome"] = self.outcome
        results["model"] = "n/a"

        # handle zip codes
        if "sdi" not in results.columns:
            results["sdi"] = None
        if "us_zip_code" in results.columns:
            mapz_arr = np.array(prevent_data.SDI_GROUP)
            mask = results["sdi"].isna() & results["us_zip_code"].notna()
            results.loc[mask, "sdi"] = mapz_arr[results.loc[mask, "us_zip_code"].astype(int)]

        # handle correct model
        results["model_score_hba1c"] = "hba1c_percent" in results.columns and results["hba1c_percent"] > 0
        results["model_score_acr"] = "acr" in results.columns and results["acr"] > 0
        results["model_score_sdi"] = "sdi" in results.columns and results["sdi"] > 0
        results["model_score"] = \
            (results["model_score_hba1c"]).astype(int) + \
            (results["model_score_acr"]).astype(int) + \
            (results["model_score_sdi"]).astype(int)

        # Process by sex for efficiency
        for sex in ["male", "female"]:
            sex_mask = results["sex"].str.lower() == sex
            if not sex_mask.any():
                continue

            # Sub-loop by model
            for model in ["full", "base", "hba1c", "uacr", "sdi"]:
                if model == "full":
                    model_mask = sex_mask & (results["model_score"] > 1)
                elif model == "base":
                    model_mask = sex_mask & (results["model_score"] == 0)
                elif model == "hba1c":
                    model_mask = sex_mask & results["model_score_hba1c"] & (results["model_score"] == 1)
                elif model == "uacr":
                    model_mask = sex_mask & results["model_score_acr"] & (results["model_score"] == 1)
                elif model == "sdi":
                    model_mask = sex_mask & results["model_score_sdi"] & (results["model_score"] == 1)
                if not model_mask.any():
                    continue
                results.loc[model_mask, "model"] = model

                coeffs = prevent_data.COEFFS[f"{model}_{sex}_{self.outcome}"]
                data = results[model_mask]

                # Vectorized transformations (inputs are in mmol/L, converted internally)
                tf = self._transform_variables(
                    data["age"].values, data["total_cholesterol"].values,
                    data["hdl_cholesterol"].values, data["systolic_bp"].values,
                    data["egfr"].values, data["bmi"].values
                )

                # Indicator variables
                i_db = data["diabetes"].fillna(False).astype(float).values
                i_smk = data["smoking"].fillna(False).astype(float).values
                i_htn = data["antihypertensive"].fillna(False).astype(float).values
                i_statin = data["statin_use"].fillna(False).astype(float).values

                # Add SDI
                if coeffs[24] and 'sdi' in results.columns:
                    cond_tier_1 = data['sdi'].notna() & (data['sdi'] == 1)
                    cond_tier_2 = data['sdi'].notna() & (data['sdi'] == 2)
                    cond_tier_3 = data['sdi'].notna() & (data['sdi'] == 3)
                    conditions = [cond_tier_1, cond_tier_2, cond_tier_3]
                    choices = [0, coeffs[24], coeffs[25]]
                    sdi_term = np.select(conditions, choices, default=coeffs[26])
                else:
                    sdi_term = coeffs[26]

                # Add uACR
                if 'acr' in results.columns:
                    acr_array = pd.to_numeric(data["acr"], errors='coerce').to_numpy(dtype=float) # convert to a native float array (None automatically becomes np.nan)
                    uacr_term = np.full(acr_array.shape, coeffs[28], dtype=float)
                    valid_mask = acr_array > 0
                    uacr_term[valid_mask] = np.log(acr_array[valid_mask]) * coeffs[27]
                else:
                    uacr_term = coeffs[28]

                # Add HbA1c
                if 'hba1c_percent' in results.columns:
                    hba1c_array = pd.to_numeric(data["hba1c_percent"], errors='coerce').to_numpy(dtype=float) # convert to a native float array (None automatically becomes np.nan)
                    cond_no_diabetes = (hba1c_array > 0) & ~data['diabetes']
                    cond_diabetes = (hba1c_array > 0) & data['diabetes']
                    choices = [
                        (hba1c_array - 5.3) * coeffs[29],
                        (hba1c_array - 5.3) * coeffs[30]
                    ]
                    hba1c_term = np.select([cond_no_diabetes, cond_diabetes], choices, default=coeffs[31])
                else:
                    hba1c_term = coeffs[31]

                # Vectorized log-odds calculation
                log_odds = \
                    coeffs[0] + \
                    coeffs[1] * tf['t_age'] + \
                    coeffs[2] * tf['t_age'] ** 2 + \
                    coeffs[3] * tf['t_cnhdl'] + \
                    coeffs[4] * tf['t_chdl'] + \
                    coeffs[5] * tf['t_csbp'] + \
                    coeffs[6] * tf['t_csbp2'] + \
                    coeffs[7] * i_db + \
                    coeffs[8] * i_smk + \
                    coeffs[9] * tf['t_cbmi'] + \
                    coeffs[10] * tf['t_cbmi2'] + \
                    coeffs[11] * tf['t_cegfr'] + \
                    coeffs[12] * tf['t_cegfr2'] + \
                    coeffs[13] * i_htn + \
                    coeffs[14] * i_statin + \
                    coeffs[15] * (i_htn * tf['t_csbp2']) + \
                    coeffs[16] * (i_statin * tf['t_cnhdl']) + \
                    coeffs[17] * (tf['t_age'] * tf['t_cnhdl']) + \
                    coeffs[18] * (tf['t_age'] * tf['t_chdl']) + \
                    coeffs[19] * (tf['t_age'] * tf['t_csbp2']) + \
                    coeffs[20] * (tf['t_age'] * i_db) + \
                    coeffs[21] * (tf['t_age'] * i_smk) + \
                    coeffs[22] * (tf['t_age'] * tf['t_cbmi2']) + \
                    coeffs[23] * (tf['t_age'] * tf['t_cegfr']) + \
                    sdi_term + \
                    uacr_term + \
                    hba1c_term

                # Convert to risk percentage
                e = np.exp(log_odds)
                risk_percent = e/(1+e) * 100.0
                risk_percent = np.round(np.clip(risk_percent, 0.0, 100.0), 16 if self._debug else 2)

                results.loc[model_mask, "risk_score"] = risk_percent
                if self._debug:
                    results.loc[model_mask, "dbgA"] = sdi_term
                    results.loc[model_mask, "dbgB"] = uacr_term
                    results.loc[model_mask, "dbgC"] = hba1c_term

        # Vectorized categorization
        moderate_cuts, high_cuts = self._get_risk_thresholds(results["age"].values)

        conditions = [
            results["risk_score"] < moderate_cuts,
            (results["risk_score"] >= moderate_cuts) & (results["risk_score"] < high_cuts),
            results["risk_score"] >= high_cuts
        ]
        choices = ["low", "moderate", "high"]
        results["risk_category"] = np.select(conditions, choices, default="unknown")

        return results
