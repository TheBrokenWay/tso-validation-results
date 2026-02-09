"""
Simulation_Engine.py
PK Simulation Engine v1
One-compartment, first-order absorption and elimination.
Deterministic, no fabrication, exposure-only (non-clinical).
"""

from typing import Dict, Any, List
import math
import json
import os
from datetime import datetime, timezone


class SimulationEngine:
    """
    PK Simulation Engine.
    Computes concentration-time profiles for a single virtual patient
    given dose, regimen, and basic ADMET/OPE parameters.
    """

    def __init__(self, time_step_h: float = 0.5):
        if time_step_h <= 0:
            raise ValueError("time_step_h must be positive")
        self.time_step_h = time_step_h
        self.metadata = {"role": "LabSim", "version": "1.2.0-GAIP-PK"}

    def simulate_one_compartment(
        self,
        dose_mg: float,
        duration_h: float,
        dosing_interval_h: float,
        patient: Dict[str, Any],
        admet: Dict[str, Any],
    ) -> Dict[str, Any]:
        """
        Simulate one-compartment PK with first-order absorption and elimination.

        Inputs:
            dose_mg: single-dose amount in mg
            duration_h: total simulation time in hours
            dosing_interval_h: interval between doses (for multiple dosing)
            patient: dict with at least 'weight_kg'
            admet: dict with optional clearance and Vd info

        Returns:
            dict with time grid, concentration profile, and summary metrics.
        """

        if dose_mg <= 0:
            raise ValueError("dose_mg must be positive")
        if duration_h <= 0:
            raise ValueError("duration_h must be positive")
        if dosing_interval_h <= 0:
            raise ValueError("dosing_interval_h must be positive")

        weight_kg = patient.get("weight_kg", 70.0)

        # --- Derive basic PK parameters from ADMET (with safe defaults) ---

        # Volume of distribution (L/kg) – default 0.7 L/kg (approx. total body water)
        vd_L_per_kg = admet.get("distribution", {}).get("predicted_vd_L_per_kg")
        if vd_L_per_kg is None:
            vd_L_per_kg = 0.7
        vd_L_base = vd_L_per_kg * weight_kg
        
        # v2.0-PHASE2: Apply IIV vd_factor if present
        vd_factor = patient.get("vd_factor", 1.0)
        vd_L = vd_L_base * vd_factor

        # Clearance (L/h/kg) – default 0.05 L/h/kg (arbitrary but deterministic)
        cl_L_per_h_per_kg = admet.get("metabolism", {}).get("predicted_clearance_L_per_h_per_kg")
        if cl_L_per_h_per_kg is None:
            cl_L_per_h_per_kg = 0.05
        cl_L_per_h_base = cl_L_per_h_per_kg * weight_kg
        
        # v2.0-PHASE2: Apply IIV clearance_factor if present
        clearance_factor = patient.get("clearance_factor", 1.0)
        cl_L_per_h = cl_L_per_h_base * clearance_factor

        # Elimination rate constant ke = CL / Vd
        ke = cl_L_per_h / max(vd_L, 1e-6)

        # Absorption rate constant ka – fixed for now, deterministic
        ka_base = 1.0  # 1/h
        
        # v2.0-PHASE2: Apply IIV ka_factor if present
        ka_factor = patient.get("ka_factor", 1.0)
        ka = ka_base * ka_factor

        # Bioavailability – if unknown, assume 1.0 (exposure-only, not clinical)
        f = admet.get("absorption", {}).get("predicted_bioavailability")
        if f is None:
            f = 1.0

        # --- Time grid ---
        n_steps = int(duration_h / self.time_step_h) + 1
        time_grid: List[float] = [i * self.time_step_h for i in range(n_steps)]

        # --- Dosing schedule (repeated doses) ---
        dose_times = []
        t = 0.0
        while t <= duration_h + 1e-6:
            dose_times.append(t)
            t += dosing_interval_h

        # --- Concentration-time profile ---
        conc: List[float] = []
        for t in time_grid:
            c_t = 0.0
            for t_dose in dose_times:
                if t < t_dose:
                    continue
                dt = t - t_dose
                # Standard one-compartment with first-order absorption:
                # C(t) = (F * Dose * ka / (Vd * (ka - ke))) * (exp(-ke*dt) - exp(-ka*dt))
                if abs(ka - ke) < 1e-8:
                    # Avoid division by zero; use limiting form
                    term = dt * math.exp(-ke * dt)
                    c_t += (f * dose_mg / vd_L) * term
                else:
                    term = math.exp(-ke * dt) - math.exp(-ka * dt)
                    c_t += (f * dose_mg * ka / (vd_L * (ka - ke))) * term
            conc.append(max(c_t, 0.0))

        # --- Summary metrics ---
        cmax = max(conc) if conc else 0.0
        tmax = time_grid[conc.index(cmax)] if conc else 0.0

        # Trapezoidal AUC
        auc = 0.0
        for i in range(1, len(time_grid)):
            dt = time_grid[i] - time_grid[i - 1]
            auc += 0.5 * (conc[i] + conc[i - 1]) * dt

        cmin_ss = conc[-1] if conc else 0.0

        return {
            "model": "ONE_COMPARTMENT_FIRST_ORDER",
            "time_grid_h": time_grid,
            "concentration_mg_per_L": conc,
            "summary": {
                "cmax_mg_per_L": cmax,
                "tmax_h": tmax,
                "auc_mg_h_per_L": auc,
                "cmin_steady_state_mg_per_L": cmin_ss,
            },
            "parameters": {
                "dose_mg": dose_mg,
                "duration_h": duration_h,
                "dosing_interval_h": dosing_interval_h,
                "weight_kg": weight_kg,
                "vd_L_per_kg": vd_L_per_kg,
                "vd_L": vd_L,
                "clearance_L_per_h_per_kg": cl_L_per_h_per_kg,
                "clearance_L_per_h": cl_L_per_h,
                "ka_1_per_h": ka,
                "ke_1_per_h": ke,
                "bioavailability": f,
            },
            "constitutional": {
                "status": "SIMULATED",
                "engine": "PK_ONE_COMPARTMENT_V1",
                "notes": "Exposure-only PK simulation; not clinical; parameters deterministic with safe defaults."
            },
        }

    def materialize_candidate(self, worldline_id, coherence):
        """
        Materialize candidate by pulling pre-calculated toxicity from the WorldLine.
        Only materialize if the manifold has reached resonance.
        """
        if coherence < 0.80:
            return {"status": "VOID", "reason": "Insufficient Manifold Coherence"}

        # Pull actual pre-calculated toxicity from the physical_realization data in the WorldLine
        # Use worldline_id to find the file; path aligned with PX_Warehouse.WorldLine_Database
        task_id = worldline_id.replace("WL-", "")
        try:
            from PX_Warehouse.WorldLine_Database import DEFAULT_WORLDLINES_PATH
            wl_dir = DEFAULT_WORLDLINES_PATH
        except Exception:
            from pathlib import Path
            _root = Path(__file__).resolve().parents[1]
            wl_dir = str(_root / "PX_Warehouse" / "WorldLines")
        # Worldlines may be in root or in tier subdirs (Bronze, Silver, Gold, Diamond)
        wl_path = os.path.join(wl_dir, f"{task_id}.worldline")
        if not os.path.isfile(wl_path):
            for sub in ("Bronze", "Silver", "Gold", "Diamond"):
                cand = os.path.join(wl_dir, sub, f"{task_id}.worldline")
                if os.path.isfile(cand):
                    wl_path = cand
                    break
        
        try:
            with open(wl_path, "r") as f:
                data = json.load(f)
            
            phys = data.get("physical_realization", {})
            tox = phys.get("toxicity_index")
            affinity = phys.get("binding_affinity_kj")
            
            if tox is None:
                # Fail-fast if toxicity_index is missing
                raise ValueError(f"TraceabilityError: toxicity_index missing in {worldline_id}")

            return {
                "candidate_id": f"PX-REAL-{task_id[-6:]}",
                "binding_affinity_kj": affinity,
                "toxicity_index": tox,
                "status": "READY_FOR_SYNTHESIS",
                "timestamp": datetime.now(timezone.utc).isoformat(),
                "source_worldline": worldline_id
            }
        except FileNotFoundError:
            # Fallback for legacy or tests if file doesn't exist, but follow deterministic rules
            raise FileNotFoundError(f"TraceabilityError: WorldLine file not found for {worldline_id}")
        except Exception as e:
            raise e
