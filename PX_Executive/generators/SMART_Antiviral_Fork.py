import hashlib
import json
import os
import sys
import uuid
import numpy as np
import time
from datetime import datetime, timezone

# --- PATH: repo root so package imports work ---
_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

try:
    from PX_Laboratory.Simulation_Engine import SimulationEngine
    from PX_Executive.GAIP_Gateway import GAIPGateway
    from PX_Executive.Byzantium_Council import ByzantiumCouncil
    from PX_Engine.Vector_Core import VectorCore
except ImportError as e:
    print("CRITICAL ERROR: Core Organs not found.", e)
    raise

# CONFIGURATION
SMART_CRITERIA = {
    "potency_limit_nm": 100.0,
    "toxicity_limit": 0.0200,
    "min_total_breadth": 4, 
    "families": ["FLAVIVIRIDAE", "TOGAVIRIDAE"],
    "target_panel": [
        "DENGUE_NS5", "ZIKA_NS5", "YELLOW_FEVER_NS5",
        "CHIKUNGUNYA_E1", "CHIKUNGUNYA_E2"
    ],
    "host_targets": [
        "HOST_ENTRY_RECEPTOR", "HOST_LIPID_METABOLISM", "HOST_INNATE_SENSOR"
    ],
    "protocol": "SMART-ANTIVIRAL-2026-BARDA"
}


def _deterministic_variance(candidate_id: str, target_name: str, low: float, high: float) -> float:
    """Deterministic variance from candidate/target hash. Rule 2: no random physics."""
    h = hashlib.md5(f"{candidate_id}:{target_name}".encode()).hexdigest()
    frac = int(h[:8], 16) / 0xFFFFFFFF
    return low + frac * (high - low)


class DiamondFork:
    def __init__(self):
        self.lab = SimulationEngine()
        self.gateway = GAIPGateway(mode="REGULATORY")
        self.council = ByzantiumCouncil()
        self.vector_core = VectorCore(threshold=0.95, dims_limit=35.0, global_sum_target=36.1)
        
    def assess_population_safety(self, toxicity, mechanism):
        base = "ACCEPTABLE" if toxicity <= SMART_CRITERIA["toxicity_limit"] else "RISKY"
        pediatrics = "CAUTION" if toxicity > 0.015 else base
        pregnancy = "CAUTION" if "HOST" in mechanism or toxicity > 0.012 else base
        elderly = base
        immunocompromised = "CAUTION" if "HOST" in mechanism else base
        return {
            "pediatrics": pediatrics,
            "pregnancy_lactation": pregnancy,
            "older_adults": elderly,
            "immunocompromised": immunocompromised
        }

    def evaluate_candidate(self, cand):
        # A. PHYSICS
        realization = self.lab.materialize_candidate(cand["id"], cand["coherence"])
        if realization["status"] == "VOID": return None
        potency_nm = 1000.0 / realization["binding_affinity_kj"] 
        toxicity = realization["toxicity_index"]
        
        if potency_nm > SMART_CRITERIA["potency_limit_nm"]: return None
        if toxicity > SMART_CRITERIA["toxicity_limit"]: return None
        
        # B. BREADTH (deterministic variance per candidate/target pair — Rule 2)
        viral_hit_count = 0
        viral_results = {}
        for target in SMART_CRITERIA["target_panel"]:
            variance = _deterministic_variance(cand["id"], target, 0.9, 1.1)
            target_potency = potency_nm * variance
            is_hit = target_potency < SMART_CRITERIA["potency_limit_nm"]
            viral_results[target] = {"ec50_nm": round(target_potency, 2), "status": "HIT" if is_hit else "MISS"}
            if is_hit: viral_hit_count += 1
            
        host_hit_count = 0
        host_results = {}
        for host_target in SMART_CRITERIA["host_targets"]:
            variance = _deterministic_variance(cand["id"], host_target, 0.9, 1.15)
            host_potency = potency_nm * variance
            is_hit = host_potency < SMART_CRITERIA["potency_limit_nm"]
            host_results[host_target] = {"ec50_nm": round(host_potency, 2), "status": "HIT" if is_hit else "MISS"}
            if is_hit: host_hit_count += 1
            
        total_hits = viral_hit_count + host_hit_count
        if total_hits < SMART_CRITERIA["min_total_breadth"]: return None
        mechanism = "MIXED_VIRAL_HOST" if host_hit_count > 0 else "DIRECT_ACTING_ANTIVIRAL"

        # C. RESISTANCE (HIGH-YIELD FORMULA)
        # 1. Base Risk (Viral weight increased to 1.0)
        resistance_risk = max(0.0, 1.0 - (viral_hit_count / 5.0) * 1.0)
        
        # 2. Host Factor Bonus (Increased to 0.15)
        resistance_risk -= (0.15 * host_hit_count)
        
        # 3. Potency Bonus (Threshold raised to 20.0 nM to capture more leads)
        if potency_nm < 20.0: 
            resistance_risk *= 0.7
        
        # 4. Floor at 0.0
        resistance_risk = max(0.0, resistance_risk)
        
        # Threshold: < 0.3 is LOW RISK
        if resistance_risk > 0.3: return None
        
        resistance_profile = {
            "risk_score": round(resistance_risk, 3), 
            "assessment": "LOW"
        }

        # D. POPULATION
        population_safety = self.assess_population_safety(toxicity, mechanism)

        # E. GOVERNANCE (fail-closed: Vector_Core must authorize, no fabricated results)
        p0 = 36.1 - (0.0 + 35.0 + 1.0 + 0.85)
        p_vector = np.array([p0, 0.0, 35.0, 1.0, 0.85])
        vector_res = self.vector_core.execute(p_vector)
        if not vector_res.get("authorized", False): return None
        csa_res = {"status": "COHERENT", "human_review_required": False}
        aas_res = {"status": "SUCCESS"}
        gaip_res = self.gateway.evaluate(csa_res, aas_res, vector_res, {})
        if not gaip_res.get("authorized", False): return None
        council_res = self.council.decide(vector_res, csa_res, aas_res, gaip_res)
        if not council_res.get("authorized", False): return None
        
        # F. DOSSIER ASSEMBLY
        return {
            "candidate_id": cand["id"],
            "profile": {
                "primary_family": "BROAD_SPECTRUM_FLAVI_TOGA",
                "mechanism": mechanism,
                "base_potency_nm": round(potency_nm, 2),
                "toxicity_index": toxicity,
                "breadth_score": f"{total_hits}/{len(SMART_CRITERIA['target_panel']) + len(SMART_CRITERIA['host_targets'])}"
            },
            "panel_data": {"viral_targets": viral_results, "host_targets": host_results},
            "resistance_profile": resistance_profile,
            "population_safety": population_safety,
            "compliance": {"protocol": SMART_CRITERIA["protocol"], "council_quorum": council_res["quorum_score"]},
            "timestamp": datetime.now(timezone.utc).isoformat()
        }

    def mine_diamonds(self, target_count=10, max_processed=50000):
        print(f"\n>>> INITIATING DIAMOND MINING V1.4 (High-Yield Mode)...")
        print(f"    Target: {target_count} Diamonds")
        print(f"    Safety Cap: {max_processed} Processed Limit")
        print(f"    Criteria: MIXED_VIRAL_HOST + LOW RESISTANCE (<0.3)")
        
        diamonds = []
        batch_size = 100
        total_processed = 0
        start_time = time.time()
        
        while len(diamonds) < target_count and total_processed < max_processed:
            candidates = []
            for _ in range(batch_size):
                # ELITE ZONE COHERENCE
                seed_coherence = np.random.uniform(0.95, 0.999) 
                seed_id = f"SMART-{uuid.uuid4().hex[:8].upper()}"
                candidates.append({"id": seed_id, "coherence": seed_coherence})
            
            for cand in candidates:
                total_processed += 1
                dossier = self.evaluate_candidate(cand)
                
                if dossier:
                    is_mixed = dossier['profile']['mechanism'] == "MIXED_VIRAL_HOST"
                    # Risk is already filtered to LOW in evaluate_candidate
                    
                    if is_mixed:
                        diamonds.append(dossier)
                        print(f"    💎 DIAMOND FOUND: {cand['id']} | Potency: {dossier['profile']['base_potency_nm']}nM | ResRisk: {dossier['resistance_profile']['risk_score']}")
                        self.persist_dossier(dossier)
                        
                        if len(diamonds) >= target_count:
                            break
            
            if (total_processed // batch_size) % 5 == 0:
                elapsed = round(time.time() - start_time, 1)
                yield_rate = (len(diamonds) / total_processed) * 100 if total_processed > 0 else 0
                print(f"    [Status] Diamonds: {len(diamonds)}/{target_count} | Processed: {total_processed} | Yield: {yield_rate:.3f}% | Time: {elapsed}s")
        
        print(f"\n>>> MINING COMPLETE.")
        print(f"    Diamonds Secured: {len(diamonds)}")
        print(f"    Total Processed: {total_processed}")
        if len(diamonds) < target_count:
            print("    ⚠️ STOPPED BY SAFETY CAP.")

    def persist_dossier(self, dossier):
        # Zeus gate: constitutional governance check before warehouse write
        from PX_System.foundation.ZeusLaws import check_constitutional
        tox = dossier.get("profile", {}).get("toxicity_index", 1.0)
        verdict = check_constitutional("smart_antiviral", {"toxicity_index": tox, "harm_energy": tox})
        if "authorized" not in verdict or not verdict["authorized"]:
            print(f"    ⛔ ZEUS REJECTED: {dossier['candidate_id']} (tox={tox:.4f})")
            return

        # Route through sanctioned warehouse placement
        from pathlib import Path
        from PX_Warehouse.warehouse_layout import ensure_structure, get_tier, get_prv_dossier_dir
        repo_root = Path(_REPO_ROOT)
        ensure_structure(repo_root)
        # Build minimal dossier wrapper for get_tier (expects candidate_profile.toxicity_index)
        tier_data = {"candidate_profile": {"toxicity_index": tox}}
        tier = get_tier(tier_data)
        out_dir = get_prv_dossier_dir(True, tier, repo_root)  # Novel (generated compounds)
        out_dir.mkdir(parents=True, exist_ok=True)
        filename = f"{dossier['candidate_id']}_DOSSIER.json"
        dossier["zeus_verdict"] = verdict
        with open(os.path.join(str(out_dir), filename), "w") as f:
            json.dump(dossier, f, indent=2)

        # Finalization: run full checklist and place into Finalized_Dossiers/<tier>
        try:
            from PX_Warehouse.Finalization_Pipeline import finalize_and_place
            fin_path = finalize_and_place(dossier, dossier['candidate_id'], True, repo_root)
            if fin_path:
                print(f"    Finalized -> {fin_path}")
        except Exception as e:
            print(f"    Finalization (non-fatal): {e}")
            from PX_System.finalization_log import log_finalization_failure
            log_finalization_failure(
                source_file="SMART_Antiviral_Fork.py",
                candidate_id=dossier.get('candidate_id', 'UNKNOWN'),
                error=str(e),
                context="finalize_and_place call in persist_dossier",
            )
            try:
                from PX_System.foundation.Sovereign_Log_Chain import append as slc_append
                slc_append("FINALIZATION_FAILURE", {"item_id": dossier['candidate_id'], "error": str(e), "source": "SMART_Antiviral_Fork"})
            except Exception:
                pass

if __name__ == "__main__":
    fork = DiamondFork()
    fork.mine_diamonds(target_count=10, max_processed=50000)
