import os
import json
import hashlib
from datetime import datetime

# --- PATH CONFIGURATION ---
BASE_PATH = r"E:\foundation\PX_Warehouse"
COMMERCIAL_DIR = os.path.join(BASE_PATH, "CommercialAssets")
LOG_DIR = os.path.join(BASE_PATH, "logs")
ANALYTICS_FILE = os.path.join(LOG_DIR, "portfolio_analytics.json")
STABILITY_FILE = os.path.join(LOG_DIR, "portfolio_stability.json")

def calculate_integrity_hash(data):
    """Generate SHA-256 hash of the candidate data (excluding monetization)."""
    data_copy = data.copy()
    data_copy.pop("monetization", None)
    serialized = json.dumps(data_copy, sort_keys=True).encode('utf-8')
    return hashlib.sha256(serialized).hexdigest()

def generate_portfolio_analytics():
    print("Generating Portfolio Analytics & Stability Metrics...")
    
    # Load previous analytics for drift tracking
    previous_analytics = {}
    if os.path.exists(ANALYTICS_FILE):
        try:
            with open(ANALYTICS_FILE, 'r', encoding='utf-8') as f:
                previous_analytics = json.load(f)
        except Exception:
            pass

    analytics = {
        "timestamp": datetime.now().isoformat(),
        "grade_distribution": {"DIAMOND": 0, "GOLD": 0, "SILVER": 0, "BRONZE": 0, "REJECT": 0},
        "dominance_stats": {"mean": 0, "max": 0, "min": 1, "histogram": []},
        "integrity_report": {"total_checked": 0, "failures": []},
        "anomaly_detection": {"threshold_near_misses": 0, "high_volatility_candidates": []},
        "stability_metrics": {
            "grade_migration_drift": {},
            "dominance_volatility": 0,
            "threshold_sensitivity": {}
        }
    }
    
    scores = []
    current_assets = {} # {mol_id: {grade, score}}
    
    # Scan all assets in CommercialAssets
    for root, _, files in os.walk(COMMERCIAL_DIR):
        for filename in files:
            if filename.endswith(".json"):
                path = os.path.join(root, filename)
                try:
                    with open(path, 'r', encoding='utf-8') as f:
                        data = json.load(f)
                    
                    v3_data = data.get("predator_x_v3", {})
                    grade = v3_data.get("grade", "REJECT")
                    score = v3_data.get("dominance_score", 0)
                    stored_hash = v3_data.get("integrity_hash")
                    mol_id = data.get("metadata", {}).get("id") or data.get("metadata", {}).get("smiles")
                    
                    if mol_id:
                        current_assets[mol_id] = {"grade": grade, "score": score}

                    # 1. Grade Distribution
                    analytics["grade_distribution"][grade] = analytics["grade_distribution"].get(grade, 0) + 1
                    
                    # 2. Dominance Stats
                    if grade != "REJECT":
                        scores.append(score)
                        analytics["dominance_stats"]["max"] = max(analytics["dominance_stats"]["max"], score)
                        analytics["dominance_stats"]["min"] = min(analytics["dominance_stats"]["min"], score)
                    
                    # 3. Integrity Verification
                    if stored_hash:
                        analytics["integrity_report"]["total_checked"] += 1
                        current_hash = calculate_integrity_hash(data)
                        if current_hash != stored_hash:
                            analytics["integrity_report"]["failures"].append({
                                "path": path,
                                "expected": stored_hash,
                                "actual": current_hash
                            })
                            
                    # 4. Anomaly Detection: Threshold Near Misses
                    sm = v3_data.get("metrics", {}).get("sm", 0)
                    rr = v3_data.get("metrics", {}).get("rr", 0)
                    ba = v3_data.get("metrics", {}).get("ba", 9999)
                    
                    # Sensitivity Mapping: How many assets would change grade with a 5% shift in SM?
                    if 47.5 <= sm < 50: analytics["stability_metrics"]["threshold_sensitivity"]["sm_gold_diamond"] = analytics["stability_metrics"]["threshold_sensitivity"].get("sm_gold_diamond", 0) + 1
                    if 19.0 <= sm < 20: analytics["stability_metrics"]["threshold_sensitivity"]["sm_silver_gold"] = analytics["stability_metrics"]["threshold_sensitivity"].get("sm_silver_gold", 0) + 1
                    if 9.5 <= sm < 10: analytics["stability_metrics"]["threshold_sensitivity"]["sm_bronze_silver"] = analytics["stability_metrics"]["threshold_sensitivity"].get("sm_bronze_silver", 0) + 1

                except Exception as e:
                    print(f"Error analyzing {filename}: {e}")

    # 5. Stability Metrics: Grade Migration Drift
    if previous_analytics and "assets" in previous_analytics:
        prev_assets = previous_analytics["assets"]
        drift_count = 0
        volatility_sum = 0
        for mol_id, current in current_assets.items():
            if mol_id in prev_assets:
                prev = prev_assets[mol_id]
                if current["grade"] != prev["grade"]:
                    drift_key = f"{prev['grade']}->{current['grade']}"
                    analytics["stability_metrics"]["grade_migration_drift"][drift_key] = analytics["stability_metrics"]["grade_migration_drift"].get(drift_key, 0) + 1
                    drift_count += 1
                volatility_sum += abs(current["score"] - prev["score"])
        
        if len(current_assets) > 0:
            analytics["stability_metrics"]["dominance_volatility"] = volatility_sum / len(current_assets)

    if scores:
        analytics["dominance_stats"]["mean"] = sum(scores) / len(scores)
        for i in range(10):
            low = i / 10
            high = (i + 1) / 10
            count = len([s for s in scores if low <= s < high])
            analytics["dominance_stats"]["histogram"].append({"bin": f"{low}-{high}", "count": count})

    # Store current assets for next drift calculation (but keep it separate or in a hidden field to keep analytics file clean)
    analytics_to_save = analytics.copy()
    analytics_to_save["assets"] = current_assets

    with open(ANALYTICS_FILE, 'w', encoding='utf-8') as f:
        json.dump(analytics_to_save, f, indent=2)
    
    # Save a separate stability report for quick viewing
    with open(STABILITY_FILE, 'w', encoding='utf-8') as f:
        json.dump(analytics["stability_metrics"], f, indent=2)
    
    print(f"Analytics & Stability generated.")
    print(f"Integrity Check: {analytics['integrity_report']['total_checked']} files verified, {len(analytics['integrity_report']['failures'])} failures.")
    print(f"Stability: Drift events detected: {len(analytics['stability_metrics']['grade_migration_drift'])}")

if __name__ == "__main__":
    generate_portfolio_analytics()
